import argparse
import os
import pathlib
import tarfile
from math import ceil

import hgvs.assemblymapper
import hgvs.dataproviders.uta
import hgvs.parser
from biocommons.seqrepo import SeqRepo
from dotenv import load_dotenv

load_dotenv()

parser = argparse.ArgumentParser(description='Pack genome and RNA refseq data from seqrepo')

parser.add_argument('--uta_database_schema',
                    help='UTA database schema',
                    default=os.getenv('UTA_DATABASE_SCHEMA', 'uta_20240523b'))
parser.add_argument('--uta_database_url', help='UTA database URL',
                    # Use the biocommons UTA database if we don't specify a custom one.
                    default=os.getenv('UTA_DATABASE_URL', 'postgresql://anonymous:anonymous@uta.biocommons.org/uta'))
parser.add_argument('--seqrepo_dir',
                    help='Seqrepo directory',
                    default='/usr/local/share/seqrepo/latest')
parser.add_argument('--output_dir',
                    help='Output directory',
                    default='tmp')
args = parser.parse_args()

# Also, make sure the URL uses `postgresql` instead of `postgres` as schema
database_url = f"{args.uta_database_url}/{args.uta_database_schema}".replace('postgres://', 'postgresql://')

hgvsDataProvider = hgvs.dataproviders.uta.connect(db_url=database_url)

seq_repo = SeqRepo(args.seqrepo_dir)


def code_to_hex(code):
    """
    Convert sequence code to a 3 bit integer.
    Supports only 6 codes: A, C, G, T, U and a generic placeholder (designated N) for everything else.
    Details here: https://en.wikipedia.org/wiki/FASTA_format#Sequence_representation)

    Nucleic Acid Code | Meaning                          | Mnemonic
    A                 | A                                | Adenine
    C                 | C                                | Cytosine
    G                 | G                                | Guanine
    T                 | T                                | Thymine
    U                 | U                                | Uracil
    N                 | A C G T U                        | Nucleic acid
    """
    match code.upper():
        case 'A': return 0x0
        case 'C': return 0x1
        case 'G': return 0x2
        case 'T': return 0x3
        case 'U': return 0x4
        # All other codes including `N`s
        case _: return 0x7


def encode_8_codes_to_3_bytes(codes, byte_count, ba):
    # Pad the codes list with extra `N`s if its length is less than 8. This can happen at the end of a RefSeq and these
    # extra codes get discarded by the decoder which receives the original length of the RefSeq, thus knowing when to
    # stop decoding.
    input_codes_len = len(codes)
    if input_codes_len < 8:
        codes = codes + 'N' * (8 - len(codes))

    # We need 3 bits to represent each code and, since we need to pack them into bytes, some of them will end up being
    # split across 2 bytes. The pattern repeats after LCM(3,8) = 24 (3 bytes), when the next code will start at index 0
    # in the next byte.

    c1 = code_to_hex(codes[0])
    c2 = code_to_hex(codes[1]) << 3
    # Chop off the top bit of codes[2] since it will be stored in the next byte.
    c3 = (code_to_hex(codes[2]) << 6) & 0xff
    ba[byte_count + 0] = c3 | c2 | c1

    # Bail out if there are no more codes left for the next byte.
    if input_codes_len < 3:
        return

    c1 = code_to_hex(codes[2]) >> 2
    c2 = code_to_hex(codes[3]) << 1
    c3 = code_to_hex(codes[4]) << 4
    # Chop off the top 2 bits of codes[5] since they'll be stored in the next byte.
    c4 = (code_to_hex(codes[5]) << 7) & 0xff
    ba[byte_count + 1] = c4 | c3 | c2 | c1

    # Bail out if there are no more codes left for the next byte.
    if input_codes_len < 6:
        return

    c1 = code_to_hex(codes[5]) >> 1
    c2 = code_to_hex(codes[6]) << 2
    c3 = code_to_hex(codes[7]) << 5
    ba[byte_count + 2] = c3 | c2 | c1


def pack_ref_seq(ref_seq):
    # Packing the RefSeq using 3 bits per code will yield a packed RefSeq that's 3/8 times shorter.
    packed_len = ceil(len(ref_seq) * 3 / 8)
    ba = bytearray(packed_len)

    # Pack each sequence of 8 codes into 3 bytes.
    byte_count = 0
    for i in range(0, len(ref_seq), 8):
        encode_8_codes_to_3_bytes(ref_seq[i:i + 8], byte_count, ba)
        byte_count += 3

    return ba


def compress_ref_seq(path):
    with tarfile.open(f'{path}_refseq.tar.gz', 'w:gz') as tar:
        tar.add(path, arcname=os.path.basename(path))


def pack_ref_seq_to_file(acc, output_dir):
    # Coalesce the entire sequence into a string to load it all in memory.
    # Otherwise, the data packing will be very slow.
    seq = str(seq_repo[acc])

    with open(os.path.join(output_dir, f'{acc}_{len(seq)}.refseq'), 'wb') as binary_file:
        binary_file.write(pack_ref_seq(seq))


def pack_rnas(output_dir):
    ref_seq_dir = os.path.join(output_dir, 'rna')
    pathlib.Path(ref_seq_dir).mkdir(parents=True, exist_ok=True)

    uta_rnas = hgvsDataProvider._fetchall(r"""select distinct(split_part(ac, '.', 1) || '.' || regexp_replace(split_part(ac, '.', 2), '/\d+$', '')) from transcript where ac ~ '^NM_'""")
    for rna in [rna[0] for rna in uta_rnas]:
        pack_ref_seq_to_file(rna, ref_seq_dir)

    compress_ref_seq(ref_seq_dir)


def pack_chromosomes(output_dir, build):
    ref_seq_dir = os.path.join(output_dir, build)
    pathlib.Path(ref_seq_dir).mkdir(parents=True, exist_ok=True)

    for acc in [acc for acc in hgvsDataProvider.get_assembly_map(build) if acc.startswith('NC_')]:
        pack_ref_seq_to_file(acc, ref_seq_dir)

    compress_ref_seq(ref_seq_dir)


if __name__ == '__main__':
    pathlib.Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    print('Packing GRCh37 refseq...')
    pack_chromosomes(args.output_dir, 'GRCh37')

    print('Packing GRCh38 refseq...')
    pack_chromosomes(args.output_dir, 'GRCh38')

    print('Packing RNA refseq...')
    pack_rnas(args.output_dir)

    print('All done!')
