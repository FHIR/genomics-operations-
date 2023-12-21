from glob import glob
from math import floor
from os.path import isdir
from pathlib import Path
from sys import exit
from threading import Lock

from bioutils.normalize import NormalizationMode, normalize

# Make sure the refseq folder exists locally
if not isdir('refseq'):
    exit("Missing refseq folder. Please run fetch_refseq.sh!")


def hex_to_code(hex):
    """
    Convert a 3 bit integer to a sequence code. See `code_to_hex()` in pack_refseq.py for details.
    """
    match hex:
        case 0x0: return 'A'
        case 0x1: return 'C'
        case 0x2: return 'G'
        case 0x3: return 'T'
        case 0x7: return 'N'
        # Blow up if we get any unexpected hex-encoded code
        case _:
            raise NotImplementedError(f"unexpected hex-encoded code: '{hex}'")


class RefSeq:
    def __init__(self, packed_ref_seq, len):
        self.__packed_ref_seq = packed_ref_seq
        # Store the RefSeq length so we can easily tell when the index operator receives an out of bounds subscript
        # because we want to represent each code using 3 bits and we don't want to introduce an end-of-sequence code, in
        # case we'll want to leverage the unused codes for something else.
        self.__len = len

    def __len__(self):
        return self.__len

    def __getitem__(self, subscript):
        if isinstance(subscript, slice):
            codes = []
            # TODO: It should be more efficient to process 3 bytes at a time for long ranges.
            for i in range(*subscript.indices(self.__len)):
                codes.append(self[i])
            return "".join(codes)

        if subscript >= self.__len:
            raise IndexError(f"list index out of range: '{subscript}' must be less than '{self.__len}'")

        # Select the byte(s) and the offset of the code we wish to extract from the packed RefSeq data.
        packed_subscript = floor(subscript * 3 / 8)
        data = self.__packed_ref_seq[packed_subscript]
        match subscript % 8:
            case 0: offset = 0
            case 1: offset = 3
            case 2:
                offset = 6
                # This code is packed across two bytes.
                data += self.__packed_ref_seq[packed_subscript + 1] << 8
            case 3: offset = 1
            case 4: offset = 4
            case 5:
                offset = 7
                # This code is packed across two bytes.
                data += self.__packed_ref_seq[packed_subscript + 1] << 8
            case 6: offset = 2
            case 7: offset = 5

        # Move the relevant three bits to the bottom, extract them using a bit mask and convert them back to a code.
        return hex_to_code((data >> offset) & 0x7)


def get_ref_seq(acc):
    try:
        ref_seq_file_pattern = f'refseq/**/{acc}_*.refseq'
        found_files = glob(ref_seq_file_pattern)
        # TODO: Don't store NM accessions for each build since they're redundant
        if not found_files:
            raise ValueError(f'failed to find expected refseq file for accession "{acc}"')

        ref_seq_file = found_files[0]
        ref_seq_length = int(Path(ref_seq_file).stem.rpartition('_')[-1])

        with open(ref_seq_file, mode='rb') as file:
            refSeq = RefSeq(file.read(), ref_seq_length)
    except Exception as err:
        print(f"failed to read refseq file: {err=}, {type(err)=}")
        raise

    return refSeq


ref_seq_lock = Lock()


def get_ref_seq_subseq(acc, start, end):
    # Need to serialise this if we can't keep all the RefSeq data in memory
    with ref_seq_lock:
        return get_ref_seq(acc)[start:end]


def get_normalized_spdi(acc, pos, ref, alt):
    # Need to serialise this if we can't keep all the RefSeq data in memory
    with ref_seq_lock:
        ref_seq_fasta = get_ref_seq(acc)
        new_ival, new_alleles = normalize(ref_seq_fasta,
                                          interval=(pos, pos + len(ref)),
                                          alleles=(None, alt),
                                          mode=NormalizationMode.EXPAND,
                                          anchor_length=0,
                                          trim=ref != alt)
        return (f"{acc}:{new_ival[0]}:{new_alleles[0]}:{new_alleles[1]}")
