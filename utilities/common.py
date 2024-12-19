from enum import Enum
import re
import pandas as pd
import pymongo
import pyfastx
from threading import Lock

utilities_data_client_uri = "mongodb+srv://download:download@cluster0.8ianr.mongodb.net/UtilitiesData"
utilities_client = pymongo.MongoClient(utilities_data_client_uri)
utilities_db = utilities_client.UtilitiesData
transcript_data = utilities_db.Transcripts

GERMLINE = 'Germline'
SOMATIC = 'Somatic'
MIXED = 'Mixed'
SVs = {'INS', 'DEL', 'DUP', 'CNV', 'INV'}

# Dictionary of SO Codes for use in molecular consequence pipeline
codeDict = {'coding_sequence_variant ': ['http://sequenceontology.org', 'SO:0001580'],
            'chromosome': ['http://sequenceontology.org', 'SO:0000340'],
            'duplication': ['http://sequenceontology.org', 'SO:1000035'],
            'inversion': ['http://sequenceontology.org', 'SO:1000036'], 'inframe_insertion': ['http://sequenceontology.org', 'SO:0001821'], 'disruptive_inframe_insertion': ['http://sequenceontology.org', 'SO:0001824'], 'inframe_deletion': ['http://sequenceontology.org', 'SO:0001822'], 'disruptive_inframe_deletion': ['http://sequenceontology.org', 'SO:0001826'], 'downstream_gene_variant': ['http://sequenceontology.org', 'SO:0001632'], 'exon_variant': ['http://sequenceontology.org', 'SO:0001791'], 'exon_loss_variant': ['http://sequenceontology.org', 'SO:0001572'], 'frameshift_variant': ['http://sequenceontology.org', 'SO:0001589'], 'gene_variant': ['http://sequenceontology.org', 'SO:0001564'], 'feature_ablation': ['http://sequenceontology.org', 'SO:0001879'], 'gene_fusion': ['http://sequenceontology.org', 'SO:0001565'], 'bidirectional_gene_fusion': ['http://sequenceontology.org', 'SO:0002086'], 'rearranged_at_DNA_level': ['http://sequenceontology.org', 'SO:0000904'], 'intergenic_region': ['http://sequenceontology.org', 'SO:0000605'], 'conserved_intergenic_variant': ['http://sequenceontology.org', 'SO:0002017'], 'intragenic_variant': ['http://sequenceontology.org', 'SO:0002011'], 'intron_variant': ['http://sequenceontology.org', 'SO:0001627'], 'conserved_intron_variant': ['http://sequenceontology.org', 'SO:0002018'], 'miRNA': ['http://sequenceontology.org', 'SO:0000276'], 'missense_variant': ['http://sequenceontology.org', 'SO:0001583'], 'initiator_codon_variant': ['http://sequenceontology.org', 'SO:0001582'], 'stop_retained_variant': ['http://sequenceontology.org', 'SO:0001567'], 'protein_protein_contact': ['http://sequenceontology.org', 'SO:0001093'], 'structural_interaction_variant': ['http://sequenceontology.org', 'SO:0002093'], 'rare_amino_acid_variant': ['http://sequenceontology.org', 'SO:0002008'], 'splice_acceptor_variant': ['http://sequenceontology.org', 'SO:0001574'], 'splice_donor_variant': ['http://sequenceontology.org', 'SO:0001575'], 'splice_region_variant': ['http://sequenceontology.org', 'SO:0001630'], 'stop_lost': ['http://sequenceontology.org', 'SO:0001578'], '5_prime_UTR_premature_start_codon_gain_variant': ['http://sequenceontology.org', 'SO:0001988'], 'start_lost': ['http://sequenceontology.org', 'SO:0002012'], 'stop_gained': ['http://sequenceontology.org', 'SO:0001587'], 'synonymous_variant': ['http://sequenceontology.org', 'SO:0001819'], 'start_retained': ['http://sequenceontology.org', 'SO:0002019'], 'transcript_variant': ['http://sequenceontology.org', 'SO:0001576'], 'regulatory_region_variant': ['http://sequenceontology.org', 'SO:0001566'], 'upstream_gene_variant': ['http://sequenceontology.org', 'SO:0001631'], '3_prime_UTR_variant': ['http://sequenceontology.org', 'SO:0001624'], '3_prime_UTR_truncation': ['http://sequenceontology.org', 'SO:0002015'], 'exon_loss': ['http://sequenceontology.org', 'SO:0001572'], '5_prime_UTR_variant': ['http://sequenceontology.org', 'SO:0001623'], '5_prime_UTR_truncation': ['http://sequenceontology.org', 'SO:0002013'], 'sequence_feature': ['http://sequenceontology.org', 'SO:0000110'], 'non_coding_transcript_exon_variant': ['http://sequenceontology.org', 'SO:0001792'], 'conservative_inframe_insertion': ['http://sequenceontology.org', 'SO:0001823'], 'non_coding_transcript_variant': ['http://sequenceontology.org', 'SO:0001619'], 'conservative_inframe_deletion': ['http://sequenceontology.org', 'SO:0001825'],
            'start_retained_variant': ['http://sequenceontology.org', 'SO:0002019']}


class Genomic_Source_Class(Enum):

    @classmethod
    def set_(cls):
        return set(map(lambda c: c.value, cls))

    GERMLINE = GERMLINE
    SOMATIC = SOMATIC
    MIXED = MIXED


def validate_chrom_identifier(chrom):
    chrom = extract_chrom_identifier(chrom)
    pattern = '^[1-9]$|^1[0-9]$|^2[0-2]$|^[XYM]$'
    result = re.match(pattern, chrom)
    return bool(result)


def _get_chrom(chrom_index):
    switcher = {
        23: 'X',
        24: 'Y',
        25: 'M'
    }
    return switcher.get(chrom_index, str(chrom_index))


def extract_chrom_identifier(chrom):
    chrom = chrom.upper().replace("CHR", "")
    if chrom == "MT":
        chrom = "M"
    return chrom


def get_allelic_state(record, ratio_ad_dp, sample_position):
    allelic_state = ''
    allelic_code = ''
    try:
        allelic_frequency = float(record.aaf[0])
    except ZeroDivisionError:
        allelic_frequency = 0
    sample = record.samples[sample_position]
    alleles = sample.gt_alleles
    if record.CHROM != 'M':
        # gt_type: hom_ref = 0; het = 1; hom_alt = 2; uncalled = None
        if len(alleles) >= 2 and sample.gt_type == 1:
            allelic_state = 'heterozygous'
            allelic_code = 'LA6706-1'
        elif len(alleles) >= 2:
            allelic_state = 'homozygous'
            allelic_code = 'LA6705-3'
        elif sample.gt_type is not None and len(alleles) == 1:
            allelic_state = 'hemizygous'
            allelic_code = 'LA6707-9'
        else:
            _error_log_allelicstate(record)
    elif (sample.gt_type is not None and
          len(alleles) == 1 and
          alleles[0] == '1'):
        if allelic_frequency > ratio_ad_dp:
            allelic_state = "homoplasmic"
            allelic_code = "LA6704-6"
        else:
            allelic_state = "heteroplasmic"
            allelic_code = "LA6703-8"
    else:
        _error_log_allelicstate(record)
    return {
        'ALLELE': allelic_state,
        'CODE': allelic_code,
        'FREQUENCY': allelic_frequency
    }


# got this function from https://github.com/elimuinformatics/vcf2fhir/blob/master/vcf2fhir/common.py
def get_sequence_relation(phased_rec_map):
    Relation_table = pd.DataFrame(columns=['POS1', 'POS2', 'Relation'])
    for key in phased_rec_map:
        prev_record = None
        for record in phased_rec_map[key]:
            if prev_record is None:
                prev_record = record
                continue
            prev_data = prev_record.samples[0].data
            record_data = record.samples[0].data
            if (prev_data.PS == record_data.PS):
                if prev_data.GT == record_data.GT:
                    Relation_table = Relation_table.append(
                        {
                            'POS1': prev_record.POS,
                            'POS2': record.POS,
                            'Relation': 'Cis'
                        },
                        ignore_index=True
                    )
                else:
                    Relation_table = Relation_table.append(
                        {
                            'POS1': prev_record.POS,
                            'POS2': record.POS,
                            'Relation': 'Trans'
                        },
                        ignore_index=True
                    )
            prev_record = record
    return Relation_table


def query_genes(transcript_map):

    query_string = [{'$match': {}},  # Match all documents
                    {'$project': {'_id': 0, 'ncbiGeneSymbol': 0, 'featureType': 0, 'build37RefSeq': 0, 'build37Start': 0,
                                  'build37End': 0, 'build38RefSeq': 0, 'build38Start': 0, 'build38End': 0}}]

    result = list(transcript_data.aggregate(query_string))
    for res in result:
        transcript_map[res['transcriptRefSeq']] = res['MANE']


def _error_log_allelicstate(record):
    pass



# Fasta file handles cache
fasta_cache = {}
fasta_lock = Lock()

# File paths
BUILD37_FILE = 'FASTA/GCF_000001405.25_GRCh37.p13_genomic.fna'
BUILD38_FILE = 'FASTA/GCF_000001405.40_GRCh38.p14_genomic.fna'

# Dictionaries for GRCh37 and GRCh38 chromosomes
grch37_chromosomes = {
    "1": "NC_000001.10", "2": "NC_000002.11", "3": "NC_000003.11", "4": "NC_000004.11", "5": "NC_000005.9",
    "6": "NC_000006.11", "7": "NC_000007.13", "8": "NC_000008.10", "9": "NC_000009.11", "10": "NC_000010.10",
    "11": "NC_000011.9", "12": "NC_000012.11", "13": "NC_000013.10", "14": "NC_000014.8", "15": "NC_000015.9",
    "16": "NC_000016.9", "17": "NC_000017.10", "18": "NC_000018.9", "19": "NC_000019.9", "20": "NC_000020.10",
    "21": "NC_000021.8", "22": "NC_000022.10", "X": "NC_000023.10", "Y": "NC_000024.9", "MIT": "NC_012920.1"
}
grch38_chromosomes = {
    "1": "NC_000001.11", "2": "NC_000002.12", "3": "NC_000003.12", "4": "NC_000004.12", "5": "NC_000005.10",
    "6": "NC_000006.12", "7": "NC_000007.14", "8": "NC_000008.11", "9": "NC_000009.12", "10": "NC_000010.11",
    "11": "NC_000011.10", "12": "NC_000012.12", "13": "NC_000013.11", "14": "NC_000014.9", "15": "NC_000015.10",
    "16": "NC_000016.10", "17": "NC_000017.11", "18": "NC_000018.10", "19": "NC_000019.10", "20": "NC_000020.11",
    "21": "NC_000021.9", "22": "NC_000022.11", "X": "NC_000023.11", "Y": "NC_000024.10", "MIT": "NC_012920.1"
}

def get_fasta(file):
    """
    Retrieves the pyfastx.Fasta object for a given file, using a caching mechanism.

    Args:
        file (str): Path to the FASTA file.

    Returns:
        pyfastx.Fasta: An indexed FASTA object.
    """
    with fasta_lock:
        if file not in fasta_cache:
            try:
                fasta = pyfastx.Fasta(file)
            except Exception as err:
                print(f"Unexpected {err=}, {type(err)=}")
                raise
            fasta_cache[file] = fasta
        return fasta_cache[file]

def extract_chromosome_sequence(chromosome_tag, file_version):
    """
    Extracts the sequence of the specified chromosome using a cached FASTA file.

    Args:
        chromosome_tag (str): Chromosome tag (1-22, X, Y, or MIT).
        file_version (int): Version of the genome build (37 for GRCh37, 38 for GRCh38).

    Returns:
        str: The full sequence of the specified chromosome.
    """
    if file_version == 37:
        file_path = BUILD37_FILE
        chromosome = grch37_chromosomes.get(chromosome_tag)
    elif file_version == 38:
        file_path = BUILD38_FILE
        chromosome = grch38_chromosomes.get(chromosome_tag)
    else:
        raise ValueError("Invalid file version. Use 37 for GRCh37 or 38 for GRCh38.")

    if not chromosome:
        raise ValueError(f"Invalid chromosome tag: {chromosome_tag}. Must be 1-22, X, Y, or MIT.")

    fasta = get_fasta(file_path)
    if chromosome in fasta:
        return fasta[chromosome].seq
    else:
        raise ValueError(f"Chromosome {chromosome} not found in {file_path}")






# Function to validate chromosome identifiers
def validate_chrom_identifier(chrom):
    """
    Validates if the given chromosome identifier is valid (1-22, X, Y, M).

    Args:
        chrom (str): Chromosome identifier (e.g., 'chr1', 'X', 'M').

    Returns:
        bool: True if the identifier is valid, False otherwise.
    """
    chrom = extract_chrom_identifier(chrom)
    pattern = r'^(?:[1-9]|1[0-9]|2[0-2]|X|Y|M)$'
    return bool(re.match(pattern, chrom))


# Function to normalize chromosome identifiers
def extract_chrom_identifier(chrom):
    """
    Normalizes the chromosome identifier by removing prefixes and handling special cases.

    Args:
        chrom (str): Chromosome identifier (e.g., 'chr1', 'MT').

    Returns:
        str: Normalized chromosome identifier (e.g., '1', 'X', 'M').
    """
    chrom = chrom.upper().replace("CHR", "")
    if chrom == "MT":
        chrom = "M"
    return chrom


# Function to generate reference sequence mapping for all chromosomes
def generate_chrom_to_refseq(build):
    """
    Generates a mapping of chromosomes to reference sequence identifiers for the given build.

    Args:
        build (str): Reference genome build ('GRCh37' or 'GRCh38').

    Returns:
        dict: Mapping of chromosome identifiers to reference sequence identifiers.
    """
    if build == "GRCh37":
        base = "NC_0000"
        mitochondrial = "NC_012920.1"
    elif build == "GRCh38":
        base = "NC_0000"
        mitochondrial = "NC_012920.2"
    else:
        raise ValueError("Unsupported reference build. Use 'GRCh37' or 'GRCh38'.")

    # Autosomes 1-22
    chrom_to_refseq = {str(i): f"{base}{i:02}.12" for i in range(1, 23)}
    # Sex chromosomes and mitochondrial
    chrom_to_refseq["X"] = f"{base}23.12"
    chrom_to_refseq["Y"] = f"{base}24.12"
    chrom_to_refseq["M"] = mitochondrial

    return chrom_to_refseq


# Function to fetch reference sequence by chromosome
def _get_ref_seq_by_chrom(ref_build, chrom):
    """
    Returns the reference sequence identifier for the given chromosome and build.

    Args:
        ref_build (str): Reference genome build (e.g., 'GRCh37' or 'GRCh38').
        chrom (str): Chromosome identifier (e.g., 'chr1', 'X').

    Returns:
        str: Reference sequence identifier.

    Raises:
        ValueError: If the chromosome identifier is invalid or the build is unsupported.
    """
    chrom = extract_chrom_identifier(chrom)

    if not validate_chrom_identifier(chrom):
        raise ValueError(f"Invalid chromosome identifier: {chrom}")

    chrom_to_refseq = generate_chrom_to_refseq(ref_build)
    return chrom_to_refseq.get(chrom, "UNKNOWN")
