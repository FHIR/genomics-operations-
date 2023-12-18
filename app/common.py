from collections import OrderedDict
from threading import Lock
from uuid import uuid4
import pyliftover
from datetime import datetime
import pymongo
from flask import abort
from itertools import groupby
import re
from .input_normalization import normalize

# MongoDB Client URIs
FHIR_genomics_data_client_uri = "mongodb+srv://download:download@cluster0.8ianr.mongodb.net/FHIRGenomicsData"
utilities_data_client_uri = "mongodb+srv://download:download@cluster0.8ianr.mongodb.net/UtilitiesData"

# MongoDB Clients
client = pymongo.MongoClient(FHIR_genomics_data_client_uri)
utilities_client = pymongo.MongoClient(utilities_data_client_uri)

# Databases
db = client.FHIRGenomicsData
utilities_db = utilities_client.UtilitiesData

# Collections
patients_db = db.Patients
variants_db = db.Variants
tests_db = db.Tests
genotypes_db = db.Genotypes
dxImplication_db = db.dxImplication
txImplication_db = db.txImplication

beds = db.BEDs
phase_data = db.PhaseData

chromosomes_data = utilities_db.Chromosomes
genes_data = utilities_db.Genes
transcripts_data = utilities_db.Transcripts
exons_data = utilities_db.Exons
proteins_data = utilities_db.Proteins

# Liftover utilities
liftover_cache = {}
liftover_lock = Lock()


def get_liftover(from_db, to_db):
    key = from_db + "-" + to_db
    with liftover_lock:
        if key not in liftover_cache:
            liftover = pyliftover.LiftOver(from_db, to_db)
            liftover_cache[key] = liftover
        return liftover_cache[key]


# Constants
SVTYPE_TO_DNA_CHANGE_TYPE = {
    'CNV': ['SO:0001019', 'copy_number_variation'],
    'DUP': ['SO:1000035', 'duplication'],
    'INV': ['SO:1000036', 'inversion'],
    'DEL': ['SO:0000159', 'deletion'],
    'INS': ['SO:0000667', 'insertion']
}

DNA_CHANGE_TYPE_TO_CODE = {
    'SNV': 'SO:0001483',
    'MNV': 'SO:0002007',
    'INS': 'SO:0000667',
    'DEL': 'SO:0000159',
    'INV': 'SO:1000036',
    'DUP': 'SO:0001019',
    'CNV': 'SO:0001019',
    'delins': 'SO:1000032',
    'copy_number_variation': 'SO:0001019',
    'transcription_variant': 'SO:0001549',
    'haplotype': 'SO:0001024'
}

GENOMIC_SOURCE_CLASS_TO_CODE = {
    'GERMLINE': 'LA6683-2',
    'SOMATIC': 'LA6684-0'
}

GENOMIC_BUILD_TO_CODE = {
    'GRCh37': 'LA14029-5',
    'GRCh38': ' LA26806-2'
}

ALLELIC_STATE_TO_CODE = {
    'HETEROZYGOUS': 'LA6706-1',
    'HOMOZYGOUS': 'LA6705-3',
    'HEMIZYGOUS': 'LA6707-9'
}

CLIN_SIG_TO_CODE = {
    'Pathogenic': 'LA6668-3',
    'Likely pathogenic': 'LA26332-9',
    'Uncertain significance': 'LA26333-7',
    'Likely benign': 'LA26334-5',
    'Benign': 'LA6675-8'
}

OPERATOR_CODE_TO_OPERATOR = {
    'eq': '$eq',
    'ne': '$ne',
    'lt': '$lt',
    'gt': '$gt',
    'ge': '$gte',
    'le': '$lte',
    'sa': '$gte',
    'eb': '$lte'
}

SUPPORTED_GENE_SYSTEM_URLS = r'^https?:\/\/www\.genenames.org\/geneId$'
# SUPPORTED_FEATURE_CONSEQUENCE_SYSTEM_URLS = r'^http?:\/\/www\.sequenceontology.org\/$'
SUPPORTED_FEATURE_CONSEQUENCE_SYSTEM_URLS = 'http://sequenceontology.org'

SUPPORTED_DATE_FORMAT = '%Y-%m-%d'

SUPPORTED_GENOMIC_SOURCE_CLASSES = ['germline', 'somatic']

# Utility Functions


def key_func(k):
    return k['CHROMOSOME']['RefSeq']


def is_overlaping(a, b):
    if b['RANGE']['L'] >= a['RANGE']['L'] and b['RANGE']['L'] <= a['RANGE']['H']:
        return True
    else:
        return False


def merge(arr):
    arr.sort(key=lambda x: x['RANGE']['L'])
    merged_list = []
    merged_list.append(arr[0])

    for i in range(1, len(arr)):
        pop_element = merged_list.pop()

        if is_overlaping(pop_element, arr[i]):
            pop_element['RANGE']['H'] = max(pop_element['RANGE']['H'], arr[i]['RANGE']['H'])
            merged_list.append(pop_element)
        else:
            merged_list.append(pop_element)
            merged_list.append(arr[i])

    return merged_list


def merge_ranges(ranges):
    ranges = sorted(ranges, key=key_func)

    merged_ranges = []

    for key, value in groupby(ranges, key_func):
        merged_ranges.extend(merge(list(value)))

    return merged_ranges


def validate_subject(patient_id):
    if not patients_db.find_one({"patientID": patient_id}):
        abort(400, f"Patient ({patient_id}) not found.")


def get_variant(variant):
    variant = variant.strip()
    if variant.count(":") not in [1, 3]:
        abort(400, f'variant ({variant}) is not in the correct format(SPDI|HGVS)')

    variant = variant.lstrip()

    if variant.count(":") == 1:  # HGVS expression
        SPDIs = normalize(variant)
        if not SPDIs:
            abort(400, f'Cannot normalize variant: {variant}')
        elif not SPDIs["GRCh37SPDI"] and not SPDIs["GRCh38SPDI"]:
            abort(400, f'Cannot normalize variant: {variant}')
        else:
            normalized_variant = {"variant": variant, "GRCh37": SPDIs["GRCh37SPDI"], "GRCh38": SPDIs["GRCh38SPDI"]}

    elif variant.count(":") == 3:  # SPDI expression
        SPDIs = normalize(variant)
        if not SPDIs:
            abort(400, f'Cannot normalize variant: {variant}')
        elif not SPDIs["GRCh37SPDI"] and not SPDIs["GRCh38SPDI"]:
            abort(400, f'Cannot normalize variant: {variant}')
        else:
            normalized_variant = {"variant": variant, "GRCh37": SPDIs["GRCh37SPDI"], "GRCh38": SPDIs["GRCh38SPDI"]}
    else:
        abort(400, f'variant ({variant}) is not in the correct format(SPDI|HGVS)')

    return normalized_variant


def get_gene(gene):
    gene = gene.strip()
    gene_return = {'isSystem': False, 'gene': gene, 'system': None}
    if "|" in gene:
        if gene.count("|") == 1 and (re.match(SUPPORTED_GENE_SYSTEM_URLS, gene.rsplit('|')[0])):
            gene_system_url = gene.rsplit("|")[0]
            gene = gene.rsplit("|")[1]
            gene_return['isSystem'] = True
            gene_return['gene'] = gene
            gene_return['system'] = gene_system_url
        else:
            abort(400, f'gene ({gene}) is not in the correct format(codesystem|code)')

    return gene_return


def get_feature_consequence(feature_consequence):
    feature_consequence = feature_consequence.strip()
    feature_consequence_return = {'isSystem': False, 'feature_consequence': feature_consequence, 'system': None}
    if "|" in feature_consequence:
        # if feature_consequence.count("|") == 1 and (re.match(SUPPORTED_FEATURE_CONSEQUENCE_SYSTEM_URLS, feature_consequence.rsplit('|')[0])):
        if feature_consequence.count("|") == 1 and feature_consequence.rsplit('|')[0] == 'http://sequenceontology.org':
            feature_consequence_system_url = feature_consequence.rsplit("|")[0]
            feature_consequence = feature_consequence.rsplit("|")[1]
            feature_consequence_return['isSystem'] = True
            feature_consequence_return['feature_consequence'] = feature_consequence
            feature_consequence_return['system'] = feature_consequence_system_url
        else:
            abort(400, f'feature_consequence ({feature_consequence}) is not in the correct format(codesystem|code)')

    return feature_consequence_return


def get_haplotype(haplotype):
    haplotype = haplotype.strip()
    haplotype_return = {'isSystem': False, 'haplotype': haplotype, 'system': None}
    if "|" in haplotype:
        if haplotype.count("|") == 1:
            haplotype_system_url = haplotype.rsplit("|")[0]
            haplotype = haplotype.rsplit("|")[1]
            haplotype_return['isSystem'] = True
            haplotype_return['haplotype'] = haplotype
            haplotype_return['system'] = haplotype_system_url
        else:
            abort(400, f'haplotype ({haplotype}) is not in the correct format(codesystem|code)')

    return haplotype_return


def get_treatment(treatment):
    treatment = treatment.strip()
    treatment_return = {'isSystem': False, 'treatment': (treatment if not is_int(treatment) else int(treatment)), 'system': None}
    if "|" in treatment:
        if treatment.count("|") == 1:
            treatment_system_url = treatment.rsplit("|")[0]
            treatment = int(treatment.rsplit("|")[1])
            treatment_return['isSystem'] = True
            treatment_return['treatment'] = treatment
            treatment_return['system'] = treatment_system_url
        else:
            abort(400, f'treatment ({treatment}) is not in the correct format(codesystem|code)')

    return treatment_return


def get_condition(condition):
    condition = condition.strip()
    condition_return = {'isSystem': False, 'condition': condition, 'system': None}
    if "|" in condition:
        if condition.count("|") == 1:
            condition_system_url = condition.rsplit("|")[0]
            condition = condition.rsplit("|")[1]
            condition_return['isSystem'] = True
            condition_return['condition'] = condition
            condition_return['system'] = condition_system_url
        else:
            abort(400, f'condition ({condition}) is not in the correct format(codesystem|code)')

    return condition_return


def get_range(range):
    range = range.strip()

    ref_seq = range.split(':')[0]
    chrom = get_build_and_chrom_by_ref_seq(ref_seq)

    if chrom is None:
        abort(400, f'RefSeq({ref_seq}) is not valid')

    _range = list(map(int, range.split(':')[1].split('-')))

    if _range[0] > _range[1]:
        abort(400, f'Range start({_range[0]}) should be less than or equal to Range end({_range[1]})')

    chromosome = {'CHROM': chrom["chrom"], 'RefSeq': ref_seq}
    _range = {'L': _range[0], 'H': _range[1]}

    return {'CHROMOSOME': chromosome, 'BUILD': chrom["build"], 'RANGE': _range}


def get_lift_over_range(ranges):
    ranges_to_add = []
    for range in ranges:
        rse_other_build = lift_over(range['CHROMOSOME']['RefSeq'], range['RANGE']['L'], range['RANGE']['H'])
        if rse_other_build is not None:
            other_genomic_build = get_other_build(range["BUILD"])
            other_ref_seq = get_ref_seq_by_chrom_and_build(other_genomic_build, range['CHROMOSOME']['CHROM'])
            chromosome = {'CHROM': range['CHROMOSOME']['CHROM'], 'RefSeq': other_ref_seq}
            _range = {'L': rse_other_build["start"], 'H': rse_other_build["end"]}
            ranges_to_add.append({'CHROMOSOME': chromosome, 'BUILD': other_genomic_build, 'RANGE': _range})
    ranges.extend(ranges_to_add)


def get_variants(ranges, query):
    variants = []
    chromosome_to_ranges = []

    for range_ in ranges:
        temp_dict = {}

        chromosome = range_['CHROMOSOME']
        provided_genomic_build = range_['BUILD']
        _range = range_['RANGE']

        temp_dict["CHROM"] = chromosome['CHROM']
        temp_dict["RefSeq"] = chromosome['RefSeq']
        temp_dict["PGB"] = {"RefSeq": chromosome['RefSeq'], "BUILD": provided_genomic_build, "L": _range["L"], "H": _range["H"]}

        chromosome_to_ranges.append(temp_dict)

    for chrom in chromosome_to_ranges:
        variant_q = []
        genomic_builds = [chrom["PGB"]["BUILD"]]

        query["$and"] = []
        query["$and"].append({"SVTYPE": {"$exists": False}})
        query["$and"].append({"$expr": {"$gte": [{"$subtract": [{"$add": [{"$strLenCP": "$REF"}, "$POS"]}, "$POS"]}, 1]}})
        query["$and"].append({"$or": [
            {
                "$and": [
                    {"POS": {"$lte": chrom["PGB"]["L"]}},
                    {"$expr": {"$gt": [{"$add": [{"$strLenCP": "$REF"}, "$POS"]}, chrom["PGB"]["L"]]}}
                ]
            },
            {
                "$and": [
                    {"POS": {"$gte": chrom["PGB"]["L"]}},
                    {"POS": {"$lt": chrom["PGB"]["H"]}}
                ]
            }
        ]})
        query["$and"].append({"CHROM": {"$eq": chrom["CHROM"]}})

        query["genomicBuild"] = {"$in": genomic_builds}

        try:
            variant_q = variants_db.aggregate([{"$match": query}])
            variant_q = list(variant_q)
        except Exception as e:
            print(f"DEBUG: Error{e} under find_subject_variants query={query}")
            variant_q = []

        for variant in variant_q:
            if "SPDI" in variant:
                variants.append({'BUILD': chrom["PGB"]["BUILD"], 'SPDI': variant["SPDI"]})

    del query["genomicBuild"]
    del query["$and"]

    return variants


def get_date(dateRange):
    dateRange = dateRange.strip()

    operator_code = dateRange[:2]

    operator = OPERATOR_CODE_TO_OPERATOR[operator_code]
    date_ = datetime.strptime(dateRange[2:], SUPPORTED_DATE_FORMAT)

    return {'OPERATOR': operator, 'DATE': date_}


def get_genomic_source_class(genomic_source_class):
    genomic_source_class = genomic_source_class.strip()

    if (genomic_source_class.lower() not in SUPPORTED_GENOMIC_SOURCE_CLASSES):
        abort(400, "Genomic Source Class must be either 'germline' or 'somatic'")

    return genomic_source_class.lower()


def is_int(string):
    try:
        int(string)
        return True
    except ValueError:
        return False


def get_genomics_build_presence(query):
    genomics_build_presence = {"GRCh37": False, "GRCh38": False}

    query["genomicBuild"] = "GRCh37"
    if variants_db.find_one(query):
        genomics_build_presence["GRCh37"] = True

    query["genomicBuild"] = "GRCh38"
    if variants_db.find_one(query):
        genomics_build_presence["GRCh38"] = True

    query.pop("genomicBuild")

    return genomics_build_presence


def get_genomics_build_presence_tests_db(query):
    genomics_build_presence = {"GRCh37": False, "GRCh38": False}

    query["genomicBuild"] = "GRCh37"
    if tests_db.find_one(query):
        genomics_build_presence["GRCh37"] = True

    query["genomicBuild"] = "GRCh38"
    if tests_db.find_one(query):
        genomics_build_presence["GRCh38"] = True

    query.pop("genomicBuild")

    return genomics_build_presence


def get_chromosome_to_ranges(ranges, genomics_build_presence):
    chromosome_to_ranges = []

    if ranges:
        for range_ in ranges:
            temp_dict = {}

            chromosome = range_['CHROMOSOME']
            _range = range_['RANGE']

            provided_genomic_build = get_build_and_chrom_by_ref_seq(chromosome['RefSeq'])["build"]
            other_genomic_build = get_other_build(provided_genomic_build)

            temp_dict["CHROM"] = chromosome['CHROM']
            temp_dict["RefSeq"] = chromosome['RefSeq']
            temp_dict["PGB"] = {"RefSeq": chromosome['RefSeq'], "BUILD": provided_genomic_build, "L": _range["L"], "H": _range["H"]}

            if (genomics_build_presence is None or genomics_build_presence[other_genomic_build]):
                rse_other_build = lift_over(chromosome['RefSeq'], _range["L"], _range["H"])

                if rse_other_build is None:
                    abort(422, f"Failed LiftOver ({chromosome['RefSeq']}:{_range['L']}-{_range['H']})")
                else:
                    other_ref_seq = get_ref_seq_by_chrom_and_build(other_genomic_build, chromosome['CHROM'])
                    temp_dict["OGB"] = {"RefSeq": other_ref_seq, "BUILD": other_genomic_build, "L": rse_other_build["start"], "H": rse_other_build["end"]}

            else:
                temp_dict["OGB"] = None

            chromosome_to_ranges.append(temp_dict)

    return chromosome_to_ranges


def get_dna_chg(svtype):
    dna_chg = SVTYPE_TO_DNA_CHANGE_TYPE.get(svtype.upper())
    return {"CODE": dna_chg[0], "DISPLAY": dna_chg[1]}


def create_fhir_variant_resource(record, ref_seq, subject):
    vid = f"dv-{str(record['_id'])}"

    resource = OrderedDict()
    resource["resourceType"] = "Observation"
    resource["id"] = vid
    resource["meta"] = {"profile": [
                        "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/variant"]}
    resource["status"] = "final"
    resource["category"] = [{"coding": [{"system": "http://terminology.hl7.org/CodeSystem/observation-category",
                                         "code": "laboratory"}]}]
    resource["code"] = {"coding": [{"system": "http://loinc.org",
                                   "code": "69548-6",
                                    "display": "Genetic variant assessment"}]}
    resource["subject"] = {"reference": f"Patient/{subject}"}
    resource["valueCodeableConcept"] = {"coding": [{"system": "http://loinc.org",
                                                    "code": "LA9633-4",
                                                    "display": "present"}]}

    resource["component"] = []

    if 'SVTYPE' in record:
        dna_chg = get_dna_chg(record['SVTYPE'])
        resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                           "code": "48019-4",
                                                           "display": "DNA Change Type"}]},
                                      "valueCodeableConcept": {"coding": [{"system": "http://sequenceontology.org",
                                                                           "code": f"{dna_chg['CODE']}",
                                                                           "display": f"{dna_chg['DISPLAY']}"}]}})

    # Genomic Source Class
    if 'genomicSourceClass' in record and record['genomicSourceClass'].upper() in ['GERMLINE', 'SOMATIC']:
        resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                           "code": "48002-0",
                                                           "display": "Genomic Source Class"}]},
                                      "valueCodeableConcept": {"coding": [{"system": "http://loinc.org",
                                                                           "code": f"{GENOMIC_SOURCE_CLASS_TO_CODE[record['genomicSourceClass'].upper()]}",
                                                                           "display": f"{record['genomicSourceClass']}"}]}})

    # Genomic reference sequence ID
    resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                       "code": "48013-7",
                                                       "display": "Genomic reference sequence ID"}]},
                                  "valueCodeableConcept": {"coding": [{"system": "http://www.ncbi.nlm.nih.gov/nuccore",
                                                                       "code": ref_seq}]}})

    # Allelic State
    if "allelicState" in record and (('SVTYPE' not in record) or ('SVTYPE' in record and 'genomicSourceClass' in record and record['SVTYPE'] in ['DUP', 'DEL', 'INV', 'INS'] and record['genomicSourceClass'].lower() == 'germline')):
        resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                           "code": "53034-5",
                                                           "display": "Allelic state"}]},
                                      "valueCodeableConcept": {"coding": [{"system": "http://loinc.org",
                                                                           "code": f"{ALLELIC_STATE_TO_CODE[record['allelicState'].upper()]}",
                                                                           "display": record["allelicState"]}]}})

    # Variation Code
    if "SPDI" in record:
        resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                           "code": "81252-9",
                                                           "display": "Discrete genetic variant"}]},
                                      "valueCodeableConcept": {"coding": [{"system": "https://api.ncbi.nlm.nih.gov/variation/v0/",
                                                                           "code": record["SPDI"],
                                                                           "display": record["SPDI"]}]}})

    # Allelic Frequency
    if ('SVTYPE' not in record):
        try:
            allelic_frequency = record["allelicFrequency"]
            resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                               "code": "81258-6",
                                                               "display": "Sample VAF"}]},
                                          "valueQuantity": {"value": allelic_frequency,
                                                            "unit": "relative frequency of a particular allele in the specimen",
                                                            "system": "http://unitsofmeasure.org",
                                                            "code": "1"}})
        except Exception as e:
            print(f"DEBUG: Error calculating Allelic Frequency - {e}")

    # Copy Number
    if (('SVTYPE' in record) and (record['SVTYPE'] in ['CNV', 'DUP', 'DEL']) and 'CN' in record):
        resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                           "code": "82155-3",
                                                           "display": "Genomic Structural Variant copy Number"}]},
                                      "valueQuantity": {"value": int(f"{record['CN']}"),
                                                        "system": "http://unitsofmeasure.org",
                                                        "code": "1"}})

    # Genomic Ref Allele ID
    resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                           "code": "69547-8",
                                                       "display": "Genomic Ref allele [ID]"}]},
                                  "valueString": f"{record['REF']}"})

    # Genomic Alt Allele ID
    if ((('SVTYPE' not in record) or ('SVTYPE' in record and record['SVTYPE'] in ['INS'])) and ('ALT' in record) and (record['ALT'].isalpha())):
        alt = record['ALT']
        resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                           "code": "69551-0",
                                                           "display": "Genomic Alt allele [ID]"}]},
                                      "valueString": f"{alt}"})

    # Genomic Coord System
    resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                       "code": "92822-6",
                                                       "display": "Genomic coord system"}]},
                                  "valueCodeableConcept": {"coding": [{"system": "http://loinc.org",
                                                                       "code": "LA30100-4",
                                                                       "display": "0-based interval counting"}]}})

    # Variant Exact Start and End
    if 'SVTYPE' not in record:
        exact_start = record["POS"]
        resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                           "code": "81254-5",
                                                           "display": "Variant exact start-end"}]},
                                      "valueRange": {"low": {"value": exact_start}}})

    # Variant Outer/Inner Start and End
    if 'SVTYPE' in record:
        inner_start = outer_start = record['POS']
        inner_end = outer_end = record['END']
        if 'CIPOS' in record:
            inner_start = record['POS'] + abs(record['CIPOS'][1])
            outer_start = record['POS'] - abs(record['CIPOS'][0])
        if 'CIEND' in record:
            inner_end = record['END'] - abs(record['CIEND'][0])
            outer_end = record['END'] + abs(record['CIEND'][1])

        resource["component"].append({"code": {"coding": [{"system": "http://loinc.org", "code": "81301-4", "display": "Variant outer start-end"}]}, "valueRange": {"low": {"value": outer_start}, "high": {"value": outer_end}}})

        resource["component"].append({"code": {"coding": [{"system": "http://loinc.org", "code": "81302-2", "display": "Variant inner start-end"}]}, "valueRange": {"low": {"value": inner_start}, "high": {"value": inner_end}}})

    # Variant population allele frequency
    if 'popAlleleFreq' in record:
        resource["component"].append({"code": {"coding": [{
            "system": "http://loinc.org",
            "code": "92821-8",
            "display": "Population allele frequency"
        }]},
            "valueQuantity": {"value": record['popAlleleFreq'],
                              "system": "http://unitsofmeasure.org",
                              "code": "1"}})
    return resource


def create_dx_implication_profile(implication, subject, vids):
    resource = OrderedDict()
    resource["resourceType"] = "Observation"
    resource["id"] = "dv-" + str(implication['_id'])
    resource["meta"] = {"profile": [
                        "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/diagnostic-implication"]}
    resource["status"] = "final"
    resource["category"] = [{"coding": [{"system": "http://terminology.hl7.org/CodeSystem/observation-category",
                                         "code": "laboratory"}]}]
    resource["code"] = {"coding": [{"system": "http://hl7.org/fhir/uv/genomics-reporting/CodeSystem/tbd-codes-cs",
                                   "code": "diagnostic-implication"}]}
    resource["subject"] = {"reference": f"Patient/{subject}"}
    if len(vids) > 0:
        resource["derivedFrom"] = []
        for vid in vids:
            resource["derivedFrom"].append({"reference": f"Observation/dv-{vid}"})

    resource["component"] = []

    if implication['clinicalSignificance'] in CLIN_SIG_TO_CODE.keys():
        resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                           "code": "53037-8",
                                                           "display": "Genetic variation clinical significance"}]},
                                      "valueCodeableConcept": {"coding": [{"system": "http://loinc.org",
                                                                           "code": f"{CLIN_SIG_TO_CODE[implication['clinicalSignificance']]}",
                                                                           "display": f"{implication['clinicalSignificance']}"}]}})
    else:
        resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                           "code": "53037-8",
                                                           "display": "Genetic variation clinical significance"}]},
                                      "valueCodeableConcept": {"text": f"{implication['clinicalSignificance']}"}})

    for predicted_phenotype in implication["predictedPhenotype"]:
        resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                           "code": "81259-4",
                                                           "display": "predicted phenotype"}]},
                                      "valueCodeableConcept": {"coding": [{"system": f"{predicted_phenotype['system']}",
                                                                           "code": f"{predicted_phenotype['code']}",
                                                                           "display": f"{predicted_phenotype['display']}"}]}})

    resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                       "code": "93044-6",
                                                       "display": "Level of evidence"}]},
                                  "valueCodeableConcept": {"text": f"{implication['evidenceLevel']}"}})

    return resource


def create_molecular_consequence_profile(molecular_consequence, subject, vids):
    resource = OrderedDict()
    resource["resourceType"] = "Observation"
    resource["id"] = "dv-" + str(molecular_consequence['_id'])
    resource["meta"] = {"profile": [
                        "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/molecular-consequence"]}
    resource["status"] = "final"
    resource["category"] = [{"coding": [{"system": "http://terminology.hl7.org/CodeSystem/observation-category",
                                         "code": "laboratory"}, {"system": "http://terminology.hl7.org/CodeSystem/v2-0074", "code": "GE"}]}]
    resource["code"] = {"coding": [{"system": "http://hl7.org/fhir/uv/genomics-reporting/CodeSystem/tbd-codes-cs",
                                   "code": "molecular-consequences"}]}
    resource["subject"] = {"reference": f"Patient/{subject}"}
    if len(vids) > 0:
        resource["derivedFrom"] = []
        for vid in vids:
            resource["derivedFrom"].append({"reference": f"Observation/dv-{vid}"})

    resource["component"] = []

    if 'cHGVS' in molecular_consequence:
        resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                           "code": "48004-6",
                                                           "display": "DNA change (c.HGVS)"}]},
                                      "valueCodeableConcept": {"text": f"{molecular_consequence['cHGVS']}"}})

    if 'transcriptRefSeq' in molecular_consequence:
        resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                           "code": "51958-7",
                                                           "display": "Reference Transcript"}]},
                                      "valueCodeableConcept": {"text": f"{molecular_consequence['transcriptRefSeq']}"}})

    if 'pHGVS' in molecular_consequence:
        resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                           "code": "48005-3",
                                                           "display": "Protein (Amino Acid) Change - pHGVS"}]},
                                      "valueCodeableConcept": {"text": f"{molecular_consequence['pHGVS']}"}})

    for feature_consequence in molecular_consequence["featureConsequence"]:
        resource["component"].append({"code": {"coding": [{"system": "http://hl7.org/fhir/uv/genomics-reporting/CodeSystem/tbd-codes-cs",
                                                           "code": "feature-consequence",
                                                           "display": "Feature Consequence"}]},
                                      "valueCodeableConcept": {"coding": [{"system": f"{feature_consequence['system']}",
                                                                           "code": f"{feature_consequence['code']}",
                                                                           "display": f"{feature_consequence['display']}"}]}})

    if 'impact' in molecular_consequence:
        resource["component"].append({"code": {"coding": [{"system": "http://hl7.org/fhir/uv/genomics-reporting/CodeSystem/tbd-codes-cs",
                                                           "code": "functional-effect",
                                                           "display": "Functional Effect"}]},
                                      "valueCodeableConcept": {"text": f"{molecular_consequence['impact']}"}})

    return resource


def create_tx_implication_profile_civic(implication, subject, vids):
    resource = OrderedDict()
    resource["resourceType"] = "Observation"
    resource["id"] = "dv-" + str(implication['_id'])
    resource["meta"] = {"profile": [
                        "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/therapeutic-implication"]}
    resource["status"] = "final"
    resource["category"] = [{"coding": [{"system": "http://terminology.hl7.org/CodeSystem/observation-category",
                                         "code": "laboratory"}]}]
    resource["code"] = {"coding": [{"system": "http://hl7.org/fhir/uv/genomics-reporting/CodeSystem/tbd-codes-cs",
                                   "code": "therapeutic-implication"}]}
    resource["subject"] = {"reference": f"Patient/{subject}"}
    if len(vids) > 0:
        resource["derivedFrom"] = []
        for vid in vids:
            resource["derivedFrom"].append({"reference": f"Observation/dv-{vid}"})

    resource["component"] = []

    for ptc in implication['phenotypicTreatmentContext']:
        resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                           "code": "81259-4",
                                                           "display": "phenotypic treatment context"}]},
                                      "valueCodeableConcept": {"coding": [{"system": f"{ptc['system']}",
                                                                           "code": f"{ptc['code']}",
                                                                           "display": f"{ptc['display']}"}]}})

    if 'evidenceLevel' in implication and implication['evidenceLevel']:
        resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                           "code": "93044-6",
                                                           "display": "Level of evidence"}]},
                                      "valueCodeableConcept": {"text": f"{implication['evidenceLevel']}"}})

    for med in implication['medicationAssessed']:
        resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                           "code": "51963-7",
                                                           "display": "medication-assessed"}]},
                                      "valueCodeableConcept": {"coding": [{"system": f"{med['system']}",
                                                                           "code": f"{med['code']}",
                                                                           "display": f"{med['display']}"}]}})

    resource["component"].append({"code": {"coding": [{"system": "http://hl7.org/fhir/uv/genomics-reporting/CodeSystem/tbd-codes-cs",
                                                       "code": "predicted-therapeutic-implication",
                                                       "display": "predicted-therapeutic-implication"}]},
                                  "valueCodeableConcept": {"text": implication['predictedImplication']}})

    return resource


def create_tx_implication_profile_pharmgkb(implication, subject, gids):
    resource = OrderedDict()
    resource["resourceType"] = "Observation"
    resource["id"] = "dv-" + str(implication['_id'])
    resource["meta"] = {"profile": [
                        "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/therapeutic-implication"]}
    resource["status"] = "final"
    resource["category"] = [{"coding": [{"system": "http://terminology.hl7.org/CodeSystem/observation-category",
                                         "code": "laboratory"}]}]
    resource["code"] = {"coding": [{"system": "http://hl7.org/fhir/uv/genomics-reporting/CodeSystem/tbd-codes-cs",
                                   "code": "therapeutic-implication"}]}
    resource["subject"] = {"reference": f"Patient/{subject}"}
    if len(gids) > 0:
        resource["derivedFrom"] = []
        for gid in gids:
            resource["derivedFrom"].append({"reference": f"Observation/dv-{gid}"})

    resource["component"] = []

    resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                       "code": "93044-6",
                                                       "display": "Level of evidence"}]},
                                  "valueCodeableConcept": {"text": f"{implication['evidenceLevel']}"}})

    for med in implication['medicationAssessed']:
        resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                           "code": "51963-7",
                                                           "display": "medication-assessed"}]},
                                      "valueCodeableConcept": {"coding": [{"system": f"{med['system']}",
                                                                           "code": f"{med['code']}",
                                                                           "display": f"{med['display']}"}]}})

    resource["component"].append({"code": {"coding": [{"system": "http://hl7.org/fhir/uv/genomics-reporting/CodeSystem/tbd-codes-cs",
                                                       "code": "predicted-therapeutic-implication",
                                                       "display": "predicted-therapeutic-implication"}]},
                                  "valueCodeableConcept": {"text": implication['predictedImplication']}})

    return resource


def create_haplotype_profile(haplotype, subject, gid):
    resource = OrderedDict()
    resource["resourceType"] = "Observation"
    resource["id"] = "dv-" + str(haplotype['_id'])
    resource["meta"] = {"profile": [
                        "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/haplotype"]}
    resource["status"] = "final"
    resource["category"] = [{"coding": [{"system": "http://terminology.hl7.org/CodeSystem/observation-category",
                                         "code": "laboratory"}]}]
    resource["code"] = {"coding": [{"system": "http://hl7.org/fhir/uv/genomics-reporting/CodeSystem/tbd-codes-cs",
                                   "code": "diagnostic-implication"}]}
    resource["subject"] = {"reference": f"Patient/{subject}"}
    if gid:
        resource["derivedFrom"] = [{"reference": f"Observation/{gid}"}]

    return resource


def create_genotype_profile(genotype, subject, gids):
    resource = OrderedDict()
    resource["resourceType"] = "Observation"
    resource["id"] = "dv-" + str(genotype['_id'])
    resource["meta"] = {"profile": [
                        "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/genotype"]}
    resource["status"] = "final"
    resource["category"] = [{"coding": [{"system": "http://terminology.hl7.org/CodeSystem/observation-category",
                                         "code": "laboratory"}]}]
    resource["code"] = {"coding": [{"system": "http://loinc.org",
                                    "code": "84413-4",
                                    "display": "Genotype display name"}]}

    resource["valueCodeableConcept"] = {"coding": [{"system": f"{genotype['genotypeCodeSystem']}",
                                                    "code": f"{genotype['genotypeCode']}",
                                                    "display": f"{genotype['genotypeDesc']}"}]}

    resource["subject"] = {
        "reference": f"Patient/{subject}"
    }

    resource["component"] = []
    resource["component"].append({"code": {"coding": [{"system": "http://loinc.org",
                                                       "code": "48018-6",
                                                       "display": "Gene studied"}]},

                                  "valueCodeableConcept": {"coding": [{"system": "http://www.genenames.org/geneId",
                                                                       "code": f"{genotype['geneCode']}",
                                                                       "display": f"{genotype['geneDesc']}"}]}})

    return resource


def add_variation_id(resource, variation_id):
    spdi_index = next((i for i, item in enumerate(resource["component"]) if item["code"]["coding"][0]["code"] == "81252-9"), None)

    for var_id in variation_id:
        if spdi_index is not None:
            resource["component"][spdi_index]["valueCodeableConcept"]["coding"].append(
                {"system": f'{var_id["system"]}',
                 "code": f'{var_id["code"]}'})
        else:
            resource["component"].append(
                {"code": {"coding": [{"system": "http://loinc.org",
                                      "code": "81252-9",
                                      "display": "Discrete genetic variant"}]},
                 "valueCodeableConcept": {"coding": [{"system": f'{var_id["system"]}',
                                                      "code": f'{var_id["code"]}'}]}})


def create_sequence_phase_relationship(subject, sequence_phase_data):
    resource = OrderedDict()
    resource["resourceType"] = "Observation"
    resource["id"] = f"sid-{str(sequence_phase_data['_id'])}"
    resource["meta"] = {"profile": [
                        "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/sequence-phase-relationship"]}
    resource["status"] = "final"
    resource["category"] = [{"coding": [{"system": "http://terminology.hl7.org/CodeSystem/observation-category",
                                         "code": "laboratory"}]}]
    resource["code"] = {"coding": [{"system": "http://loinc.org",
                                    "code": "82120-7",
                                    "display": "Allelic phase"}]}
    resource["subject"] = {"reference": f"Patient/{subject}"}
    resource["valueCodeableConcept"] = {"coding": [{"system": "http://hl7.org/fhir/uv/genomics-reporting/CodeSystem/sequence-phase-relationship-cs",
                                        "code": sequence_phase_data["phase"],
                                                    "display": sequence_phase_data["phase"]}]}
    resource["derivedFrom"] = [{"reference": f"Observation/dv-{sequence_phase_data['variantID1']}"},
                               {"reference": f"Observation/dv-{sequence_phase_data['variantID2']}"}]

    return resource


def get_sequence_phase_data(patient_id):
    try:
        results = phase_data.aggregate([{'$match': {"patientID": {"$eq": patient_id}}}])
        results = list(results)
    except Exception as e:
        print(f"DEBUG: Error({e}) under get_sequence_phase_data(patient_id={patient_id})")
        results = []

    return results


def lift_over(ref_seq, start, end):
    try:
        start = int(start)
        end = int(end)
    except ValueError:
        return None

    tempDict = get_build_and_chrom_by_ref_seq(ref_seq)
    if tempDict:
        chrom = tempDict['chrom']
        build = tempDict['build']
    else:
        return None

    # convert build37 to build38
    if build == 'GRCh37':
        new_ref_seq = get_ref_seq_by_chrom_and_build('GRCh38', chrom)
        lo_b37_to_b38 = get_liftover('hg19', 'hg38')
        new_coor = lo_b37_to_b38.convert_coordinate(chrom, start, '+')
        if len(new_coor) != 1 or chrom != new_coor[0][0]:
            return None
        else:
            new_start = new_coor[0][1]
        new_coor = lo_b37_to_b38.convert_coordinate(chrom, end, '+')
        if len(new_coor) != 1 or chrom != new_coor[0][0]:
            return None
        else:
            new_end = new_coor[0][1]
        if new_start > new_end:
            return None
        else:
            return {"refSeq": new_ref_seq, "start": new_start, "end": new_end}

    # convert build38 to build37
    elif build == 'GRCh38':
        new_ref_seq = get_ref_seq_by_chrom_and_build('GRCh37', chrom)
        lo_b38_to_b37 = get_liftover('hg38', 'hg19')
        new_coor = lo_b38_to_b37.convert_coordinate(chrom, start, '+')
        if len(new_coor) != 1 or chrom != new_coor[0][0]:
            return None
        else:
            new_start = new_coor[0][1]
        new_coor = lo_b38_to_b37.convert_coordinate(chrom, end, '+')
        if len(new_coor) != 1 or chrom != new_coor[0][0]:
            return None
        else:
            new_end = new_coor[0][1]
        if new_start > new_end:
            return None
        else:
            return {"refSeq": new_ref_seq, "start": new_start, "end": new_end}
    else:
        return None


def get_build_and_chrom_by_ref_seq(ref_seq):
    chrom_data = chromosomes_data.find_one({'$or': [{'build37RefSeq': {'$eq': ref_seq}}, {'build38RefSeq': {'$eq': ref_seq}}]})

    if chrom_data:
        if chrom_data['build37RefSeq'] == ref_seq:
            return {"chrom": chrom_data['chr'], "build": 'GRCh37'}
        else:
            return {"chrom": chrom_data['chr'], "build": 'GRCh38'}

    return None


def get_ref_seq_by_chrom_and_build(build, chrom):
    chrom_data = chromosomes_data.find_one({'chr': {'$eq': chrom}})

    if chrom_data:
        if build == 'GRCh37':
            return chrom_data['build37RefSeq']
        if build == 'GRCh38':
            return chrom_data['build38RefSeq']

    return None


def get_other_build(build):
    return "GRCh37" if build == 'GRCh38' else "GRCh38"


def get_intersected_regions(bed_id, build, chrom, start, end, intersected_regions):
    try:
        result = beds.aggregate([
            {"$match": {"$and": [{"BedID": {"$eq": bed_id}},
                                 {"BED.Chromosome": {"$eq": chrom}},
                                 {'$or': [
                                     {"$and": [
                                         {"BED.Start": {"$gte": start}},
                                         {"BED.Start": {"$lte": end}}
                                     ]},
                                     {"$and": [
                                         {"BED.Start": {"$lte": start}},
                                         {"BED.End": {"$gte": start}}
                                     ]}]}]}},
            {"$unwind": "$BED"},
            {"$match": {"$and": [{"BED.Chromosome": {"$eq": chrom}},
                                 {'$or': [
                                     {"$and": [
                                         {"BED.Start": {"$gte": start}},
                                         {"BED.Start": {"$lte": end}}
                                     ]},
                                     {"$and": [
                                         {"BED.Start": {"$lte": start}},
                                         {"BED.End": {"$gte": start}}
                                     ]}]}]}},
            {"$group": {"_id": "$BedID", "BED": {"$push": "$$ROOT.BED"}}}])

        result = list(result)
    except Exception as e:
        print(f"DEBUG: Error({e}) under get_intersected_regions(bed_id={bed_id}, build={build}, chrom={chrom}, start={start}, end={end}, intersected_regions={intersected_regions})")
        result = result

    if result:
        result = result[0]["BED"]

        for csePair in result:
            ref_seq = get_ref_seq_by_chrom_and_build(build, csePair["Chromosome"])
            intersected_regions.append(f'{ref_seq}:{max(start, csePair["Start"])}-{min(end, csePair["End"])}')


def query_clinvar_by_variants(normalized_variant_list, code_list, query, population=False):
    variant_list = []
    for item in normalized_variant_list:
        if "GRCh37" in item:
            variant_list.append(item["GRCh37"])
        if "GRCh38" in item:
            variant_list.append(item["GRCh38"])

    pipeline_part = [{'$match': {'$expr': {'$and': [{'$or': [{'$eq': ['$b37SPDI', '$$mySPDI']},
                                                             {'$eq': ['$b38SPDI', '$$mySPDI']}]}]}}},
                     {'$addFields': {}}]
    if code_list != []:
        pipeline_part.append({'$match': {'$or': []}})
        or_query = []
        for condition in code_list:
            if condition['isSystem']:
                or_query.append({'$and': [{'predictedPhenotype.code': {'$eq': condition['condition']}}, {'predictedPhenotype.system': {'$eq': condition['system']}}]})
            else:
                or_query.append({'$or': [
                    {'predictedPhenotype.code': {'$regex': ".*"+str(condition['condition']).replace('*', r'\*')+".*"}},
                    {'predictedPhenotype.display': {'$regex': ".*"+str(condition['condition']).replace('*', r'\*')+".*"}}
                ]})
        pipeline_part[-1]['$match']['$or'] = or_query
        pipeline_part.append({"$unwind": "$predictedPhenotype"})
        pipeline_part.append({'$match': {'$or': or_query}})
        pipeline_part.append({"$group": {"_id": "$_id",
                                         "variationID": {
                                             "$first": "$$ROOT.variationID"
                                         },
                                         "b37SPDI": {
                                             "$first": "$$ROOT.b37SPDI"
                                         },
                                         "b38SPDI": {
                                             "$first": "$$ROOT.b38SPDI"
                                         },
                                         "evidenceLevel": {
                                             "$first": "$$ROOT.evidenceLevel"
                                         },
                                         "clinicalSignificance": {
                                             "$first": "$$ROOT.clinicalSignificance"
                                         },
                                         "predictedPhenotype": {
                                             "$push": "$$ROOT.predictedPhenotype"
                                         }}})

    query['SPDI'] = {'$in': variant_list}

    query_string = [{'$match': query},
                    {'$lookup': {'from': 'dxImplication', 'let': {'mySPDI': '$SPDI'}, 'pipeline': pipeline_part,
                                 'as': 'dxImplicationMatches'}},
                    {'$addFields': {}},
                    {'$match': {'dxImplicationMatches': {'$exists': True, '$not': {'$size': 0}}}}]

    if population:
        query_string.append({'$group': {'_id': '$patientID'}})

    try:
        results = variants_db.aggregate(query_string)
        results = list(results)
    except Exception as e:
        print(f"DEBUG: Error({e}) under query_clinvar_by_variants(normalized_variant_list={normalized_variant_list}, code_list={code_list}, query={query}, population={population})")
        results = []

    query_results = []

    for item in results:
        item["UUID"] = str(uuid4())
        query_results.append(item)

    return query_results


def query_clinvar_by_condition(code_list, query):
    condition_query = {'$or': []}

    or_query = []
    for condition in code_list:
        if condition['isSystem']:
            or_query.append({'$and': [{'predictedPhenotype.code': {'$eq': condition['condition']}}, {'predictedPhenotype.system': {'$eq': condition['system']}}]})
        else:
            or_query.append({'$or': [
                {'predictedPhenotype.code': {'$regex': ".*"+str(condition['condition']).replace('*', r'\*')+".*"}},
                {'predictedPhenotype.display': {'$regex': ".*"+str(condition['condition']).replace('*', r'\*')+".*"}}
            ]})

    condition_query['$or'] = or_query

    query_string = [
        {'$match': condition_query},
        {"$unwind": "$predictedPhenotype"},
        {'$match': {'$or': or_query}},
        {"$group": {"_id": "$_id",
                    "variationID": {
                        "$first": "$$ROOT.variationID"
                    },
                    "b37SPDI": {
                        "$first": "$$ROOT.b37SPDI"
                    },
                    "b38SPDI": {
                        "$first": "$$ROOT.b38SPDI"
                    },
                    "evidenceLevel": {
                        "$first": "$$ROOT.evidenceLevel"
                    },
                    "clinicalSignificance": {
                        "$first": "$$ROOT.clinicalSignificance"
                    },
                    "predictedPhenotype": {
                        "$push": "$$ROOT.predictedPhenotype"
                    }}},
        {'$lookup': {'from': 'Variants', 'let': {'b37SPDI': '$b37SPDI', 'b38SPDI': '$b38SPDI'},
                     'pipeline': [{'$match': query},
                                  {'$addFields': {}},
                                  {'$match': {'$expr': {'$and': [{'$or': [
                                      {'$eq': ['$SPDI', '$$b37SPDI']},
                                      {'$eq': ['$SPDI', '$$b38SPDI']}
                                  ]}]}}},
                                  {'$addFields': {}}],
                     'as': 'patientMatches'}},
        {'$addFields': {}},
        {'$match': {'patientMatches': {'$exists': True, '$not': {'$size': 0}}}}
    ]

    try:
        results = dxImplication_db.aggregate(query_string)
        results = list(results)
    except Exception as e:
        print(f"DEBUG: Error({e}) under query_clinvar_by_condition(code_list={code_list}, query={query})")
        results = []

    query_results = []

    for item in results:
        for match in item["patientMatches"]:
            match["UUID"] = str(uuid4())
        query_results.append(item)

    return query_results


def query_CIVIC_by_variants(normalized_variant_list, code_list, treatment_list, query, population=False):
    variant_list = []
    for item in normalized_variant_list:
        if "GRCh37" in item:
            variant_list.append(item["GRCh37"])
        if "GRCh38" in item:
            variant_list.append(item["GRCh38"])

    pipeline_part = [{'$match': {'$expr': {'$and': [{'$or': [{'$eq': ['$b37SPDI', '$$mySPDI']},
                                                             {'$eq': ['$b38SPDI', '$$mySPDI']}]}]}}},
                     {'$addFields': {}}]

    if code_list != []:
        pipeline_part.append({'$match': {'$or': []}})
        or_query = []
        for condition in code_list:
            if condition['isSystem']:
                or_query.append({'$and': [{'phenotypicTreatmentContext.code': {'$eq': condition['condition']}}, {'phenotypicTreatmentContext.system': {'$eq': condition['system']}}]})
            else:
                or_query.append({'$or': [
                    {'phenotypicTreatmentContext.code': {'$regex': ".*"+str(condition['condition']).replace('*', r'\*')+".*"}},
                    {'phenotypicTreatmentContext.display': {'$regex': ".*"+str(condition['condition']).replace('*', r'\*')+".*"}}
                ]})
        pipeline_part[-1]['$match']['$or'] = or_query
        pipeline_part.append({"$unwind": "$phenotypicTreatmentContext"})
        pipeline_part.append({'$match': {'$or': or_query}})
        pipeline_part.append({"$group": {"_id": "$_id",
                                         "variationID": {
                                             "$first": "$$ROOT.variationID"
                                         },
                                         "b37SPDI": {
                                             "$first": "$$ROOT.b37SPDI"
                                         },
                                         "b38SPDI": {
                                             "$first": "$$ROOT.b38SPDI"
                                         },
                                         "evidenceLevel": {
                                             "$first": "$$ROOT.evidenceLevel"
                                         },
                                         "predictedImplication": {
                                             "$first": "$$ROOT.predictedImplication"
                                         },
                                         "medicationAssessed": {
                                             "$first": "$$ROOT.medicationAssessed"
                                         },
                                         "phenotypicTreatmentContext": {
                                             "$push": "$$ROOT.phenotypicTreatmentContext"
                                         }}})

    if treatment_list != []:
        pipeline_part.append({'$match': {'$or': []}})
        or_query = []
        for treatment in treatment_list:
            if treatment['isSystem']:
                or_query.append({'$and': [{'medicationAssessed.code': {'$eq': treatment['treatment']}}, {'medicationAssessed.system': {'$eq': treatment['system']}}]})
            else:
                or_query.append({'$or': [
                    {'medicationAssessed.code': {'$regex': ".*"+str(treatment['treatment']).replace('*', r'\*')+".*"}},
                    {'medicationAssessed.display': {'$regex': ".*"+str(treatment['treatment']).replace('*', r'\*')+".*"}}
                ]})
        pipeline_part[-1]['$match']['$or'] = or_query
        pipeline_part.append({"$unwind": "$medicationAssessed"})
        pipeline_part.append({'$match': {'$or': or_query}})
        pipeline_part.append({"$group": {"_id": "$_id",
                                         "variationID": {
                                             "$first": "$$ROOT.variationID"
                                         },
                                         "b37SPDI": {
                                             "$first": "$$ROOT.b37SPDI"
                                         },
                                         "b38SPDI": {
                                             "$first": "$$ROOT.b38SPDI"
                                         },
                                         "evidenceLevel": {
                                             "$first": "$$ROOT.evidenceLevel"
                                         },
                                         "predictedImplication": {
                                             "$first": "$$ROOT.predictedImplication"
                                         },
                                         "phenotypicTreatmentContext": {
                                             "$first": "$$ROOT.phenotypicTreatmentContext"
                                         },
                                         "medicationAssessed": {
                                             "$push": "$$ROOT.medicationAssessed"
                                         }}})

    query['SPDI'] = {'$in': variant_list}

    query_string = [{'$match': query},
                    {'$lookup': {'from': 'txImplication', 'let': {'mySPDI': '$SPDI'}, 'pipeline': pipeline_part,
                                 'as': 'txImplicationMatches'}},
                    {'$addFields': {}},
                    {'$match': {'txImplicationMatches': {'$exists': True, '$not': {'$size': 0}}}}]

    if population:
        query_string.append({'$group': {'_id': '$patientID'}})

    try:
        results = variants_db.aggregate(query_string)
        results = list(results)
    except Exception as e:
        print(f"DEBUG: Error({e}) under query_CIVIC_by_variants(normalized_variant_list={normalized_variant_list}, code_list={code_list}, treatment_list={treatment_list}, query={query}, population={population})")
        results = []

    query_results = []

    for item in results:
        item["UUID"] = str(uuid4())
        query_results.append(item)

    return query_results


def query_CIVIC_by_condition(code_list, treatment_list, query):
    condition_query = {"$or": []}
    if code_list != []:
        condition_or_query = []
        for condition in code_list:
            if condition['isSystem']:
                condition_or_query.append({'$and': [{'phenotypicTreatmentContext.code': {'$eq': condition['condition']}}, {'phenotypicTreatmentContext.system': {'$eq': condition['system']}}]})
            else:
                condition_or_query.append({'$or': [
                    {'phenotypicTreatmentContext.code': {'$regex': ".*"+str(condition['condition']).replace('*', r'\*')+".*"}},
                    {'phenotypicTreatmentContext.display': {'$regex': ".*"+str(condition['condition']).replace('*', r'\*')+".*"}}
                ]})

        condition_query['$or'] = condition_or_query

    treatment_query = {"$or": []}
    if treatment_list != []:
        treatment_or_query = []
        for treatment in treatment_list:
            if treatment['isSystem']:
                treatment_or_query.append({'$and': [{'medicationAssessed.code': {'$eq': str(treatment['treatment'])}}, {'medicationAssessed.system': {'$eq': treatment['system']}}]})
            else:
                treatment_or_query.append({'$or': [
                    {'medicationAssessed.code': {'$regex': ".*"+str(treatment['treatment']).replace('*', r'\*')+".*"}},
                    {'medicationAssessed.display': {'$regex': ".*"+str(treatment['treatment']).replace('*', r'\*')+".*"}}
                ]})

        treatment_query['$or'] = treatment_or_query

    condition_query_string = [{'$match': condition_query},
                              {"$unwind": "$phenotypicTreatmentContext"},
                              {'$match': condition_query},
                              {"$group": {"_id": "$_id",
                                          "variationID": {
                                              "$first": "$$ROOT.variationID"
                                          },
                                          "b37SPDI": {
                                              "$first": "$$ROOT.b37SPDI"
                                          },
                                          "b38SPDI": {
                                              "$first": "$$ROOT.b38SPDI"
                                          },
                                          "evidenceLevel": {
                                              "$first": "$$ROOT.evidenceLevel"
                                          },
                                          "predictedImplication": {
                                              "$first": "$$ROOT.predictedImplication"
                                          },
                                          "medicationAssessed": {
                                              "$first": "$$ROOT.medicationAssessed"
                                          },
                                          "phenotypicTreatmentContext": {
                                              "$push": "$$ROOT.phenotypicTreatmentContext"
                                          }}}]

    treatment_query_string = [{'$match': treatment_query},
                              {"$unwind": "$medicationAssessed"},
                              {'$match': treatment_query},
                              {"$group": {"_id": "$_id",
                                          "variationID": {
                                              "$first": "$$ROOT.variationID"
                                          },
                                          "b37SPDI": {
                                              "$first": "$$ROOT.b37SPDI"
                                          },
                                          "b38SPDI": {
                                              "$first": "$$ROOT.b38SPDI"
                                          },
                                          "evidenceLevel": {
                                              "$first": "$$ROOT.evidenceLevel"
                                          },
                                          "predictedImplication": {
                                              "$first": "$$ROOT.predictedImplication"
                                          },
                                          "medicationAssessed": {
                                              "$push": "$$ROOT.medicationAssessed"
                                          },
                                          "phenotypicTreatmentContext": {
                                              "$first": "$$ROOT.phenotypicTreatmentContext"
                                          }}}]

    query_string = []
    if condition_query['$or']:
        query_string.extend(condition_query_string)

    if treatment_query['$or']:
        query_string.extend(treatment_query_string)

    query_string.extend([{'$lookup': {'from': 'Variants', 'let': {'b37SPDI': '$b37SPDI', 'b38SPDI': '$b38SPDI'},
                                      'pipeline': [{'$match': query},
                                                   {'$match': {'$expr': {'$and': [{'$or': [
                                                       {'$eq': ['$SPDI', '$$b37SPDI']},
                                                       {'$eq': ['$SPDI', '$$b38SPDI']}
                                                   ]}]}}},
                                                   {'$addFields': {}}],
                                      'as': 'patientMatches'}},
                         {'$addFields': {}},
                         {'$match': {'patientMatches': {'$exists': True, '$not': {'$size': 0}}}}])

    try:
        results = txImplication_db.aggregate(query_string)
        results = list(results)
    except Exception as e:
        print(f"DEBUG: Error({e}) under query_CIVIC_by_condition(code_list={code_list}, treatment_list={treatment_list}, query={query})")
        results = []

    query_results = []

    for item in results:
        for match in item["patientMatches"]:
            match["UUID"] = str(uuid4())
        query_results.append(item)

    return query_results


def query_PharmGKB_by_haplotypes(normalizedHaplotypeList, treatmentCodeList, query, population=False):
    pipeline_part = [{'$match': {'$expr': {'$and': [{'$eq': ['$genotype.code', '$$mygenotypeCode']}]}}},
                     {'$addFields': {}}]

    if treatmentCodeList != []:
        pipeline_part.append({'$match': {'$or': []}})
        or_query = []
        for treatment in treatmentCodeList:
            if treatment['isSystem']:
                or_query.append({'$and': [{'medicationAssessed.code': {'$eq': treatment['treatment']}}, {'medicationAssessed.system': {'$eq': treatment['system']}}]})
            else:
                or_query.append({'$or': [
                    {'medicationAssessed.code': {'$regex': ".*"+str(treatment['treatment']).replace('*', r'\*')+".*"}},
                    {'medicationAssessed.display': {'$regex': ".*"+str(treatment['treatment']).replace('*', r'\*')+".*"}}
                ]})
        pipeline_part[-1]['$match']['$or'] = or_query
        pipeline_part.append({"$unwind": "$medicationAssessed"})
        pipeline_part.append({'$match': {'$or': or_query}})
        pipeline_part.append({"$group": {"_id": "$_id",
                                         "variationID": {
                                             "$first": "$$ROOT.variationID"
                                         },
                                         "b37SPDI": {
                                             "$first": "$$ROOT.b37SPDI"
                                         },
                                         "b38SPDI": {
                                             "$first": "$$ROOT.b38SPDI"
                                         },
                                         "evidenceLevel": {
                                             "$first": "$$ROOT.evidenceLevel"
                                         },
                                         "predictedImplication": {
                                             "$first": "$$ROOT.predictedImplication"
                                         },
                                         "genotype_code": {
                                             "$first": {"$first": "$$ROOT.genotype.code"}
                                         },
                                         "genotype_system": {
                                             "$first": {"$first": "$$ROOT.genotype.system"}
                                         },
                                         "genotype_display": {
                                             "$first": {"$first": "$$ROOT.genotype.display"}
                                         },
                                         "medicationAssessed": {
                                             "$push": "$$ROOT.medicationAssessed"
                                         }}})

    query['$or'] = []

    for haplotype in normalizedHaplotypeList:

        if haplotype['isSystem']:
            query['$or'].append({'genotypeCode': {"$eq": haplotype['haplotype']}})
        else:
            query['$or'].append({'$or': [
                {'genotypeCode': {'$regex': ".*"+str(haplotype['haplotype']).replace('*', r'\*')+".*"}},
                {'genotypeDesc': {'$regex': ".*"+str(haplotype['haplotype']).replace('*', r'\*')+".*"}}
            ]})

    query_string = [{'$match': query},
                    {'$lookup': {'from': 'txImplication',
                                 'localField': 'genotypeCode',
                                 'foreignField': 'genotype.code',
                                 'as': 'txImplicationMatches'}},
                    {'$addFields': {}},
                    {'$match': {'txImplicationMatches': {'$exists': True, '$not': {'$size': 0}}}}]

    if population:
        query_string.append({'$group': {'_id': '$patientID'}})

    try:
        results = genotypes_db.aggregate(query_string)
        results = list(results)
    except Exception as e:
        print(f"DEBUG: Error({e}) under query_PharmGKB_by_haplotypes(normalizedHaplotypeList={normalizedHaplotypeList}, treatmentCodeList={treatmentCodeList}, query={query}, population={population})")
        results = []

    query_results = []

    for item in results:
        item["UUID"] = str(uuid4())
        query_results.append(item)

    return query_results


def query_PharmGKB_by_treatments(code_list, treatment_list, query):
    condition_query = {"$or": []}
    if code_list != []:
        condition_or_query = []
        for condition in code_list:
            if condition['isSystem']:
                condition_or_query.append({'$and': [{'phenotypicTreatmentContext.code': {'$eq': condition['condition']}}, {'phenotypicTreatmentContext.system': {'$eq': condition['system']}}]})
            else:
                condition_or_query.append({'$or': [
                    {'phenotypicTreatmentContext.code': {'$regex': ".*"+str(condition['condition']).replace('*', r'\*')+".*"}},
                    {'phenotypicTreatmentContext.display': {'$regex': ".*"+str(condition['condition']).replace('*', r'\*')+".*"}}
                ]})

        condition_query['$or'] = condition_or_query

    treatment_query = {"$or": []}
    if treatment_list != []:
        treatment_or_query = []
        for treatment in treatment_list:
            if treatment['isSystem']:
                treatment_or_query.append({'$and': [{'medicationAssessed.code': {'$eq': str(treatment['treatment'])}}, {'medicationAssessed.system': {'$eq': treatment['system']}}]})
            else:
                treatment_or_query.append({'$or': [
                    {'medicationAssessed.code': {'$regex': ".*"+str(treatment['treatment']).replace('*', r'\*')+".*"}},
                    {'medicationAssessed.display': {'$regex': ".*"+str(treatment['treatment']).replace('*', r'\*')+".*"}}
                ]})

        treatment_query['$or'] = treatment_or_query

    condition_query_string = [{'$match': condition_query},
                              {"$unwind": "$phenotypicTreatmentContext"},
                              {'$match': condition_query},
                              {"$group": {"_id": "$_id",
                                          "variationID": {
                                              "$first": "$$ROOT.variationID"
                                          },
                                          "evidenceLevel": {
                                              "$first": "$$ROOT.evidenceLevel"
                                          },
                                          "predictedImplication": {
                                              "$first": "$$ROOT.predictedImplication"
                                          },
                                          "genotype_code": {
                                              "$first": {"$first": "$$ROOT.genotype.code"}
                                          },
                                          "genotype_system": {
                                              "$first": {"$first": "$$ROOT.genotype.system"}
                                          },
                                          "genotype_display": {
                                              "$first": {"$first": "$$ROOT.genotype.display"}
                                          },
                                          "medicationAssessed": {
                                              "$first": "$$ROOT.medicationAssessed"
                                          },
                                          "phenotypicTreatmentContext": {
                                              "$push": "$$ROOT.phenotypicTreatmentContext"
                                          }}}]

    treatment_query_string = [{'$match': treatment_query},
                              {"$unwind": "$medicationAssessed"},
                              {'$match': treatment_query},
                              {"$group": {"_id": "$_id",
                                          "variationID": {
                                              "$first": "$$ROOT.variationID"
                                          },
                                          "evidenceLevel": {
                                              "$first": "$$ROOT.evidenceLevel"
                                          },
                                          "predictedImplication": {
                                              "$first": "$$ROOT.predictedImplication"
                                          },
                                          "genotype_code": {
                                              "$first": {"$first": "$$ROOT.genotype.code"}
                                          },
                                          "genotype_system": {
                                              "$first": {"$first": "$$ROOT.genotype.system"}
                                          },
                                          "genotype_display": {
                                              "$first": {"$first": "$$ROOT.genotype.display"}
                                          },
                                          "medicationAssessed": {
                                              "$push": "$$ROOT.medicationAssessed"
                                          },
                                          "phenotypicTreatmentContext": {
                                              "$first": "$$ROOT.phenotypicTreatmentContext"
                                          }}}
                              ]

    query_string = []
    if condition_query['$or']:
        query_string.extend(condition_query_string)

    if treatment_query['$or']:
        query_string.extend(treatment_query_string)

    if 'genomicSourceClass' in query:
        query.pop('genomicSourceClass')

    query_string.extend([
        {'$lookup': {'from': 'Genotypes', 'let': {'genotype_code_v': '$genotype_code'},
                     'pipeline': [{'$match': query},
                                  {'$match': {'$expr': {'$and': [{'$eq': ['$genotypeCode', '$$genotype_code_v']}]}}},
                                  {'$addFields': {}}],
                     'as': 'patientMatches'}},
        {'$addFields': {}},
        {'$match': {'patientMatches': {'$exists': True, '$not': {'$size': 0}}}}])

    try:
        results = txImplication_db.aggregate(query_string)
        results = list(results)
    except Exception as e:
        print(f"DEBUG: Error({e}) under query_PharmGKB_by_treatments(code_list={code_list}, treatment_list={treatment_list}, query={query})")
        results = []

    query_results = []

    for item in results:
        for match in item["patientMatches"]:
            match["UUID"] = str(uuid4())
        query_results.append(item)

    return query_results


def query_genes(gene):
    query = {}

    pipeline_part = [{'$match': {'$expr': {'$eq': ['$ncbiGeneSymbol', '$$geneSymbol']}}},
                     {'$addFields': {}}]

    query['$or'] = [{"HGNCgeneId": {"$eq": gene}},
                    {"HGNCgeneSymbol": {"$eq": gene}}]

    query_string = [{'$match': query},
                    {'$lookup': {'from': 'Transcripts', 'let': {'geneSymbol': '$NCBIgeneSymbol'}, 'pipeline': pipeline_part,
                                 'as': 'transcriptMatches'}},
                    {'$addFields': {}}]

    try:
        results = genes_data.aggregate(query_string)
        results = list(results)
    except Exception as e:
        print(f"DEBUG: Error({e}) under get_feature_coordinates(gene={gene})")
        results = []

    return results


def query_genes_range(given_range):
    pipeline_part = [{'$match': {'$expr': {'$eq': ['$ncbiGeneSymbol', '$$geneSymbol']}}},
                     {'$addFields': {}}]

    ref_seq = given_range['CHROMOSOME']['RefSeq']
    build = get_build_and_chrom_by_ref_seq(ref_seq)['build']
    start = given_range['RANGE']['L']
    end = given_range['RANGE']['H']

    if build == 'GRCh37':
        query = {"$and": [{"build37RefSeq": {"$eq": ref_seq}},
                          {'$or': [
                              {"$and": [
                                  {"build37Start": {"$gte": start}},
                                  {"build37Start": {"$lte": end}}
                              ]},
                              {"$and": [
                                  {"build37Start": {"$lte": start}},
                                  {"build37End": {"$gte": start}}
                              ]}]}]}
    else:
        query = {"$and": [{"build38RefSeq": {"$eq": ref_seq}},
                          {'$or': [
                              {"$and": [
                                  {"build38Start": {"$gte": start}},
                                  {"build38Start": {"$lte": end}}
                              ]},
                              {"$and": [
                                  {"build38Start": {"$lte": start}},
                                  {"build38End": {"$gte": start}}
                              ]}]}]}

    query_string = [{'$match': query},
                    {'$lookup': {'from': 'Transcripts', 'let': {'geneSymbol': '$NCBIgeneSymbol'}, 'pipeline': pipeline_part,
                                 'as': 'transcriptMatches'}},
                    {'$addFields': {}}]

    try:
        results = genes_data.aggregate(query_string)
        results = list(results)
    except Exception as e:
        print(f"DEBUG: Error({e}) under find_the_gene(given_range={given_range})")
        results = []

    return results


def query_transcript(transcript):
    query = {}

    gene_pipeline_part = [{'$match': {'$expr': {'$eq': ['$NCBIgeneSymbol', '$$geneSymbol']}}},
                          {'$addFields': {}}]

    exon_pipeline_part = [{'$match': {'$expr': {'$eq': ['$transcript', '$$refSeq']}}},
                          {'$addFields': {}}]

    cds_pipeline_part = [{'$match': {'$expr': {'$eq': ['$transcript', '$$refSeq']}}},
                         {'$addFields': {}}]

    query['transcriptRefSeq'] = {"$regex": ".*"+str(transcript).replace('*', r'\*')+".*"}

    query_string = [{'$match': query},
                    {'$lookup': {'from': 'Genes', 'let': {'geneSymbol': '$ncbiGeneSymbol'}, 'pipeline': gene_pipeline_part,
                                 'as': 'geneMatches'}},
                    {'$addFields': {}},
                    {'$lookup': {'from': 'Exons', 'let': {'refSeq': '$transcriptRefSeq'}, 'pipeline': exon_pipeline_part,
                                 'as': 'exonMatches'}},
                    {'$addFields': {}},
                    {'$lookup': {'from': 'CDS', 'let': {'refSeq': '$transcriptRefSeq'}, 'pipeline': cds_pipeline_part,
                                 'as': 'cdsMatches'}},
                    {'$addFields': {}}]

    try:
        results = transcripts_data.aggregate(query_string)
        results = list(results)
    except Exception as e:
        print(f"DEBUG: Error({e}) under get_feature_coordinates(transcript={transcript})")
        results = []

    if len(results) > 1:
        abort(400, "Unable to provide information on this transcript at this time")

    return results


def query_molecular_consequences_by_variants(normalized_variant_list, feature_consequence_list, query):
    variant_list = []
    for item in normalized_variant_list:
        if "GRCh37" in item:
            variant_list.append(item["GRCh37"])
        if "GRCh38" in item:
            variant_list.append(item["GRCh38"])

    pipeline_part = [{'$match': {'$expr': {'$and': [{'$or': [{'$eq': ['$variantID', '$$myvariant_id']}]}]}}}]

    if feature_consequence_list != []:
        pipeline_part.append({'$match': {'$or': []}})
        or_query = []

        for feature_consequence in feature_consequence_list:
            if feature_consequence['isSystem']:
                or_query.append({'$and': [{'featureConsequence.code': {'$eq': feature_consequence['feature_consequence']}}, {'featureConsequence.system': {'$eq': feature_consequence['system']}}]})
            else:
                or_query.append({'$or': [
                    {'featureConsequence.code': {'$regex': ".*"+str(feature_consequence['feature_consequence']).replace('*', r'\*')+".*"}},
                    {'featureConsequence.display': {'$regex': ".*"+str(feature_consequence['feature_consequence']).replace('*', r'\*')+".*"}}
                ]})
        pipeline_part[-1]['$match']['$or'] = or_query

    query['SPDI'] = {'$in': variant_list}

    query_string = [{'$match': query},
                    {'$lookup': {'from': 'MolecConseq', 'let': {'myvariant_id': '$_id'}, 'pipeline': pipeline_part,
                                 'as': 'molecularConsequenceMatches'}},
                    {'$match': {'molecularConsequenceMatches': {'$exists': True, '$not': {'$size': 0}}}}]
    try:
        results = variants_db.aggregate(query_string)
        results = list(results)
    except Exception as e:
        print(f"DEBUG: Error {e} under query_molecular_consequences_by_variants query={query}")
        results = []

    query_results = []

    for item in results:
        query_results.append(item)

    return query_results
