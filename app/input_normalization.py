import os

import hgvs.assemblymapper
import hgvs.dataproviders.uta
import hgvs.parser

from utilities.SPDI_Normalization import get_normalized_spdi

# Set the HGVS_SEQREPO_URL env var so the hgvs library will use the local `utilities/seqfetcher` endpoint instead of
# making NCBI API calls.
port = os.getenv('PORT', 5000)  # The localhost debugger starts the app on port 5000
os.environ['HGVS_SEQREPO_URL'] = f"http://localhost:{port}/utilities/seqfetcher"

database_schema = os.getenv('UTA_DATABASE_SCHEMA', 'uta_20210129b')
# Use the biocommons UTA database if we don't specify a custom one.
# Also, make sure the URL uses `postgresql` instead of `postgres` as schema
database_url = f"{os.getenv('UTA_DATABASE_URL', 'postgresql://anonymous:anonymous@uta.biocommons.org/uta')}/{database_schema}".replace('postgres://', 'postgresql://')

hgvsParser = hgvs.parser.Parser()
hgvsDataProvider = hgvs.dataproviders.uta.connect(db_url=database_url)

# Note: One can set `replace_reference=False` and `prevalidation_level=None` to skip data validation
# Also, until https://github.com/biocommons/hgvs/issues/704 is addressed, the following config settings need to also be set.
# hgvs.config.global_config.normalizer.validate = False
# hgvs.global_config.mapping.prevalidation_level = None  # TODO: Open issue
# hgvs.global_config.mapping.replace_reference = False  # TODO: Open issue
# Until https://github.com/biocommons/hgvs/issues/705 is addressed, validation cannot be disabled completely.

b37hgvsAssemblyMapper = hgvs.assemblymapper.AssemblyMapper(
    hgvsDataProvider, assembly_name='GRCh37', alt_aln_method='splign')
b38hgvsAssemblyMapper = hgvs.assemblymapper.AssemblyMapper(
    hgvsDataProvider, assembly_name='GRCh38', alt_aln_method='splign')
# ------------- point to latest data source ------------------------
# at unix command line: export UTA_DB_URL=postgresql://anonymous:anonymous@uta.biocommons.org/uta/uta_20210129b

# ------------------ PARSE -------------


def parse_variant(variant):
    parsed_variant_dict = dict()
    parsed_variant_dict['parsed'] = hgvsParser.parse_hgvs_variant(variant)
    return parsed_variant_dict

# ------------------ PROJECT -------------


def project_variant(parsed_variant):
    projected_variant_dict = dict()
    projected_variant_dict['b37projected'] = b37hgvsAssemblyMapper.c_to_g(
        parsed_variant)
    projected_variant_dict['b38projected'] = b38hgvsAssemblyMapper.c_to_g(
        parsed_variant)
    return projected_variant_dict

# ---------------- NORMALIZE to canonical SPDIs ---------------


def normalize_variant(parsed_variant, build):
    pos = parsed_variant.posedit.pos.start.base-1
    if parsed_variant.posedit.edit.ref:
        ref = parsed_variant.posedit.edit.ref
    else:  # ref is blank for insertions
        ref = ''
        pos = pos+1
    if str(parsed_variant.posedit.edit) == 'dup':
        alt = ref+ref
    elif parsed_variant.posedit.edit.alt:
        alt = parsed_variant.posedit.edit.alt
    else:  # alt is blank for deletions
        alt = ''
    return get_normalized_spdi(parsed_variant.ac, pos, ref, alt)

# ---------------- CONVERT NM_HGVS to canonical SPDIs ---------------


def process_NM_HGVS(NM_HGVS):
    parsed_variant_dict = parse_variant(NM_HGVS)
    print(f"parsed: {parsed_variant_dict['parsed']}")

    projected_variant_dict = project_variant(parsed_variant_dict['parsed'])
    print(
        f"b37projected: {projected_variant_dict['b37projected']}; b38projected: {projected_variant_dict['b38projected']}")

    b37SPDI = normalize_variant(
        projected_variant_dict['b37projected'], 'GRCh37')
    b38SPDI = normalize_variant(
        projected_variant_dict['b38projected'], 'GRCh38')
    print(f"b37normalized: {b37SPDI}; b38normalized: {b38SPDI}")

    return {'GRCh37SPDI': b37SPDI, 'GRCh38SPDI': b38SPDI, 'GRCh37HGVS': str(projected_variant_dict['b37projected']), 'GRCh38HGVS': str(projected_variant_dict['b38projected'])}


# ---------------- CONVERT NC_HGVS to canonical SPDIs ---------------


def process_NC_HGVS(NC_HGVS):
    parsed_variant_dict = parse_variant(NC_HGVS)
    parsed_variant = parsed_variant_dict['parsed']
    print(f"parsed: {parsed_variant_dict['parsed']}")

    transcripts = b38hgvsAssemblyMapper.relevant_transcripts(parsed_variant)
    # Since we might not have access to all the transcripts that UTA contains, we select the "Longest Compatible
    # Remaining" transcript as described here: https://github.com/GenomicMedLab/cool-seq-tool/blob/main/docs/TranscriptSelectionPriority.md
    # - If there is a tie, choose the first-published transcript (lowest-numbered accession for RefSeq/Ensembl) among
    #   those transcripts meeting this criterion.
    #   Note: We want the most recent version of a transcript associated with an assembly.
    # Eventually, hgvs should include pyliftover for this purpose: https://github.com/biocommons/hgvs/issues/711
    nm_transcripts = [t for t in transcripts if 'NM_' in t]
    # Since a transcript looks like 'NM_006015.6', we sort them based on the number following 'NM_'
    ordered_transcripts = sorted(nm_transcripts, key=lambda t: int(t[3:t.find('.')]))
    # Pick the first transcript accession
    transcript_acc = ordered_transcripts[0].split('.')[0]
    # Search for its most recent version among nm_transcripts
    relevant_transcript = max((t for t in nm_transcripts if t.startswith(transcript_acc)), key=lambda t: int(t.split('.')[-1]))

    var_c = b38hgvsAssemblyMapper.g_to_c(
        parsed_variant, relevant_transcript)

    projected_variant_dict = project_variant(var_c)
    print(
        f"b37projected: {projected_variant_dict['b37projected']}; b38projected: {projected_variant_dict['b38projected']}")

    b37SPDI = normalize_variant(
        projected_variant_dict['b37projected'], 'GRCh37')
    b38SPDI = normalize_variant(
        projected_variant_dict['b38projected'], 'GRCh38')
    print(f"b37normalized: {b37SPDI}; b38normalized: {b38SPDI}")

    return {'GRCh37SPDI': b37SPDI, 'GRCh38SPDI': b38SPDI, 'GRCh37HGVS': str(projected_variant_dict['b37projected']), 'GRCh38HGVS': str(projected_variant_dict['b38projected'])}

# ---------------- CONVERT NM_SPDI to canonical SPDIs ---------------


def process_NM_SPDI(NM_SPDI):
    # convert SPDI into NM_HGVS then use NM_HGVS pipeline
    refSeq = NM_SPDI.split(":")[0]
    pos = int(NM_SPDI.split(":")[1])+1
    ref = NM_SPDI.split(":")[2]
    alt = NM_SPDI.split(":")[3]

    if len(ref) == len(alt) == 1:  # SNV
        var_n = hgvsParser.parse_hgvs_variant(
            refSeq+":n."+str(pos)+ref+">"+alt)
    elif len(ref) == 0:  # INS (e.g. NM_007294.3:c.5533_5534insG)
        start = pos-1
        end = start+1
        var_n = hgvsParser.parse_hgvs_variant(
            refSeq+":n."+str(start)+"_"+str(end)+'ins'+alt)
    elif len(alt) == 0:  # DEL (e.g. NM_000527.5:c.1350_1355del)
        start = pos
        end = start+len(ref)-1
        var_n = hgvsParser.parse_hgvs_variant(
            refSeq+":n."+str(start)+"_"+str(end)+'del')
    elif len(alt) != 0 and len(ref) != 0:  # DELINS (e.g. NM_007294.3:c.5359_5363delinsAGTGA)
        start = pos
        end = start+len(ref)-1
        var_n = hgvsParser.parse_hgvs_variant(
            refSeq+":n."+str(start)+"_"+str(end)+'delins'+alt)
    NM_HGVS = b38hgvsAssemblyMapper.n_to_c(var_n)

    return process_NM_HGVS(str(NM_HGVS))

# ---------------- CONVERT NC_SPDI to canonical SPDIs ---------------


def process_NC_SPDI(NC_SPDI):
    # convert SPDI into NC_HGVS then use NC_HGVS pipeline
    refSeq = NC_SPDI.split(":")[0]
    pos = int(NC_SPDI.split(":")[1])+1
    ref = NC_SPDI.split(":")[2]
    alt = NC_SPDI.split(":")[3]

    if len(ref) == len(alt) == 1:  # SNV
        NC_HGVS = (refSeq+":g."+str(pos)+ref+">"+alt)
    elif len(ref) == 0:  # INS (e.g. NM_007294.3:c.5533_5534insG)
        start = pos-1
        end = start+1
        NC_HGVS = (refSeq+":g."+str(start)+"_"+str(end)+'ins'+alt)
    elif len(alt) == 0:  # DEL (e.g. NM_000527.5:c.1350_1355del)
        start = pos
        end = start+len(ref)-1
        NC_HGVS = (refSeq+":g."+str(start)+"_"+str(end)+'del')
    elif len(alt) != 0 and len(ref) != 0:  # DELINS (e.g. NM_007294.3:c.5359_5363delinsAGTGA)
        start = pos
        end = start+len(ref)-1
        NC_HGVS = (refSeq+":g."+str(start)+"_"+str(end)+'delins'+alt)

    return process_NC_HGVS(str(NC_HGVS))


def normalize(variant):
    print(f"submitted: {variant}")
    if variant.upper().startswith('NM'):
        if variant.count(':') == 3:
            return process_NM_SPDI(variant)
        else:
            return process_NM_HGVS(variant)
    elif variant.upper().startswith('NC'):
        if variant.count(':') == 3:
            return process_NC_SPDI(variant)
        else:
            return process_NC_HGVS(variant)
