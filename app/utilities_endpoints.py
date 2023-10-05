from flask import abort, jsonify
from collections import OrderedDict
from app import common
from utilities import SPDI_Normalization


def get_feature_coordinates(
        chromosome=None, gene=None, transcript=None, protein=None):

    if not chromosome and not gene and not transcript and not protein:
        abort(400, "You must provide one parameter.")

    if chromosome:
        if gene or transcript or protein:
            abort(400, "Only one argument can be suuplied")

        chromosome = chromosome.strip()

        try:
            result = common.chromosomes_data.aggregate([{"$match": {"chr": {"$eq": f"chr{chromosome[3:].upper()}"}}}])
            result = list(result)
        except Exception as e:
            print(f"DEBUG: Error({e}) under get_feature_coordinates(chromosome={chromosome})")
            result = []

        if not result:
            return jsonify([])

        result = result[0]

        output = OrderedDict()
        output["chromosome"] = result["chr"]
        output["chromosomeLink"] = f'https://medlineplus.gov/genetics/chromosome/{chromosome[3:].lower() if chromosome[3:].lower() != "m" else "mitochondrial-dna"}/'
        output["build37RefSeq"] = result["build37RefSeq"]
        output["build38RefSeq"] = result["build38RefSeq"]

        return jsonify(output)

    if gene:
        if transcript or protein:
            abort(400, "Only one argument can be suuplied")

        gene = gene.strip().upper()

        if common.is_int(gene):
            gene = f"HGNC:{gene}"

        result = common.query_genes(gene)

        output = []

        for res in result:
            ord_dict = OrderedDict()

            ord_dict['geneId'] = res['HGNCgeneId']
            ord_dict['geneSymbol'] = res['HGNCgeneSymbol']
            ord_dict['geneLink'] = f"https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/{res['HGNCgeneId']}"
            ord_dict['build37Coordinates'] = f"{res['build37RefSeq']}:{res['build37Start']}-{res['build37End']}"
            ord_dict['build38Coordinates'] = f"{res['build38RefSeq']}:{res['build38Start']}-{res['build38End']}"
            ord_dict['transcripts'] = []
            ord_dict['MANE'] = []
            for transcript in res['transcriptMatches']:
                ord_dict['transcripts'].append(transcript['transcriptRefSeq'])
                if transcript['MANE'] == 1:
                    ord_dict['MANE'].append(transcript['transcriptRefSeq'])

            if not ord_dict['transcripts']:
                ord_dict.pop('transcripts')
            if not ord_dict['MANE']:
                ord_dict.pop('MANE')

            output.append(ord_dict)

        return (jsonify(output))

    if transcript:
        if protein:
            abort(400, "Only one argument can be suuplied")

        transcript = transcript.strip().upper()

        if '.' in transcript:
            transcript = transcript.split('.')[0]

        result = common.query_transcript(transcript)

        output = []

        for res in result:
            ord_dict = OrderedDict()

            ord_dict['transcript'] = res['transcriptRefSeq']
            ord_dict['transcriptLink'] = f"https://www.ncbi.nlm.nih.gov/nuccore/{res['transcriptRefSeq']}"

            if res['geneMatches']:
                ord_dict['geneSymbol'] = res['geneMatches'][0]['HGNCgeneSymbol']

            ord_dict['build37Coordinates'] = f"{res['build37RefSeq']}:{res['build37Start']}-{res['build37End']}"
            ord_dict['build38Coordinates'] = f"{res['build38RefSeq']}:{res['build38Start']}-{res['build38End']}"
            ord_dict['MANE'] = True if res['MANE'] == 1 else False

            if res['cdsMatches']:
                ord_dict['build38CDSStart'] = res['cdsMatches'][0]['build38CDSStart']
                ord_dict['build38CDSEnd'] = res['cdsMatches'][0]['build38CDSEnd']
                ord_dict['build38Strand'] = res['cdsMatches'][0]['build38Strand']
                ord_dict['build37CDSStart'] = res['cdsMatches'][0]['build37CDSStart']
                ord_dict['build37CDSEnd'] = res['cdsMatches'][0]['build37CDSEnd']
                ord_dict['build37Strand'] = res['cdsMatches'][0]['build37Strand']

            ord_dict['exons'] = []
            for exon in res['exonMatches']:
                exon_dict = OrderedDict()

                exon_dict['exonNumber'] = exon['exonNumber']
                exon_dict['build37Coordinates'] = f"{exon['build37RefSeq']}:{exon['build37Start']}-{exon['build37End']}"
                exon_dict['build38Coordinates'] = f"{exon['build38RefSeq']}:{exon['build38Start']}-{exon['build38End']}"

                ord_dict['exons'].append(exon_dict)

            if not ord_dict['exons']:
                ord_dict.pop('exons')
            else:
                ord_dict['exons'] = sorted(ord_dict['exons'], key=lambda d: d['exonNumber'])

            output.append(ord_dict)

        return (jsonify(output))

    if protein:

        protein = protein.strip().upper()

        if '.' in protein:
            protein = protein.split('.')[0]

        try:
            result = common.proteins_data.aggregate(
                [{"$match": {"proteinRefSeq": {'$regex': ".*"+str(protein).replace('*', r'\*')+".*"}}}])
            result = list(result)
        except Exception as e:
            print(f"DEBUG: Error({e}) under get_feature_coordinates(protein={protein})")
            result = []

        if len(result) > 1:
            abort(400, "Unable to provide information on this transcript at this time")

        if not result:
            return jsonify([])

        result = result[0]

        output = OrderedDict()
        output["protein"] = result["proteinRefSeq"]
        output["proteinLink"] = f'https://www.ncbi.nlm.nih.gov/protein/{output["protein"]}/'
        output["proteinName"] = result["proteinName"]
        output["transcript"] = result["transcript"]

        return (jsonify(output))


def find_the_gene(range=None):
    if not range:
        abort(400, "Range is required.")

    given_range = common.get_range(range)

    result = common.query_genes_range(given_range)

    output = []

    for res in result:
        ord_dict = OrderedDict()

        ord_dict['geneId'] = res['HGNCgeneId']
        ord_dict['geneSymbol'] = res['HGNCgeneSymbol']
        ord_dict['geneLink'] = f"https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/{res['HGNCgeneId']}"
        ord_dict['build37Coordinates'] = f"{res['build37RefSeq']}:{res['build37Start']}-{res['build37End']}"
        ord_dict['build38Coordinates'] = f"{res['build38RefSeq']}:{res['build38Start']}-{res['build38End']}"
        ord_dict['transcripts'] = []
        ord_dict['MANE'] = []
        for transcript in res['transcriptMatches']:
            ord_dict['transcripts'].append(transcript['transcriptRefSeq'])
            if transcript['MANE'] == 1:
                ord_dict['MANE'].append(transcript['transcriptRefSeq'])

        if not ord_dict['transcripts']:
            ord_dict.pop('transcripts')
        if not ord_dict['MANE']:
            ord_dict.pop('MANE')

        output.append(ord_dict)

    return (jsonify(output))


def seqfetcher(ref_seq, start, end):
    try:
        subseq = SPDI_Normalization.get_ref_seq_subseq('GRCh37', ref_seq, start, end)
    except Exception:
        subseq = SPDI_Normalization.get_ref_seq_subseq('GRCh38', ref_seq, start, end)
    return f'>{ref_seq}:{start}-{end} Homo sapiens chromosome 1, GRCh37.p13 Primary Assembly\n{subseq}\n\n'
