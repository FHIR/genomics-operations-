import vcf2json
import json
import os
import common
import pandas as pd
from tqdm import tqdm
from collections import OrderedDict
import parse_non_variant
import uuid
from common import _get_ref_seq_by_chrom, extract_chrom_identifier

def extractData(csv_file_path, variants_data, molecular_output, phase_output):
    try:
        # Fetching all Transcripts data and storing it in transcript_map
        transcript_map = {}
        common.query_genes(transcript_map)

        df = pd.read_csv(csv_file_path)
        for _, row in tqdm(df.iterrows(), total=len(df)-1, desc="Processing Data"):
            phased_rec_map = {}
            vcf2json.vcf2json(row['vcf_filename'],
                              row['ref_build'],
                              row['patient_id'],
                              row['test_date'],
                              row['test_id'],
                              row['specimen_id'],
                              row['genomic_source_class'],
                              row['ratio_ad_dp'],
                              row['sample_position'], transcript_map, variants_data, molecular_output, phased_rec_map)

            file = pd.read_json("convertedVCF.json", orient=str)

            vcf2json.add_phased_relationship_obv(row['patient_id'],
                                                 row['test_id'],
                                                 row['specimen_id'],
                                                 row['ref_build'], phase_output, file, phased_rec_map)

    except Exception as e:
        print("Error with csv file:", e)


def extract_non_variant_data(csv_file_path):
    """
    Extracts non-variant data from a CSV file and generates combined JSON outputs.

    Args:
        csv_file_path (str): Path to the CSV file containing GVCF metadata.

    Returns:
        tuple: Lists containing non-variant JSON data and BED JSON entries.
    """
    df = pd.read_csv(csv_file_path)
    processed_files = set()  # To track processed GVCF files
    bed_entries = []  # Collect BED entries for all files
    non_variant_entries = []  # Collect non-variant entries for all files

    for _, row in tqdm(df.iterrows(), total=len(df), desc="Processing GVCF Data"):
        gvcf_filename = row['vcf_filename']
        ref_build = row['ref_build']  # Dynamically set build version

        # Skip already processed files
        if gvcf_filename in processed_files:
            continue

        # Extract BED entries and detailed records from the GVCF file
        if gvcf_filename.lower().endswith(".g.vcf"):

            studied_regions, detailed_blocks = parse_non_variant.parse_gvcf(
                    gvcf_filename, int(ref_build[-2:])
                )  # GRCh37 or GRCh38

            # Update studied regions (BED data) with metadata
            for region in studied_regions:
                region.update({
                    "REF_BUILD": ref_build,
                    "PATIENT_ID": row['patient_id'],
                    "TEST_DATE": row['test_date'],
                    "TEST_ID": row['test_id'],
                    "SPECIMEN_ID": row['specimen_id'],
                    "GENOMIC_SOURCE_CLASS": row['genomic_source_class'],
                    "RATIO_AD_DP": row['ratio_ad_dp'],
                    "SAMPLE_POSITION": row['sample_position'],
                })
                bed_entries.append(region)

            # Update detailed blocks with metadata
            for block in detailed_blocks:
                block.update({
                    "REF_BUILD": ref_build,
                    "PATIENT_ID": row['patient_id'],
                    "TEST_DATE": row['test_date'],
                    "TEST_ID": row['test_id'],
                    "SPECIMEN_ID": row['specimen_id'],
                })
                non_variant_entries.append(block)

            # Mark file as processed
            processed_files.add(gvcf_filename)

        if gvcf_filename.lower().endswith(".vcf") and gvcf_filename not in processed_files:
            vcf_reader = vcf.Reader(filename=gvcf_filename)

            for record in tqdm(vcf_reader, desc=f"Parsing {gvcf_filename}"):
                chrom = f"chr{record.CHROM}"
                pos = record.POS
                end = record.INFO.get("END", pos)  # Use END from INFO if available, else use POS
                ref_allele = record.REF
                alt_allele = record.ALT[0] if record.ALT else ref_allele
                filter_status = record.FILTER[0] if record.FILTER else None

                # Extract genotype (GT) from the first sample in the VCF
                gt = "0" if chrom == "chrM" else record.samples[0].data.GT

                # Determine allelic state
                allelic_state = determine_allelic_state(chrom, gt)

                # Collect metadata
                non_variant_entry = {
                    "CHROM": chrom,
                    "POS": pos,
                    "END": end,
                    "REF_ALLELE": ref_allele.upper(),  # Capitalized REF_ALLELE
                    "ALT": alt_allele.upper(),
                    "FILTER": filter_status,
                    "GT": gt,
                    "allelicState": allelic_state,
                    "REF_BUILD": ref_build,
                    "PATIENT_ID": row["patient_id"],
                    "TEST_DATE": row["test_date"],
                    "TEST_ID": row["test_id"],
                    "SPECIMEN_ID": row["specimen_id"],
                }

                # Add the entry to the non-variant entries list
                non_variant_entries.append(non_variant_entry)

    # Mark file as processed
    processed_files.add(gvcf_filename)



    return non_variant_entries, bed_entries

def transform_non_variant_output(non_variant_data, bed_file_output):
    """
    Transforms non-variant data and BED file data to create the required output format,
    including SPDI for non-variant data.

    Args:
        non_variant_data (list): List of non-variant entries with START and REF_ALLELE.
        bed_file_output (list or None): List of BED file entries with Chromosome, START, END, and PATIENT_ID.

    Returns:
        tuple: Transformed non-variant data and transformed BED data as a list.
    """
    transformed_data = []
    transformed_bed_data = []

    # Check if bed_file_output is empty or None
    if not bed_file_output:
        print("No BED file output provided. Skipping BED transformation.")
        return transformed_data, transformed_bed_data  # Return empty BED data

    # Create a mapping of (CHROM, START, END) regions from bed_file_output
    bed_region_map = []
    for record in bed_file_output:
        bed_region_map.append({
            "CHROM": record.get("CHROM"),
            "START": record.get("START"),
            "END": record.get("END"),
            "PATIENT_ID": record.get("PATIENT_ID")
        })

    # Transform Non-Variant Data
    for record in tqdm(non_variant_data, desc="Transforming Non-Variant Data"):
        chrom = record.get("CHROM")  # Safely get CHROM to avoid KeyError
        ref_allele = record.get("REF_ALLELE", "")
        FILTER = record.get("FILTER", "")
        ref_build = record.get("REF_BUILD", "")
        patient_id = record.get("PATIENT_ID", "")
        test_date = record.get("TEST_DATE", "")
        test_id = record.get("TEST_ID", "")
        specimen_id = record.get("SPECIMEN_ID", "")
        gt = record.get("GT", "")
        allelicState = record.get("allelicState", "")
        pos = record.get("POS", 0)  # Default to 0 if POS is missing

        # Fetch reference sequence identifier
        ref_seq = _get_ref_seq_by_chrom(ref_build, extract_chrom_identifier(chrom))

        # Find the matching END from bed_file_output
        matched_end = None
        for bed_region in bed_region_map:
            if (
                bed_region["CHROM"] == chrom
                and bed_region["PATIENT_ID"] == patient_id
                and bed_region["START"] <= pos <= bed_region["END"]
            ):
                matched_end = bed_region["END"]
                break

        # If no match is found, use the default END from the non-variant record
        end = matched_end if matched_end is not None else record.get("END", 0)

        # Create a record for each position using REF_ALLELE length
        for idx, allele in enumerate(ref_allele):
            pos_adjusted = pos + idx  # Adjust position for each allele
            alt = allele  # For non-variants, ALT matches REF

            # Calculate SPDI
            spdi = f"{ref_seq}:{pos_adjusted - 1}:{allele}:{alt}"

            transformed_data.append({
                "_id": uuid.uuid4().hex,  # Assign unique identifier
                "CHROM": chrom,
                "REF": allele,  # Single allele
                "ALT": alt,  # Same as REF for non-variants
                "POS": pos_adjusted,  # Increment position
                "END": end,  # Use the matched END
                "GT": gt,
                "FILTER": FILTER,
                "allelicState": allelicState,
                "genomicBuild": ref_build,
                "SPDI": spdi,  # Add SPDI
                "patientID": patient_id,
                "testDate": test_date,
                "testID": test_id,
                "specimenID": specimen_id,
            })

    # Transform BED File Data
    patient_regions = {}

    for record in tqdm(bed_file_output, desc="Transforming BED File Data"):
        chrom = record.get("CHROM")  # Correct case-sensitive key
        start = record.get("START", 0)
        end = record.get("END", 0)
        patient_id = record.get("PATIENT_ID", "")

        # Initialize a new list for the patient if not already present
        if patient_id not in patient_regions:
            patient_regions[patient_id] = []

        # Add the region to the patient's list
        patient_regions[patient_id].append({
            "Chromosome": chrom,
            "Start": start,
            "End": end,
        })

    # Convert grouped patient regions into the final list format
    for patient_id, regions in patient_regions.items():
        transformed_bed_data.append({
            "BedID": f"BED_{patient_id}",
            "BED": regions
        })

    return transformed_data, transformed_bed_data

def run_vcf2json():
    variants_data = []
    molecular_output = []
    phase_output = []

    # Set to True to include non-variant JSON output
    allow_non_variants = True
    files_to_write = []

    if allow_non_variants:
        non_variant_output, bed_file_output = extract_non_variant_data("vcfData.csv")

        # Transform non-variant data into the desired format
        transformed_data, transformed_data2 = transform_non_variant_output(non_variant_output, bed_file_output)

        # Add non-variant files to write list
        files_to_write.append(("non_variant_output.json", transformed_data))  # Always add transformed data

        # Only add transformed_data2 if it's not empty
        if transformed_data2 and transformed_data2 != [] and transformed_data2 != {}:
            files_to_write.append(("bed_file_output.json", transformed_data2))

    extractData("vcfData.csv", variants_data, molecular_output, phase_output)

    # Add variant data files to the write list
    files_to_write.extend([
        ("variantsData.json", variants_data),
        ("molecularConsequences.json", molecular_output),
        ("phaseData.json", phase_output),
    ])

    # Write files
    for filename, data in tqdm(files_to_write, desc="Writing files"):
        with open(filename, "w") as f:
            f.write(json.dumps(data, indent=4))

    # Clean up temporary files
    try:
        os.remove("convertedVCF.json")
    except FileNotFoundError:
        pass

    print("Data Generated.")


if __name__ == "__main__":
    run_vcf2json()
