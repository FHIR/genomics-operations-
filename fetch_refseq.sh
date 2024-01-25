#!/bin/sh

if [ -d ./refseq ]; then
    echo "refseq files already fetched."
    exit 0
fi

mkdir -p ./refseq
cd ./refseq

echo "Downloading refseq files..."

curl -sLO https://github.com/FHIR/genomics-operations/releases/download/113c119/GRCh37_refseq.tar.gz
curl -sLO https://github.com/FHIR/genomics-operations/releases/download/113c119/GRCh38_refseq.tar.gz
curl -sLO https://github.com/FHIR/genomics-operations/releases/download/113c119/rna_refseq.tar.gz

echo "Extracting refseq files..."

tar -xzf GRCh37_refseq.tar.gz
tar -xzf GRCh38_refseq.tar.gz
tar -xzf rna_refseq.tar.gz

echo "Finished extracting refseq files."

# TODO: REFACTOR THIS!
cd ..
mkdir -p ./pyard
cd ./pyard
curl -sLO https://github.com/FHIR/genomics-operations/releases/download/113c119/pyard.sqlite3.tar.gz

echo "Extracting pyard database..."
tar -xzf pyard.sqlite3.tar.gz
echo "Finished extracting pyard database."
