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

echo "Extracting refseq files..."

tar -xzf GRCh37_refseq.tar.gz
tar -xzf GRCh38_refseq.tar.gz

echo "Finished extracting refseq files."
