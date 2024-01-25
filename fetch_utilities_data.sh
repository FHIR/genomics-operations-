#!/bin/sh

set -eu

. ./.env

if [ -d ./data ]; then
    echo "files already fetched."
    exit 0
fi

mkdir -p ./data/refseq
(
    cd ./data/refseq

    echo "Downloading refseq files..."

    curl -sLO https://github.com/FHIR/genomics-operations/releases/download/${UTILITIES_DATA_VERSION}/GRCh37_refseq.tar.gz
    curl -sLO https://github.com/FHIR/genomics-operations/releases/download/${UTILITIES_DATA_VERSION}/GRCh38_refseq.tar.gz
    curl -sLO https://github.com/FHIR/genomics-operations/releases/download/${UTILITIES_DATA_VERSION}/rna_refseq.tar.gz

    echo "Extracting refseq files..."

    tar -xzf GRCh37_refseq.tar.gz
    tar -xzf GRCh38_refseq.tar.gz
    tar -xzf rna_refseq.tar.gz

    echo "Finished extracting refseq files."
)

mkdir -p ./data/pyard
(
    cd ./data/pyard

    echo "Downloading pyard database..."

    curl -sLO https://github.com/FHIR/genomics-operations/releases/download/${UTILITIES_DATA_VERSION}/pyard.sqlite3.tar.gz

    echo "Extracting pyard database..."

    tar -xzf pyard.sqlite3.tar.gz

    echo "Finished extracting pyard database."
)
