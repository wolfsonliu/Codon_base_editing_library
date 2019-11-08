#! /usr/bin/bash
set -e -u -o pipefail

BASE_DIR=$(pwd)
CODE_DIR=${BASE_DIR}/bin
DATA_DIR=${BASE_DIR}/data
REF_DIR=${DATA_DIR}/reference

REF_FA=${REF_DIR}/genome.fa
REF_GTF=${REF_DIR}/genes.gtf
REF_UNIPROT=${REF_DIR}/UP000005640.fa

####################

# AA protein position
python3 ${CODE_DIR}/protein_aa_pos.py \
    -f ${REF_UNIPROT} -a K \
    -o ${DATA_DIR}/k_protein_pos.txt

# AA protein position mapped to genome position

python3 ${CODE_DIR}/protein_aa_pos_2_genome_pos.py \
    -f ${REF_FA} -g ${REF_GTF} \
    -p ${DATA_DIR}/k_protein_pos.txt \
    -o ${DATA_DIR}/k_relationship_protein_genome_pos_no_splicing.txt \
    -s ${DATA_DIR}/k_relationship_protein_genome_pos_splicing.txt

####################
