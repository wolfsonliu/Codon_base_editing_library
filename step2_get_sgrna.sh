#! /usr/bin/bash
set -e -u -o pipefail

BASE_DIR=$(pwd)
CODE_DIR=${BASE_DIR}/bin
DATA_DIR=${BASE_DIR}/data
REF_DIR=${DATA_DIR}/reference

REF_FA=${REF_DIR}/genome.fa
REF_GTF=${REF_DIR}/genes.gtf
REF_UNIPROT=${REF_DIR}/UP000005640.fa

POS_DIR=${DATA_DIR}/pos
SGRNA_DIR=${DATA_DIR}/sgrna

####################

# no splicing
#  1 seqname
#  7 strand
# 29 aa_genome_pos_start
# 30 aa_genome_pos_end
# 31 aa_genome_seq

awk -F "\t" '{print $1 FS $7 FS $29 FS $30 FS $31}' \
    ${DATA_DIR}/k_relationship_protein_genome_pos_no_splicing.txt | \
    grep -v seqname | sort | uniq > \
    ${POS_DIR}/k_genome_pos_no_splicing.txt

# splicing 12 link
awk -F "\t" '$2 ~ /+/ {
    print $1 FS $2 FS $14 FS $15 FS $17$18;
}
$2 ~ /-/ {
    print $1 FS $2 FS $15 FS $14 FS $18$17;
}' ${DATA_DIR}/k_relationship_protein_genome_pos_splicing_12link.txt \
    | grep -v seqname | sort | uniq > \
    ${POS_DIR}/k_genome_pos_splicing_12link.txt

for x in no_splicing splicing_12link; do
    python3 ${CODE_DIR}/get_sgrna.py \
        -r ${REF_FA} \
        -i ${POS_DIR}/k_genome_pos_${x}.txt \
        -a 19 -l 21 \
        -o ${SGRNA_DIR}/k_sgrna_${x}_19_21.txt &

    python3 ${CODE_DIR}/get_sgrna.py \
        -r ${REF_FA} \
        -i ${POS_DIR}/k_genome_pos_${x}.txt \
        -a 13 14 18 -l 20 \
        -o ${SGRNA_DIR}/k_sgrna_${x}_13-14-18_20.txt &

    python3 ${CODE_DIR}/get_sgrna.py \
        -r ${REF_FA} \
        -i ${POS_DIR}/k_genome_pos_${x}.txt \
        -a 15 16 17 -l 19 \
        -o ${SGRNA_DIR}/k_sgrna_${x}_15-16-17_19.txt
done

# 23link
awk -F "\t" '{
    print $1 FS $2 FS $14 FS $14 FS $17;
}' ${DATA_DIR}/k_relationship_protein_genome_pos_splicing_23link.txt \
    | grep -v seqname | sort | uniq > \
    ${POS_DIR}/k_genome_pos_splicing_23link_1.txt

awk -F "\t" '{
    print $1 FS $2 FS $15 FS $15 FS $18;
}' ${DATA_DIR}/k_relationship_protein_genome_pos_splicing_23link.txt \
    | grep -v seqname | sort | uniq > \
    ${POS_DIR}/k_genome_pos_splicing_23link_2.txt

for x in 1 2; do
    python3 ${CODE_DIR}/get_sgrna.py \
        -r ${REF_FA} \
        -i ${POS_DIR}/k_genome_pos_splicing_23link_${x}.txt \
        -a 18 -l 21 \
        -o ${SGRNA_DIR}/k_sgrna_splicing_23link_${x}_18_21.txt &

    python3 ${CODE_DIR}/get_sgrna.py \
        -r ${REF_FA} \
        -i ${POS_DIR}/k_genome_pos_splicing_23link_${x}.txt \
        -a 13 17 -l 20 \
        -o ${SGRNA_DIR}/k_sgrna_splicing_23link_${x}_13-17_20.txt &

    python3 ${CODE_DIR}/get_sgrna.py \
        -r ${REF_FA} \
        -i ${POS_DIR}/k_genome_pos_splicing_23link_${x}.txt \
        -a 14 15 16 -l 19 \
        -o ${SGRNA_DIR}/k_sgrna_splicing_23link_${x}_14-15-16_19.txt
done

####################
