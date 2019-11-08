#! /usr/bin/bash
set -e -u -o pipefail

BASE_DIR=$(pwd)
CODE_DIR=${BASE_DIR}/bin
DATA_DIR=${BASE_DIR}/data
REF_DIR=${DATA_DIR}/reference
REF_BOWTIE_INDEX_DIR=${REF_DIR}

REF_FA=${REF_DIR}/genome.fa
REF_GTF=${REF_DIR}/genes.gtf
REF_UNIPROT=${REF_DIR}/UP000005640.fa

POS_DIR=${DATA_DIR}/pos
SGRNA_DIR=${DATA_DIR}/sgrna
SELECTED_SGRNA_DIR=${DATA_DIR}/selected

OFF_DIR=${DATA_DIR}/offtarget

####################

filelabels=(k_sgrna_no_splicing
    k_sgrna_splicing_12link
    k_sgrna_splicing_23link_1
    k_sgrna_splicing_23link_2)

lenlabels21=(19_21
    19_21
    18_21
    18_21)

lenlabels20=(13-14-18_20
    13-14-18_20
    13-17_20
    13-17_20)

lenlabels19=(15-16-17_19
    15-16-17_19
    14-15-16_19
    14-15-16_19)

for idx in `seq 1 4`; do
    filelabel=${filelabels[$idx]}
    len21=${lenlabels21[$idx]}
    len20=${lenlabels20[$idx]}
    len19=${lenlabels19[$idx]}

    Rscript ${CODE_DIR}/select_sgrna.r ${filelabel} \
        ${SELECTED_SGRNA_DIR} \
        ${SGRNA_DIR}/${filelabel}_${len19}.txt \
        ${OFF_DIR}/${filelabel}_${len19}.off \
        ${SGRNA_DIR}/${filelabel}_${len20}.txt \
        ${OFF_DIR}/${filelabel}_${len20}.off \
        ${SGRNA_DIR}/${filelabel}_${len21}.txt \
        ${OFF_DIR}/${filelabel}_${len21}.off

done
####################
