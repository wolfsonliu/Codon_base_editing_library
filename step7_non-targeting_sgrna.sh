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

NON_TARGET_DIR=${DATA_DIR}/non-target

####################

# script to find the offtarget sites of sgRNAs using blat software

THREAD=1

# setting variables
# scripts
chromosomes=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11
chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX
chrY)

SGLENGTHS=(19 20 21)

for sglength in ${SGLENGTHS[*]}; do
    # Data and ouput
    OUT_DIR=${NON_TARGET_DIR}/result_${sglength}

    OUT_FILE=${NON_TARGET_DIR}/non-target_${sglength}nt.txt

    if [ ! -e ${OUT_DIR} ]; then
        mkdir -p ${OUT_DIR}
    fi

    NUMBERS=0

    while [ "${NUMBERS}" -le 3000 ]; do
        ${CODE_DIR}/randomseq.py ${sglength} 1000000 | \
            awk '{print ">"$0"\n"$0;}' \
            > ${OUT_DIR}/random_sequence_${sglength}nt.fa

        bowtie -k 1 -v 3 -p ${THREAD} -f \
            --un ${OUT_DIR}/unalign_${sglength}.txt \
            ${REF_BOWTIE_INDEX_DIR}/genome \
            ${OUT_DIR}/random_sequence_${sglength}nt.fa \
            ${OUT_DIR}/random_${sglength}nt.out

        awk 'FNR % 2 == 0 {print $0;}' \
            ${OUT_DIR}/unalign_${sglength}.txt >> ${OUT_FILE}

        rm -f ${OUT_DIR}/unalign_${sglength}.txt
        rm -f ${OUT_DIR}/random_${sglength}nt.out

        NUMBERS=`wc -l ${OUT_FILE} | cut -d " " -f 1`
    done

    rm -f ${OUT_DIR}/random_sequence_${sglength}nt.fa

done
####################
