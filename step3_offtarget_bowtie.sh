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

OFF_DIR=${DATA_DIR}/offtarget

####################

# functions
function generate_fa {
    awk -F "\t" 'FNR > 1 {print $10;}' ${1} | \
        sort | uniq | awk '{print ">"$0"\n"$0;}'
}

function run_bowtie {
    # ${1} should be the reference basename
    # ${2} should be the input fa file path
    # ${3} should be the output file path
    bowtie -a -v 1  -p ${THREAD} -f ${1} ${2} ${3}
}

####################

THREAD=4

DATA_DIR=${BASE_DIR}/sgrna

# reference

export chromosomes=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10
chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21
chr22 chrX chrY)

filenames=(k_sgrna_no_splicing_13-14-18_20.txt
k_sgrna_no_splicing_15-16-17_19.txt
k_sgrna_no_splicing_19_21.txt
k_sgrna_splicing_12link_13-14-18_20.txt
k_sgrna_splicing_12link_15-16-17_19.txt
k_sgrna_splicing_12link_19_21.txt
k_sgrna_splicing_23link_1_13-17_20.txt
k_sgrna_splicing_23link_1_14-15-16_19.txt
k_sgrna_splicing_23link_1_18_21.txt
k_sgrna_splicing_23link_2_13-17_20.txt
k_sgrna_splicing_23link_2_14-15-16_19.txt
k_sgrna_splicing_23link_2_18_21.txt)

####################
# sgRNA design info
# start site of nt in sgrna to check the offtarget
#   -- for sgRNA with length 21, SS = 4, to generate 18nt mapping sequence
#   -- for sgRNA with length 20, SS = 3, to generate 18nt mapping sequence
#   -- for sgRNA with length 19, SS = 2, to generate 18nt mapping sequence

for file in ${filenames[*]}; do

    # generate fa for bowtie
    generate_fa ${SGRNA_DIR}/${file} > ${OFF_DIR}/${file/txt/fa}

    # running bowtie
    echo "==== Bowtie started - ["`date -Iminutes`"]"
    run_bowtie ${REF_BOWTIE_INDEX_DIR}/genome \
        ${OFF_DIR}/${file/txt/fa} \
        ${OFF_DIR}/${file/txt/out}
    echo "==== Bowtie finished - ["`date -Iminutes`"]"
    # get
    python3 ${CODE_DIR}/get_pam_by_coord_bowtie.py \
        -r ${REF_FA} \
        -i ${OFF_DIR}/${file/txt/out} \
        -o ${OFF_DIR}/${file/txt/off}
done

####################
