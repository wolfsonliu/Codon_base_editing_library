# Codon base editing library

This pipeline is for condon targeting base editing sgRNA library
design.

There are seven standard steps in total. But some steps can be changed
or omitted for different purposes. To run the codes, open each file,
and make sure all requird files and directories exist and all the
variables well set. The processing scripts stored in the `bin` directory.

Some reference files are required:
* genome fasta file
* protein fasta file (from uniprot proteomes)
* gene gtf file


1. `step1_aa_genome_location.sh`: get genomic locations for all codons from protein FASTA reference files.
2. `step2_get_sgrna.sh`: acquire all the possible sgRNA sequences.
3. `step3_offtarget_bowtie.sh`: evaluate the offtargets for the sgRNAs.
4. `step4_select_sgrna.sh`: filter sgRNAs by offtargets, GC contents and so on.
5. `step5_merge_library.r`: merge all the sgRNA files.
6. `step6_check.r`: quality check.
7. `step7_non-targeting_sgrna.sh`: generate non-targeting controls.
