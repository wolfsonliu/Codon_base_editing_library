args <- commandArgs(trailingOnly=TRUE)


splicing <- read.table(
    args[1], 
    header=TRUE, stringsAsFactors=FALSE, sep='\t'
)

splicing$codon <- paste0(
    splicing$aa_genome_seq_1,
    splicing$aa_genome_seq_2,
    splicing$aa_genome_seq_3
)

write.table(
    splicing[
        abs(splicing$aa_genome_pos_2 - splicing$aa_genome_pos_1) == 1,
        c(
            "seqname", "strand", "gene_id", "gene_name", "p_id", "transcript_id",
            "uniprot_1", "uniprot_2", "aa", "aa_prot_pos",
            "aa_cds_pos_1", "aa_cds_pos_2", "aa_cds_pos_3",
            "aa_genome_pos_1", "aa_genome_pos_2", "aa_genome_pos_3",
            "aa_genome_seq_1", "aa_genome_seq_2", "aa_genome_seq_3"
        )
    ],
    file=args[2],
    row.names=FALSE, quote=FALSE, sep='\t'
)

write.table(
    splicing[
        abs(splicing$aa_genome_pos_2 - splicing$aa_genome_pos_1) != 1,
        c(
            "seqname", "strand", "gene_id", "gene_name", "p_id", "transcript_id",
            "uniprot_1", "uniprot_2", "aa", "aa_prot_pos",
            "aa_cds_pos_1", "aa_cds_pos_2", "aa_cds_pos_3",
            "aa_genome_pos_1", "aa_genome_pos_2", "aa_genome_pos_3",
            "aa_genome_seq_1", "aa_genome_seq_2", "aa_genome_seq_3"
        )
    ],
    file=args[3],
    row.names=FALSE, quote=FALSE, sep='\t'
)
