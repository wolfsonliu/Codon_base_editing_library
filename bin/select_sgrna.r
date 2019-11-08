#! /usr/bin/env Rscript
library(dplyr)
library(tidyr)

reverse_complement <- function(x) {
    compvect <- c('A'='T', 'T'='A', 'G'='C', 'C'='G')
    splitx <- unlist(strsplit(x, ''))
    paste0(rev(compvect[splitx]), collapse='')
}

####################
args <- commandArgs(trailingOnly = TRUE)

filelabel <- args[1]
outdir <- args[2]

sgrnafiles <- list()
offtargetfiles <- list()

sgrnafiles[[1]] <- args[3]
offtargetfiles[[1]] <- args[4]
sgrnafiles[[2]] <- args[5]
offtargetfiles[[2]] <- args[6]
sgrnafiles[[3]] <- args[7]
offtargetfiles[[3]] <- args[8]

rawdata <- list()

for (i in seq(length(sgrnafiles))) {

    off <- read.table(
        offtargetfiles[[i]], header=FALSE, stringsAsFactors=FALSE, sep='\t'
    )
    colnames(off) <- c(
        'name', 'strand', 'seqname', 'start', 'seq',
        'quality', 'n', 'mismatch', 'pstart', 'pend', 'pam'
    )
    off$name <- toupper(off$name)

    off$pam <- toupper(off$pam)
    off <- off[grepl('GG$', off$pam),]
    off$mismatch.num <- unlist(
        lapply(off$mismatch, function(z) {sum(grepl(':',z))})
    )

    offstat <- spread(
        off %>% group_by(name, mismatch.num) %>% summarize(n=n()),
        'mismatch.num', 'n'
    )
    colnames(offstat) <- c('sgrna_seq', 'm0', 'm1')
    ## seqname strand aa_genome_pos_start aa_genome_pos_end
    ## aa_genome_seq sgrna_genome_start sgrna_genome_end
    ## sgrna_genome_seq genome_pam sgrna_seq pam codon_1stnt_pos gc
    readdata <- read.table(
        sgrnafiles[[i]], header=TRUE, stringsAsFactors=FALSE, sep='\t'
    )

    rawdata[[i]] <- replace_na(
        left_join(
            readdata, offstat, by='sgrna_seq'
        ),
        list('m0'=0, 'm1'=0)
    )

    rm(off)
    rm(offstat)
    gc()
}

####################
## sgRNA Selection
data <- Reduce(rbind, rawdata)
rm(rawdata)

## Remove TTTT
data <- data[!grepl('TTTT', data$sgrna_seq),]

## GC content
data <- data[data$gc >= 0.2 & data$gc <= 0.8,]

data$genome_browser_search <- paste0(
    data$seqname, ':', data$sgrna_genome_start, '-', data$sgrna_genome_end
)

write.table(
    data, file.path(outdir, paste(filelabel, 'library.txt', sep='-')),
    row.names=FALSE,quote=FALSE, sep='\t'
)

data <- data[
    order(
        data$seqname, data$aa_genome_pos_start,
        data$sgrna_genome_start, nchar(data$sgrna_seq),
        decreasing=FALSE
    ),
]

## sgRNA per AA

data.sitensg <- data %>% group_by(
                        seqname, aa_genome_pos_start
                    ) %>% mutate(n=n())

data.selected <- data.sitensg[data.sitensg$n <= 3,]

data.moresg <-  data.sitensg[data.sitensg$n > 3,]

if (dim(data.moresg)[1] > 0) {
    data.moresg$select.score <- 0

    ## score with 4 continue same nt
    data.moresg$select.score <- data.moresg$select.score + ifelse(
                                                               grepl(
                                                                   '(AAAA|GGGG|CCCC)',
                                                                   data.moresg$sgrna_seq
                                                               ),
                                                               2, 1
                                                           )

    ## score with offtarget
    data.moresg$select.score <- data.moresg$select.score + (
        4 * unlist(lapply(data.moresg$m0 - 1, max, 0)) + data.moresg$m1
    )

    select_guide <- function(x) {
        x[
            order(
                nchar(x$sgrna_seq),
                x$select.score,
                rnorm(dim(x)[1]),
                decreasing=FALSE
            )[1:3],
            ]
    }

    data.3sg <- data.moresg %>% group_by(
                                    seqname, aa_genome_pos_start
                                ) %>% do(
                                          select_guide(.data)
                                      ) %>% ungroup()

    data.3sg$select.score <- NULL

    result <- rbind(as.data.frame(data.selected), as.data.frame(data.3sg))

    result$n <- NULL
} else {
    result <- data
}

## output synthesis file without duplicated sequences.
output <- result[
    c('sgrna_seq', 'gc', 'm0', 'm1')
] %>% group_by(
          sgrna_seq, gc, m0, m1
      ) %>% summarize(number_in_library=n())

write.table(
    output,
    file.path(outdir, paste(filelabel, 'library_synthesis.txt', sep='-')),
    row.names=FALSE, quote=FALSE, sep='\t'
)

################################################################################
