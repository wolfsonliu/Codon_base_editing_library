args <- commandArgs(trailingOnly = TRUE)

dirs <- list(
    'selected'=file.path(getwd(), 'data', 'selected'),
    'checked'=file.path(getwd(), 'data', 'checked'),
)

####################

data <- list()

data[['no_splicing']] <- read.table(
    file.path(dirs[['selected']], 'k_sgrna_no_splicing-library.txt'),
    header=TRUE, stringsAsFactors=FALSE
)

data[['splicing_12link']] <- read.table(
    file.path(dirs[['selected']],
              'k_sgrna_splicing_12link-library.txt'),
    header=TRUE, stringsAsFactors=FALSE
)

data[['splicing_23link_1']] <- read.table(
    file.path(dirs[['selected']],
              'k_sgrna_splicing_23link_1-library.txt'),
    header=TRUE, stringsAsFactors=FALSE
)

data[['splicing_23link_2']] <- read.table(
    file.path(dirs[['selected']],
              'k_sgrna_splicing_23link_2-library.txt'),
    header=TRUE, stringsAsFactors=FALSE
)

for (x in names(data)) {
    for (i in seq(10)) {
        write.table(
            data[[x]][
                c(
                    sample(which(nchar(data[[x]]$sgrna_seq) == 19),
                           40, replace=FALSE),
                    sample(which(nchar(data[[x]]$sgrna_seq) == 20),
                           40, replace=FALSE),
                    sample(which(nchar(data[[x]]$sgrna_seq) == 21),
                           20, replace=FALSE)
                ),
                ],
            file.path(dirs[['check']], paste(x, i, 'txt', sep='.')),
            row.names=FALSE, quote=FALSE, sep='\t'
        )
    }
}

####################
