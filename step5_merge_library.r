library(dplyr)

dirs <- list(
    'selected'=file.path(getwd(), 'data', 'selected')
)

####################

synfilenames <- c(
    'k_sgrna_no_splicing-library_synthesis.txt',
    'k_sgrna_splicing_12link-library_synthesis.txt',
    'k_sgrna_splicing_23link_1-library_synthesis.txt',
    'k_sgrna_splicing_23link_2-library_synthesis.txt'
)

syn <- list()
for (x in synfilenames) {
    syn[[x]] <- read.table(
        file.path(dirs[['selected']], x), header=TRUE, stringsAsFactors=FALSE
    )
    syn[[x]] <- syn[[x]][syn[[x]]$m0 < 100 | syn[[x]]$number_in_library > 20,]
}

synthesis <- Reduce(
    rbind, syn
) %>% group_by(
          sgrna_seq, gc, m0, m1
      ) %>% summarize(number_in_library=sum(number_in_library))

write.table(
    synthesis,
    file.path(dirs[['selected']], 'k_library_synthesis.txt'),
    row.names=FALSE, quote=FALSE, sep='\t'
)

####################

filenames <- c(
    'k_sgrna_no_splicing-library.txt',
    'k_sgrna_splicing_12link-library.txt',
    'k_sgrna_splicing_23link_1-library.txt',
    'k_sgrna_splicing_23link_2-library.txt'
)

data <- list()

for (x in filenames) {
    data[[x]] <- read.table(
        file.path(dirs[['selected']], x),
        header=TRUE, stringsAsFactors=FALSE
    )
    data[[x]] <- data[[x]][data[[x]]$sgrna_seq %in% synthesis$sgrna_seq,]
}

fdata <- Reduce(rbind, data)

write.table(
    fdata,
    file.path(dirs[['selected']], 'k_library.txt'),
    row.names=FALSE, quote=FALSE, sep='\t'
)
####################
