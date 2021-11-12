library(dplyr)
library(readr)
library(purrr)

colnames <- c("chr", "start", "end", "strand", "motif", "annot", "n_unique", "n_multi", "max_overhang") 
coltypes = "cdddddddd"
cat.sjs <- map_dfr(snakemake@input, ~ readr::read_tsv(.x, col_names = colnames, col_types = coltypes))
sum.sjs <- cat.sjs %>% 
    group_by(chr, start, end, strand, motif, annot) %>%
    summarise(
        n_unique = sum(n_unique),
        n_multi = sum(n_multi), 
        max_overhang = max(max_overhang)
        ) %>% 
        arrange(chr, start)

readr::write_tsv(sum.sjs, snakemake@output[[1]], col_names = F)
