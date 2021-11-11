library(dplyr)
library(readr)
library(purrr)

colnames <- c("chr", "start", "end", "strand", "motif", "annot", "n_unique", "n_multi", "max_overhang") 
cat.sjs <- map_dfr(snakemake@input, readr::read_tsv(col_names = colnames))
sum.sjs <- cat.sjs %>% 
    group_by(chr, start, end, strand, motif, annot) %>%
    summarise(
        n_unique = sum(n_unique),
        n_multi = sum(n_multi), 
        max_overhang = max(max_overhang)
        ) %>% 
        arrange(chr, start)

readr::write(sum.sjs, snakemake@output[[1]])
