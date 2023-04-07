library(data.table)
setDTthreads(snakemake@threads)

dat <- fread(snakemake@input[['sumstats']], header = T, sep = '\t')
bedfile <- fread(snakemake@input[['lifted']], header = F, sep = '\t')

names(bedfile) <- c('CHR38', 'BP38', 'BP38+1', 'SNPID')

bedfile <- bedfile[, .(SNPID, CHR38, BP38)]
# TODO not sure how to use the SNPID column in the merge

