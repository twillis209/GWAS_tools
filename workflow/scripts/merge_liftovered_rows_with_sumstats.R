library(data.table)
setDTthreads(snakemake@threads)

dat <- fread(snakemake@input[['sumstats']], header = T, sep = '\t')
bedfile <- fread(snakemake@input[['lifted']], header = F, sep = '\t')
build <- fread(snakemake@input[['build']], header = T, sep = '\t')

assembly <- names(which.max(build[1]))

names(bedfile) <- c('CHR38', 'BP38', 'BP2', 'SNPID')

bedfile <- bedfile[, .(SNPID, CHR38, BP38)]

dat <- merge(dat, bedfile, by = 'SNPID')

dat[, CHR38 := stringr::str_remove(CHR38, 'chr')]

if(assembly == 'hg38') {
  dat[, c('CHR', 'BP') := NULL]
}

fwrite(dat, file = snakemake@output[[1]], sep = '\t')
