library(data.table)
setDTthreads(snakemake@threads)

dat <- fread(snakemake@input[['sumstats']], header = T, sep = '\t')

if('CHR38' %in% names(dat) & 'BP38' %in% names(dat)) {
  dat[, CHR38 := stringr::str_remove(CHR38, 'chr')]

  fwrite(dat, file = snakemake@output[[1]], sep = '\t')
} else {
  bedfile <- fread(snakemake@input[['lifted']], header = F, sep = '\t')

  names(bedfile) <- c('CHR38', 'BP38', 'BP2', 'SNPID')

  bedfile <- bedfile[, .(SNPID, CHR38, BP38)]

  dat[, CHR19 := stringr::str_remove(CHR19, 'chr')]

  dat <- dat[SNPID != '.']

  bedfile <- bedfile[SNPID != '.']

  dat <- merge(dat, bedfile, by = 'SNPID')

  dat[, CHR38 := stringr::str_remove(CHR38, 'chr')]

  dat[, SNPID := paste(CHR38, BP38, ALT, REF, sep = ':')]

  fwrite(dat, file = snakemake@output[[1]], sep = '\t')
}
