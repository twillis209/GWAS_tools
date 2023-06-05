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

  dat <- merge(dat, bedfile, by = 'SNPID')

  dat[, CHR38 := stringr::str_remove(CHR38, 'chr')]

  # Update chr:bp:alt:ref SNPID to hg38 coordinates
  if(dat[SNPID %like% "\\d+:\\d+:\\w+:\\w:", .N] > 100) {
    dat[, SNPID := paste(CHR38, BP38, ALT, REF, sep = '_')]
  }

  fwrite(dat, file = snakemake@output[[1]], sep = '\t')
}
