library(data.table)
setDTthreads(snakemake@threads)

dat <- fread(snakemake@input[[1]], header = T, sep = '\t')
dat <- dat[, .(CHR, BP, BP2 = BP+1, SNPID)]
dat[, `:=` (CHR = paste0('chr', CHR), BP = format(BP, scientific = F), BP2 = format(BP2, scientific = F))]

fwrite(dat[, .(CHR, BP, BP2, SNPID)], file = snakemake@output[[1]], sep = '\t', col.names = F)
