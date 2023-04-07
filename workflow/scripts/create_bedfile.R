library(data.table)
setDTthreads(snakemake@threads)

dat <- fread(snakemake@input[[1]], header = T, sep = '\t')

fwrite(dat[, .(CHR, BP, BP+1, SNPID)], file = snakemake@output[[1]], sep = '\t', col.names = F)
