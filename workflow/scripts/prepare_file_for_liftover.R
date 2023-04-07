library(data.table)
setDTthreads(snakemake@threads)

library(stringr)

dat <- fread(snakemake@input[['sumstats']], header = T, sep = '\t')
build <- fread(snakemake@input[['build']], header = T, sep = '\t')

dat[, BP := format(BP, scientific = F)]

colnames <- names(dat)

setcolorder(dat, c('SNPID', colnames[colnames != 'SNPID']))

assembly <- names(which.max(build[1]))

setnames(dat, c('CHR', 'BP'), c(paste0('CHR', str_remove(assembly, 'hg')), paste0('BP', str_remove(assembly, 'hg'))))

fwrite(dat, file = snakemake@output[['prepared']], sep = '\t')
fwrite(data.table(build = assembly), file = snakemake@output[['build']], sep = ' ', col.names = F)
