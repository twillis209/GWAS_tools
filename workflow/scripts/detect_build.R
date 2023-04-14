library(data.table)
setDTthreads(snakemake@threads)

dat <- fread(snakemake@input[['sumstats']], sep = '\t', header = T)

manifest <- fread(snakemake@input[['manifest']], sep = '\t', header = T)

dat[, CHR := as.character(CHR)]
manifest[, `:=` (CHR18 = as.character(CHR18), CHR19 = as.character(CHR19), CHR38 = as.character(CHR38))]

hg18 <- merge(dat[, .(CHR, BP)], manifest[, .(CHR18, BP18)], by.x = c('CHR', 'BP'), by.y = c('CHR18', 'BP18'))
hg19 <- merge(dat[, .(CHR, BP)], manifest[, .(CHR19, BP19)], by.x = c('CHR', 'BP'), by.y = c('CHR19', 'BP19'))
hg38 <- merge(dat[, .(CHR, BP)], manifest[, .(CHR38, BP38)], by.x = c('CHR', 'BP'), by.y = c('CHR38', 'BP38'))

fwrite(data.table(hg18 = hg18[, .N], hg19 = hg19[, .N], hg38 = hg38[, .N]), file = snakemake@output[[1]], sep = '\t')
