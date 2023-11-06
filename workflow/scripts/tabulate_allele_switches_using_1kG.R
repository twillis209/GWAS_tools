library(data.table)
setDTthreads(snakemake@threads)
library(pidPipelineCode)

leg_dat <- fread(snakemake@input[['legend']], sep = ' ', header = F)
names(leg_dat) <- c('ID', 'BP', 'ALT.1kG', 'REF.1kG')

leg_dat[, 'CHR' := tstrsplit(ID, split = ':', keep = 1)]

dat <- fread(snakemake@input[['sum_stats']], sep = '\t', header = T)

dat[, snakemake@params[['chr_col']] := as.character(get(snakemake@params[['chr_col']]))]

dat <- merge(dat, leg_dat[, .(CHR, BP, REF.1kG, ALT.1kG)], by.x = c(snakemake@params[['chr_col']], snakemake@params[['bp_col']]), by.y = c('CHR', 'BP'), all.x = T)

dat[, A.CODE := paste(REF, ALT, sep = '/')]
dat[, B.CODE := paste(REF.1kG, ALT.1kG, sep = '/')]
dat[, code := g.class(A.CODE, B.CODE)]

fwrite(dat, file = snakemake@output[['merged_stats_with_code']], sep = '\t')

fwrite(data.table(code = g.class(dat[,A.CODE], dat[, B.CODE]))[, .N, by = code], file = snakemake@output[['table']], sep = '\t')
