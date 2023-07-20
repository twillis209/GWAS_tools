library(data.table)
setDTthreads(snakemake@threads)

dat <- fread(snakemake@input[[1]], sep = '\t')

if(snakemake@params[['strand_policy']] == 'nonswitched') {
  dat <- dat[code %in% c('nochange', 'rev', 'ambig')]

  dat[code == 'rev', `:=` (REF = REF.1kG, ALT = ALT.1kG, BETA = -1 * BETA, Z = -1 * Z)]
  dat[(code == 'ambig' & (REF == ALT.1kG & ALT == REF.1kG)), `:=` (REF = REF.1kG, ALT = ALT.1kG, BETA = -1 * BETA, Z = -1 * Z)]
} else {
  stop('Behaviour not yet implemented for this strand policy')
}

dat[, c('code', 'REF.1kG', 'ALT.1kG') := NULL]

fwrite(dat, file = snakemake@output[['aligned_stats']], sep = '\t')
