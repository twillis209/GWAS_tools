library(data.table)
setDTthreads(snakemake@threads)

dat <- fread(snakemake@input[[1]], sep = '\t')

dat <- dat[code %in% snakemake@params[['strand_policy']][['codes_to_retain']]]

if('ambig' %in% snakemake@params[['strand_policy']][['flip']]) {
  # Drop ambig variants which neither need flipping nor can be fixed with a flip
  dat <- dat[!(code == 'ambig' & !((REF == REF.1kG & ALT == ALT.1kG) | (REF == ALT.1kG & ALT == REF.1kG)))]
  dat[code == 'ambig' & REF == ALT.1kG & ALT == REF.1kG, `:=` (REF = REF.1kG, ALT = ALT.1kG, BETA = -1 * BETA)]
  # NB: don't need to change other class of ambig variants as alleles already match
}

if('rev' %in%  snakemake@params[['strand_policy']][['flip']]) {
  dat[code == 'rev', `:=` (REF = REF.1kG, ALT = ALT.1kG, BETA = -1 * BETA)]
}

if('comp' %in%  snakemake@params[['strand_policy']][['flip']]) {
  dat[code == 'rev', `:=` (REF = REF.1kG, ALT = ALT.1kG)]
}

if('revcomp' %in%  snakemake@params[['strand_policy']][['flip']]) {
  dat[code == 'revcomp', `:=` (REF = REF.1kG, ALT = ALT.1kG, BETA = -1 * BETA)]
}

dat[, c('code', 'REF.1kG', 'ALT.1kG', 'A.CODE', 'B.CODE') := NULL]

fwrite(dat, file = snakemake@output[[1]], sep = '\t')
