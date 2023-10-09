library(data.table)
setDTthreads(snakemake@threads)

dat <- fread(snakemake@input[[1]], sep = '\t')

if(snakemake@params[['strand_policy']] == 'nonswitched') {
  dat <- dat[code %in% c('nochange', 'rev', 'ambig')]

  dat[code == 'rev', `:=` (REF = REF.1kG, ALT = ALT.1kG, BETA = -1 * BETA)]

  if('OR' %in% names(dat)) {
    dat[code == 'rev', OR := 1/OR]
  }

  if('Z' %in% names(dat)) {
    dat[code == 'rev', Z := -1 * Z]
  }

  if('OR' %in% names(dat)) {
    dat[(code == 'ambig' & (REF == ALT.1kG & ALT == REF.1kG)), OR := 1/OR]
  }

  if('Z' %in% names(dat)) {
    dat[(code == 'ambig' & (REF == ALT.1kG & ALT == REF.1kG)), Z := -1 * Z]
  }

  dat[(code == 'ambig' & (REF == ALT.1kG & ALT == REF.1kG)), `:=` (REF = REF.1kG, ALT = ALT.1kG, BETA = -1 * BETA)]

} else {
  stop('Behaviour not yet implemented for this strand policy')
}

dat[, c('code', 'REF.1kG', 'ALT.1kG', 'A.CODE', 'B.CODE') := NULL]

fwrite(dat, file = snakemake@output[[1]], sep = '\t')
