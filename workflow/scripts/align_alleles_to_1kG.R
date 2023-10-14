library(data.table)
setDTthreads(snakemake@threads)

dat <- fread(snakemake@input[[1]], sep = '\t')

dat <- dat[code %in% snakemake@params[['strand_policy']][['codes_to_retain']]]

if(snakemake@params[['strand_policy']][['status']] == 'nonswitched') {
  dat[code %in% snakemake@params[['strand_policy']][['flip']], `:=` (REF = REF.1kG, ALT = ALT.1kG, BETA = -1 * BETA)]

  if('OR' %in% names(dat)) {
    dat[code %in% snakemake@params[['strand_policy']][['flip']], OR := 1/OR]
  }

  if('Z' %in% names(dat)) {
    dat[code %in% snakemake@params[['strand_policy']][['flip']], Z := -1 * Z]
  }
} else {
  stop('Behaviour not yet implemented for this strand policy')
}

dat[, c('code', 'REF.1kG', 'ALT.1kG', 'A.CODE', 'B.CODE') := NULL]

fwrite(dat, file = snakemake@output[[1]], sep = '\t')
