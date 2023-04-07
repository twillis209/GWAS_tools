library(data.table)
setDTthreads(snakemake@threads)

dat <- fread(snakemake@input[[1]], header = T, sep = '\t')

if('OR' %in% names(dat) & !('BETA' %in% names(dat))) {
  dat[, BETA := log(OR)]
}

if(!('SE' %in% names(dat))) {
  if(('BETA' %in% names(dat) & 'P' %in% names(dat))) {
    dat[, Z := sign(BETA)* abs(qnorm(P/2))]
    dat[, SE := BETA/Z]
  } else if(('BETA' %in% names(dat) & 'Z' %in% names(dat))) {
    dat[, SE := BETA/Z]
  } else {
    stop("No SE column and need BETA and P or BETA and Z to calculate it")
  }
}

fwrite(dat, file = snakemake@output[[1]], sep = '\t')
