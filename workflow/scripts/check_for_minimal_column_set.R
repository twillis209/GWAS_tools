library(data.table)
setDTthreads(snakemake@threads)

dat <- fread(snakemake@input[[1]], header = T)

if(!('CHR' %in% names(dat))) {
  stop("No chromosome column")
} else if(!('BP' %in% names(dat))) {
  stop("No basepair column")
} else if(!('REF' %in% names(dat))) {
  stop("No REF column")
} else if(!('ALT' %in% names(dat))) {
  stop("No ALT column")
} else if(!('P' %in% names(dat))) {
  stop("No P column")
} else if(!('OR' %in% names(dat)) & !('BETA' %in% names(dat))) {
  stop("Neither OR nor BETA column")
} else {
  fwrite(dat, file = snakemake@output[[1]], sep = '\t')
}


