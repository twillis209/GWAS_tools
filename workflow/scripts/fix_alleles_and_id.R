library(data.table)
setDTthreads(snakemake@threads)

library(stringr)

dat <- fread(snakemake@input[[1]], header = T, sep = '\t')

dat[, REF := toupper(REF)]
dat[, ALT := toupper(ALT)]

dat[, CHR := str_remove(CHR, 'chr|CHR')]

dat[CHR == '23', CHR := 'X']

# Sometimes picked up as character column if in scientific notation
dat[, P := as.numeric(P)]

if(!('SNPID' %in% names(dat))) {
  dat[is.na(REF) | is.na(ALT), SNPID := paste(CHR, BP, sep = '_')]
  dat[!is.na(REF) & !is.na(ALT), SNPID := paste(CHR, BP, REF, ALT, sep = '_')]
} else {
  dat[(SNPID == "" | is.na(SNPID)) & (is.na(REF) | is.na(ALT)), SNPID := paste(CHR, BP, sep = '_')]
  dat[(SNPID == "" | is.na(SNPID)) & !is.na(REF) & !is.na(ALT), SNPID := paste(CHR, BP, REF, ALT, sep = '_')]
}

fwrite(dat, file = snakemake@output[[1]], sep = '\t')
