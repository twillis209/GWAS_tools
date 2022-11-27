library(data.table)
setDTthreads(snakemake@threads)

dat <- fread(snakemake@input[[1]])

fwrite(dat, file = snakemake@output[[1]], sep = snakemake@params[['output_separator']])
