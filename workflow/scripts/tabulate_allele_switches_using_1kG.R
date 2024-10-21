library(data.table)
setDTthreads(snakemake@threads)

# NB: Some code written by Dr Chris Wallace for the annotSnpStats package. Pasting this in here avoids having to create and upload the conda version of the package to anaconda just for the sake of these functions.
##' Reverse alleles in a genotype
##'
##' ie A/G -> G/A
##'
##' @param x character vector of genotypes
##' @param sep character with which to separate alleles. Default is "/".
##' @importFrom data.table tstrsplit
##' @export
##' @return character vector of reversed genotypes 
##' @examples
##' g.rev(c("A/G","A/T"))
##' @author Chris Wallace
g.rev <- function(x,sep="/") {
  tmp=tstrsplit(x,sep)
  paste(tmp[[2]],tmp[[1]],sep=sep) 
}

##' Complement genotypes
##'
##' ie A/G -> T/C
##' 
##' @param x character vector of genotypes
##' @export
##' @return character vector of genotypes on the alternative strand
##' @examples
##' g.complement(c("A/G","A/T"))
##' @author Chris Wallace
g.complement <- function(x) {
  chartr("ATCG","TAGC",toupper(x))
}

##' define possible allele switching classes
##'
##' @title g.class
##' @param x vector of allele codes from dataset X
##' @param y vector of allele codes from dataset Y, same length as x
##' @return character vector of allele switching classes
##' @export
##' @examples
##' alleles.X <- c(snp1="A/G",snp2="A/G",snp3="A/G",snp4="A/G",snp5="A/T",snp6="A/T")
##' alleles.Y <- c(snp1="A/G",snp2="G/A",snp3="T/C",snp4="C/T",snp5="A/T",snp6="T/A")
##' classes <- g.class(x=alleles.X,y=alleles.Y)
##' cbind(alleles.X,alleles.Y,classes)
##' @author Chris Wallace
g.class <- function(x,y) {
  if(!identical(names(x),names(y)))
    stop("x and y must relate to same SNPs")
  mat <- matrix(FALSE,length(x),4,dimnames=list(names(x),c("nochange","rev","comp","revcomp")))
  ## nochange
  mat[ , "nochange" ] <- x==y
  mat[, "rev"] <- x==g.rev(y)
  mat[,"comp"] <- x==g.complement(y)
  mat[,"revcomp"] <- x==g.rev(g.complement(y))
  indels <- x %in% c("I/D","D/I")
  if(any(indels))
    mat[indels,c("comp","revcomp")] <- FALSE
  ret <- character(nrow(mat))
  rs <- rowSums(mat)
  if(length(wh <- which(rs>1))) # ambiguity first
    ret[wh] <- "ambig"  
  if(length(wh <- which(rs==0))) # impossible
    ret[wh] <- "impossible"
  if(length(wh <- which(rs==1))) # impossible
    ret[wh] <- colnames(mat)[ apply(mat[wh,,drop=FALSE],1,which) ]
  return(ret)
}

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

