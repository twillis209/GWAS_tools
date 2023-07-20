##' Reverse alleles in a genotype
##'
##' ie A/G -> G/A
##'
##' @param x character vector of genotypes
##' @param sep character with which to separate alleles. Default is "/".
##' @export
##' @return character vector of reversed genotypes 
##' @examples
##' g.rev(c("A/G","A/T"))
##' @author Chris Wallace
g.rev <- function(x,sep="/") {
  tmp=tstrsplit(x,sep)
  paste(tmp[[2]],tmp[[1]],sep=sep) 
  ## sapply(strsplit(x,sep),function(g) paste(rev(g),collapse="/"))
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
  ## x <- toupper(x)
  ## switches <- c("A"="t","T"="a","C"="g","G"="c")
  ## for(i in seq_along(switches))
  ##   x <- sub(names(switches)[i],switches[i],x)
  ## toupper(x)
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
