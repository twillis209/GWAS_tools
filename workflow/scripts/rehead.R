library(data.table)
setDTthreads(snakemake@threads)

library(magrittr)
library(stringr)

dat <- fread(snakemake@input[[1]])

col_names <- names(dat)

cols_to_drop <- col_names[col_names %in% snakemake@params$columns_to_drop]

if(!is.null(snakemake@params$study_specific_params$extraneous_columns)) {
  cols_to_drop <- c(cols_to_drop, col_names[col_names %in% snakemake@params$study_specific_params$extraneous_columns])
}

if(length(cols_to_drop) > 0) {
  dat[, (cols_to_drop) := NULL]
}

# chromosome column
str_replace(col_names, "^Chr$|^chromosome$|^Chromosome$|^chr$|^Chr_ID$|^hg18chr$|^CHROMOSOME$|^#chrom$|^#CHROM$|^chrom$|^#CHR$", "CHR") %>%
 	str_replace("^Pos$|^base_pair_location$|^BP$|^BP\\(hg19\\)$|^Position$|^POS$|^pos$|^Chr_Position$|^bp$|^position$|^Position\\(hg19\\)$|^POSITION$|^bp_hg19$|^BP_hg19$|^Coordinate$|^chrloc$", "BP") %>%
 	str_replace("^íd$|^id$|^ID$|^variant_id$|^MarkerName$|^SNP$|^rsid$|^rsids$|^SNP_Name$|^snp$|^snpid$|^SNP_ID$|^rsID$|^#SNPID$|^rs_number$|^RSID$|^rs$|^db_SNP_RS_IDMarker$|^dbSNP_RS_ID$|^Variant$","SNPID") %>%
 	str_replace("^Beta$|^beta$|^Effect$|^effect$|^EFFECT$|^beta_SNP_add$|^EFFECT_ALT$|^effB$|^beta_EUR$|^all_inv_var_meta_beta$|^frequentist_add_beta_1$","BETA") %>%
 	str_replace("^standard_error$|^StdErr$|^stderr$|^sebeta_SNP_add$|^se$|^STDERR$|^sebeta$|^se_effB$|^se_EUR$|^all_inv_var_meta_sebeta$|^LOG\\(OR\\)_SE$|^frequentist_add_se_1$","SE") %>%
 	str_replace("^odds_ratio$|^Odds_ratio$|^or$|^OddsRatio$|^OR\\(A1\\)$|^ORX$","OR") %>%
 	str_replace("^p_value$|^P.value$|^pvalue$|^P-value$|^pval$|^p.value$|^Pval$|^PVALUE$|^PVAL$|^Pvalue$|^P_VALUE$|^P-val$|^p$|^All.p.value$|^P_value$|^p-value$|^GC-adjusted_P_$|^Chi-Squared__P$|^P1df$|^all_inv_var_meta_p$","P") %>%
 	str_replace("Log10p","LOG10P") %>%
 	str_replace("-log10_p-value","-LOG10P") %>%
 	str_replace("^effect_allele_frequency$|^<maf$|^<MAF$","ALT_FREQ") %>%
 	str_replace("^EMP_Beta$","EMP_BETA") %>%
 	str_replace("^EMP1$","EMP_P") %>%
 	str_replace("^EMP_se$","EMP_SE") %>%
 	str_replace("^hm_effect_allele$","hm_ALT") %>%
 	str_replace("^hm_beta$","hm_BETA") %>%
 	str_replace("^hm_pos$","hm_BP") %>%
 	str_replace("^hm_chrom$","hm_CHR") %>%
 	str_replace("^hm_odds_ratio$","hm_OR") %>%
 	str_replace("^hm_other_allele$","hm_REF") %>%
 	str_replace("^hm_rsid$","hm_SNPID") %>%
 	str_replace("^n$","N") %>%
 	str_replace("^Rsq$","RSQ") %>%
 	str_replace("^MARKER$|^íd$|^Chr:Position$","CHR:BP") %>%
  str_replace("^Zscore$|^ZSCORE$|^Z_STAT$","Z") %>%
  str_replace(sprintf("^%s$", snakemake@params$pan_ukb_beta_column), "BETA") %>%
  str_replace(sprintf("^%s$", snakemake@params$pan_ukb_se_column), "SE") -> updated_col_names

# Caution! Sometimes "other_allele" means effect allele", e.g. the LADA data set, check papers prior to running the script, and pre-rename accordingly.
if(!is.null(snakemake@params$study_specific_params$ref)) {
  updated_col_names <- str_replace(updated_col_names, snakemake@params$study_specific_params$ref, 'REF')
} else {
  updated_col_names <- str_replace(updated_col_names, "^OtherAllele$|^reference_allele$|^Ref_Allele$|^OTHER_ALLELE$|^other_allele$|^A2_other$|^NEA$|^Ref_Allele$|^Ref$|^ref$|^Allele1$|^A2$", "REF")
}

if(!is.null(snakemake@params$study_specific_params$alt)) {
  updated_col_names <- str_replace(updated_col_names, snakemake@params$study_specific_params$alt, 'ALT')
} else {
  updated_col_names <- str_replace(updated_col_names, "^effect_allele$|^Effect_Allele$|^EffectAllele$|^A1_effect$|^RISK_ALLELE$|^EA$|^Risk_Allele$|^EFFECT_ALLELE$|^Alt$|^alt$|^Allele2$|^A1$", "ALT")
}

names(dat) <- updated_col_names

if(!('P' %in% updated_col_names) & snakemake@params$pan_ukb_p_column %in% names(dat)) {
  setnames(dat, snakemake@params$pan_ukb_p_column, "P")
  dat[, P := exp(P)]
} else if(!('P' %in% updated_col_names)) {
  stop("Could not find p-value column")
}

if(dat[, class(P)] != 'numeric') {
  dat[, P := as.numeric(P)]
}

if(dat[P < 0 | P > 1, .N] > 0) {
  stop("One or more p-values falls outside [0, 1]")
}

if('CHR' %in% names(dat) & 'CHR38' %in% names(dat)) {
  dat[, CHR := NULL]
  setnames(dat, 'CHR38', 'CHR')
}

if('CHR38' %in% names(dat) & 'CHR19' %in% names(dat)) {
  # Keep newer assembly's coordinates
  dat[, CHR19 := NULL]
  setnames(dat, 'CHR38', 'CHR')
} else if('CHR38' %in% names(dat)) {
  setnames(dat, 'CHR38', 'CHR')
} else if('CHR19' %in% names(dat)) {
  setnames(dat, 'CHR19', 'CHR')
}


if('BP' %in% names(dat) & 'BP38' %in% names(dat)) {
  dat[, BP := NULL]
  setnames(dat, 'BP38', 'BP')
}


if('BP38' %in% names(dat) & 'BP19' %in% names(dat)) {
  # Keep newer assembly's coordinates
  dat[, BP19 := NULL]
  setnames(dat, 'BP38', 'BP')
} else if('BP38' %in% names(dat)) {
  setnames(dat, 'BP38', 'BP')
} else if('BP19' %in% names(dat)) {
  setnames(dat, 'BP19', 'BP')
}

fwrite(dat, file = snakemake@output[[1]], sep = '\t')
