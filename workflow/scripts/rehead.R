library(data.table)
setDTthreads(snakemake@threads)

library(magrittr)
library(stringr)

dat <- fread(snakemake@input[[1]], header = T, sep = '\t')

col_names <- names(dat)

# chromosome column
str_replace(col_names, "^Chr$|^chromosome$|^Chromosome$|^chr$|^Chr_ID$|^hg18chr$|^CHROMOSOME$|^#chrom$|^#CHROM$|^chrom$|^#CHR$", "CHR") %>%
 	str_replace("^Pos$|^base_pair_location$|^BP$|^BP(hg19)$|^Position$|^POS$|^pos$|^Chr_Position$|^bp$|^position$|^Position(hg19)$|^POSITION$|^bp_hg19$|^Coordinate$|^chrloc$", "BP") %>%
 	str_replace("^íd$|^id$|^ID$|^variant_id$|^MarkerName$|^SNP$|^rsid$|^rsids$|^SNP_Name$|^snp$|^snpid$|^SNP_ID$|^rsID$|^#SNPID$|^rs_number$|^RSID$|^rs$|^db_SNP_RS_IDMarker$|^dbSNP_RS_ID$|^Variant$","SNPID") %>%
 	str_replace("^OtherAllele$|^reference_allele$|^Ref_Allele$|^OTHER_ALLELE$|^other_allele$|^A2_other$|^NEA$|^Ref_Allele$|^Ref$|^ref$|^Allele1$","REF") %>%
 	str_replace("^effect_allele$|^Effect_Allele$|^EffectAllele$|^A1_effect$|^RISK_ALLELE$|^EA$|^Risk_Allele$|^EFFECT_ALLELE$|^Alt$|^alt$|^Allele2$","ALT") %>%
 	str_replace("^Beta$|^beta$|^Effect$|^effect$|^EFFECT$|^beta_SNP_add$|^EFFECT_ALT$|^effB$|^beta_EUR$|^all_inv_var_meta_beta$","BETA") %>%
 	str_replace("^standard_error$|^StdErr$|^stderr$|^sebeta_SNP_add$|^se$|^STDERR$|^sebeta$|^se_effB$|^se_EUR$|^all_inv_var_meta_sebeta$|^LOG(OR)_SE$","SE") %>%
 	str_replace("^odds_ratio$|^Odds_ratio$|^or$|^OddsRatio$|^OR(A1)$|^ORX$","OR") %>%
 	str_replace("^p_value$|^P.value$|^pvalue$|^P-value$|^pval$|^p.value$|^Pval$|^PVALUE$|^Pvalue$|^P_VALUE$|^P-val$|^p$|^All.p.value$|^P_value$|^p-value$|^GC-adjusted_P_$|^Chi-Squared__P$|^P1df$|^pval_EUR$|^all_inv_var_meta_p$","P") %>%
 	str_replace("Log10p","LOG10P") %>%
 	str_replace("-log10_p-value","-LOG10P") %>%
 	str_replace("^effect_allele_frequency$|^<maf$|^<MAF$","ALT_FREQ") %>%
 # Caution! Sometimes "other_allele" means effect allele", check papers prior to run the script, and pre-rename accordingly.
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
  str_replace("^Zscore$|^ZSCORE$|^Z_STAT$","Z") -> updated_col_names

names(dat) <- updated_col_names

fwrite(dat, file = snakemake@output[[1]], sep = '\t')
