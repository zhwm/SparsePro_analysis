setwd("~/scratch/SparsePro/UKB-GWAS")
library(data.table)
INTfun <- function(x) {qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}

pheno <- fread("pheno.retain.txt")
phenoEUR <- pheno[pheno$ethnicity==1001,]
removeID <- read.table("ukb45551_rel_s488264.remove.related.idx")$V1
phenoEUR <- phenoEUR[!phenoEUR$ID %in% removeID,]
phenoEUR <- as.data.frame(phenoEUR)
phenoEUR$centre <- as.factor(phenoEUR$centre)
phenoEUR$array <- ifelse(phenoEUR$array > 0, 1, 0) #two types of array
phenoEUR <- phenoEUR[!is.na(phenoEUR$PC1),] #no genotype
phenoEUR$age2 <- phenoEUR$age^2

phenoEUR$pulse_rate_resid <- NA
phenoEUR$pulse_rate_resid[!is.na(phenoEUR$pulse_rate)] <- summary(lm(pulse_rate ~ age + age2 + sex + array + centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, data = phenoEUR))$residuals
phenoEUR$pulse_rate_resid <- INTfun(phenoEUR$pulse_rate_resid)
write_output <- phenoEUR[,c("ID","ID","pulse_rate_resid")]
write.table(write_output,"GWAS/pulse_rate_INToutcome.txt",sep=" ",col.names=F,row.names=F,quote=F)
phenoEUR$ALT_resid <- NA
phenoEUR$ALT_resid[!is.na(phenoEUR$ALT)] <- summary(lm(ALT ~ age + age2 + sex + array + centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, data = phenoEUR))$residuals
phenoEUR$ALT_resid <- INTfun(phenoEUR$ALT_resid)
write_output <- phenoEUR[,c("ID","ID","ALT_resid")]
write.table(write_output,"GWAS/ALT_INToutcome.txt",sep=" ",col.names=F,row.names=F,quote=F)
phenoEUR$AST_resid <- NA
phenoEUR$AST_resid[!is.na(phenoEUR$AST)] <- summary(lm(AST ~ age + age2 + sex + array + centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, data = phenoEUR))$residuals
phenoEUR$AST_resid <- INTfun(phenoEUR$AST_resid)
write_output <- phenoEUR[,c("ID","ID","AST_resid")]
write.table(write_output,"GWAS/AST_INToutcome.txt",sep=" ",col.names=F,row.names=F,quote=F)
phenoEUR$gammaGT_resid <- NA
phenoEUR$gammaGT_resid[!is.na(phenoEUR$gammaGT)] <- summary(lm(gammaGT ~ age + age2 + sex + array + centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, data = phenoEUR))$residuals
phenoEUR$gammaGT_resid <- INTfun(phenoEUR$gammaGT_resid)
write_output <- phenoEUR[,c("ID","ID","gammaGT_resid")]
write.table(write_output,"GWAS/gammaGT_INToutcome.txt",sep=" ",col.names=F,row.names=F,quote=F)
phenoEUR$glucose_resid <- NA
phenoEUR$glucose_resid[!is.na(phenoEUR$glucose)] <- summary(lm(glucose ~ age + age2 + sex + array + centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, data = phenoEUR))$residuals
phenoEUR$glucose_resid <- INTfun(phenoEUR$glucose_resid)
write_output <- phenoEUR[,c("ID","ID","glucose_resid")]
write.table(write_output,"GWAS/glucose_INToutcome.txt",sep=" ",col.names=F,row.names=F,quote=F)
phenoEUR$tbilirubin_resid <- NA
phenoEUR$tbilirubin_resid[!is.na(phenoEUR$tbilirubin)] <- summary(lm(tbilirubin ~ age + age2 + sex + array + centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, data = phenoEUR))$residuals
phenoEUR$tbilirubin_resid <- INTfun(phenoEUR$tbilirubin_resid)
write_output <- phenoEUR[,c("ID","ID","tbilirubin_resid")]
write.table(write_output,"GWAS/tbilirubin_INToutcome.txt",sep=" ",col.names=F,row.names=F,quote=F)
phenoEUR$tprotein_resid <- NA
phenoEUR$tprotein_resid[!is.na(phenoEUR$tprotein)] <- summary(lm(tprotein ~ age + age2 + sex + array + centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, data = phenoEUR))$residuals
phenoEUR$tprotein_resid <- INTfun(phenoEUR$tprotein_resid)
write_output <- phenoEUR[,c("ID","ID","tprotein_resid")]
write.table(write_output,"GWAS/tprotein_INToutcome.txt",sep=" ",col.names=F,row.names=F,quote=F)
phenoEUR$FFR_resid <- NA
phenoEUR$FFR_resid[!is.na(phenoEUR$FFR)] <- summary(lm(FFR ~ age + age2 + sex + array + centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, data = phenoEUR))$residuals
phenoEUR$FFR_resid <- INTfun(phenoEUR$FFR_resid)
write_output <- phenoEUR[,c("ID","ID","FFR_resid")]
write.table(write_output,"GWAS/FFR_INToutcome.txt",sep=" ",col.names=F,row.names=F,quote=F)

# https://www.niddk.nih.gov/health-information/professionals/clinical-tools-patient-management/kidney-disease/laboratory-evaluation/glomerular-filtration-rate/estimating
phenoEUR$creatinine <- phenoEUR$creatinine * 0.0113
phenoEUR <- phenoEUR[,c(1,2,6,7,8:28,31,36)]
phenoEUR <- na.omit(phenoEUR)
phenoEUR$eGFR <- NA
phenoEUR$eGFR[phenoEUR$creatinine <= 0.7 & phenoEUR$sex == 0] <- 144 * ((phenoEUR$creatinine[phenoEUR$creatinine <= 0.7 & phenoEUR$sex == 0] / 0.7) ^ -0.329) * (0.993 ^ phenoEUR$age[phenoEUR$creatinine <= 0.7 & phenoEUR$sex == 0])
phenoEUR$eGFR[phenoEUR$creatinine > 0.7 & phenoEUR$sex == 0] <- 144 * ((phenoEUR$creatinine[phenoEUR$creatinine > 0.7 & phenoEUR$sex == 0] / 0.7) ^ -1.209) * (0.993 ^ phenoEUR$age[phenoEUR$creatinine > 0.7 & phenoEUR$sex == 0])
phenoEUR$eGFR[phenoEUR$creatinine <= 0.7 & phenoEUR$sex == 1] <- 141 * ((phenoEUR$creatinine[phenoEUR$creatinine <= 0.7 & phenoEUR$sex == 1] / 0.9) ^ -0.411) * (0.993 ^ phenoEUR$age[phenoEUR$creatinine <= 0.7 & phenoEUR$sex == 1])
phenoEUR$eGFR[phenoEUR$creatinine > 0.7 & phenoEUR$sex == 1] <- 141 * ((phenoEUR$creatinine[phenoEUR$creatinine > 0.7 & phenoEUR$sex == 1] / 0.9) ^ -1.209) * (0.993 ^ phenoEUR$age[phenoEUR$creatinine > 0.7 & phenoEUR$sex == 1])

phenoEUR$eGFR_resid <- NA
phenoEUR$eGFR_resid[!is.na(phenoEUR$eGFR)] <- summary(lm(eGFR ~ age + age2 + sex + array + centre + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, data = phenoEUR))$residuals
phenoEUR$eGFR_resid <- INTfun(phenoEUR$eGFR_resid)
write_output <- phenoEUR[,c("ID","ID","eGFR_resid")]
write.table(write_output,"GWAS/eGFR_INToutcome.txt",sep=" ",col.names=F,row.names=F,quote=F)
