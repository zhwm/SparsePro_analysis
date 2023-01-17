setwd("~/Desktop/SparsePro/Overleaf/Tables")

result <- read.table("../../Zenodo/doc/gw_precision_recall_size_time.csv",header = T)
sumdat <- data.frame(Method = rep(unique(result$Method),each = 3),
                     W = 0:2,
                     meanPrecision = NA,
                     medianPrecision = NA,
                     meanRecall = NA,
                     medianRecall = NA,
                     meanSize = NA,
                     medianSize = NA)
for (j in 1:nrow(sumdat)) {
  sumdat$meanPrecision[j] <- mean(result$Precision[result$Method==sumdat$Method[j] &
                                                     result$W==sumdat$W[j]])
  sumdat$medianPrecision[j] <- median(result$Precision[result$Method==sumdat$Method[j] &
                                                         result$W==sumdat$W[j]])
  sumdat$meanRecall[j] <- mean(result$Recall[result$Method==sumdat$Method[j] &
                                               result$W==sumdat$W[j]])
  sumdat$medianRecall[j] <- median(result$Recall[result$Method==sumdat$Method[j] &
                                                   result$W==sumdat$W[j]])
  sumdat$meanSize[j] <- mean(result$Size[result$Method==sumdat$Method[j] &
                                           result$W==sumdat$W[j]])
  sumdat$medianSize[j] <- median(result$Size[result$Method==sumdat$Method[j] &
                                               result$W==sumdat$W[j]])
}
write.table(sumdat,"SuppTab-gw-simulation-statistics.txt",sep = "\t",quote = F,col.names = T,row.names = F)


result <- read.table("../../Zenodo/doc/loc_cs.csv",header = T)
sumdat <- data.frame(Method = rep(unique(result$Method),each = 36),
                     K = rep(c(1,2,5,10),each = 9),
                     W = rep(0:2,each = 3),
                     Chr = 20:22,
                     meanPrecision = NA,
                     medianPrecision = NA,
                     meanRecall = NA,
                     medianRecall = NA,
                     meanSize = NA,
                     medianSize = NA)
for (j in 1:nrow(sumdat)) {
  sumdat$meanPrecision[j] <- mean(result$Precision[result$Method==sumdat$Method[j] &
                                                     result$K==sumdat$K[j] &
                                                     result$W==sumdat$W[j] &
                                                     result$CHR==sumdat$Chr[j]])
  sumdat$medianPrecision[j] <- median(result$Precision[result$Method==sumdat$Method[j] &
                                                     result$K==sumdat$K[j] &
                                                     result$W==sumdat$W[j] &
                                                     result$CHR==sumdat$Chr[j]])
  sumdat$meanRecall[j] <- mean(result$Recall[result$Method==sumdat$Method[j] &
                                                     result$K==sumdat$K[j] &
                                                     result$W==sumdat$W[j] &
                                                     result$CHR==sumdat$Chr[j]])
  sumdat$medianRecall[j] <- median(result$Recall[result$Method==sumdat$Method[j] &
                                                         result$K==sumdat$K[j] &
                                                         result$W==sumdat$W[j] &
                                                         result$CHR==sumdat$Chr[j]])
  sumdat$meanSize[j] <- mean(result$Size[result$Method==sumdat$Method[j] &
                                                     result$K==sumdat$K[j] &
                                                     result$W==sumdat$W[j] &
                                                     result$CHR==sumdat$Chr[j]])
  sumdat$medianSize[j] <- median(result$Size[result$Method==sumdat$Method[j] &
                                                         result$K==sumdat$K[j] &
                                                         result$W==sumdat$W[j] &
                                                         result$CHR==sumdat$Chr[j]])
}
sumdat$Chr <- c("Locus 1","Locus 2","Locus 3")
colnames(sumdat)[4] <- "Locus"
write.table(sumdat,"SuppTab-locus-simulation-statistics.txt",sep = "\t",quote = F,col.names = T,row.names = F)

result <- read.table("../../Zenodo/doc/gw_sparsepro_wsep.csv",sep = "\t",header = T)
result$Chromosome <- rep(1:22, each = 10)
colnames(result) <- c("Annotation","Estimated W","se","p","Real W","Chromosome")
write.table(result,"SuppTab-gw-sparsepro-enrichment.txt",sep = "\t",quote = F,col.names = T,row.names = F)

result <- read.table("../../Zenodo/doc/gw_polyfun_annot.csv",sep = "\t",header = F)
result$Chromosome <- rep(1:22, each = 20)
colnames(result) <- c("Annotation","Coefficient","Real W","Chromosome")
write.table(result,"SuppTab-gw-polyfun-enrichment.txt",sep = "\t",quote = F,col.names = T,row.names = F)

result <- rbind.data.frame(read.table("../../Zenodo/doc/loc_sparsepro_anno_20.csv",sep = "\t",header = F),
                           read.table("../../Zenodo/doc/loc_sparsepro_anno_21.csv",sep = "\t",header = F),
                           read.table("../../Zenodo/doc/loc_sparsepro_anno_22.csv",sep = "\t",header = F))
result$Locus <- rep(c("Locus 1","Locus 2","Locus 3"),each = 120)
colnames(result)[1:6] <- c("Annotation","Estimated W","se","p","K","true W")
write.table(result,"SuppTab-loc-sparsepro-enrichment.txt",sep = "\t",quote = F,col.names = T,row.names = F)

result <- rbind.data.frame(read.table("../../Zenodo/doc/loc_paintor_anno_20.csv",sep = "\t",header = F),
                           read.table("../../Zenodo/doc/loc_paintor_anno_21.csv",sep = "\t",header = F),
                           read.table("../../Zenodo/doc/loc_paintor_anno_22.csv",sep = "\t",header = F))
result <- result[,-1]
colnames(result) <- c("baseline","Conserved_LindbladToh","DHS_Trynka","H3K27ac_Hnisz",
                      "H3K4me3_Trynka","Transcr_Hoffman","TSS_Hoffman","UTR_3_UCSC","UTR_5_UCSC",
                      "non_synonymous","Human_Promoter_Villar_ExAC","K","Real W")
write.table(result,"SuppTab-loc-paintor-enrichment.txt",sep = "\t",quote = F,col.names = T,row.names = F)

result <- read.table("../../Zenodo/doc/eGFR.cgene.csv",sep = "\t")
colnames(result) <- c("Region","Chr","Start","End","Nearest_gene","Distance_to_nearest_gene",
                      "Causal_set","PIP","Causal_effect_estimate")
result$Method <- "SparsePro-"
newres <- read.table("../../Zenodo/doc/eGFR.anno.cgene.csv",sep = "\t")
colnames(newres) <- c("Region","Chr","Start","End","Nearest_gene","Distance_to_nearest_gene",
                      "Causal_set","PIP","Causal_effect_estimate")
newres$Method <- "SparsePro+"
result <- rbind.data.frame(result,newres)
write.table(result,"SuppTab-realdata-eGFR-sparsepro.txt",sep = "\t",quote = F,col.names = T,row.names = F)

result <- read.table("../../Zenodo/doc/eGFR.susie.cgene.csv",sep = "\t")
colnames(result) <- c("Region","Chr","Start","End","Nearest_gene","Distance_to_nearest_gene",
                      "Causal_set","PIP")
result$Method <- "SuSiE"
newres <- read.table("../../Zenodo/doc/eGFR.susiepoly.cgene.csv",sep = "\t")
colnames(newres) <- c("Region","Chr","Start","End","Nearest_gene","Distance_to_nearest_gene",
                      "Causal_set","PIP")
newres$Method <- "SuSiE+PolyFun"
result <- rbind.data.frame(result,newres)
write.table(result,"SuppTab-realdata-eGFR-susie.txt",sep = "\t",quote = F,col.names = T,row.names = F)

result <- read.table("../../Zenodo/doc/FFR.cgene.csv",sep = "\t")
colnames(result) <- c("Region","Chr","Start","End","Nearest_gene","Distance_to_nearest_gene",
                      "Causal_set","PIP","Causal_effect_estimate")
result$Method <- "SparsePro-"
newres <- read.table("../../Zenodo/doc/FFR.anno.cgene.csv",sep = "\t")
colnames(newres) <- c("Region","Chr","Start","End","Nearest_gene","Distance_to_nearest_gene",
                      "Causal_set","PIP","Causal_effect_estimate")
newres$Method <- "SparsePro+"
result <- rbind.data.frame(result,newres)
write.table(result,"SuppTab-realdata-FFR-sparsepro.txt",sep = "\t",quote = F,col.names = T,row.names = F)

result <- read.table("../../Zenodo/doc/FFR.susie.cgene.csv",sep = "\t")
colnames(result) <- c("Region","Chr","Start","End","Nearest_gene","Distance_to_nearest_gene",
                      "Causal_set","PIP")
result$Method <- "SuSiE"
newres <- read.table("../../Zenodo/doc/FFR.susiepoly.cgene.csv",sep = "\t")
colnames(newres) <- c("Region","Chr","Start","End","Nearest_gene","Distance_to_nearest_gene",
                      "Causal_set","PIP")
newres$Method <- "SuSiE+PolyFun"
result <- rbind.data.frame(result,newres)
write.table(result,"SuppTab-realdata-FFR-susie.txt",sep = "\t",quote = F,col.names = T,row.names = F)

result <- read.table("../../Zenodo/doc/glucose.cgene.csv",sep = "\t")
colnames(result) <- c("Region","Chr","Start","End","Nearest_gene","Distance_to_nearest_gene",
                      "Causal_set","PIP","Causal_effect_estimate")
result$Method <- "SparsePro-"
newres <- read.table("../../Zenodo/doc/glucose.anno.cgene.csv",sep = "\t")
colnames(newres) <- c("Region","Chr","Start","End","Nearest_gene","Distance_to_nearest_gene",
                      "Causal_set","PIP","Causal_effect_estimate")
newres$Method <- "SparsePro+"
result <- rbind.data.frame(result,newres)
write.table(result,"SuppTab-realdata-glucose-sparsepro.txt",sep = "\t",quote = F,col.names = T,row.names = F)

result <- read.table("../../Zenodo/doc/glucose.susie.cgene.csv",sep = "\t")
colnames(result) <- c("Region","Chr","Start","End","Nearest_gene","Distance_to_nearest_gene",
                      "Causal_set","PIP")
result$Method <- "SuSiE"
newres <- read.table("../../Zenodo/doc/glucose.susiepoly.cgene.csv",sep = "\t")
colnames(newres) <- c("Region","Chr","Start","End","Nearest_gene","Distance_to_nearest_gene",
                      "Causal_set","PIP")
newres$Method <- "SuSiE+PolyFun"
result <- rbind.data.frame(result,newres)
write.table(result,"SuppTab-realdata-glucose-susie.txt",sep = "\t",quote = F,col.names = T,row.names = F)

result <- read.table("../../Zenodo/doc/gammaGT.cgene.csv",sep = "\t")
colnames(result) <- c("Region","Chr","Start","End","Nearest_gene","Distance_to_nearest_gene",
                      "Causal_set","PIP","Causal_effect_estimate")
result$Method <- "SparsePro-"
newres <- read.table("../../Zenodo/doc/gammaGT.anno.cgene.csv",sep = "\t")
colnames(newres) <- c("Region","Chr","Start","End","Nearest_gene","Distance_to_nearest_gene",
                      "Causal_set","PIP","Causal_effect_estimate")
newres$Method <- "SparsePro+"
result <- rbind.data.frame(result,newres)
write.table(result,"SuppTab-realdata-gammaGT-sparsepro.txt",sep = "\t",quote = F,col.names = T,row.names = F)

result <- read.table("../../Zenodo/doc/gammaGT.susie.cgene.csv",sep = "\t")
colnames(result) <- c("Region","Chr","Start","End","Nearest_gene","Distance_to_nearest_gene",
                      "Causal_set","PIP")
result$Method <- "SuSiE"
newres <- read.table("../../Zenodo/doc/gammaGT.susiepoly.cgene.csv",sep = "\t")
colnames(newres) <- c("Region","Chr","Start","End","Nearest_gene","Distance_to_nearest_gene",
                      "Causal_set","PIP")
newres$Method <- "SuSiE+PolyFun"
result <- rbind.data.frame(result,newres)
write.table(result,"SuppTab-realdata-gammaGT-susie.txt",sep = "\t",quote = F,col.names = T,row.names = F)

result <- read.table("../../Zenodo/doc/prate.cgene.csv",sep = "\t")
colnames(result) <- c("Region","Chr","Start","End","Nearest_gene","Distance_to_nearest_gene",
                      "Causal_set","PIP","Causal_effect_estimate")
result$Method <- "SparsePro-"
newres <- read.table("../../Zenodo/doc/prate.anno.cgene.csv",sep = "\t")
colnames(newres) <- c("Region","Chr","Start","End","Nearest_gene","Distance_to_nearest_gene",
                      "Causal_set","PIP","Causal_effect_estimate")
newres$Method <- "SparsePro+"
result <- rbind.data.frame(result,newres)
write.table(result,"SuppTab-realdata-prate-sparsepro.txt",sep = "\t",quote = F,col.names = T,row.names = F)

result <- read.table("../../Zenodo/doc/prate.susie.cgene.csv",sep = "\t")
colnames(result) <- c("Region","Chr","Start","End","Nearest_gene","Distance_to_nearest_gene",
                      "Causal_set","PIP")
result$Method <- "SuSiE"
newres <- read.table("../../Zenodo/doc/prate.susiepoly.cgene.csv",sep = "\t")
colnames(newres) <- c("Region","Chr","Start","End","Nearest_gene","Distance_to_nearest_gene",
                      "Causal_set","PIP")
newres$Method <- "SuSiE+PolyFun"
result <- rbind.data.frame(result,newres)
write.table(result,"SuppTab-realdata-prate-susie.txt",sep = "\t",quote = F,col.names = T,row.names = F)

