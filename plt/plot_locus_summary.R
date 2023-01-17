library(ggplot2)
auc = read.table("loc_auprc.csv",header = T)
auc22 = auc[auc$CHR==22,]

plotdat = data.frame(AUPRC = c(auc22$FINEMAP,auc22$SuSiE,auc22$SparsePro.,auc22$SparsePro..1,auc22$PAINTOR.,auc22$PAINTOR..1),
                      K = paste0("K = ",auc22$K),
                      W = paste0("W = ",auc22$W),
                      Method = rep(c("FINEMAP","SuSiE","SparsePro-","SparsePro+","PAINTOR-","PAINTOR+"),each = nrow(auc22)))

plotdat$Method = factor(plotdat$Method, levels = c("FINEMAP","SuSiE","PAINTOR-","PAINTOR+","SparsePro-","SparsePro+"))

plotdat$K = factor(plotdat$K, levels = paste0("K = ",c(1,2,5,10)))

plotdat = plotdat[!(plotdat$Method=="SparsePro+" & plotdat$W=="W = 0"),]
ggplot(plotdat, aes(x = K, y = AUPRC, fill = Method)) +
  facet_grid(~W)+
  geom_bar(stat = "identity",position = "dodge") +
  theme_classic() +
  theme_bw() +
  scale_fill_manual(values = c("green","royalblue","orchid","darkmagenta","orange","red")) +
  theme(legend.position = "top",
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        strip.text = element_text(size = 13),
        legend.title = element_blank(),
        legend.text = element_text(size = 13)) +
  xlab("") +
  coord_cartesian(ylim = c(0.2,1)) -> plt
ggsave("loc_auprc.pdf",plt,height = 4,width = 10)

auc21 = auc[auc$CHR==21,]

plotdat = data.frame(AUPRC = c(auc21$FINEMAP,auc21$SuSiE,auc21$SparsePro.,auc21$SparsePro..1,auc21$PAINTOR.,auc21$PAINTOR..1),
                      K = paste0("K = ",auc21$K),
                      W = paste0("W = ",auc21$W),
                      Method = rep(c("FINEMAP","SuSiE","SparsePro-","SparsePro+","PAINTOR-","PAINTOR+"),each = nrow(auc21)))

plotdat$Method = factor(plotdat$Method, levels = c("FINEMAP","SuSiE","PAINTOR-","PAINTOR+","SparsePro-","SparsePro+"))

plotdat$K = factor(plotdat$K, levels = paste0("K = ",c(1,2,5,10)))

plotdat = plotdat[!(plotdat$Method=="SparsePro+" & plotdat$W=="W = 0"),]

ggplot(plotdat, aes(x = K, y = AUPRC, fill = Method)) +
  facet_grid(~W)+
  geom_bar(stat = "identity",position = "dodge") +
  theme_classic() +
  theme_bw() +
  scale_fill_manual(values = c("green","royalblue","orchid","darkmagenta","orange","red")) +
  theme(legend.position = "top",
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        strip.text = element_text(size = 13),
        legend.title = element_blank(),
        legend.text = element_text(size = 13)) +
  xlab("") +
  coord_cartesian(ylim = c(0.2,1)) -> plt
ggsave("loc_auprc_chr21.pdf",plt,height = 4,width = 10)

auc20 = auc[auc$CHR==20,]

plotdat = data.frame(AUPRC = c(auc20$FINEMAP,auc20$SuSiE,auc20$SparsePro.,auc20$SparsePro..1,auc20$PAINTOR.,auc20$PAINTOR..1),
                      K = paste0("K = ",auc20$K),
                      W = paste0("W = ",auc20$W),
                      Method = rep(c("FINEMAP","SuSiE","SparsePro-","SparsePro+","PAINTOR-","PAINTOR+"),each = nrow(auc20)))

plotdat$Method = factor(plotdat$Method, levels = c("FINEMAP","SuSiE","PAINTOR-","PAINTOR+","SparsePro-","SparsePro+"))

plotdat$K = factor(plotdat$K, levels = paste0("K = ",c(1,2,5,10)))

plotdat = plotdat[!(plotdat$Method=="SparsePro+" & plotdat$W=="W = 0"),]

ggplot(plotdat, aes(x = K, y = AUPRC, fill = Method)) +
  facet_grid(~W)+
  geom_bar(stat = "identity",position = "dodge") +
  theme_classic() +
  theme_bw() +
  scale_fill_manual(values = c("green","royalblue","orchid","darkmagenta","orange","red")) +
  theme(legend.position = "top",
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        strip.text = element_text(size = 13),
        legend.title = element_blank(),
        legend.text = element_text(size = 13)) +
  xlab("") +
  coord_cartesian(ylim = c(0.2,1)) -> plt

ggsave("loc_auprc_chr20.pdf",plt,height = 4,width = 10)



cs_res = read.table("loc_cs.csv",header = T)
plotdat = data.frame(Value = cs_res$Precision,
                      Method = cs_res$Method,
                      W = paste0("W = ",cs_res$W),
                      K = paste0("K = ",cs_res$K))
plotdat$K = factor(plotdat$K, levels = paste0("K = ",c(1,2,5,10)))
plotdat$Method = factor(plotdat$Method, levels = c("FINEMAP","SuSiE","PAINTOR-","PAINTOR+","SparsePro-","SparsePro+"))
ggplot(plotdat, aes(x = Method, y = Value, fill = Method)) +
#  geom_violin() +
  geom_boxplot(width = 0.2) +
  facet_grid(K ~ W) +
  theme_bw() +
  scale_fill_manual(values = c("green","royalblue","orchid","darkmagenta","orange","red")) +
  theme(legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text.x = element_text(angle = 45,size = 13,vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 13),
        strip.text = element_text(size = 13)) +
  ylab("Precision") +
  xlab("") -> plt
ggsave("loc_cs_precision.pdf",plt,height = 8,width = 10)

plotdat = data.frame(Value = cs_res$Recall,
                      Method = cs_res$Method,
                      W = paste0("W = ",cs_res$W),
                      K = paste0("K = ",cs_res$K))
plotdat$K = factor(plotdat$K, levels = paste0("K = ",c(1,2,5,10)))
plotdat$Method = factor(plotdat$Method, levels = c("FINEMAP","SuSiE","PAINTOR-","PAINTOR+","SparsePro-","SparsePro+"))
ggplot(plotdat, aes(x = Method, y = Value, fill = Method)) +
#  geom_violin() +
  geom_boxplot(width = 0.2) +
  facet_grid(K ~ W) +
  theme_bw() +
  scale_fill_manual(values = c("green","royalblue","orchid","darkmagenta","orange","red")) +
  theme(legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text.x = element_text(angle = 45,size = 13,vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 13),
        strip.text = element_text(size = 13)) +
  ylim(0,1) +
  ylab("Recall") +
  xlab("") -> plt
ggsave("loc_cs_recall.pdf",plt,height = 8,width = 10)

plotdat = data.frame(Value = log10(cs_res$Size+1),
                      Method = cs_res$Method,
                      W = paste0("W = ",cs_res$W),
                      K = paste0("K = ",cs_res$K))
plotdat$K = factor(plotdat$K, levels = paste0("K = ",c(1,2,5,10)))
plotdat$Method = factor(plotdat$Method, levels = c("FINEMAP","SuSiE","PAINTOR-","PAINTOR+","SparsePro-","SparsePro+"))
ggplot(plotdat, aes(x = Method, y = Value, fill = Method)) +
#  geom_violin() +
  geom_boxplot(width = 0.2) +
  facet_grid(K ~ W) +
  theme_bw() +
  scale_fill_manual(values = c("green","royalblue","orchid","darkmagenta","orange","red")) +
  theme(legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text.x = element_text(angle = 45,size = 13,vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 13),
        strip.text = element_text(size = 13)) +
  ylab(expression(log[10](Size~of~credible~set+1))) +
  xlab("") -> plt
ggsave("loc_cs_size.pdf",plt,height = 8,width = 10)


df_time = read.table("loc_time.csv",header = T)
plotdat = data.frame(T = as.vector(tapply(df_time$time,list(df_time$K,df_time$Method),mean)),
                     SD = as.vector(tapply(df_time$time,list(df_time$K,df_time$Method),sd)),
                     Method = rep(c('FINEMAP','PAINTOR+','PAINTOR-','SparsePro+','SparsePro-','SuSiE'),each=4),
                     K = rep(c(1,2,5,10),6))
plotdat$Method = factor(plotdat$Method, levels=c('FINEMAP','SuSiE','PAINTOR-','PAINTOR+','SparsePro-','SparsePro+'))
ggplot(plotdat, aes(x=K, y=log(T), colour=Method)) + 
    geom_point() +
    geom_line() +
    scale_color_manual(values = c("green","royalblue","orchid","darkmagenta","orange","red")) +
    geom_errorbar(aes(ymin=log(T-SD), ymax=log(T+SD)),width = 0) +
    theme_bw() +
    theme(legend.position = "top",
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        strip.text = element_text(size = 13),
        legend.title = element_blank(),
        legend.text = element_text(size = 13)) +
    xlab("Number of causal variant (K)") +
#    ylim(2,5.1) +
    scale_x_continuous(breaks = c(1,2,5,10)) +
    ylab(expression(log[10](Total~computation~time))) -> plt
ggsave("loc_time.pdf",plt,height = 4,width = 4)

spenrich = read.table("loc_sparsepro_anno_22.csv")
spenrich$K = factor(paste0("K = ",spenrich$V5),levels = paste0("K = ",c(1,2,5,10)))
spenrich$W = paste0("True W = ",spenrich$V6)
spenrich$V1 = factor(spenrich$V1, levels = rev(unique(spenrich$V1)[c(1,2,3,4,9,5,6,7,8,10)]))

ggplot(spenrich, aes(x = V2, y = V1, size=-log10(V4))) +
  geom_vline(xintercept = 0, lty =2 , color = "grey") +
  geom_point() +
  geom_errorbarh(aes(xmin = V2 - 1.96 * V3,
                     xmax = V2 + 1.96 * V3,size=1),height = 0,
                 position = position_dodge(width = 0.3),show.legend = F) +
  facet_grid(K~W) +
  ylab("") +
  xlab("Estimated W") +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 13),
        legend.position = "bottom",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14)) +
  labs(color = "", size = expression(-log[10](p-value))) -> plt

ggsave("loc_sparsepro_anno_22.pdf",plt,height = 12,width = 10)

spenrich = read.table("loc_sparsepro_anno_21.csv")
spenrich$K = factor(paste0("K = ",spenrich$V5),levels = paste0("K = ",c(1,2,5,10)))
spenrich$W = paste0("True W = ",spenrich$V6)
spenrich$V1 = factor(spenrich$V1, levels = rev(unique(spenrich$V1)[c(1,2,3,4,9,5,6,7,8,10)]))

ggplot(spenrich, aes(x = V2, y = V1, size=-log10(V4))) +
  geom_vline(xintercept = 0, lty =2 , color = "grey") +
  geom_point() +
  geom_errorbarh(aes(xmin = V2 - 1.96 * V3,
                     xmax = V2 + 1.96 * V3,size=1),height = 0,
                 position = position_dodge(width = 0.3),show.legend = F) +
  facet_grid(K~W) +
  ylab("") +
  xlab("Estimated W") +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 13),
        legend.position = "bottom",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14)) +
  labs(color = "", size = expression(-log[10](p-value))) -> plt

ggsave("loc_sparsepro_anno_21.pdf",plt,height = 12,width = 10)

spenrich = read.table("loc_sparsepro_anno_20.csv")
spenrich$K = factor(paste0("K = ",spenrich$V5),levels = paste0("K = ",c(1,2,5,10)))
spenrich$W = paste0("True W = ",spenrich$V6)
spenrich$V1 = factor(spenrich$V1, levels = rev(unique(spenrich$V1)[c(1,2,3,4,9,5,6,7,8,10)]))

ggplot(spenrich, aes(x = V2, y = V1, size=-log10(V4))) +
  geom_vline(xintercept = 0, lty =2 , color = "grey") +
  geom_point() +
  geom_errorbarh(aes(xmin = V2 - 1.96 * V3,
                     xmax = V2 + 1.96 * V3,size=1),height = 0,
                 position = position_dodge(width = 0.3),show.legend = F) +
  facet_grid(K~W) +
  ylab("") +
  xlab("Estimated W") +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 13),
        legend.position = "bottom",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14)) +
  labs(color = "", size = expression(-log[10](p-value))) -> plt

ggsave("loc_sparsepro_anno_20.pdf",plt,height = 12,width = 10)

paintor = read.table("loc_paintor_anno_22.csv",header = F)

for (j in 1:nrow(paintor)) {
  paintor[j,3:12] = log((exp(paintor[j,2]) + 1)/(exp(paintor[j,3:12] + paintor[j,2]) + 1))
}

plotdat = data.frame(Value = unlist(paintor[,2:11]),
                      K = rep(paste0("K = ",paintor$V13),10),
                      W = rep(paste0("W = ",paintor$V14),10),
                      annot = rep(c("Conserved_LindbladToh",
                                    "DHS_Trynka",
                                    "H3K27ac_Hnisz",
                                    "H3K4me3_Trynka",
                                    "Transcr_Hoffman",
                                    "TSS_Hoffman",
                                    "UTR_3_UCSC",
                                    "UTR_5_UCSC",
                                    "non_synonymous",
                                    "Human_Promoter_Villar_ExAC"),each = 12))

plotdat$K <- factor(plotdat$K, levels = paste0("K = ",c(1,2,5,10)))
plotdat$annot <- factor(plotdat$annot, levels = rev(c("Conserved_LindbladToh",
                                    "DHS_Trynka",
                                    "H3K27ac_Hnisz",
                                    "H3K4me3_Trynka",
                                    "non_synonymous",
                                    "Transcr_Hoffman",
                                    "TSS_Hoffman",
                                    "UTR_3_UCSC",
                                    "UTR_5_UCSC",
                                    "Human_Promoter_Villar_ExAC")))

ggplot(plotdat, aes(x = Value, y = annot)) +
  geom_vline(xintercept = 0, lty =2 , color = "grey") +
  geom_point() +
  facet_grid(K~W) +
  ylab("") +
  xlab("Coefficient") +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 13)) +
  xlim(-5,5) -> plt

paintor = read.table("loc_paintor_anno_21.csv",header = F)

for (j in 1:nrow(paintor)) {
  paintor[j,3:12] = log((exp(paintor[j,2]) + 1)/(exp(paintor[j,3:12] + paintor[j,2]) + 1))
}

plotdat = data.frame(Value = unlist(paintor[,2:11]),
                      K = rep(paste0("K = ",paintor$V13),10),
                      W = rep(paste0("W = ",paintor$V14),10),
                      annot = rep(c("Conserved_LindbladToh",
                                    "DHS_Trynka",
                                    "H3K27ac_Hnisz",
                                    "H3K4me3_Trynka",
                                    "Transcr_Hoffman",
                                    "TSS_Hoffman",
                                    "UTR_3_UCSC",
                                    "UTR_5_UCSC",
                                    "non_synonymous",
                                    "Human_Promoter_Villar_ExAC"),each = 12))

plotdat$K <- factor(plotdat$K, levels = paste0("K = ",c(1,2,5,10)))
plotdat$annot <- factor(plotdat$annot, levels = c("Conserved_LindbladToh",
                                    "DHS_Trynka",
                                    "H3K27ac_Hnisz",
                                    "H3K4me3_Trynka",
                                    "non_synonymous",
                                    "Transcr_Hoffman",
                                    "TSS_Hoffman",
                                    "UTR_3_UCSC",
                                    "UTR_5_UCSC",
                                    "Human_Promoter_Villar_ExAC"))

ggplot(plotdat, aes(x = Value, y = annot)) +
  geom_vline(xintercept = 0, lty =2 , color = "grey") +
  geom_point() +
  facet_grid(K~W) +
  ylab("") +
  xlab("Coefficient") +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 13)) +
  xlim(-5,5) -> plt

ggsave("loc_paintor_anno_21.pdf",plt,height = 12,width = 10)

