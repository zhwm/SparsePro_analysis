library(ggplot2)
traits = c('eGFR','FFR','glucose','gammaGT','prate')
methods = c('sparsepro','sparsepro_anno','susie','susie_poly')
W_cts = se_cts = p_cts = W_anno = se_anno = p_anno = c()
for (trait in traits){
    for (method in methods){
        wsep_cts=read.table(paste0("../dat/UKB-GWAS/GWAS/",trait,'/',method,'/',trait,'_cts.wsep'))
        wsep_anno=read.table(paste0("../dat/UKB-GWAS/GWAS/",trait,'/',method,'/',trait,'_anno.wsep'))
        W_cts= c(W_cts,wsep_cts$W)
        se_cts = c(se_cts,wsep_cts$se)
        p_cts = c(p_cts,wsep_cts$p)
        W_anno= c(W_anno,wsep_anno$W)
        se_anno = c(se_anno,wsep_anno$se)
        p_anno = c(p_anno,wsep_anno$p)
    }
}
df_mt=expand.grid(c('SparsePro-','SparsePro+','SuSiE','SuSiE+PolyFun'),c('eGFR','FFR','Glucose','gamma-GT','Pulse rate'))
df_cts = data.frame(W = W_cts,
                    se = se_cts,
                    p = p_cts,
                    trait = rep(df_mt$Var2,each=10),
                    method = rep(df_mt$Var1,each=10),
                    anno = rep(rownames(wsep_cts),20))
df_anno = data.frame(W = W_anno,
                    se = se_anno,
                    p = p_anno,
                    trait = rep(df_mt$Var2,each=10),
                    method = rep(df_mt$Var1,each=10),
                    anno = rep(rownames(wsep_anno),20))
df_anno$anno = factor(df_anno$anno,levels=rev(rownames(wsep_anno)[c(1,2,3,4,9,5,6,7,8,10)]))
ggplot(df_anno,aes(x = exp(W), y = anno, color = method, size = -log10(p))) +
  geom_point(position = position_dodge(width = 0.5)) +
  facet_grid(~trait) +
  geom_vline(xintercept = 1, color = "grey", lty = 2) +
  geom_errorbarh(aes(xmin = exp(W - 1.96 * se),
                     xmax = exp(W + 1.96 * se),size = 1),height = 0,
                 position = position_dodge(width = 0.5),show.legend = F) +
  scale_color_manual(values = c("orange","red","royalblue","darkblue")) +
  theme_classic() +
  ylab("") +
  xlab("Relative enrichment") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 13),
        legend.position = "bottom",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14)) +
  labs(color = "", size = expression(-log[10](p-value))) -> plt
ggsave("enrichment.anno.pdf",height = 7,width = 12)
df_cts$anno[df_cts$anno=="Adrenal_Pancreas"] <- "Adrenal/Pancreatic"
df_cts$anno[df_cts$anno=="Connective_Bone"] <- "Connective"
df_cts$anno[df_cts$anno=="SkeletalMuscle"] <- "Musculoskeletal"
df_cts$anno[df_cts$anno=="GI"] <- "Gastrointestinal"
ggplot(df_cts,aes(x = exp(W), y = anno, color = method, size = -log10(p))) +
  geom_point(position = position_dodge(width = 0.5)) +
  facet_grid(~trait) +
  geom_vline(xintercept = 1, color = "grey", lty = 2) +
  geom_errorbarh(aes(xmin = exp(W - 1.96 * se),
                     xmax = exp(W + 1.96 * se),size = 1),height = 0,
                 position = position_dodge(width = 0.5),show.legend = F) +
  scale_color_manual(values = c("orange","red","royalblue","darkblue")) +
  theme_classic() +
  ylab("") +
  xlab("Relative enrichment") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 13),
        legend.position = "bottom",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14)) +
  labs(color = "", size = expression(-log[10](p-value))) -> plt
ggsave("enrichment.cts.pdf",height = 7,width = 12)


h2 = read.table("eFggp.h2")
h2$trait = rep(c('eGFR','FFR','gamma-GT','Glucose','Pulse rate'),each=4)
h2$trait = factor(h2$trait,levels=c('eGFR','FFR','Glucose','gamma-GT','Pulse rate'))
h2$method = rep(c('SparsePro+','SparsePro-','SuSiE','SuSiE+PolyFun'),5)
ggplot(h2,aes(x = method, y = V5, fill = method)) +
  geom_bar(stat = "identity") +
  facet_grid(~trait) + 
  scale_fill_manual(values = c("orange","red",'royalblue','darkblue')) +
  theme_classic() +
  ylab("Proportion of variance explained by causal set") +
  ylim(0,0.15) +
#  xlab("Method") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 13),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 14)) -> plt
ggsave("h2.pdf",height = 7,width = 9)
