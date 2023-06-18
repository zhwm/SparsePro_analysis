library(ggplot2)

paintor = read.table("../doc/loc_paintor_anno.txt",header = T)

for (j in 1:nrow(paintor)) {
  paintor[j,2:11] = log((exp(paintor[j,1]) + 1)/(exp(paintor[j,2:11] + paintor[j,1]) + 1))
}

plotdat = data.frame(Value = c(unlist(paintor[,2:11]),
                               rep(c(1,1,1,1,0,0,0,0,1,0), each=12) * rep(paintor$W,10)),
                     K = rep(paste0("K = ",paintor$K),20),
                     W = rep(paste0("W = ",paintor$W),20),
                     annot = rep(rep(c("Conserved_LindbladToh",
                                   "DHS_Trynka",
                                   "H3K27ac_Hnisz",
                                   "H3K4me3_Trynka",
                                   "Transcr_Hoffman",
                                   "TSS_Hoffman",
                                   "UTR_3_UCSC",
                                   "UTR_5_UCSC",
                                   "non_synonymous",
                                   "Human_Promoter_Villar_ExAC"),each = 12),2),
                     color = rep(c('Estimated', 'Simulated'), each=12*10))

plotdat$K <- factor(plotdat$K, levels = paste0("K = ",c(1,2,5,10)))
plotdat$annot <- factor(plotdat$annot, 
                        levels = rev(c("Conserved_LindbladToh",
                                                      "DHS_Trynka",
                                                      "H3K27ac_Hnisz",
                                                      "H3K4me3_Trynka",
                                                      "non_synonymous",
                                                      "Transcr_Hoffman",
                                                      "TSS_Hoffman",
                                                      "UTR_3_UCSC",
                                                      "UTR_5_UCSC",
                                                      "Human_Promoter_Villar_ExAC")))

ggplot(plotdat, aes(x = Value, y = annot, color=color, shape=color)) +
  geom_vline(xintercept = 0, lty =2 , color = "grey") +
  geom_point(position = position_dodge(width = 0.5)) +
  facet_grid(K~W) +
  scale_color_manual(values = c("darkblue","gold")) +
  ylab("") +
  xlab("Enrichment weight") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 13),
        legend.position = "top",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14)) +
  labs(color = "", shape="") +
  xlim(-5,5) -> plt
ggsave("../Fig/loc_paintor_anno.pdf", plt, width = 10, height = 12)


spenrich = read.table("../doc/loc_sparsepro_anno.txt",header = T)
spenrich$K = factor(paste0("K = ",spenrich$K),levels = paste0("K = ",c(1,2,5,10)))
spenrich$W = paste0("True W = ",spenrich$W)
spenrich$Annotation = factor(spenrich$Annotation, levels = rev(unique(spenrich$Annotation)[c(1,2,3,4,9,5,6,7,8,10)]))

spplot <- data.frame(Value = c(spenrich$Simulated_W, spenrich$Estimated_W),
                     sd = c(rep(0, nrow(spenrich)), spenrich$Standard_error),
                     K = rep(spenrich$K, 2),
                     W = rep(spenrich$W, 2),
                     anno = rep(spenrich$Annotation, 2),
                     color = rep(c('Simulated', 'Estimated'), each=nrow(spenrich)))

ggplot(spplot, aes(x = Value, y = anno, color=color, shape=color)) +
  geom_vline(xintercept = 0, lty =2 , color = "grey") +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbarh(aes(xmin = Value - 1.96 * sd,
                     xmax = Value + 1.96 * sd), height = 0,
                 position = position_dodge(width = 0.5),show.legend = F) +
  facet_grid(K~W) +
  ylab("") +
  xlab("Enrichment weight") +
  scale_color_manual(values = c("darkblue","gold")) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 13),
        legend.position = "top",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14)) +
  labs(color = "", shape= "") -> plt
ggsave("../Fig/loc_sparsepro_anno.pdf", plt, height = 12, width = 10)


