library(ggplot2)

df_anno = read.table("../doc/dat_trait_anno.txt", sep='\t', header=T)

ggplot(df_anno,aes(x = exp(W), y = Annotation)) +
  geom_point(position = position_dodge(width = 0.5)) +
  facet_grid(~Trait) +
  geom_vline(xintercept = 1, color = "grey", lty = 2) +
  geom_errorbarh(aes(xmin = exp(W - 1.96 * W_se),
                     xmax = exp(W + 1.96 * W_se)), height = 0,
                 position = position_dodge(width = 0.5), show.legend = F) +
  theme_classic() +
  ylab("") +
  xlab("Enrichment fold") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 13),
        legend.position = "bottom",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14)) -> plt
ggsave("../Fig/dat_trait_anno.pdf", plt, height = 12, width = 10)


df_cts = read.table("../doc/dat_trait_cts.txt", sep='\t', header=T)
ggplot(df_cts,aes(x = exp(W), y = Tissue, color = Method)) +
  geom_point(position = position_dodge(width = 0.5)) +
  facet_grid(~Trait) +
  geom_vline(xintercept = 1, color = "grey", lty = 2) +
  geom_errorbarh(aes(xmin = exp(W - 1.96 * W_se),
                     xmax = exp(W + 1.96 * W_se)),height = 0,
                 position = position_dodge(width = 0.5),show.legend = F) +
  scale_color_manual(values = c("orange","red","royalblue","darkblue")) +
  theme_classic() +
  ylab("") +
  xlab("Enrichment fold") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 13),
        legend.position = "bottom",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14)) +
  labs(color = "", size = expression(-log[10](p-value))) -> plt
ggsave("../Fig/dat_trait_cts.pdf",height = 7,width = 12)


## LocusZoom 

loci = 'GCKR'
tr = c('eGFR', 'glucose', 'prate', 'gammaGT')
traitname = c('eGFR', 'Glucose', 'Pulse rate', 'gamma-GT')
lead = '2.27730940.T.C'
rsid = 'rs1260326'
trait_mapping <- as.list(setNames(traitname, tr))
for (trait in tr){
  print(trait)
  sumstats = read.table(paste0('../doc/LocusZoom/', loci, '_', trait, '_', lead,'.txt'), header=T)
  sumstats$lead <- as.factor(ifelse(sumstats$lead=='True',"Lead causal variant","Other variant"))
  ggplot(sumstats,aes(x = pos / 1000000, y = -log10(p), color = ld^2, shape = lead, size = lead)) +
    geom_point() +
    theme_classic() +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 13),
          title = element_text(size = 14)) +
    ylab(expression(-log[10](p-value))) +
    labs(color = expression(r^2)) +
    scale_color_stepsn(n.breaks = 6, colours = c("darkblue","blue","green","orange","red")) +
    xlab(paste0("Chr",sumstats$V2[1]," (Mb)")) +
    scale_shape_manual(values = c(18,16)) +
    scale_size_manual(values = c(4,2)) +
    labs(shape = "", size = "") +
    ggtitle(paste0(trait_mapping[[trait]],": ",loci," (",rsid,")")) -> plt
  ggsave(paste0("../Fig/",trait,"_",loci,"_",lead,".pdf"),plt,height = 6,width = 7)
  ggplot(sumstats,aes(x = pos / 1000000, y = PIP, color = ld^2, shape = lead, size = lead)) +
    geom_point() +
    theme_classic() +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 13),
          title = element_text(size = 14)) +
    ylab('PIP') +
    ylim(0,1) +
    labs(color = expression(r^2)) +
    scale_color_stepsn(n.breaks = 6, colours = c("darkblue","blue","green","orange","red")) +
    xlab(paste0("Chr",sumstats$V2[1]," (Mb)")) +
    scale_shape_manual(values = c(18,16)) +
    scale_size_manual(values = c(4,2)) +
    labs(shape = "", size = "") +
    ggtitle(paste0(trait_mapping[[trait]],": ",loci," (",rsid,")")) -> plt
  ggsave(paste0("../Fig/",trait,"_",loci,"_",lead,".pip.pdf"),plt,height = 6,width = 7)
}
