library(ggplot2)

results = read.table("gw_precision_recall_size_time.csv",header = T)

plotdat = data.frame(Mean = c(as.vector(tapply(results$Precision, list(results$Method,results$W), mean)),
                                 as.vector(tapply(results$Recall, list(results$Method,results$W), mean)),
                                 as.vector(tapply(results$Size, list(results$Method,results$W), mean)),
                                 as.vector(tapply(results$Time, list(results$Method,results$W), mean))),
                      SD = c(as.vector(tapply(results$Precision, list(results$Method,results$W), sd)),
                                 as.vector(tapply(results$Recall, list(results$Method,results$W), sd)),
                                 as.vector(tapply(results$Size, list(results$Method,results$W), sd)),
                                 as.vector(tapply(results$Time, list(results$Method,results$W), sd))),
                      method = c("SparsePro-","SparsePro+","SuSiE","SuSiE+PolyFun"),
                      W = rep(c(0,1,2),each = 4),
                      metric = rep(c("Precision","Recall","Size","Time"),each = 12))
plotdat <- data.frame(Value = c(results$Precision,
                                results$Recall,
                                results$Size,
                                results$Time / 3600),
                      Method = results$Method,
                      W = paste0("W = ",results$W),
                      metric = rep(c("Precision","Recall","Size of credible set","Time (hour)"),each = nrow(results)))
ggplot(plotdat, aes(x = Method, y = Value, fill = Method)) +
  geom_violin() +
  geom_boxplot(width = 0.2) +
  facet_grid(metric ~ W, scales = "free") +
  theme_bw() +
  scale_fill_manual(values = c("orange","red","royalblue","darkblue")) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45,size = 13,vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 13),
        strip.text = element_text(size = 13)) +
  ylab("") +
  xlab("") -> plt
ggsave("summary_gw_precision_recall_size_time.pdf",plt,height = 10,width = 10)

sp = read.table("gw_sparsepro_wsep.csv",header = T,sep = "\t")
poly = read.table("gw_polyfun_annot.csv")
plotdat = data.frame(W = as.vector(tapply(sp$estimated.W,list(sp$annotation,sp$real.W),mean)),
                      SD = as.vector(tapply(sp$estimated.W,list(sp$annotation,sp$real.W),sd)),
                      trueW = rep(paste0("True W = ",c(0,1,2)),each = 10),
                      annotation =rownames(tapply(sp$estimated.W,list(sp$annotation,sp$real.W),sd)) )
plotdat$annotation = factor(plotdat$annotation, levels = rev(rownames(tapply(sp$estimated.W,list(sp$annotation,sp$real.W),sd))[c(1,2,3,4,6,5,7,8,9,10)]))
ggplot(plotdat,aes(x = W, y = annotation)) +
  facet_grid(~trueW) +
  geom_vline(xintercept = 0, lty = 2, color = "grey") +
  geom_point() +
  geom_errorbarh(aes(xmin = W - 1.96 * SD,
                     xmax = W + 1.96 * SD),height = 0) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 13)) +
  xlab("Estimated W") +
  ylab("") -> plt
ggsave("gw_sparsepro_wsep.pdf",plt,height = 5,width = 12)

plotdat = data.frame(W = as.vector(tapply(poly$V2,list(poly$V1,poly$V3),mean)),
                      SD = as.vector(tapply(poly$V2,list(poly$V1,poly$V3),sd)),
                      trueW = rep(paste0("True W = ",c(0,1,2)),each = 20),
                      annotation =rownames(tapply(poly$V2,list(poly$V1,poly$V3),mean)) )
plotdat$group = c("Common","Low-frequency")
plotdat$annotation = gsub("_lowfreq","",gsub("_common","",plotdat$annotation))
plotdat$annotation = factor(plotdat$annotation, levels = rev(unique(plotdat$annotation)[c(1,2,3,4,6,5,7,8,9,10)]))
ggplot(plotdat,aes(x = W, y = annotation, shape = group)) +
  facet_grid(~trueW, scales = "free") +
  geom_vline(xintercept = 0, lty = 2, color = "grey") +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbarh(aes(xmin = W - 1.96 * SD,
                     xmax = W + 1.96 * SD),height = 0,
                 position = position_dodge(width = 0.5)) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 13),
        axis.text.x = element_text(size = 8),
        strip.text = element_text(size = 13),
        legend.title = element_blank(),
        legend.text = element_text(size = 13),
        legend.position = "bottom") +
  xlab("Coefficient") +
  ylab("") +
  scale_shape_manual(values = c(24,25)) -> plt
ggsave("gw_polyfun_annot.pdf",plt,height = 5,width = 12)
