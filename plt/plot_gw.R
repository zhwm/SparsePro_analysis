library(ggplot2)

spenrich = read.table("../doc/gw_sparsepro_anno.txt",header = T)
spenrich$Annotation = factor(spenrich$Annotation, levels = rev(unique(spenrich$Annotation)[c(1,2,3,4,9,5,6,7,8,10)]))
spenrich$W = paste0("True W = ",spenrich$W)

spplot <- data.frame(Value = c(spenrich$Simulated_W, spenrich$Estimated_W),
                     sd = c(rep(0, nrow(spenrich)), spenrich$Standard_error),
                     W = rep(spenrich$W, 2),
                     anno = rep(spenrich$Annotation, 2),
                     color = rep(c('Simulated', 'Estimated'), each=nrow(spenrich)))

ggplot(spplot, aes(x = Value, y = anno, color=color, shape=color)) +
  geom_vline(xintercept = 0, lty =2 , color = "grey") +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbarh(aes(xmin = Value - 1.96 * sd,
                     xmax = Value + 1.96 * sd), height = 0,
                 position = position_dodge(width = 0.5),show.legend = F) +
  facet_grid(W~.) +
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
ggsave("../Fig/gw_sparsepro_anno.pdf", plt, height = 12, width = 10)


spenrich = read.table("../doc/gw_sparsepro_w1.txt",header = T)
spenrich$Annotation = factor(spenrich$Annotation, levels = rev(unique(spenrich$Annotation)[c(1,2,3,4,9,5,6,7,8,10)]))
spenrich$W = paste0("True W = ",spenrich$W)

spplot <- data.frame(Value = c(spenrich$Simulated_W, spenrich$Estimated_W),
                     sd = c(rep(0, nrow(spenrich)), spenrich$Standard_error),
                     W = rep(spenrich$W, 2),
                     anno = rep(spenrich$Annotation, 2),
                     color = rep(c('Simulated', 'Estimated'), each=nrow(spenrich)))

ggplot(spplot, aes(x = Value, y = anno, color=color, shape=color)) +
  geom_vline(xintercept = 0, lty =2 , color = "grey") +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbarh(aes(xmin = Value - 1.96 * sd,
                     xmax = Value + 1.96 * sd), height = 0,
                 position = position_dodge(width = 0.5),show.legend = F) +
  facet_grid(W~.) +
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
ggsave("../Fig/gw_sparsepro_w1.pdf", plt, height = 12, width = 10)

spenrich = read.table("../doc/gw_polyfun_anno.txt",header = T)
spenrich$ANNOTATION = factor(spenrich$ANNOTATION, levels = rev(unique(spenrich$ANNOTATION)[c(1,2,3,4,9,5,6,7,8,10)]))
spenrich$W = paste0("True W = ",spenrich$W)
spplot <- data.frame(Value = c(spenrich$Simulated_W, spenrich$ANNOTATION_COEFFICIENT),
                     W = rep(spenrich$W, 2),
                     anno = rep(spenrich$ANNOTATION, 2),
                     color = rep(c('Simulated', 'Estimated'), each=nrow(spenrich)))

ggplot(spplot, aes(x = Value, y = anno, color=color, shape=color)) +
  geom_vline(xintercept = 0, lty =2 , color = "grey") +
  geom_point(position = position_dodge(width = 0.5)) +
#  geom_errorbarh(aes(xmin = Value - 1.96 * sd,
#                     xmax = Value + 1.96 * sd), height = 0,
#                 position = position_dodge(width = 0.5),show.legend = F) +
  facet_grid(W~.) +
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
ggsave("../Fig/gw_polyfun_anno.pdf", plt, height = 12, width = 10)
