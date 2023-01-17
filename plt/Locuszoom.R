setwd("~/Desktop/SparsePro/Zenodo/dat/UKB-GWAS/GWAS/LocusZoom/")
library(ggplot2)

for (ldname in list.files()[grep(".ld",list.files())]) {
    traitname <- strsplit(ldname,"_")[[1]][1]
    genename <- strsplit(ldname,"_")[[1]][2]
    snpname <- gsub(".ld","",strsplit(ldname,"_")[[1]][3])
    sumstatsname <- paste0(genename,"_",traitname,".locus.sumstats")
    traitdisplay <- ifelse(traitname == "gammaGT","gamma-GT",
                          ifelse(traitname == "prate","Pulse rate",
                                ifelse(traitname == "glucose","Glucose",traitname)))
    sumstats <- read.table(sumstatsname)
    sumstats <- sumstats[!duplicated(sumstats$V1),]
    rownames(sumstats) <- sumstats$V1
    ld <- read.table(ldname,header = T)
    ld <- ld[!duplicated(ld$SNP_B),]
    rownames(ld) <- ld$SNP_B
    retainID <- intersect(rownames(sumstats),rownames(ld))
    sumstats <- sumstats[retainID,]
    ld <- ld[retainID,]
    sumstats$ld <- ld$R2
    sumstats$lead <- as.factor(ifelse(sumstats$V1==snpname,"Lead causal variant","Other variant"))
    pip <- read.table(paste0("~/Desktop/SparsePro/Zenodo/dat/UKB-GWAS/GWAS/",traitname,"/sparsepro_anno/",traitname,"_",sumstats$V2[1],".pip"))
    pip$ID <- unlist(lapply(strsplit(pip$V1,"[.]"), function(x) x[2]))
    pip <- pip[pip$ID %in% sumstats$V3,]
    sumstats$pip <- pip$V3                   
    ggplot(sumstats,aes(x = V3 / 1000000, y = -log10(V6), color = ld, shape = lead, size = lead)) +
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
      ggtitle(paste0(traitdisplay,": ",genename," (",snpname,")")) -> plt
    ggsave(paste0("~/Desktop/SparsePro/Zenodo/doc/",traitname,"_",genename,"_",snpname,".pdf"),plt,height = 6,width = 7)
    ggplot(sumstats,aes(x = V3 / 1000000, y = pip, color = ld, shape = lead, size = lead)) +
      geom_point() +
      theme_classic() +
      theme(axis.title = element_text(size = 14),
           axis.text = element_text(size = 13),
           title = element_text(size = 14)) +
      ylab("PIP") +
      labs(color = expression(r^2)) +
      scale_color_stepsn(n.breaks = 6, colours = c("darkblue","blue","green","orange","red")) +
      xlab(paste0("Chr",sumstats$V2[1]," (Mb)")) +
      scale_shape_manual(values = c(18,16)) +
      scale_size_manual(values = c(4,2)) +
      labs(shape = "", size = "") +
      ggtitle(paste0(traitdisplay,": ",genename," (",snpname,")")) -> plt
    ggsave(paste0("~/Desktop/SparsePro/Zenodo/doc/",traitname,"_",genename,"_",snpname,".pip.pdf"),plt,height = 6,width = 7)
}


