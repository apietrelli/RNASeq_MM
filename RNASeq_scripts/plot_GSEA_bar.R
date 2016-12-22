# Load library
library(ggplot2)
library(data.table)
library(reshape2)
library(dplyr)
# Load file
HET_path = "../GSEA/my_analysis.Gsea.1461755156354.Hallmark/gsea_report_for_HET_1461755156354.xls"
WT_path = "../GSEA/my_analysis.Gsea.1461755156354.Hallmark/gsea_report_for_WT_1461755156354.xls"
HET_report = fread(HET_path)
WT_report = fread(WT_path)
# Reshape data and merge DT from WT and HET
HET_report = HET_report %>% 
  filter(`FDR q-val` <= 0.1) %>%
  select(NAME, SIZE, NES, ES, `FDR q-val`, `FWER p-val`)
HET_report[,class:="HET"]

WT_report = WT_report %>% 
  filter(`FDR q-val` <= 0.1) %>%
  select(NAME, SIZE, NES, ES, `FDR q-val`, `FWER p-val`) 
WT_report[,class:="WT"]

gsea_report = rbindlist(list(HET_report,WT_report))

setorder(gsea_report, NES)
gsea_report$NAME = factor(gsea_report$NAME, gsea_report$NAME)
gsea_report$hjust <- ifelse(gsea_report$NES > 0, 1.3, -0.4)
#Alpha trasparency based on FDR
gsea_report$alpha = ifelse(gsea_report$`FDR q-val` < 0.05, 1, 0)

palette = c("firebrick1","steelblue3")

plot = ggplot(gsea_report, aes(NAME, NES,
                               label = gsub("_"," ",gsub("hallmark_", "", tolower(NAME))),
                hjust = hjust)) +
  geom_bar(stat = "identity", width = 0.5, aes(fill = class,
                                               alpha = alpha)) +
  scale_alpha(range=c(0.5,1)) +
  theme_bw() +
  geom_text(aes(family="Helvetica", y=0), size=5) +
  coord_flip() +
  xlab("") +
  theme(legend.position = "none",
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid = element_blank()) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = palette)

plot

pdf("../GSEA/HSC_Hallmark_R_0.1.pdf", paper ="a4")
print(plot)
dev.off()
