#!/usr/bin/env Rscript

rm(list=ls());
library(myfunc);
source("~/mybiotools/r/myfunc.R");
options("scipen"=1, "digits"=4, stringsAsFactors=FALSE);
args<-commandArgs(T)

meta.d = "/gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/hiv_meta_residV2"
gse.d = "/gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/GSE217633/TCA/hiv_ewas_preART"
celltypes_paper_name = c("PBMC", "CD4", "CD8", "B", "NK", "M")
celltypes = c("PBMC", "CD4T", "CD8T", "Bcell", "NK", "Mono")
(gse.ewas.f = paste0(gse.d, "/", celltypes, "_methylome_hiv_ewas/top5000.xls"))
(meta.ewas.f = paste0(meta.d, "/", celltypes, "/top5000.xls"))

for (i in 1:length(celltypes)) {
	ct = celltypes_paper_name[i]
	outf = paste0(ct, "_diff_venn.pdf")
	comm_out = paste0(ct, "_comm_out.txt")
	CMD = paste("plotvenn.pl -i", gse.ewas.f[i], meta.ewas.f[i], "-row 'UCSC_REFGENE_NAME' 'Nearest gene (s)' -s ';' -d T -r gg -si 0.05 -sc BH_fdr_adj -k GSE217633 Meta -o", outf, "-comm_out", comm_out, "&")
	print (CMD)
	system (CMD)
}

