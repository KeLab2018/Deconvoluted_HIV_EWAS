#!/usr/bin/env Rscript

rm(list=ls());
library(myfunc);
source("~/mybiotools/r/myfunc.R");
options("scipen"=1, "digits"=4, stringsAsFactors=FALSE);
args<-commandArgs(T)

GSE.d = "/gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/GSE217633/TCA/hiv_ewas_preART"
vacs.d = "/gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/hiv_450k_comm_allVACS/CTmethylome_6CT"
celltypes_paper_name = c("PBMC", "CD4", "CD8", "B", "NK", "M", "Gran")
celltypes = c("PBMC", "CD4T", "CD8T", "Bcell", "NK", "Mono", "Gran")
(vacs.ewas.f = paste0(vacs.d, "/", celltypes, "_methylome_hiv_ewas/glm_pv_adj_anno.xls"))
(gse.ewas.f = paste0(GSE.d, "/", celltypes, "_methylome_hiv_ewas/glm_pv_adj_anno.xls"))

for (i in 1:length(celltypes)) {
	ct = celltypes_paper_name[i]
	outf = paste0(ct, "_diff_venn.pdf")
	comm_out = paste0(ct, "_comm_out.txt")
	#fdr & pv 0.05
	CMD = paste("plotvenn.pl -norb -i", vacs.ewas.f[i], gse.ewas.f[i], "-row 'UCSC_REFGENE_NAME' 'UCSC_REFGENE_NAME' -s ';' -d T -r gg -si 0.05 -sc BH_fdr_adj Pvalue -k Cohort_1 Cohort_3 -o", outf, "-comm_out", comm_out, "&")
	#fdr 0.05
	#CMD = paste("plotvenn.pl -norb -i", vacs.ewas.f[i], gse.ewas.f[i], "-row 'UCSC_REFGENE_NAME' 'UCSC_REFGENE_NAME' -s ';' -d T -r gg -si 0.05 -sc BH_fdr_adj -k Cohort_1 Cohort_3 -o", outf, "-comm_out", comm_out, "&")
	#pv 0.05
	#CMD = paste("plotvenn.pl -norb -i", vacs.ewas.f[i], gse.ewas.f[i], "-row 'UCSC_REFGENE_NAME' 'UCSC_REFGENE_NAME' -s ';' -d T -r gg -si 0.05 -sc Pvalue -k Cohort_1 Cohort_3 -o", outf, "-comm_out", comm_out, "&")
	#probe based
	#CMD = paste("plotvenn.pl -norb -i", vacs.ewas.f[i], gse.ewas.f[i], "-row 'probe' 'probe' -d T -r gg -si 0.05 -sc Pvalue -k Cohort_1 Cohort_3 -o", outf, "-comm_out", comm_out, "&")
	print (CMD)
	system (CMD)
}

