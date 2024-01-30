#!/usr/bin/env Rscript

rm(list=ls());
library(myfunc);
source("~/mybiotools/r/myfunc.R");
options("scipen"=1, "digits"=4, stringsAsFactors=FALSE);
args<-commandArgs(T)

outd="/gpfs/ycga/project/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/hiv_meta_3datasets"
vacs="/gpfs/ycga/project/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/hiv_450k_comm_allVACS/CTmethylome"
brad="/gpfs/ycga/project/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/brad_hiv_comm_residV2/CTmethylome"
gse="/gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/GSE217633/TCA/hiv_ewas_preART"

setwd (outd)

celltype=qw("CD4T CD8T Mono Bcell NK Gran PBMC")

for (ct in celltype)
{ 
	system(paste0("mkdir ", ct))
	setwd (paste0(outd, "/", ct))
	cat ('sample.size = c(718, 436, 229)', "\n\n", file="ewas_meta_anal_mult_para.R", append=F)
	cat (paste0("glm.out = c(", 
				paste0('"', vacs, "/", ct, "_methylome_hiv_ewas/glm_pv_adj_anno.xls", '"', ", "), 
				paste0('"', brad, "/", ct, "_methylome_hiv_ewas/glm_pv_adj_anno.xls", '"', ", "),
				paste0('"', gse, "/", ct, "_methylome_hiv_ewas/glm_pv_adj_anno.xls", '"'), 
				")"), "\n\n", file="ewas_meta_anal_mult_para.R", append=T)
	cat ('platform = c("450K", "EPIC", "EPIC")', "\n\n", file="ewas_meta_anal_mult_para.R", append=T)
	system("Rscript ~/mybiotools/r/ewas_meta_anal_mult_v2.R")
	system("META_post_metal.sh 0.1 &")
	setwd(outd)
}

