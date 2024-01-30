#!/usr/bin/env Rscript

rm(list=ls());
library(myfunc);
source("~/mybiotools/r/myfunc.R");
options("echo"=T, "scipen"=1, "digits"=4, stringsAsFactors=FALSE);
args<-commandArgs(T)

tss.dat = r.table("out.tss.xls", T)

anno = r.table("/gpfs/ycga/project/xu_ke/xz345/work/VACS/Methylation/illumina_870k_0317/annotation/anno_850k_brief.txt", T)
alltss= gene_symbol_trim(anno$UCSC_REFGENE_GROUP)
alltss2 = data.frame(table(alltss))
colnames(alltss2) = c("TSS", "Count")
alltss2

promoter.name = c("1stExon", "5'UTR", "TSS1500", "TSS200")
promoter.name = c("5'UTR", "TSS1500", "TSS200")
promoter.name = c("TSS1500", "TSS200")

calc_promoter_ind <- function (x) {
	anno = tss.dat[x,]
	total.count = sum(anno$Count)
	promoter.ind = which(anno$TSS %in% promoter.name)
	promoter.count = sum(anno$Count[promoter.ind])
	return (c(promoter.count, total.count))
}
calc_promoter <- function (anno) {
	total.count = sum(anno$Count)
	promoter.ind = which(anno$TSS %in% promoter.name)
	promoter.count = sum(anno$Count[promoter.ind])
	return (c(promoter.count, total.count))
}

ct.count = tapply(1:nrow(tss.dat), tss.dat$Celltype, FUN=calc_promoter_ind)

allprobe.count = calc_promoter(alltss2)

#C1 <- # 显著的 CpG 总数
#C2 <- # 显著的 CpG 中在启动子的数量
#T1 <- # 所有 CpG 的总数
#T2 <- # 所有 CpG 中在启动子的数量
prom.test <- function (c, t) {
	C2 = c[1]
	C1 = c[2]
	T2 = t[1]
	T1 = t[2]
	table <- matrix(c(C2, C1 - C2, T2, T1 - T2), nrow = 2)
	test <- chisq.test(table)
	out = c(test$p.value, test$statistic)
	names(out) = c("Pvalue", "X-squared")
	return (out)
}

prom.test.results = lapply(ct.count, function(x){
		   prom.test(x, allprobe.count)
})
prom.test.results



