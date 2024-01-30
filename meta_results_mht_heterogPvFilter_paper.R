#!/usr/bin/env Rscript
# to draw UpSetR & overall MHT plot
rm(list=ls());
source("~/mybiotools/r/myfunc.R");
options("scipen"=1, "digits"=4, stringsAsFactors=FALSE);
args<-commandArgs(T)
outkw = args[1] # CT methylome ewas folds key word
original_pbmc_glmF.fold = args[2] 
fdr.cutoff = as.numeric(args[3]) # for UpSetR fdr cutoff
legend.po = args[4]
legend.cex = as.numeric(args[5])
glmf = args[6]

if (is.na(outkw)) outkw = ""
if (is.na(original_pbmc_glmF.fold)) original_pbmc_glmF.fold = "PBMC"
if (is.na(fdr.cutoff)) fdr.cutoff = 0.05
if (is.na(legend.po)) legend.po = "top"
if (is.na(legend.cex)) legend.cex = 2.2
if (is.na(glmf)) glmf = "EWAS_Meta_analysis_sampleSize_adj.xls"

library(ggplot2)
library(UpSetR)

###debug
#outkw = "_methylome_hiv_ewas"
#outkw = "_InHIVP_Rmethy_Csmkcigs_Ccelltype_Cadhmed_C30contPC_C30resiPC"
###debug

cols <- rainbow(10, alpha=0.7);
library(RColorBrewer)
my_cols2 <- brewer.pal(8, "Set1")[c(1,2,3,4,5,8)]
my_cols <- rainbow(6, alpha = 0.9)
my_cols[2] = my_cols2[5]
my_cols = add.alpha(c("red", "orange", "green", "purple", "blue"), alpha = 0.9)
my_cols_pbmc = add.alpha("black", alpha = 0.9)

limit_hight_num = 20
premht.PBMC = NULL

celltypes = c("CD4T", "CD8T", "Mono", "Bcell", "NK", "Gran")
(ewas_fold = paste0(celltypes, outkw))

celltypes_paper = c("CD4T", "CD8T", "Bcell", "NK", "Mono")
celltypes_paper_name = c("CD4", "CD8", "B", "NK", "M")
celltypes.i = match(celltypes_paper, celltypes)
(celltypes = celltypes[celltypes.i])
(ewas_fold = ewas_fold[celltypes.i])


listInput = vector("list", length=(length(celltypes)))
names(listInput) = c(celltypes_paper_name)
fdr.listInput = vector("list", length=(length(celltypes)))
names(fdr.listInput) = c(celltypes_paper_name)
original_pbmc_glmF = paste0(original_pbmc_glmF.fold, "/", glmf)

#450k
#anno.450k.f = "~/work/VACS/Methylation/rawdata/methy_probe_annotion_brief.xls"
#pr.450k = r.table(anno.450k.f,T)[,1]

pbmc.p.cutoff = c()
if (file.exists(original_pbmc_glmF))
{
	glm = r.table(original_pbmc_glmF, T)
	glm_top500 = glm[,1][1:500]
	head(glm_top500)
	listInput[["PBMC"]] = glm_top500

#to find PBMC pv cutoff
	#glm.pv = glm[, c("Chromosome", "Position", "Meta P")]
	glm.pv = glm[, c("Chromosome", "Position", "Meta P")]
	(pbmc.p.cutoff = determine_bonf_fdr_cutoff(glm.pv, 3, 0.05))

	glm = glm[which(glm$"Hat P" > 0.05), ]
	#for fdr cutoff plot
	fdr.i = which(glm$BH_fdr_adj <= fdr.cutoff)
	if (length(fdr.i) > 0) {
		fdr.listInput[["PBMC"]] = glm[fdr.i, 1]
	} else {
		fdr.listInput[["PBMC"]] = NULL
	}

	#for draw PBMC mht
	ct = "PBMC"
	fold = "PBMC_hiv_ewas"
	ctewas = glm
	#premht = data.frame(ctewas[, c("probe", "Chromosome", "Position", "Meta P")])
	premht = data.frame(ctewas[, c("probe", "Chromosome", "Position", "Meta P")])
	premht$hight.col = my_cols_pbmc
	premht$CT = ct
	h.i = which(premht[,4] <= pbmc.p.cutoff[1])
	symbols = gene_symbol_trim(ctewas[h.i, "Nearest gene(s)"])
	#if (length(symbols) > 0) symbols = paste0(ct, "_", symbols)
	gene.i = which(symbols != "")
	premht[h.i[gene.i], 1] = symbols[gene.i]
	premht.PBMC = premht
}

#to draw UpSetR
i = 1
pb <- txtProgressBar(min = 0, max = length(celltypes), style = 3)
for (i in 1:length(celltypes)) {
	ct = celltypes_paper_name[i]
	fold = ewas_fold[i]
	glm.500.f = paste0(fold,"/top5000.xls")
	if (file.exists(glm.500.f)) {
		ctewas_top500 = r.table(glm.500.f, T)[1:500, 1]
	} else {
		next;
	}
	#ctewas_top500 = r.table(paste0(fold,"/top500.xls"), T)[,1]
	listInput[[ct]] = ctewas_top500
	setTxtProgressBar(pb, i)
}
close(pb)
listInput = listInput[!sapply(listInput, is.null)]
source("~/mybiotools/r/myupset.R")
#ct.n.i = match(names(listInput), c(celltypes_paper_name, "PBMC"))
ct.n.i = match(names(listInput)[order(sapply(listInput, length))], c(celltypes_paper_name, "PBMC"))
upset.col = c(my_cols, my_cols_pbmc)[ct.n.i]

creat_upset_mat = function(data, nsets = ncol(data), order.by = "freq", group.by = "degree", main.bar.color = "gray23") {
	startend <- UpSetR:::FindStartEnd(data)
	first.col <- startend[1]
	last.col <- startend[2]
	Set_names <- UpSetR:::FindMostFreq(data, first.col, last.col, nsets)
	Sets_to_remove <- UpSetR:::Remove(data, first.col, last.col, Set_names)
	New_data <- UpSetR:::Wanted(data, Sets_to_remove)
	Num_of_set <- UpSetR:::Number_of_sets(Set_names)
	Set_names <- UpSetR:::order_sets(New_data, Set_names)
	All_Freqs <- UpSetR:::Counter(New_data, Num_of_set, first.col, 
								  Set_names, 40, main.bar.color, order.by,
								  group.by, NULL, NULL, c(T, F))
	Matrix_setup <- UpSetR:::Create_matrix(All_Freqs)
	Matrix_setup[nrow(Matrix_setup):1, ]
}

#l.len = length(listInput)
my.mat = creat_upset_mat(fromList(listInput), order.by = "freq", main.bar.color = my_cols2[3])
my.mat
stopifnot(nrow(my.mat)==length(upset.col))
mtxcol = data.frame(x=rep(1:ncol(my.mat), each=nrow(my.mat)),
						  color=rep(upset.col, times=ncol(my.mat)))
#mtxcol <- data.frame(x=rep(1:12,each=3), 
#					 color=rep(c("maroon","blue","orange"),each=12))
p2 = myupset(fromList(listInput), 
		   order.by = "freq", 
		   nsets=length(listInput), 
		   #matrix.color = "black",
		   main.bar.color = my_cols2[3],
		   sets.bar.color = rev(upset.col),
		   point.size=3,
		   mat_col=mtxcol,
		   line.size=0.5)
pdf(paste0("UpSetR_top500", ".pdf"), 6, 4.5, family="ArialMT")
print (p2)
dev.off()

###############
# only used for calculate Pvalue bonffonie/ FDR cutoff
dat.allct = c()
i = 1
pb <- txtProgressBar(min = 0, max = length(celltypes), style = 3)
for (i in 1:length(celltypes)) {
	ct = celltypes_paper_name[i]
	fold = ewas_fold[i]
	glm.f = paste0(fold, "/", glmf)
	if (file.exists(glm.f)) {
		ctewas = r.table(glm.f, T, 1)
	} else {
		next;
	}
	premht = data.frame(ctewas[, c("Chromosome", "Position", "Meta P")])
	dat.allct = rbind(dat.allct, premht)
	setTxtProgressBar(pb, i)
}
###############
head(dat.allct)
#bonferroni_adj = signif(p.adjust(dat.allct[, 3],method = "bonferroni"), 3)
#BH_fdr_adj = signif(p.adjust(dat.allct[, 3],method = "BH"), 3)
(p.cutoff = determine_bonf_fdr_cutoff(dat.allct, 3))

#to draw merged MHT
premht.allct = data.frame()
i = 1
pb <- txtProgressBar(min = 0, max = length(celltypes), style = 3)
for (i in 1:length(celltypes)) {
	ct = celltypes_paper_name[i]
	fold = ewas_fold[i]
	glm.f = paste0(fold,"/", glmf)
	if (file.exists(glm.f)) {
		ctewas = r.table(glm.f, T)
	} else {
		next;
	}

	#filter meta with hetero pv
	ctewas = ctewas[which(ctewas$"Hat P" > 0.05), ]
	#for fdr cutoff plot
	fdr.i = which(ctewas$BH_fdr_adj <= fdr.cutoff)
	if (length(fdr.i) > 0) {
		fdr.listInput[[ct]] = ctewas[fdr.i, 1]
	} else {
		fdr.listInput[[ct]] = NULL
	}

	premht = data.frame(ctewas[, c("probe", "Chromosome", "Position", "Meta P")])
	premht$hight.col = my_cols[i]
	premht$CT = ct
	h.i = which(premht[,4] <= p.cutoff[1])
	symbols = gene_symbol_trim(ctewas[h.i, "Nearest gene(s)"])
	#if (length(symbols) > 0) symbols = paste0(ct, "_", symbols)
	gene.i = which(symbols != "")
	premht[h.i[gene.i], 1] = symbols[gene.i]
	premht.allct = rbind(premht.allct, premht)
	setTxtProgressBar(pb, i)
}
close(pb)
str(premht.allct)
str(fdr.listInput)

fdr.listInput = fdr.listInput[!sapply(fdr.listInput, is.null)]
if (F) { #for UpSetR plot without Gran
	fdr.listInput = fdr.listInput[names(fdr.listInput) != "Gran"]
	names(fdr.listInput)
}
ct.n.i = match(names(fdr.listInput)[order(sapply(fdr.listInput, length))], c(celltypes_paper_name, "PBMC"))
upset.col = c(my_cols, my_cols_pbmc)[ct.n.i]
my.mat = creat_upset_mat(fromList(fdr.listInput), order.by = "freq", main.bar.color = my_cols2[3])
my.mat
stopifnot(nrow(my.mat)==length(upset.col))
mtxcol = data.frame(x=rep(1:ncol(my.mat), each=nrow(my.mat)),
                          color=rep(upset.col, times=ncol(my.mat)))

p1 = upset(fromList(fdr.listInput), 
		   order.by = "freq", 
		   nsets=length(fdr.listInput), 
		   #matrix.color = mtxcol,
		   matrix.color = "black",
		   main.bar.color = my_cols2[3],
		   sets.bar.color = rev(upset.col),
		   point.size=3,
		   line.size=0.5)
pdf(paste0("UpSetR_fdr", fdr.cutoff, ".pdf"), 6, 4.5, family="ArialMT")
print (p1)
dev.off()

pre_multiple_mht_premht <- function (premht.df, mht.pv.cutoff) {
	na.i = which(is.na(premht.df[,4]))
	if (length(na.i)>0) {
		premht.df = premht.df[-na.i, ]
	}
	premht.df = premht.df[which(premht.df$Chromosome != "" & premht.df$Chromosome != "X" & premht.df$Chromosome != "Y"),]
	premht.df <- premht.df[order(as.numeric(premht.df$Chromosome), as.numeric(premht.df$Position), premht.df$CT),]
	h.i = which(premht.df[,4] <= mht.pv.cutoff)
	cg.i = grep("^cg", premht.df[h.i, 1])
	#if (length(cg.i)>0) h.i = h.i[-cg.i]
	if (length(cg.i)>0) {
		premht.df[h.i[cg.i], 1] = "SIGGENES"
	}
	highlight = premht.df[h.i, 1]
	#to change gene name 
	allg.i = which(premht.df$probe %in% highlight)
	if(length(allg.i) > length(h.i)) {
		diff.i = setdiff(allg.i, h.i)
		premht.df = premht.df[-diff.i, ]
	}
	hight.col = premht.df$hight.col[h.i]
	return (list(premht.df = premht.df, highlight=highlight, hight.col=hight.col))
}

premht.out = pre_multiple_mht_premht(premht.allct,  p.cutoff[1])
premht.allct2 = premht.out[["premht.df"]]
highlight = premht.out[["highlight"]]
hight.col = premht.out[["hight.col"]]
l.i = match(names(table(hight.col)), my_cols)
#legend.name = celltypes[l.i]
#legend.cols  = my_cols[l.i]
(legend.name = celltypes_paper_name)
(legend.cols  = my_cols)
cutoff = -log10(p.cutoff)
dataF = premht.allct2[, 1:4]
pic_name = paste0("meta_CellType_mhtplot_paper.png")
qq_pic_name = paste0("meta_CellType_qqplot_paper.png")
annotatePval = p.cutoff[1]
#annotatePval = NULL
cex.label = 3
cex.text = 1.5

premht.allct2[which(premht.allct2$probe %in% highlight), ]
highlight
hight.col
#d[which(d$SNP %in% highlight), ]
legend.name
legend.cols
legend_order = NULL
if (legend.po == "top") legend_order=matrix(1:6, ncol=2, byrow=T)
source("~/mybiotools/r/gwas_mht_qq_plot_multiple.R")
MHT_QQ_plot_multy(dataF = dataF, 
				  pic_name = pic_name, 
				  qq_pic_name = qq_pic_name, 
				  highlight=highlight, 
				  cutoff = cutoff,
				  hight.col = hight.col,
				  legend.name = legend.name,
				  legend.cols = legend.cols,
				  legend.po = legend.po,
				  legend.cex = legend.cex,
				  cex.label = cex.label,
				  cex.text = cex.text,
				  legend_order = legend_order,
				  wordcloud = T,
				  plotrix = F,
				  textxy = F,
				  annotatePval=annotatePval
				  )

# to draw PBMC
if (!is.null(premht.PBMC)) {
	str(premht.PBMC)
	premht.out = pre_multiple_mht_premht(premht.PBMC,  pbmc.p.cutoff[1])
	premht.PBMC2 = premht.out[["premht.df"]]
	highlight = premht.out[["highlight"]]
	hight.col = premht.out[["hight.col"]]
	legend.name = "PBMC"
	legend.cols  = my_cols_pbmc
	dataF = premht.PBMC2[, 1:4]
	pic_name = paste0("mht_PBMC_v3_paper.png")
	qq_pic_name = paste0("qq_PBMC_v3_paper.png")

	source("~/mybiotools/r/gwas_mht_qq_plot_multiple.R")
	MHT_QQ_plot_multy(dataF = dataF,
					  pic_name = pic_name, 
					  qq_pic_name = qq_pic_name, 
					  highlight=highlight,
					  cutoff = -log10(pbmc.p.cutoff),
					  hight.col = hight.col,
					  legend.name = legend.name,
					  legend.cols = legend.cols,
					  legend.po = legend.po,
					  legend.cex = legend.cex,
					  cex.label = cex.label,
					  cex.text = cex.text,
					  wordcloud = T,
					  #show.lines = T,
					  plotrix = F,
					  textxy = F,
					  annotatePval=pbmc.p.cutoff[1]
					  )

}


