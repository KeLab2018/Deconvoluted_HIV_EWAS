#!/usr/bin/env Rscript

rm(list=ls());
library(myfunc);
source("~/mybiotools/r/myfunc.R");
options("scipen"=1, "digits"=4, stringsAsFactors=FALSE);
args<-commandArgs(T)

library(reshape2)

celltypes_paper_name = c("CD4", "CD8", "B", "NK", "M", "G")
celltypes = c("CD4T", "CD8T", "Bcell", "NK", "Mono", "Gran")


########################################### VSCS ################################################

arrayid = "Meth450Kmicroarrayid"
dataf = "/gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/hiv_450k_comm_allVACS/splitfile/beta0109_X_000000.txt"

X = r.table(dataf, T, 1)
############for phenotype################
phefile = "/gpfs/ycga/project/xu_ke/xz345/work/VACS/phenotype/10212015_new_cross_link/allmethy.phe"
phe1 <- read.delim(phefile, header=T, row.names=NULL, check.names =F, sep = "\t");

phe1 = phe1[!is.na(phe1[,arrayid]),]
phe1 = phe1[!is.na(phe1[,"hiv"]),]
phe1 = phe1[which(phe1[,"dmgsex"]==1),]

dup.i = which(duplicated(phe1[, arrayid]));
if (length(dup.i>0)) phe1 = phe1[-dup.i,]

###### for 450K ######
rownames(phe1) = phe1[, arrayid]
#

dat1 = X
comm = colnames(dat1)[which(colnames(dat1) %in% rownames(phe1))]
phe=phe1[comm,]
dat=dat1[,comm]

X = dat[,rownames(phe)]
dim(dat)
dim(phe)
stopifnot(all(colnames(dat) == rownames(phe)))

#attach(phe)
model.data = data.frame(probe=rownames(phe));
model.data[,"AGEBL"] <- as.numeric(phe[,"AGEBL"])
model.data[,"WBC_new"] <- as.numeric(phe[,"WBC_new"])
model.data[,"CD8T"] <- as.numeric(phe[,"CD8T"])
model.data[,"CD4T"] <- as.numeric(phe[,"CD4T"])
model.data[,"Gran"] <- as.numeric(phe[,"Gran"])
model.data[,"NK"] <- as.numeric(phe[,"NK"])
model.data[,"Bcell"] <- as.numeric(phe[,"Bcell"])
model.data[,"Mono"] <- as.numeric(phe[,"Mono"])
model.data[,"hiv"]  <- as.factor(phe[,"hiv"])
model.data[,"smkcigs"]  <- as.factor(phe[,"smkcigs"])
model.data[,"RACECOMG"]  <- as.factor(phe[,"RACECOMG"])

##### phenotype #####
W = as.matrix(model.data[, 4:9])
W[which(W<0)] = 0
sum1 = function(x) {
	x.sum = sum(x,na.rm=T)
	newx = c()
	for (i in x) {
		i = i * 1 / x.sum
		newx = c(newx, i)
	}
	newx
}

if (F) {
	W2 = apply(W, 1, sum1)
	W2 = t(W2)
	colnames(W2) = colnames(W)
	W = W2
	apply(W, 1, sum)
	rownames(W) = model.data[, "probe"]
	#rm Gran
	W2 = W[, -3]

	#### plot cell type percentage jitter_plot
	pdffile = "hiv_celltypes_jitterplot.pdf"
	cols = 3
	plot_celltype(W, pdffile, cols = cols)
	pdffile = "hiv_5celltypes_jitterplot.pdf"
	plot_celltype(W2, pdffile, cols = cols)
	w.table(W, file="hiv_celltypes_jitterplot.xls", T, T, "sample")
}

head(W)
W = W[, celltypes]
head(W)
colnames(W) = celltypes_paper_name
head(W)
W2 = melt(W)
states.ct.vacs = as.factor(W2[, "Var2"])
ct.vacs = data.matrix(rbind(VACS=W2[,"value"]))
plot_celltype(W, pdffile = "vacs_celltypes_jitterplot_6ct.pdf", cols = 6)


########################################### WIHS ################################################

###### SETUP ######
arrayid = "microarrayID"
###### SETUP ######


dataf = "/gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/brad_hiv_comm_residV2/splitfile/beta0109_X_000000.txt"
X = r.table(dataf, T, 1)

############for phenotype################
phefile = "/gpfs/ycga/project/xu_ke/xl535/other/Brad_022020_Methy/phenotype/Brad_phe_v3.txt"
#phefile = "/gpfs/ycga/project/xu_ke/xl535/other/Brad_022020_Methy/phenotype/Brad_phe_v2_NDRNKWK.Average.txt"
phe1 <- read.delim(phefile, header=T, row.names=NULL, check.names =F, sep = "\t");

phefile_add = "/gpfs/ycga/project/xu_ke/xl535/other/Brad_022020_Methy/PCA/PositiveControlPCA/cont_PC_Brad.txt";
phe_add <- read.delim(phefile_add, header=T, row.names=NULL, check.names =F, sep = "\t");
phe_add[phe_add==""] = NA
phe_add = phe_add[!is.na(phe_add[,"microarrayID"]),];
dup.i = which(duplicated(phe_add[, "microarrayID"]));
if (length(dup.i>0)) phe_add = phe_add[-dup.i,];
phe1.2 = unique(merge(phe1, phe_add, by="microarrayID", all.x = T, all.y = F, sort = FALSE));
phe1 = merge_all_col(phe1.2)

phefile_add = "/gpfs/ycga/project/xu_ke/xl535/other/Brad_022020_Methy/Probe450K/smoker_HIV/residuals/residualPCs.txt";
phe_add <- read.delim(phefile_add, header=T, row.names=NULL, check.names =F, sep = "\t");
phe_add[phe_add==""] = NA
phe_add = phe_add[!is.na(phe_add[,"microarrayID"]),];
dup.i = which(duplicated(phe_add[, "microarrayID"]));
if (length(dup.i>0)) phe_add = phe_add[-dup.i,];
phe1.2 = unique(merge(phe1, phe_add, by="microarrayID", all.x = T, all.y = F, sort = FALSE));
phe1 = merge_all_col(phe1.2)

phe1 = phe1[!is.na(phe1[,arrayid]),]
phe1 = phe1[!is.na(phe1[,"STATUS"]),]

dup.i = which(duplicated(phe1[, arrayid]));
if (length(dup.i>0)) phe1 = phe1[-dup.i,]

rownames(phe1) = phe1[, arrayid]

dat1 = X
comm = colnames(dat1)[which(colnames(dat1) %in% rownames(phe1))]
phe1=phe1[comm,]
dat1=dat1[,comm]

# consider all individuals

phe1$race=phe1$RACECAT
phe1$race[grepl(pattern="African-American",phe1$race)]=1
phe1$race[grepl(pattern="White",phe1$race)]=2
phe1$race=as.factor(phe1$race)

phe1$CUMUSE.binomial=phe1$COMBINED_SUBSTANCE_CUMULATIVE_USE
phe1$CUMUSE.binomial[which(phe1$CUMUSE.binomial!="NEVER USE")]="1"
phe1$CUMUSE.binomial[which(phe1$CUMUSE.binomial=="NEVER USE")]="0"
phe1$CUMUSE.binomial=as.factor(phe1$CUMUSE.binomial)

phe1$smoker=phe1$SMOKING
phe1$smoker[phe1$smoker=="Current smoker"]="1"
phe1$smoker[phe1$smoker=="Former smoker"]="1"
phe1$smoker[phe1$smoker=="Never Smoker"]="0"

phe1$drinker=phe1$NDRNKWK
phe1$drinker[phe1$drinker=="HIGH_RISK"]="1"
phe1$drinker[phe1$drinker=="LOW_RISK"]="1"
phe1$drinker[phe1$drinker=="NONDRINKER"]="0"

phe1$HIV=phe1$STATUS
phe1$HIV[phe1$HIV=="HIVpos"]="1"
phe1$HIV[phe1$HIV=="HIVneg"]="0"

phe1$log10VL=log10(phe1$VLOAD)

phe=phe1
dat=dat1[,rownames(phe)]

X = dat
dim(dat)
dim(phe)
stopifnot(all(colnames(dat) == rownames(phe)))


#attach(phe)
model.data = data.frame(probe=rownames(phe));
model.data[,"AGEATVIS"] <- as.numeric(phe[,"AGEATVIS"])
model.data[,"CD8T"] <- as.numeric(phe[,"CD8T"])
model.data[,"CD4T"] <- as.numeric(phe[,"CD4T"])
model.data[,"Gran"] <- as.numeric(phe[,"Gran"])
model.data[,"NK"] <- as.numeric(phe[,"NK"])
model.data[,"Bcell"] <- as.numeric(phe[,"Bcell"])
model.data[,"Mono"] <- as.numeric(phe[,"Mono"])
model.data[,"Cont_Pr_Brad_PC1"] <- as.numeric(phe[,"Cont_Pr_Brad_PC1"])
model.data[,"Cont_Pr_Brad_PC2"] <- as.numeric(phe[,"Cont_Pr_Brad_PC2"])
model.data[,"Cont_Pr_Brad_PC3"] <- as.numeric(phe[,"Cont_Pr_Brad_PC3"])
model.data[,"Cont_Pr_Brad_PC4"] <- as.numeric(phe[,"Cont_Pr_Brad_PC4"])
model.data[,"Cont_Pr_Brad_PC5"] <- as.numeric(phe[,"Cont_Pr_Brad_PC5"])
model.data[,"Cont_Pr_Brad_PC6"] <- as.numeric(phe[,"Cont_Pr_Brad_PC6"])
model.data[,"Cont_Pr_Brad_PC7"] <- as.numeric(phe[,"Cont_Pr_Brad_PC7"])
model.data[,"Cont_Pr_Brad_PC8"] <- as.numeric(phe[,"Cont_Pr_Brad_PC8"])
model.data[,"Cont_Pr_Brad_PC9"] <- as.numeric(phe[,"Cont_Pr_Brad_PC9"])
model.data[,"Cont_Pr_Brad_PC10"] <- as.numeric(phe[,"Cont_Pr_Brad_PC10"])
model.data[,"Cont_Pr_Brad_PC11"] <- as.numeric(phe[,"Cont_Pr_Brad_PC11"])
model.data[,"Cont_Pr_Brad_PC12"] <- as.numeric(phe[,"Cont_Pr_Brad_PC12"])
model.data[,"Cont_Pr_Brad_PC13"] <- as.numeric(phe[,"Cont_Pr_Brad_PC13"])
model.data[,"Cont_Pr_Brad_PC14"] <- as.numeric(phe[,"Cont_Pr_Brad_PC14"])
model.data[,"Cont_Pr_Brad_PC15"] <- as.numeric(phe[,"Cont_Pr_Brad_PC15"])
model.data[,"Cont_Pr_Brad_PC16"] <- as.numeric(phe[,"Cont_Pr_Brad_PC16"])
model.data[,"Cont_Pr_Brad_PC17"] <- as.numeric(phe[,"Cont_Pr_Brad_PC17"])
model.data[,"Cont_Pr_Brad_PC18"] <- as.numeric(phe[,"Cont_Pr_Brad_PC18"])
model.data[,"Cont_Pr_Brad_PC19"] <- as.numeric(phe[,"Cont_Pr_Brad_PC19"])
model.data[,"Cont_Pr_Brad_PC20"] <- as.numeric(phe[,"Cont_Pr_Brad_PC20"])
model.data[,"Cont_Pr_Brad_PC21"] <- as.numeric(phe[,"Cont_Pr_Brad_PC21"])
model.data[,"Cont_Pr_Brad_PC22"] <- as.numeric(phe[,"Cont_Pr_Brad_PC22"])
model.data[,"Cont_Pr_Brad_PC23"] <- as.numeric(phe[,"Cont_Pr_Brad_PC23"])
model.data[,"Cont_Pr_Brad_PC24"] <- as.numeric(phe[,"Cont_Pr_Brad_PC24"])
model.data[,"Cont_Pr_Brad_PC25"] <- as.numeric(phe[,"Cont_Pr_Brad_PC25"])
model.data[,"Cont_Pr_Brad_PC26"] <- as.numeric(phe[,"Cont_Pr_Brad_PC26"])
model.data[,"Cont_Pr_Brad_PC27"] <- as.numeric(phe[,"Cont_Pr_Brad_PC27"])
model.data[,"Cont_Pr_Brad_PC28"] <- as.numeric(phe[,"Cont_Pr_Brad_PC28"])
model.data[,"Cont_Pr_Brad_PC29"] <- as.numeric(phe[,"Cont_Pr_Brad_PC29"])
model.data[,"Cont_Pr_Brad_PC30"] <- as.numeric(phe[,"Cont_Pr_Brad_PC30"])
model.data[,"PC1"] <- as.numeric(phe[,"PC1"])
model.data[,"PC2"] <- as.numeric(phe[,"PC2"])
model.data[,"PC3"] <- as.numeric(phe[,"PC3"])
model.data[,"PC4"] <- as.numeric(phe[,"PC4"])
model.data[,"PC5"] <- as.numeric(phe[,"PC5"])
model.data[,"PC6"] <- as.numeric(phe[,"PC6"])
model.data[,"PC7"] <- as.numeric(phe[,"PC7"])
model.data[,"PC8"] <- as.numeric(phe[,"PC8"])
model.data[,"PC9"] <- as.numeric(phe[,"PC9"])
model.data[,"PC10"] <- as.numeric(phe[,"PC10"])
model.data[,"PC11"] <- as.numeric(phe[,"PC11"])
model.data[,"PC12"] <- as.numeric(phe[,"PC12"])
model.data[,"PC13"] <- as.numeric(phe[,"PC13"])
model.data[,"PC14"] <- as.numeric(phe[,"PC14"])
model.data[,"PC15"] <- as.numeric(phe[,"PC15"])
model.data[,"PC16"] <- as.numeric(phe[,"PC16"])
model.data[,"PC17"] <- as.numeric(phe[,"PC17"])
model.data[,"PC18"] <- as.numeric(phe[,"PC18"])
model.data[,"PC19"] <- as.numeric(phe[,"PC19"])
model.data[,"PC20"] <- as.numeric(phe[,"PC20"])
model.data[,"HIV"]  <- as.factor(phe[,"HIV"])
model.data[,"drinker"]  <- as.factor(phe[,"drinker"])
model.data[,"race"]  <- as.factor(phe[,"race"])
model.data[,"smoker"]  <- as.factor(phe[,"smoker"])
names(model.data)

ct.i = 3:8 # cell types
W = as.matrix(model.data[, ct.i])
W[which(W<0)] = 0
sum1 = function(x) {
	x.sum = sum(x,na.rm=T)
	newx = c()
	for (i in x) {
		i = i * 1 / x.sum
		newx = c(newx, i)
	}
	newx
}
W2 = apply(W, 1, sum1)
W2 = t(W2)
colnames(W2) = colnames(W)
W = W2
apply(W, 1, sum)
sum(is.na(W))
rownames(W) = phe[, arrayid]

head(W)
W = W[, celltypes]
head(W)
colnames(W) = celltypes_paper_name
head(W)
W2 = melt(W)
states.ct.wihs = as.factor(W2[, "Var2"])
ct.wihs = data.matrix(rbind(VACS=W2[,"value"]))
plot_celltype(W, pdffile = "wihs_celltypes_jitterplot_6ct.pdf", cols = 6)


########################################### GSE ################################################

W = r.table("/gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/GSE217633/TCA/PBMC_beta_ct.txt", T, 1)
W = as.matrix(W)

head(W)
W = W[, celltypes]
head(W)
colnames(W) = celltypes_paper_name
head(W)
W2 = melt(W)
states.ct.gse = as.factor(W2[, "Var2"])
ct.gse = data.matrix(rbind(VACS=W2[,"value"]))


plot_celltype(W, pdffile = "GSE_celltypes_jitterplot_6ct.pdf", cols = 6)
