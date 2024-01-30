#!/usr/bin/env Rscript

rm(list=ls());
library(myfunc);
source("~/mybiotools/r/myfunc.R");
options("echo"=T, "scipen"=1, "digits"=4, stringsAsFactors=FALSE);
args<-commandArgs(T)

topnum = 500
cohort.num = 3

metafiles = list.files(".", pattern="EWAS_Meta_analysis_sampleSize_adj.xls", recursive=T, full.names=T)
celltypes = do.call(rbind, strsplit(metafiles, '/'))[,2]

path.cohort = c("/gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/hiv_450k_comm_allVACS",
			   "/gpfs/ycga/project/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/brad_hiv_comm_residV2",
			   "/gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/GSE217633/TCA")

phe.list = list()

############################
#for cohort 1 phe
### read data
phefile = "/home/xz345/phenotype/10212015_new_cross_link/allmethy.phe"
phe1 <- read.delim(phefile, header=T, row.names=NULL, check.names =F, sep = "\t");

phefile_add = "/gpfs/ycga/project/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/hiv_450k_comm_allVACS/residual/residualPCs_noHIVnoVL.txt";
phe_add <- read.delim(phefile_add, header=T, row.names=NULL, check.names =F, sep = "\t");
#emp.i = which(phe_add[,"hiv"] == "");
#if (length(emp.i)>0) phe_add[emp.i,"hiv"] = NA
phe_add[phe_add==""] = NA
phe_add = phe_add[!is.na(phe_add[,"Meth450Kmicroarrayid"]),];
dup.i = which(duplicated(phe_add[, "Meth450Kmicroarrayid"]));
if (length(dup.i>0)) phe_add = phe_add[-dup.i,];
phe1.2 = unique(merge(phe1, phe_add, by="Meth450Kmicroarrayid", all.x = T, all.y = F, sort = FALSE));
phe1 = merge_all_col(phe1.2)

#########################
#added 11/05/2015 to avoid na/duplicated id in row col
phe1 = phe1[!is.na(phe1[,"Meth450Kmicroarrayid"]),]
phe1 = phe1[!is.na(phe1[,"hiv"]),]
dup.i = which(duplicated(phe1[, "Meth450Kmicroarrayid"]));
if (length(dup.i>0)) phe1 = phe1[-dup.i,]
#########################
rownames(phe1) = phe1[, "Meth450Kmicroarrayid"]
phe1 = phe1[which(phe1[,"dmgsex"]==1),]
phe1$HIV = phe1$hiv
phe.list[[1]] = phe1

############################
#for cohort 2 phe
### read data
phefile = "/gpfs/ycga/project/xu_ke/xl535/other/Brad_022020_Methy/phenotype/Brad_phe_v3.txt"
phe1 <- read.delim(phefile, header=T, row.names=NULL, check.names =F, sep = "\t");

phefile_add = "/gpfs/ycga/project/xu_ke/xl535/other/Brad_022020_Methy/PCA/PositiveControlPCA/cont_PC_Brad.txt";
phe_add <- read.delim(phefile_add, header=T, row.names=NULL, check.names =F, sep = "\t");
#emp.i = which(phe_add[,"HIV"] == "");
#if (length(emp.i)>0) phe_add[emp.i,"HIV"] = NA
phe_add[phe_add==""] = NA
phe_add = phe_add[!is.na(phe_add[,"microarrayID"]),];
dup.i = which(duplicated(phe_add[, "microarrayID"]));
if (length(dup.i>0)) phe_add = phe_add[-dup.i,];
phe1.2 = unique(merge(phe1, phe_add, by="microarrayID", all.x = T, all.y = F, sort = FALSE));
phe1 = merge_all_col(phe1.2)

phefile_add = "/gpfs/ycga/project/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/brad_hiv_comm_residV2/residual/residualPC_noVL.txt";
phe_add <- read.delim(phefile_add, header=T, row.names=NULL, check.names =F, sep = "\t");
#emp.i = which(phe_add[,"HIV"] == "");
#if (length(emp.i)>0) phe_add[emp.i,"HIV"] = NA
phe_add[phe_add==""] = NA
phe_add = phe_add[!is.na(phe_add[,"microarrayID"]),];
dup.i = which(duplicated(phe_add[, "microarrayID"]));
if (length(dup.i>0)) phe_add = phe_add[-dup.i,];
phe1.2 = unique(merge(phe1, phe_add, by="microarrayID", all.x = T, all.y = F, sort = FALSE));
phe1 = merge_all_col(phe1.2)

#########################
#added 11/05/2015 to avoid na/duplicated id in row col
phe1 = phe1[!is.na(phe1[,"microarrayID"]),]
dup.i = which(duplicated(phe1[, "microarrayID"]));
if (length(dup.i>0)) phe1 = phe1[-dup.i,]
rownames(phe1) = phe1[, "microarrayID"]
phe1$HIV=phe1$STATUS
phe1$HIV[phe1$HIV=="HIVpos"]="1"
phe1$HIV[phe1$HIV=="HIVneg"]="0"
phe1$HIV=as.factor(phe1$HIV)
phe.list[[2]] = phe1

############################
#for cohort 3 phe
### read data
phefile = "/gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/GSE217633/geo_phenotype.txt"
trueResp = "HIV";
phe1 <- read.delim(phefile, header=T, row.names=NULL, check.names =F, sep = "\t");

phefile_add = "/gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/GSE217633/PCA/PositiveControlPCA/cont_PC_850k.txt";
phe_add <- read.delim(phefile_add, header=T, row.names=NULL, check.names =F, sep = "\t");
#emp.i = which(phe_add[,"HIV"] == "");
#if (length(emp.i)>0) phe_add[emp.i,"HIV"] = NA
phe_add[phe_add==""] = NA
phe_add = phe_add[!is.na(phe_add[,"id"]),];
dup.i = which(duplicated(phe_add[, "id"]));
if (length(dup.i>0)) phe_add = phe_add[-dup.i,];
phe1.2 = unique(merge(phe1, phe_add, by.x="id", by.y="id", all.x = T, all.y = F, sort = FALSE));
phe1 = merge_all_col(phe1.2)

phefile_add = "/gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/GSE217633/PCA/residual_C30contPCA/resiPCA.txt";
phe_add <- read.delim(phefile_add, header=T, row.names=NULL, check.names =F, sep = "\t");
#emp.i = which(phe_add[,"HIV"] == "");
#if (length(emp.i)>0) phe_add[emp.i,"HIV"] = NA
phe_add[phe_add==""] = NA
phe_add = phe_add[!is.na(phe_add[,"id"]),];
dup.i = which(duplicated(phe_add[, "id"]));
if (length(dup.i>0)) phe_add = phe_add[-dup.i,];
phe1.2 = unique(merge(phe1, phe_add, by.x="id", by.y="id", all.x = T, all.y = F, sort = FALSE));
phe1 = merge_all_col(phe1.2)

phefile_add = "/gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/GSE217633/TCA/PBMC_beta_ct.txt";
phe_add <- read.delim(phefile_add, header=T, row.names=NULL, check.names =F, sep = "\t");
#emp.i = which(phe_add[,"HIV"] == "");
#if (length(emp.i)>0) phe_add[emp.i,"HIV"] = NA
phe_add[phe_add==""] = NA
phe_add[phe_add==""] = NA
phe_add = phe_add[!is.na(phe_add[,"id"]),];
dup.i = which(duplicated(phe_add[, "id"]));
if (length(dup.i>0)) phe_add = phe_add[-dup.i,];
phe1.2 = unique(merge(phe1, phe_add, by.x="id", by.y="id", all.x = T, all.y = F, sort = FALSE));
phe1 = merge_all_col(phe1.2)

#########################
#added 11/05/2015 to avoid na/duplicated id in row col
phe1 = phe1[!is.na(phe1[,"id"]),]
dup.i = which(duplicated(phe1[, "id"]));
if (length(dup.i>0)) phe1 = phe1[-dup.i,]
#########################
rownames(phe1) = phe1[, "id"]
phe1$HIV = phe1$time
phe1 = phe1[-which(phe1$HIV == "sample post ART"),]
phe1$HIV[phe1$HIV=="sample pre ART"]="1"
phe1$HIV[phe1$HIV=="control sample"]="0"
phe1$HIV=as.factor(phe1$HIV)
phe.list[[3]] = phe1

lapply(phe.list, function(x) table(x[, "HIV"]))

alldata.l = list()
i = 1
j = 1
for (i in 1:length(celltypes))
{
	ct = celltypes[i]
	metadat = r.table(metafiles[i])
	topprobes = metadat$probe[1:topnum]
	for (j in 1:cohort.num) {
		path.c = path.cohort[j]
		phe.c = phe.list[[j]]
		pv.cohort = metadat[1:topnum, paste0("Cohort ", j," P")]
		rawdatafile = ifelse (ct == "PBMC", "beta0109_X.txt", paste0(ct, "_tcaMethy.txt"))
		dat1 = r.table(paste0(path.c, "/", rawdatafile), T, 1)
		phe1 = phe.c
		comm = colnames(dat1)[which(colnames(dat1) %in% rownames(phe1))]
		phe1=phe1[comm,]
		dat1=dat1[topprobes, comm]
		states = phe1$HIV
		states1 = as.factor(paste0("Cohort", j, "_", ct, "_HIV", phe1$HIV))
		for (k in 1:topnum) {
			cpg = topprobes[k]
			boxplot_df1 = data.frame(value=as.numeric(dat1[cpg, ]), 
									group=states1)
			pval.anno1  = data.frame(group1=levels(states1)[1],
									 group2=levels(states1)[2],
									 pval=pv.cohort[k])
			alldata.l[[ct]][[cpg]][["boxplot_df"]] = rbind(alldata.l[[ct]][[cpg]][["boxplot_df"]], boxplot_df1)
			alldata.l[[ct]][[cpg]][["pval.anno"]] = rbind(alldata.l[[ct]][[cpg]][["pval.anno"]], pval.anno1)
		}

	}
}

source("~/mybiotools/r/myfunc.R");
outfolder = "sigCpG_BoxPlot"
if (!dir.exists(outfolder)) dir.create(outfolder)
for (i in 1:length(alldata.l))
{
	ct = names(alldata.l)[i]
	message (paste0("Processing ", ct))
	ct.d = alldata.l[[ct]]
	pdf(paste0(outfolder, "/", ct, "_top", topnum, "_sigCpG_BoxPlot.pdf"), 5, 5, family="ArialMT")
	for (j in 1:length(ct.d))
	{
		cpg = names(ct.d)[j]
		myv = ct.d[[cpg]][["boxplot_df"]]$value 
		myv[myv<0] = NA
		myv[myv>1] = NA
		ct.d[[cpg]][["boxplot_df"]]$value = myv
		
		myg = ct.d[[cpg]][["boxplot_df"]]$group
		mean.g.v = tapply(myv, myg, mean, na.rm=T)
		mean.g.v.ratio = c(abs(mean.g.v[1] - mean.g.v[2])/mean.g.v[1], 
						   abs(mean.g.v[3] - mean.g.v[4])/mean.g.v[3],
						   abs(mean.g.v[5] - mean.g.v[6])/mean.g.v[5])
		if (sum(mean.g.v.ratio >= 0.05, na.rm=T) > 1) {
			p = gg.scatter.plot.rainbow(ct.d[[cpg]][["boxplot_df"]], ct.d[[cpg]][["pval.anno"]], 
										color.mode = "paired",
										title = cpg, angle=90)
			print (p)
		}
	}
	dev.off()
}



