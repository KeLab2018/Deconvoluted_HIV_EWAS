#!/usr/bin/env Rscript
# to draw effect size plot for each ct
rm(list=ls());
source("~/mybiotools/r/myfunc.R");
options("scipen"=1, "digits"=5, stringsAsFactors=FALSE);
args<-commandArgs(T)
outkw = args[1] # CT methylome ewas folds key word
root.d1 = args[2]
root.d2 = args[3]
root.d3 = args[4]
glmf = args[6]
ifMC = as.logical(args[7])

print (args)

if (is.na(outkw)) outkw = "_methylome_hiv_ewas"
if (is.na(root.d1)) root.d1 = "/gpfs/ycga/project/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/hiv_450k_comm_allVACS/CTmethylome_6CT"
if (is.na(root.d2)) root.d2 = "/gpfs/ycga/project/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/brad_hiv_comm_residV2/CTmethylome"
if (is.na(root.d3)) root.d3 = "/gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/GSE217633/TCA/hiv_ewas_preART"
#if (is.na(glmf)) glmf = "top5000.xls"
if (is.na(glmf)) glmf = "glm_pv_adj_anno.xls"
if (is.na(ifMC)) ifMC = F

comp.value = "t.value"

root.meta = "/gpfs/ycga/work/xu_ke/xz345/work/FlowSorted_methy/FlowSorted450k_Public_deconvolution_TCA/hiv_meta_3datasets"
metaf = "EWAS_Meta_analysis_sampleSize_adj.xls"


#library(ggplot2)
library(lattice)
library("plot3D")
library("car")
library(RColorBrewer)
library("scatterplot3d")

pv.cutoff= 0.05

kw = "^cg"
if (ifMC) {
	kw = "|_|"
}


#celltypes = c("CD4T", "CD8T", "Mono", "Bcell", "NK", "Gran")
#if (!ifGran) celltypes = c("CD4T", "CD8T", "Mono", "Bcell", "NK")

#(ewas_fold = paste0(celltypes, outkw))
(ewas_fold = list.files(root.d1, pattern= paste0(outkw, '$'), all.files = F, full.names = T))
(celltypes = gsub (outkw, "", ewas_fold))
(celltypes = gsub (paste0(root.d1, "/"), "", celltypes))
c.i = which(!grepl("\\.|logs", celltypes))
(celltypes = celltypes[c.i])
(ewas_fold = ewas_fold[c.i])
(ewas_fold2 = list.files(root.d2, pattern= paste0(outkw, '$'), all.files = F, full.names = T))
(ewas_fold2 = ewas_fold2[c.i])
(ewas_fold3 = list.files(root.d3, pattern= paste0(outkw, '$'), all.files = F, full.names = T))
(ewas_fold3 = ewas_fold3[c.i])

(celltypes = celltypes[which(celltypes != "Gran")])

#celltypes_paper = c("CD4T", "CD8T", "Bcell", "NK", "Mono")
#celltypes.i = match(celltypes_paper, celltypes)
#(celltypes = celltypes[celltypes.i])
#(ewas_fold = ewas_fold[celltypes.i])

#renew
my_cols = add.alpha(c("black", "red", "orange", "green", "purple", "blue"), alpha = 0.5)
celltypes = c("PBMC", "CD4T", "CD8T", "Bcell", "NK", "Mono")
celltypes_paper_name = c("PBMC", "CD4", "CD8", "B", "NK", "M")
(ewas_fold = paste0(root.d1, "/", celltypes, outkw))
(ewas_fold2 = paste0(root.d2, "/", celltypes, outkw))
(ewas_fold3 = paste0(root.d3, "/", celltypes, outkw))
(meta_fold = paste0(root.meta, "/", celltypes_paper_name))

listInput = vector("list", length=(length(celltypes)))
names(listInput) = c(celltypes_paper_name)
premht.allct = data.frame()
sig.probe = list()

require(plyr)
glm2mc = function(glmdata) {
	# transfer MC header name to glm header name
	glmout = rename(glmdata, c("Chr"="CHR", 
							   "Start"="MAPINFO", 
							   "Gene Name"="UCSC_REFGENE_NAME", 
							   "PeakID"="probe" 
							   ))
	return (glmout)
}

my.names = c("probe", "t.value", "Pvalue", "Coefficients", "bonferroni_adj", "BH_fdr_adj", "BY_adj", "CHR", "MAPINFO", "UCSC_REFGENE_NAME", "UCSC_REFGENE_GROUP", "RELATION_TO_UCSC_CPG_ISLAND", "CT")

# meta
premht.allct.meta = data.frame()
i = 1
pb <- txtProgressBar(min = 0, max = length(celltypes), style = 3)
for (i in 1:length(celltypes)) {
	ct = celltypes_paper_name[i]
	fold = meta_fold[i]
	#glm.f = paste0(fold,"/glm_pv_adj_anno.xls")
    glm.f = paste0(fold,"/", metaf)
    if (file.exists(glm.f)) {
        ctewas = r.table(glm.f, T)
        #sig.i = which(ctewas$BH_fdr_adj <= 0.05)
        sig.i = which(ctewas$BH_fdr_adj <= 1) #v3
        if (length(sig.i) > 0) {
            ctewas = ctewas[sig.i, ]
        } else {
            next
        }
        if (ifMC) ctewas = glm2mc(ctewas)
        ctewas$CT = ct
    } else {
        next;
    }
    premht.allct.meta = rbind(premht.allct.meta, ctewas[,])
    sig.probe[[length(sig.probe)+1]] = ctewas$probe
    setTxtProgressBar(pb, i)
}
close(pb)
str(premht.allct.meta)

# fold vacs
i = 1
pb <- txtProgressBar(min = 0, max = length(celltypes), style = 3)
for (i in 1:length(celltypes)) {
	ct = celltypes_paper_name[i]
	fold = ewas_fold[i]
	#glm.f = paste0(fold,"/glm_pv_adj_anno.xls")
    glm.f = paste0(fold,"/", glmf)
    if (file.exists(glm.f)) {
        ctewas = r.table(glm.f, T)
        #sig.i = which(ctewas$BH_fdr_adj <= 0.05)
        sig.i = which(ctewas$Pvalue <= pv.cutoff) #v3
        if (length(sig.i) > 0) {
            ctewas = ctewas[sig.i, ]
        } else {
            next
        }
        if (ifMC) ctewas = glm2mc(ctewas)
        ctewas$CT = ct
    } else {
        next;
    }
    premht.allct = rbind(premht.allct, ctewas[, my.names])
    sig.probe[[length(sig.probe)+1]] = ctewas$probe
    setTxtProgressBar(pb, i)
}
close(pb)
str(premht.allct)
names(sig.probe) = celltypes_paper_name
str(sig.probe)

# fold wihs
premht.allct2 = data.frame()
i = 1
pb <- txtProgressBar(min = 0, max = length(celltypes), style = 3)
for (i in 1:length(celltypes)) {
	ct = celltypes_paper_name[i]
	fold = ewas_fold2[i]
	#glm.f = paste0(fold,"/glm_pv_adj_anno.xls")
	glm.f = paste0(fold,"/", glmf)
	if (file.exists(glm.f)) {
		ctewas = r.table(glm.f, T)
		#sig.i = which(ctewas$BH_fdr_adj <= 0.05)
		sig1.i = match(sig.probe[[ct]], ctewas$probe)
		#sig.i = which(ctewas$Pvalue[sig1.i] <= 0.05)
		sig.i = which(ctewas$Pvalue[sig1.i] <= pv.cutoff) #v3
		if (length(sig.i) > 0) {
			ctewas = ctewas[sig1.i[sig.i], ]
		} else {
			next
		}
		if (ifMC) ctewas = glm2mc(ctewas)
		ctewas$CT = ct
	} else {
		next;
	}
	premht.allct2 = rbind(premht.allct2, ctewas[, my.names])
	setTxtProgressBar(pb, i)
}
close(pb)
str(premht.allct2)

# fold gse
premht.allct3 = data.frame()
i = 1
pb <- txtProgressBar(min = 0, max = length(celltypes), style = 3)
for (i in 1:length(celltypes)) {
	ct = celltypes_paper_name[i]
	fold = ewas_fold3[i]
	#glm.f = paste0(fold,"/glm_pv_adj_anno.xls")
	glm.f = paste0(fold,"/", glmf)
	if (file.exists(glm.f)) {
		ctewas = r.table(glm.f, T)
		#sig.i = which(ctewas$BH_fdr_adj <= 0.05)
		sig1.i = match(sig.probe[[ct]], ctewas$probe)
		#sig.i = which(ctewas$Pvalue[sig1.i] <= 0.05)
		sig.i = which(ctewas$Pvalue[sig1.i] <= pv.cutoff) #v3
		if (length(sig.i) > 0) {
			ctewas = ctewas[sig1.i[sig.i], ]
		} else {
			next
		}
		if (ifMC) ctewas = glm2mc(ctewas)
		ctewas$CT = ct
	} else {
		next;
	}
	premht.allct3 = rbind(premht.allct3, ctewas[, my.names])
	setTxtProgressBar(pb, i)
}
close(pb)

premht.allct$ctprobe = paste0(premht.allct$probe, "_", premht.allct$CT)
premht.allct2$ctprobe = paste0(premht.allct2$probe, "_", premht.allct2$CT)
premht.allct3$ctprobe = paste0(premht.allct3$probe, "_", premht.allct3$CT)
premht.allct.meta$ctprobe = paste0(premht.allct.meta$probe, "_", premht.allct.meta$CT)
str(premht.allct)
str(premht.allct2)
str(premht.allct3)
str(sig.probe)

comm.probe = unique(intersect(premht.allct3$ctprobe, intersect(premht.allct$ctprobe, premht.allct2$ctprobe)))
length(comm.probe)
# 80

#m.i = match(premht.allct2$probe, premht.allct$probe)
m.i = match(comm.probe, premht.allct$ctprobe)
premht.allct_ = premht.allct[m.i,]
m.i = match(comm.probe, premht.allct2$ctprobe)
premht.allct2_ = premht.allct2[m.i,]
m.i = match(comm.probe, premht.allct3$ctprobe)
premht.allct3_ = premht.allct3[m.i,]

table(premht.allct_$CT)

output = cbind(premht.allct_[,1:3], premht.allct2_[,2:3], premht.allct3_[,2:3], premht.allct_[, 8:14])
colnames(output)[2:7] = paste0(rep(c("VACS", "WIHS", "GSE217633"), each=2), rep(c(".tvalue", ".Pvalue"), 3))
Directioni1 = ifelse(output[,2]<0, "-", "+")
Directioni2 = ifelse(output[,4]<0, "-", "+")
Directioni3 = ifelse(output[,6]<0, "-", "+")
output$Direction = paste0(Directioni1, Directioni2, Directioni3)
head(output)
w.table(output, "Three_cohort_tvalue_consistency.xls")

stopifnot(all(premht.allct_$probe == premht.allct2_$probe))
data = data.frame(VACS=premht.allct_[, comp.value], WIHS=premht.allct2_[, comp.value], GSE217633=premht.allct3_[, comp.value], CellType=premht.allct_$CT)
data$CellType = factor(data$CellType, levels=celltypes_paper_name)
str(data)
table(data$CellType)

ct.cor.l = tapply(1:nrow(data), data$CellType, function(i) {
	   x = data[i, 1:3]
	   cor.fit = psych::corr.test(x)
	   #r = cor.fit$estimate
	   #pv = signif(cor.fit$p.value, 2)
	   #cor.res = signif(c(r, pv), 2)
	   cor.fit
})

sink("corr.test.log")
lapply(ct.cor.l, function(x) {list(R=x$r, P=x$p)})
sink()

ct.cor = lapply(ct.cor.l, function(x) {
					if (is.null(x)) return (NULL)
					R = x$r
					P = x$p
					non.diag.mean = function(A) {
						mean(A[row(R) != col(R)], na.rm=T)
					}
					r = non.diag.mean(A=R)
					pv = non.diag.mean(A=P)
					cor.res = signif(c(r, pv), 2)
					cor.res
})


cor1 = do.call(rbind, ct.cor)
cor2 = cbind(rownames(cor1), cor1)
ct2 = paste0(cor2[, 1], ", r=", cor2[,2], ", p=", cor2[, 3])
levels(data$CellType) <- ct2

myct = names(table(data$CellType))
my_cols2 <- brewer.pal(8, "Set1")[c(1,2,3,4,5,8)][1:length(myct)]
cols = add.alpha(my_cols2, 0.7)
cols1 = add.alpha(my_cols2, 0.7)
cols2 = add.alpha(my_cols2, 0.5)

#cols = c()
mycol = c()
for (k in 1:nrow(data)) {
	if (data[k, 1] > 0 && data[k, 2] > 0 && data[k, 3] > 0)
	{
		mycol = c(mycol, "Consistent")
	} else if (data[k, 1] < 0 && data[k, 2] < 0 && data[k, 3] < 0) {
		mycol = c(mycol, "Consistent")
	} else if ((data[k, 1] * data[k, 3] ) > 0 ) {
		mycol = c(mycol, "PartialConsistent")
	} else {
		mycol = c(mycol, "Inconsistent")
	}
}
mycol2 = mycol3 = mycol
mycol2[mycol == "PartialConsistent"] = "Consistent"
mycol3[mycol == "PartialConsistent"] = "Inconsistent"
data$Consistency = as.factor(mycol2)
data$PartialConsistency = as.factor(mycol3)
table(data$Consistency)
table(data$PartialConsistency)
data2 = data
data2$PartialConsistent = mycol
(data3 = table(data2$CellType, data2$PartialConsistent))
data3 = as.data.frame(cbind(CellType=rownames(data3), data3))
w.table.excel(data3, "PartialConsistency.xlsx")

shapes1 = c(16, 17) 
#colors1 <- c("#E69F00", "#56B4E9")
colors1 <- add.alpha(my_cols2[1:2], 0.7)

pdf("tvalue_celltype_plot_3cohort_paper.pdf", 9, 6, family="ArialMT") #v3
opar <- par(no.readonly=TRUE)
par(mfrow = c(2, 3))
par(mar = c(5, 5, 8, 2), mgp = c(3, 1, 0))
for (i in 1:length(myct)) 
{
	ct = myct[i]
	dat = data[which(data$CellType == ct),]
	shapes <- shapes1[as.numeric(dat$PartialConsistency)]
	colors <- colors1[as.numeric(dat$Consistency)]
	scatterplot3d(dat[,1:3], pch = shapes, cex.symbols = 1.5, main = ct, color=colors)
	if (i == 1) legend("topleft", legend = levels(factor(mycol)),
					   col =  c(colors1, colors1[1]),
					   pch = c(shapes1, shapes1[2]),
					   cex = 0.8,
					   inset = -0.06,
					   bty = "n",
					   xpd = TRUE, horiz = TRUE
					   )
}
par(opar)
dev.off()


if (F) {
pdf("tvalue_celltype_plot_3cohort_paper_v4.pdf", 9, 6, family="ArialMT") #v3
opar <- par(no.readonly=TRUE)
#par(mfrow = c(2, 3))
layout(matrix(c(2, 2), 1))
par(mar = c(5, 5, 8, 2), mgp = c(3, 1, 0))
for (i in 1:length(myct)) 
{
	ct = myct[i]
	dat = data[which(data$CellType == ct),]
	shapes <- shapes1[as.numeric(dat$PartialConsistency)]
	colors <- colors1[as.numeric(dat$Consistency)]
	scatter3d(x=dat[,1], y=dat[,2], z=dat[,3], bty = "g", 
			  xlab = colnames(dat)[1],
			  ylab = colnames(dat)[2],
			  zlab = colnames(dat)[3],
			  #pch = shapes,
			  point.col = colors,
			  size = 2)
}
par(opar)
dev.off()


pdf("tvalue_celltype_plot_3cohort_paper_v2.pdf", 9, 6, family="ArialMT") #v3
opar <- par(no.readonly=TRUE)
par(mfrow = c(2, 3))
par(mar = c(5, 5, 8, 2), mgp = c(3, 1, 0))
for (i in 1:length(myct)) 
{
	ct = myct[i]
	dat = data[which(data$CellType == ct),]
	shapes <- shapes1[as.numeric(dat$Consistency)]
	colors <- colors1[as.numeric(dat$Consistency)]
	scatter3D(x=dat[,1], y=dat[,2], z=dat[,3], bty = "g", 
			  pch = shapes,
			  col = colors,
			  cex = 2,
			  ticktype = "detailed")
	if (i == 1) legend("topleft", legend = levels(dat$Consistency),
		   col =  colors1,
		   pch = shapes1,
		   xpd = TRUE, horiz = TRUE)
}
par(opar)
dev.off()


pdf("tvalue_celltype_plot_3cohort_paper_v3.pdf", 9, 6, family="ArialMT") #v3
opar <- par(no.readonly=TRUE)
par(mfrow = c(2, 3))
par(mar = c(5, 5, 8, 2), mgp = c(3, 1, 0))
for (i in 1:length(myct)) 
{
	ct = myct[i]
	dat = data[which(data$CellType == ct),]
	shapes <- shapes1[as.numeric(dat$Consistency)]
	colors <- colors1[as.numeric(dat$Consistency)]
	scatter3D(x=dat[,1], y=dat[,2], z=dat[,3], bty = "g", 
			  phi = 0,
			  pch = shapes,
			  col = colors,
			  cex = 2,
			  ticktype = "detailed")
	if (i == 1) legend("topleft", legend = levels(dat$Consistency),
		   col =  colors1,
		   pch = shapes1,
		   xpd = TRUE, horiz = TRUE)
}
par(opar)
dev.off()

}


