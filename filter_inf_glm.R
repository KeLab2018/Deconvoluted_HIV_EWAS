#!/usr/bin/env Rscript

rm(list=ls());
library(myfunc);
source("~/mybiotools/r/myfunc.R");
options("scipen"=1, "digits"=4, stringsAsFactors=FALSE);
args<-commandArgs(T)

glm.f = "glm_pv_adj_anno.xls"
ctewas = r.table(glm.f, T)
ctewas = ctewas[!is.infinite(ctewas$t.value), ]
w.table(ctewas, glm.f, T, F)

