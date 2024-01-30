
phe1$HIV = phe1$time
phe1 = phe1[-which(phe1$HIV == "sample pre ART"),]
table(phe1$HIV)

phe = phe1
dat = dat1[,rownames(phe)]

