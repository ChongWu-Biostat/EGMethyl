setwd("/Users/uniquechong/Dropbox (Personal)/Chong Wu/Undergoing/TWAS/Codes/Methy4/ResAna/")

# get the common gene set



res.out = as.data.frame(matrix(NA,3,13))


file.name = "enhancer_MCF_scz1"
#file.name = "V1_scz2_ALL_hui"
data = readRDS(paste0("res",file.name,".rds"))
data = data[!duplicated(data[,c(2,3,4)]),]
dim(data)
data[,3] = data[,3] +500
data[,4] = data[,4] -500

mcf.gene.list = data[,2]
rownames(data) = paste0(data[,2],data[,3],data[,4])
cis = readRDS("summary_best_hui2_cis_scz1.rds")

rownames(cis) = paste0(cis[,4],cis[,2],cis[,3])
cis = cis[rownames(data),]

cis2 = readRDS("summary_best_hui2_cis_scz2.rds")
rownames(cis2) = paste0(cis2[,4],cis2[,2],cis2[,3])
cis2 = cis2[rownames(data),]

sum(data[,2] %in% cis[,4])
data1 = as.data.frame(matrix(NA,dim(data)[1],12))

data1[,1:8] = data
data1[,9:10] = cis[,c(6,5)]
data1[,11:12] = cis2[,c(6,5)]

data = data1
data[,3:12] = apply(data[,3:12], 2, FUN=as.numeric)

cutoff = 0.05/ dim(data)[1]
data.tmp = data[data[,6] < cutoff,]


data.tmp = data.tmp[data.tmp[,10] > 5 * 10^{-8},]
dim(data.tmp)
#data.tmp = data.tmp[!duplicated(data.tmp[,2]),]

res.out[1,3] = paste0(dim(data.tmp)[1], "/", sum(data.tmp[,12]< 5*10^{-8}))

data.tmp = data[data[,7] < cutoff,]
data.tmp = data.tmp[data.tmp[,10] > 5 * 10^{-8},]

res.out[2,3] = paste0(dim(data.tmp)[1],"/", sum(data.tmp[,12]< 5*10^{-8}))

res.out[3,3] = dim(data)[1]


file.name = "enhancer_hip_scz1"
data = readRDS(paste0("res",file.name,".rds"))
data = data[!duplicated(data[,c(2,3,4)]),]
dim(data)
data[,3] = data[,3] +500
data[,4] = data[,4] -500

rownames(data) = paste0(data[,2],data[,3],data[,4])

cis = readRDS("summary_best_cis_scz1.rds")

rownames(cis) = paste0(cis[,4],cis[,2],cis[,3])
cis = cis[rownames(data),]

cis2 = readRDS("summary_best_cis_scz2.rds")

rownames(cis2) = paste0(cis2[,4],cis2[,2],cis2[,3])
cis2 = cis2[rownames(data),]

sum(data[,2] %in% cis[,4])

hip.gene.list = data[,2]

data1 = as.data.frame(matrix(NA,dim(data)[1],12))

data1[,1:8] = data
data1[,9:10] = cis[,c(6,5)]
data1[,11:12] = cis2[,c(6,5)]

data = data1
data[,3:12] = apply(data[,3:12], 2, FUN=as.numeric)

cutoff = 0.05/dim(data)[1]

data.tmp = data[data[,6] < cutoff,]


data.tmp = data.tmp[data.tmp[,10] > 5 * 10^{-8},]
dim(data.tmp)
#data.tmp = data.tmp[!duplicated(data.tmp[,2]),]



res.out[1,4] = paste0(dim(data.tmp)[1], "/", sum(data.tmp[,12]< 5*10^{-8}))

data.tmp = data[data[,7] < cutoff,]
data.tmp = data.tmp[data.tmp[,10] > 5 * 10^{-8},]

res.out[2,4] = paste0(dim(data.tmp)[1],"/", sum(data.tmp[,12]< 5*10^{-8}))
res.out[3,4] = dim(data)[1]





file.name = "new_scz1_ALL_hui2"
data = readRDS(paste0("res",file.name,".rds"))
data = data[!duplicated(data[,c(2,3,4)]),]
dim(data)
rownames(data) = paste0(data[,2],data[,3],data[,4])

cis = readRDS("summary_best_hui2_cis_scz1.rds")
rownames(cis) = paste0(cis[,4],cis[,2],cis[,3])
cis = cis[rownames(data),]

cis2 = readRDS("summary_best_hui2_cis_scz2.rds")

rownames(cis2) = paste0(cis2[,4],cis2[,2],cis2[,3])
cis2 = cis2[rownames(data),]
sum(data[,2] %in% cis[,4])

data1 = as.data.frame(matrix(NA,dim(data)[1],18))

data1[,1:14] = data[,c(1:7,10,11,8,9,16,17)]
data1[,15:16] = cis[,c(6,5)]
data1[,17:18] = cis2[,c(6,5)]

data = data1
data[,8:18] = apply(data[,8:18], 2, FUN=as.numeric)

data = data[data[,2] %in% mcf.gene.list,] # if there is no information in enhancer, there is no point to conduct test.
gene.name = data[,2]

cutoff = 0.05/dim(data)[1]

data.tmp = data[data[,10] < cutoff,]

data.tmp = data.tmp[data.tmp[,16] > 5 * 10^{-8},]
dim(data.tmp)
#data.tmp = data.tmp[!duplicated(data.tmp[,2]),]


res.out[1,12] = paste0(sum(data.tmp[,18]< 5*10^{-8}) / dim(data.tmp)[1], "/",phyper(sum(data.tmp[,18]< 5*10^{-8})-1,sum(data[,18] <  5*10^{-8}), dim(data)[1] - dim(data.tmp)[1], dim(data.tmp)[1] , lower.tail = FALSE) )

res.out[1,1] = paste0(dim(data.tmp)[1], "/", sum(data.tmp[,18]< 5*10^{-8}))

data.tmp = data[data[,11] < cutoff,]
data.tmp = data.tmp[data.tmp[,16] > 5 * 10^{-8},]

res.out[2,1] = paste0(dim(data.tmp)[1],"/", sum(data.tmp[,18]< 5*10^{-8}))

res.out[2,12] = paste0(sum(data.tmp[,18]< 5*10^{-8}) / dim(data.tmp)[1], "/",phyper(sum(data.tmp[,18]< 5*10^{-8})-1,sum(data[,18] <  5*10^{-8}), dim(data)[1] - dim(data.tmp)[1], dim(data.tmp)[1] , lower.tail = FALSE) )

res.out[3,1] = dim(data)[1]


file.name = "new_scz1_ALL_hip"
data = readRDS(paste0("res",file.name,".rds"))
data = data[!duplicated(data[,c(2,3,4)]),]
dim(data)
rownames(data) = paste0(data[,2],data[,3],data[,4])

cis = readRDS("summary_best_cis_scz1.rds")
rownames(cis) = paste0(cis[,4],cis[,2],cis[,3])
cis = cis[rownames(data),]

cis2 = readRDS("summary_best_cis_scz2.rds")
rownames(cis2) = paste0(cis2[,4],cis2[,2],cis2[,3])
cis2 = cis2[rownames(data),]

sum(data[,2] %in% cis[,4])
data1 = as.data.frame(matrix(NA,dim(data)[1],18))


data1[,1:14] = data[,c(1:7,10,11,8,9,16,17)]
data1[,15:16] = cis[,c(6,5)]
data1[,17:18] = cis2[,c(6,5)]

gene.name = c(gene.name,data[,2])

data = data1
data[,8:18] = apply(data[,8:18], 2, FUN=as.numeric)

data = data[data[,2] %in% hip.gene.list,]

cutoff = 0.05/dim(data)[1]

data.tmp = data[data[,10] < cutoff,]


data.tmp = data.tmp[data.tmp[,16] > 5 * 10^{-8},]
dim(data.tmp)
#data.tmp = data.tmp[!duplicated(data.tmp[,2]),]

res.out[1,2] = paste0(dim(data.tmp)[1], "/", sum(data.tmp[,18]< 5*10^{-8}))


res.out[1,13] = paste0(sum(data.tmp[,18]< 5*10^{-8}) / dim(data.tmp)[1], "/",phyper(sum(data.tmp[,18]< 5*10^{-8})-1,sum(data[,18] <  5*10^{-8}), dim(data)[1] - dim(data.tmp)[1], dim(data.tmp)[1] , lower.tail = FALSE) )

data.tmp = data[data[,11] < cutoff,]
data.tmp = data.tmp[data.tmp[,16] > 5 * 10^{-8},]

res.out[2,2] = paste0(dim(data.tmp)[1],"/", sum(data.tmp[,18]< 5*10^{-8}))
res.out[3,2] = dim(data)[1]

res.out[2,13] = paste0(sum(data.tmp[,18]< 5*10^{-8}) / dim(data.tmp)[1], "/",phyper(sum(data.tmp[,18]< 5*10^{-8})-1,sum(data[,18] <  5*10^{-8}), dim(data)[1] - dim(data.tmp)[1], dim(data.tmp)[1] , lower.tail = FALSE) )







file.name = "std_scz1"

#file.name = "V1_scz2_ALL_hui"
data = readRDS(paste0("res",file.name,".rds"))
data = data[!duplicated(data[,c(2,3,4)]),]
dim(data)
data[,3] = data[,3] +500
data[,4] = data[,4] -500

rownames(data) = paste0(data[,2],data[,3],data[,4])

cis = readRDS("summary_best_hui2_cis_scz1.rds")
#cis = readRDS("summary_best_cis_scz1.rds")

rownames(cis) = paste0(cis[,4],cis[,2],cis[,3])
cis = cis[rownames(data),]

cis2 = readRDS("summary_best_hui2_cis_scz2.rds")
#cis2 = readRDS("summary_best_cis_scz2.rds")

rownames(cis2) = paste0(cis2[,4],cis2[,2],cis2[,3])
cis2 = cis2[rownames(data),]

sum(data[,2] %in% cis[,4])

data1 = as.data.frame(matrix(NA,dim(data)[1],12))

data1[,1:8] = data
data1[,9:10] = cis[,c(6,5)]
data1[,11:12] = cis2[,c(6,5)]

data = data1
data[,3:12] = apply(data[,3:12], 2, FUN=as.numeric)

cutoff = 0.05/dim(data)[1]

data.tmp = data[data[,6] < cutoff,]


data.tmp = data.tmp[data.tmp[,10] > 5 * 10^{-8},]
dim(data.tmp)
#data.tmp = data.tmp[!duplicated(data.tmp[,2]),]

res.out[1,5] = paste0(dim(data.tmp)[1], "/", sum(data.tmp[,12]< 5*10^{-8}))

data.tmp = data[data[,7] < cutoff,]
data.tmp = data.tmp[data.tmp[,10] > 5 * 10^{-8},]

res.out[2,5] = paste0(dim(data.tmp)[1],"/", sum(data.tmp[,12]< 5*10^{-8}))
res.out[3,5] = dim(data)[1]


###
file.name = "summary_scz1_YFS"
data = readRDS(paste0(file.name,".rds"))


SPU1 = data$`SPU(1)`
SPU2 = data$`SPU(2)`


fina.aSPU = aSPU
fina.SPU1 = SPU1
fina.SPU2 = SPU2

fina.aSPU = SPU1
YFS = readRDS("scz1_YFS.rds")
cutoff = 0.05/dim(YFS)[1]

fina.aSPU = fina.aSPU[fina.aSPU[,6] < cutoff,]
fina.aSPU = fina.aSPU[!duplicated(fina.aSPU[,2]),]


fina.aSPU = fina.aSPU[fina.aSPU[,17] > 5* 10^{-8},]

res.out[1,6] = paste0(dim(fina.aSPU)[1], "/",sum(fina.aSPU[,18]< 5*10^{-8}))


fina.aSPU = SPU2
fina.aSPU = fina.aSPU[fina.aSPU[,7] < cutoff,]
fina.aSPU = fina.aSPU[!duplicated(fina.aSPU[,2]),]

fina.aSPU = fina.aSPU[fina.aSPU[,17] > 5* 10^{-8},]

res.out[2,6] = paste0(dim(fina.aSPU)[1], "/",sum(fina.aSPU[,18]< 5*10^{-8}))

res.out[3,6] = dim(YFS)[1]

file.name = "summary_scz1_NTR"
data = readRDS(paste0(file.name,".rds"))

aSPU = data$aSPU
SPU1 = data$`SPU(1)`
SPU2 = data$`SPU(2)`

fina.aSPU = rbind(fina.aSPU,aSPU)
fina.SPU1 = rbind(fina.SPU1,SPU1)
fina.SPU2 = rbind(fina.SPU2,SPU2)

YFS = readRDS("scz1_NTR.rds")
cutoff = 0.05/dim(YFS)[1]


fina.aSPU = SPU1
fina.aSPU = fina.aSPU[fina.aSPU[,6] < cutoff,]
fina.aSPU = fina.aSPU[!duplicated(fina.aSPU[,2]),]
fina.aSPU = fina.aSPU[fina.aSPU[,17] > 5* 10^{-8},]

res.out[1,7] = paste0(dim(fina.aSPU)[1], "/",sum(fina.aSPU[,18]< 5*10^{-8}))

res.out[3,7] = dim(YFS)[1]

fina.aSPU = SPU2
fina.aSPU = fina.aSPU[fina.aSPU[,7] < cutoff,]
fina.aSPU = fina.aSPU[!duplicated(fina.aSPU[,2]),]

fina.aSPU = fina.aSPU[fina.aSPU[,17] > 5* 10^{-8},]

res.out[2,7] = paste0(dim(fina.aSPU)[1], "/",sum(fina.aSPU[,18]< 5*10^{-8}))

res.out[3,7] = dim(YFS)[1]


file.name = "summary_scz1_METSIM"
data = readRDS(paste0(file.name,".rds"))

aSPU = data$aSPU
SPU1 = data$`SPU(1)`
SPU2 = data$`SPU(2)`

fina.aSPU = rbind(fina.aSPU,aSPU)
fina.SPU1 = rbind(fina.SPU1,SPU1)
fina.SPU2 = rbind(fina.SPU2,SPU2)

fina.aSPU = SPU1

YFS = readRDS("scz1_METSIM.rds")
cutoff = 0.05/dim(YFS)[1]

fina.aSPU = fina.aSPU[fina.aSPU[,6] < cutoff,]
fina.aSPU = fina.aSPU[!duplicated(fina.aSPU[,2]),]
fina.aSPU = fina.aSPU[fina.aSPU[,17] > 5* 10^{-8},]


res.out[1,8] = paste0(dim(fina.aSPU)[1],"/",sum(fina.aSPU[,18]< 5*10^{-8}))


fina.aSPU = SPU2
fina.aSPU = fina.aSPU[fina.aSPU[,7] < cutoff,]
fina.aSPU = fina.aSPU[!duplicated(fina.aSPU[,2]),]

fina.aSPU = fina.aSPU[fina.aSPU[,17] > 5* 10^{-8},]

res.out[2,8] = paste0(dim(fina.aSPU)[1],"/",sum(fina.aSPU[,18]< 5*10^{-8}))

res.out[3,8] = dim(YFS)[1]

file.name = "summary_scz1_CMC2"
data = readRDS(paste0(file.name,".rds"))

aSPU = data$aSPU
SPU1 = data$`SPU(1)`
SPU2 = data$`SPU(2)`

fina.aSPU = rbind(fina.aSPU,aSPU)
fina.SPU1 = rbind(fina.SPU1,SPU1)
fina.SPU2 = rbind(fina.SPU2,SPU2)

fina.aSPU = SPU1

YFS = readRDS("scz1_CMC2.rds")
cutoff = 0.05/dim(YFS)[1]

fina.aSPU = fina.aSPU[fina.aSPU[,6] < cutoff,]
fina.aSPU = fina.aSPU[!duplicated(fina.aSPU[,2]),]
fina.aSPU = fina.aSPU[fina.aSPU[,17] > 5* 10^{-8},]

res.out[1,9] = paste0(dim(fina.aSPU)[1], "/",sum(fina.aSPU[,18]< 5*10^{-8}))


fina.aSPU = SPU2
fina.aSPU = fina.aSPU[fina.aSPU[,7] < cutoff,]
fina.aSPU = fina.aSPU[!duplicated(fina.aSPU[,2]),]
fina.aSPU = fina.aSPU[fina.aSPU[,17] > 5* 10^{-8},]

res.out[2,9] = paste0(dim(fina.aSPU)[1], "/",sum(fina.aSPU[,18]< 5*10^{-8}))

res.out[3,9] = dim(YFS)[1]

# mQTL 500 kb results, newly added

file.name = "V1_scz1_mQTL_20kb"
data = readRDS(paste0("res",file.name,".rds"))
data = data[!duplicated(data[,c(2,3,4)]),]

dim(data)
#data[,3] = data[,3] +500
#data[,4] = data[,4] -500

dim(data)
rownames(data) = paste0(data[,2],data[,3],data[,4])

cis = readRDS("summary_best_hui2_cis_scz1.rds")
rownames(cis) = paste0(cis[,4],cis[,2],cis[,3])
cis = cis[rownames(data),]

cis2 = readRDS("summary_best_hui2_cis_scz2.rds")
rownames(cis2) = paste0(cis2[,4],cis2[,2],cis2[,3])
cis2 = cis2[rownames(data),]

sum(data[,2] %in% cis[,4])

data1 = as.data.frame(matrix(NA,dim(data)[1],12))

data1[,1:8] = data[,c(1:5,7:9)]
data1[,9:10] = cis[,c(6,5)]
data1[,11:12] = cis2[,c(6,5)]

data = data1
data[,3:12] = apply(data[,3:12], 2, FUN=as.numeric)

gene.name = data[,2]
cutoff = 0.05/10000
data.tmp = data[data[,6] < cutoff,]


cutoff = 0.05/ dim(data)[1]
data.tmp = data[data[,6] < cutoff,]

data.tmp = data.tmp[data.tmp[,10] > 5 * 10^{-8},]
dim(data.tmp)
#data.tmp = data.tmp[!duplicated(data.tmp[,2]),]

res.out[1,10] = paste0(dim(data.tmp)[1], "/", sum(data.tmp[,12]< 5*10^{-8}))

data.tmp = data[data[,7] < cutoff,]
data.tmp = data.tmp[data.tmp[,10] > 5 * 10^{-8},]

res.out[2,10] = paste0(dim(data.tmp)[1],"/", sum(data.tmp[,12]< 5*10^{-8}))
res.out[3,10] = dim(data)[1]





file.name = "V1_scz1_mQTL_500kb"
data = readRDS(paste0("res",file.name,".rds"))
data = data[!duplicated(data[,c(2,3,4)]),]

dim(data)
#data[,3] = data[,3] +500
#data[,4] = data[,4] -500

dim(data)
rownames(data) = paste0(data[,2],data[,3],data[,4])

cis = readRDS("summary_best_hui2_cis_scz1.rds")
rownames(cis) = paste0(cis[,4],cis[,2],cis[,3])
cis = cis[rownames(data),]

cis2 = readRDS("summary_best_hui2_cis_scz2.rds")
rownames(cis2) = paste0(cis2[,4],cis2[,2],cis2[,3])
cis2 = cis2[rownames(data),]

sum(data[,2] %in% cis[,4])

data1 = as.data.frame(matrix(NA,dim(data)[1],12))

data1[,1:8] = data[,c(1:5,7:9)]
data1[,9:10] = cis[,c(6,5)]
data1[,11:12] = cis2[,c(6,5)]

data = data1
data[,3:12] = apply(data[,3:12], 2, FUN=as.numeric)

gene.name = data[,2]
cutoff = 0.05/10000
data.tmp = data[data[,6] < cutoff,]


cutoff = 0.05/ dim(data)[1]
data.tmp = data[data[,6] < cutoff,]

data.tmp = data.tmp[data.tmp[,10] > 5 * 10^{-8},]
dim(data.tmp)
#data.tmp = data.tmp[!duplicated(data.tmp[,2]),]

res.out[1,11] = paste0(dim(data.tmp)[1], "/", sum(data.tmp[,12]< 5*10^{-8}))

data.tmp = data[data[,7] < cutoff,]
data.tmp = data.tmp[data.tmp[,10] > 5 * 10^{-8},]

res.out[2,11] = paste0(dim(data.tmp)[1],"/", sum(data.tmp[,12]< 5*10^{-8}))
res.out[3,11] = dim(data)[1]



res.out[c(1,2,3),] = res.out[c(3,1,2),]


res.out

res.out = res.out[,1:9]
library(xtable)
#res.out = res.out[,1:6]
xtable(res.out)



