######################################################
# Multiple Imputation for cell type proportion    ####
# Version 1.00                                    ####
# Feb 2, 2015                                     ####
# Author: Chong Wu, Weihua Guan                   ####
######################################################

setwd("/Users/uniquechong/Dropbox (Personal)/Chong Wu/Undergoing/TWAS/Codes/Methy4/ResAna/")

cutoff = 0.05/10000
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

data.tmp = data[data[,6] < cutoff | data[,7] < cutoff,]

data.tmp2 = data.tmp[!duplicated(data.tmp[,2]),]

data.tmp2 = data.tmp2[data.tmp2[,10] > 5e-8,]
enhancer.only = data.tmp2
enhancer.only.gene = unique(gene.name) # for easier to reuse the codes, here, ehnacer.only stands for the results of standard gene based analysis with several kb extension and only focus on the mQTLs.




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

gene.name = data[,2]


data.MCF = data[data[,10] < cutoff | data[,11] < cutoff,]

data.MCF$weight = "MCF7"



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

data.tmp = data[data[,10] < cutoff | data[,11] < cutoff,]
data.tmp$weight = "HIP"

data.new = rbind(data.tmp,data.MCF)
data.tmp = data.new

data.tmp2 = data.tmp[!duplicated(data.tmp[,2]),]
data.tmp2 = data.tmp2[data.tmp2[,16] > 5e-8,]
enhancer.methy = data.tmp2
enhancer.methy.gene = unique(gene.name)






file.name = "enhancer_MCF_scz1"
#file.name = "V1_scz2_ALL_hui"
data = readRDS(paste0("res",file.name,".rds"))
data = data[!duplicated(data[,c(2,3,4)]),]

dim(data)
data[,3] = data[,3] +500
data[,4] = data[,4] -500

dim(data)
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
gene.name = data[,2]

data[,3:12] = apply(data[,3:12], 2, FUN=as.numeric)


#data.tmp = data[data[,6] > 0.05 & data[,7] > 0.05 & data[,8] >0.05,]
data.MCF = data[data[,6] < cutoff | data[,7] < cutoff,]


data.MCF$weight = "MCF7"
file.name = "enhancer_hip_scz1"
#file.name = "V1_scz2_ALL_hui"
data = readRDS(paste0("res",file.name,".rds"))
data = data[!duplicated(data[,c(2,3,4)]),]
dim(data)
data[,3] = data[,3] +500
data[,4] = data[,4] -500

data = data[!duplicated(data[,c(2,3,4)]),]
dim(data)
rownames(data) = paste0(data[,2],data[,3],data[,4])

cis = readRDS("summary_best_cis_scz1.rds")

rownames(cis) = paste0(cis[,4],cis[,2],cis[,3])
cis = cis[rownames(data),]
cis2 = readRDS("summary_best_cis_scz2.rds")

rownames(cis2) = paste0(cis2[,4],cis2[,2],cis2[,3])
cis2 = cis2[rownames(data),]

gene.name = c(gene.name,data[,2])


sum(data[,2] %in% cis[,4])
data1 = as.data.frame(matrix(NA,dim(data)[1],12))

data1[,1:8] = data
data1[,9:10] = cis[,c(6,5)]
data1[,11:12] = cis2[,c(6,5)]
data = data1
data[,3:12] = apply(data[,3:12], 2, FUN=as.numeric)

data.tmp = data[data[,6] < cutoff | data[,7] < cutoff,]
data.tmp$weight = "HIP"

data.new = rbind(data.tmp,data.MCF)
data.tmp = data.new
data.tmp2 = data.tmp[!duplicated(data.tmp[,2]),]

data.tmp2 = data.tmp2[data.tmp2[,10] > 5e-8,]

enhancer.coding = data.tmp2
enhancer.coding.gene = unique(gene.name)


file.name = "std_scz1"
data = readRDS(paste0("res",file.name,".rds"))
data = data[!duplicated(data[,c(2,3,4)]),]
dim(data)
data[,3] = data[,3] +500
data[,4] = data[,4] -500
rownames(data) = paste0(data[,2],data[,3],data[,4])


cis = readRDS("summary_best_hui2_cis_scz1.rds")
#cis = readRDS("summary_best_cis_scz1.rds")

rownames(cis) = paste0(cis[,4],cis[,2],cis[,3])

sum(data[,2] %in% cis[,4])

cis = cis[rownames(data),]

cis2 = readRDS("summary_best_hui2_cis_scz2.rds")
#cis2 = readRDS("summary_best_cis_scz2.rds")

rownames(cis2) = paste0(cis2[,4],cis2[,2],cis2[,3])
cis2 = cis2[rownames(data),]



data1 = as.data.frame(matrix(NA,dim(data)[1],12))

data1[,1:8] = data
data1[,9:10] = cis[,c(6,5)]
data1[,11:12] = cis2[,c(6,5)]

data = data1
data[,3:12] = apply(data[,3:12], 2, FUN=as.numeric)


standard = data[data[,6] < cutoff | data[,7] < cutoff,]

standard = standard[standard[,10] >5e-8,]
standard.gene = unique(data[,2])



## TWAS results
YFS = readRDS("scz1_YFS.rds")
CMC2 = readRDS("scz1_CMC2.rds")
NTR = readRDS("scz1_NTR.rds")
METSIM = readRDS("scz1_METSIM.rds")

gene.expression = rbind(YFS,CMC2)
gene.expression = rbind(gene.expression,NTR)
gene.expression = rbind(gene.expression,METSIM)

length(unique(gene.expression[,2]))
gene.expression[,3:16] = apply(gene.expression[,3:16], 2, FUN=as.numeric)
gene.expression.org = gene.expression
TWAS.name = gene.expression[,2]

file.name = "summary_scz1_YFS"
data = readRDS(paste0(file.name,".rds"))


SPU1 = data$`SPU(1)`
SPU2 = data$`SPU(2)`


SPU1 = SPU1[SPU1[,6] < cutoff,]
SPU1 = SPU1[!duplicated(SPU1[,2]),]
SPU1 = SPU1[SPU1[,17] > 5e-8,]


SPU2 = SPU2[SPU2[,7] < cutoff,]
SPU2 = SPU2[!duplicated(SPU2[,2]),]
SPU2 = SPU2[SPU2[,17] > 5e-8,]

fina.SPU1 = SPU1
fina.SPU2 = SPU2


file.name = "summary_scz1_NTR"
data = readRDS(paste0(file.name,".rds"))

SPU1 = data$`SPU(1)`
SPU2 = data$`SPU(2)`


SPU1 = SPU1[SPU1[,6] < cutoff,]
SPU1 = SPU1[!duplicated(SPU1[,2]),]
SPU1 = SPU1[SPU1[,17] > 5e-8,]

SPU2 = SPU2[SPU2[,7] < cutoff,]
SPU2 = SPU2[!duplicated(SPU2[,2]),]
SPU2 = SPU2[SPU2[,17] > 5e-8,]

fina.SPU1 = rbind(fina.SPU1,SPU1)
fina.SPU2 = rbind(fina.SPU2,SPU2)

file.name = "summary_scz1_METSIM"
data = readRDS(paste0(file.name,".rds"))

aSPU = data$aSPU
SPU1 = data$`SPU(1)`
SPU2 = data$`SPU(2)`
SPU1 = SPU1[SPU1[,6] < cutoff,]
SPU1 = SPU1[!duplicated(SPU1[,2]),]
SPU1 = SPU1[SPU1[,17] > 5e-8,]

SPU2 = SPU2[SPU2[,7] < cutoff,]
SPU2 = SPU2[!duplicated(SPU2[,2]),]
SPU2 = SPU2[SPU2[,17] > 5e-8,]

fina.SPU1 = rbind(fina.SPU1,SPU1)
fina.SPU2 = rbind(fina.SPU2,SPU2)

file.name = "summary_scz1_CMC2"
data = readRDS(paste0(file.name,".rds"))

aSPU = data$aSPU
SPU1 = data$`SPU(1)`
SPU2 = data$`SPU(2)`

SPU1 = SPU1[SPU1[,6] < cutoff,]
SPU1 = SPU1[!duplicated(SPU1[,2]),]
SPU1 = SPU1[SPU1[,17] > 5e-8,]

SPU2 = SPU2[SPU2[,7] < cutoff,]
SPU2 = SPU2[!duplicated(SPU2[,2]),]
SPU2 = SPU2[SPU2[,17] > 5e-8,]

TWAS.SPU1 = rbind(fina.SPU1,SPU1)
TWAS.SPU2 = rbind(fina.SPU2,SPU2)

TWAS = rbind(TWAS.SPU1,TWAS.SPU2)

TWAS.sig = TWAS[,2]

## Gene expression
YFS = readRDS("scz1_YFS.rds")
CMC2 = readRDS("scz1_CMC2.rds")
NTR = readRDS("scz1_NTR.rds")
METSIM = readRDS("scz1_METSIM.rds")

gene.expression = rbind(YFS,CMC2)
gene.expression = rbind(gene.expression,NTR)
gene.expression = rbind(gene.expression,METSIM)

length(unique(gene.expression[,2]))
gene.expression[,3:16] = apply(gene.expression[,3:16], 2, FUN=as.numeric)
gene.expression.org = gene.expression
TWAS.gene = unique(gene.expression[,2])



enhancer.only.sig = enhancer.only[,2]
#enhancer.only.gene

enhancer.methy.sig = enhancer.methy[,2]
#enhancer.methy.gene

enhancer.coding.sig = enhancer.coding[,2]
#enhancer.coding.gene

standard.sig = standard[,2]
#standard.gene

TWAS.sig = as.character(TWAS.sig)
#TWAS.gene

gene.used = intersect(TWAS.gene,enhancer.coding.gene)
gene.used = intersect(gene.used, enhancer.methy.gene)
gene.used = intersect(gene.used,standard.gene) #3521

#gene.used = intersect(gene.used, enhancer.only.gene)


standard.sig = standard.sig[standard.sig%in% gene.used]
enhancer.only.sig =enhancer.only.sig[enhancer.only.sig %in% gene.used]
enhancer.methy.sig =enhancer.methy.sig[enhancer.methy.sig %in% gene.used]
enhancer.coding.sig = enhancer.coding.sig[enhancer.coding.sig%in% gene.used]
TWAS.sig = TWAS.sig[TWAS.sig%in%gene.used ]

standard.sig = unique(standard.sig)
enhancer.only.sig = unique(enhancer.only.sig)
enhancer.methy.sig = unique(enhancer.methy.sig)
enhancer.coding.sig = unique(enhancer.coding.sig)
TWAS.sig = unique(TWAS.sig)
#TWAS.sig = enhancer.methy.sig
saveRDS(TWAS.sig,"twas.sig.gene.rds")

enhancer.coding.sig[!enhancer.coding.sig %in% TWAS.sig]
library(VennDiagram)


setwd("/Users/uniquechong/Dropbox (Personal)/Chong Wu/Undergoing/TWAS/Codes/Methy4/ResAna/Res_manuscript_final/EGM_results_final")


pdf("Venn_scz1_2.pdf")

#standard.sig = enhancer.only.sig

enhancer.only.sig = enhancer.methy.sig
venn.plot = draw.quad.venn(
area1 = length(enhancer.coding.sig),
area2 = length(enhancer.only.sig),
area3 = length(TWAS.sig),
area4 = length(standard.sig),
n12 = sum(enhancer.coding.sig %in% enhancer.only.sig),
n13 = sum(enhancer.coding.sig %in% TWAS.sig),
n14 = sum(enhancer.coding.sig %in% standard.sig),
n23 = sum(enhancer.only.sig%in% TWAS.sig),
n24 = sum(enhancer.only.sig %in% standard.sig),
n34 = sum(TWAS.sig %in% standard.sig),
n123 = length(intersect(intersect(enhancer.coding.sig, enhancer.only.sig),TWAS.sig)),
n124 = length(intersect(intersect(enhancer.coding.sig, enhancer.only.sig),standard.sig)),
n134 = length(intersect(intersect(enhancer.coding.sig, TWAS.sig),standard.sig)),
n234 = length(intersect(intersect(enhancer.only.sig, TWAS.sig),standard.sig)),
n1234 = length(intersect(intersect(intersect(enhancer.coding.sig, enhancer.only.sig), TWAS.sig),standard.sig)),
category = c("E+G", "E+G+Methyl", "TWAS","STD"),
lty = "blank",
cex = 2,
margin = 0.1,
fill = c("skyblue", "pink1",
"mediumorchid", "orange"));
dev.off()



pdf("Venn_scz1_5.pdf")
venn.plot <- draw.quintuple.venn(
area1 = length(enhancer.coding.sig),
area2 = length(enhancer.only.sig),
area3 = length(TWAS.sig),
area4 = length(standard.sig),
area5 = length(enhancer.methy.sig),
n12 = sum(enhancer.coding.sig %in% enhancer.only.sig),
n13 = sum(enhancer.coding.sig %in% TWAS.sig),
n14 = sum(enhancer.coding.sig %in% standard.sig),
n15 = sum(enhancer.coding.sig %in% enhancer.methy.sig),
n23 = sum(enhancer.only.sig%in% TWAS.sig),
n24 = sum(enhancer.only.sig %in% standard.sig),
n25 = sum(enhancer.only.sig %in% enhancer.methy.sig),
n34 = sum(TWAS.sig %in% standard.sig),
n35 = sum(TWAS.sig %in% enhancer.methy.sig),
n45 = sum(standard.sig %in% enhancer.methy.sig),
n123 = length(intersect(intersect(enhancer.coding.sig, enhancer.only.sig),TWAS.sig)),
n124 = length(intersect(intersect(enhancer.coding.sig, enhancer.only.sig),standard.sig)),
n125 = length(intersect(intersect(enhancer.coding.sig, enhancer.only.sig),enhancer.methy.sig)),
n134 = length(intersect(intersect(enhancer.coding.sig, TWAS.sig),standard.sig)),
n135 = length(intersect(intersect(enhancer.coding.sig, TWAS.sig),enhancer.methy.sig)),
n145 = length(intersect(intersect(enhancer.coding.sig, standard.sig),enhancer.methy.sig)),
n234 = length(intersect(intersect(enhancer.only.sig, TWAS.sig),standard.sig)),
n235 = length(intersect(intersect(enhancer.only.sig, TWAS.sig),enhancer.methy.sig)),
n245 = length(intersect(intersect(enhancer.only.sig, standard.sig),enhancer.methy.sig)),
n345 = length(intersect(intersect(TWAS.sig, standard.sig),enhancer.methy.sig)),
n1234 = length(intersect(intersect(intersect(enhancer.coding.sig, enhancer.only.sig), TWAS.sig),standard.sig)),
n1235 = length(intersect(intersect(intersect(enhancer.coding.sig, enhancer.only.sig), TWAS.sig),enhancer.methy.sig)),
n1245 = length(intersect(intersect(intersect(enhancer.coding.sig, enhancer.only.sig), standard.sig),enhancer.methy.sig)),
n1345 = length(intersect(intersect(intersect(enhancer.coding.sig, TWAS.sig),standard.sig ),enhancer.methy.sig)),
n2345 = length(intersect(intersect(intersect(enhancer.only.sig, TWAS.sig), standard.sig),enhancer.methy.sig)),
n12345 = length(intersect(intersect(intersect(intersect(enhancer.coding.sig, enhancer.only.sig), TWAS.sig),standard.sig),enhancer.methy.sig)),
category = c("E+G", "G+Methyl", "TWAS", "STD", "E+G+Methyl"),
fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
cat.cex = 1,
margin = 0.1,
cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
ind = TRUE
);

dev.off()
