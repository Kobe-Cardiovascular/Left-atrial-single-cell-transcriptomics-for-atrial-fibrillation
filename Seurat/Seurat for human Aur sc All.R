# ライブラリの読み込み----------------------------------------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

# split the dataset into a list of two seurat objects (stim and CTRL)----
Aur.list<-list()

data.Aur1Endo <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur1/outs/per_sample_outs/Aur1_Endo/count/sample_feature_bc_matrix")
data.Aur1Epi <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur1/outs/per_sample_outs/Aur1_Epi/count/sample_feature_bc_matrix")
data.Aur2 <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur2_without_CellPlex/outs/filtered_feature_bc_matrix")
data.Aur3Endo <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur3/outs/per_sample_outs/Aur3_Endo/count/sample_feature_bc_matrix")
data.Aur3Epi <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur3/outs/per_sample_outs/Aur3_Epi/count/sample_feature_bc_matrix")
data.Aur4Endo <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur4/outs/per_sample_outs/Aur4_Endo/count/sample_feature_bc_matrix")
data.Aur4Epi <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur4/outs/per_sample_outs/Aur4_Epi/count/sample_feature_bc_matrix")
data.Aur3Blood <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur3/outs/per_sample_outs/Aur3_Blood/count/sample_feature_bc_matrix/")
data.Aur4Blood <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur4/outs/per_sample_outs/Aur4_Blood/count/sample_feature_bc_matrix/")
data.Aur12 <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur12/outs/filtered_feature_bc_matrix")
data.Aur13 <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur13/outs/filtered_feature_bc_matrix")

data.Aur10 <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur10/outs/filtered_feature_bc_matrix")
data.Aur16 <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur16/outs/filtered_feature_bc_matrix")
data.Aur19 <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur19/outs/filtered_feature_bc_matrix")


Aur.list$Aur1Endo<- CreateSeuratObject(counts = data.Aur1Endo$`Gene Expression`, project = "A", min.cells = 3, min.features = 200)
Aur.list$Aur1Endo <- AddMetaData(Aur.list$Aur1Endo , "Aur1Endo", col.name = "stim")
Aur.list$Aur1Epi <- CreateSeuratObject(counts = data.Aur1Epi$`Gene Expression`, project = "A",min.cells = 3, min.features = 200)
Aur.list$Aur1Epi <- AddMetaData(Aur.list$Aur1Epi , "Aur1Epi", col.name = "stim")
Aur.list$Aur2 <- CreateSeuratObject(counts = data.Aur2, project = "A",min.cells = 3, min.features = 200)
Aur.list$Aur2 <- AddMetaData(Aur.list$Aur2 , "Aur2", col.name = "stim")
Aur.list$Aur3Endo <- CreateSeuratObject(counts = data.Aur3Endo$`Gene Expression`, project = "A",min.cells = 3, min.features = 200)
Aur.list$Aur3Endo <- AddMetaData(Aur.list$Aur3Endo , "Aur3Endo", col.name = "stim")
Aur.list$Aur3Epi <- CreateSeuratObject(counts = data.Aur3Epi$`Gene Expression`, project = "A",min.cells = 3, min.features = 200)
Aur.list$Aur3Epi <- AddMetaData(Aur.list$Aur3Epi , "Aur3Epi", col.name = "stim")
Aur.list$Aur4Endo <- CreateSeuratObject(counts = data.Aur4Endo$`Gene Expression`, project = "A",min.cells = 3, min.features = 200)
Aur.list$Aur4Endo <- AddMetaData(Aur.list$Aur4Endo , "Aur4Endo", col.name = "stim")
Aur.list$Aur4Epi <- CreateSeuratObject(counts = data.Aur4Epi$`Gene Expression`, project = "A",min.cells = 3, min.features = 200)
Aur.list$Aur4Epi <- AddMetaData(Aur.list$Aur4Epi , "Aur4Epi", col.name = "stim")
Aur.list$Aur3Blood <- CreateSeuratObject(counts = data.Aur3Blood$`Gene Expression`, project = "A",min.cells = 3, min.features = 200)
Aur.list$Aur3Blood <- AddMetaData(Aur.list$Aur3Blood , "Aur3Blood", col.name = "stim")
Aur.list$Aur4Blood <- CreateSeuratObject(counts = data.Aur4Blood$`Gene Expression`, project = "A",min.cells = 3, min.features = 200)
Aur.list$Aur4Blood <- AddMetaData(Aur.list$Aur4Blood , "Aur4Blood", col.name = "stim")
Aur.list$Aur12 <- CreateSeuratObject(counts = data.Aur12, project = "A",min.cells = 3, min.features = 200)
Aur.list$Aur12 <- AddMetaData(Aur.list$Aur12 , "Aur12", col.name = "stim")
Aur.list$Aur13 <- CreateSeuratObject(counts = data.Aur13, project = "A",min.cells = 3, min.features = 200)
Aur.list$Aur13 <- AddMetaData(Aur.list$Aur13 , "Aur13", col.name = "stim")

Aur.list$Aur10 <- CreateSeuratObject(counts = data.Aur10, project = "A",min.cells = 3, min.features = 200)
Aur.list$Aur10 <- AddMetaData(Aur.list$Aur10 , "Aur10", col.name = "stim")
Aur.list$Aur16 <- CreateSeuratObject(counts = data.Aur16, project = "A",min.cells = 3, min.features = 200)
Aur.list$Aur16 <- AddMetaData(Aur.list$Aur16 , "Aur16", col.name = "stim")
Aur.list$Aur19 <- CreateSeuratObject(counts = data.Aur19, project = "A",min.cells = 3, min.features = 200)
Aur.list$Aur19 <- AddMetaData(Aur.list$Aur19 , "Aur19", col.name = "stim")


# low-quality cellの確認①------------------------------------------------------
# MT-から始まるミトコンドリアRNAを"percent.mt"列としてデータに追加

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats

Aur.list$Aur1Endo[["percent.mt"]] <- PercentageFeatureSet(Aur.list$Aur1Endo, pattern = "^MT-")
Aur.list$Aur1Epi[["percent.mt"]] <- PercentageFeatureSet(Aur.list$Aur1Epi, pattern = "^MT-")
Aur.list$Aur2[["percent.mt"]] <- PercentageFeatureSet(Aur.list$Aur2, pattern = "^MT-")
Aur.list$Aur3Endo[["percent.mt"]] <- PercentageFeatureSet(Aur.list$Aur3Endo, pattern = "^MT-")
Aur.list$Aur3Epi[["percent.mt"]] <- PercentageFeatureSet(Aur.list$Aur3Epi, pattern = "^MT-")
Aur.list$Aur4Endo[["percent.mt"]] <- PercentageFeatureSet(Aur.list$Aur4Endo, pattern = "^MT-")
Aur.list$Aur4Epi[["percent.mt"]] <- PercentageFeatureSet(Aur.list$Aur4Epi, pattern = "^MT-")
Aur.list$Aur3Blood[["percent.mt"]] <- PercentageFeatureSet(Aur.list$Aur3Blood, pattern = "^MT-")
Aur.list$Aur4Blood[["percent.mt"]] <- PercentageFeatureSet(Aur.list$Aur4Blood, pattern = "^MT-")
Aur.list$Aur12[["percent.mt"]] <- PercentageFeatureSet(Aur.list$Aur12, pattern = "^MT-")
Aur.list$Aur13[["percent.mt"]] <- PercentageFeatureSet(Aur.list$Aur13, pattern = "^MT-")

Aur.list$Aur10[["percent.mt"]] <- PercentageFeatureSet(Aur.list$Aur10, pattern = "^MT-")
Aur.list$Aur16[["percent.mt"]] <- PercentageFeatureSet(Aur.list$Aur16, pattern = "^MT-")
Aur.list$Aur19[["percent.mt"]] <- PercentageFeatureSet(Aur.list$Aur19, pattern = "^MT-")

VlnPlot(Aur.list$Aur1Endo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Aur.list$Aur1Epi, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Aur.list$Aur2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Aur.list$Aur3Endo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Aur.list$Aur3Epi, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Aur.list$Aur4Endo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Aur.list$Aur4Epi, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Aur.list$Aur3Blood, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Aur.list$Aur4Blood, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Aur.list$Aur12, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Aur.list$Aur13, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

VlnPlot(Aur.list$Aur10, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Aur.list$Aur16, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Aur.list$Aur19, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)



# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

# percent.mtが極端に大きい（10%以上）、またはnFeature_RNAが小さすぎる（200以下）・大きすぎる(2500以上)細胞は、
# 死細胞やdoubletの可能性が高い


plot1 <- FeatureScatter(Aur.list$Aur1Endo, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Aur.list$Aur1Endo, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(Aur.list$Aur1Epi, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Aur.list$Aur1Epi, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(Aur.list$Aur2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Aur.list$Aur2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(Aur.list$Aur3Endo, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Aur.list$Aur3Endo, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(Aur.list$Aur3Epi, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Aur.list$Aur3Epi, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(Aur.list$Aur4Endo, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Aur.list$Aur4Endo, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(Aur.list$Aur4Epi, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Aur.list$Aur4Epi, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(Aur.list$Aur3Blood, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Aur.list$Aur3Blood, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(Aur.list$Aur4Blood, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Aur.list$Aur4Blood, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(Aur.list$Aur12, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Aur.list$Aur12, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(Aur.list$Aur13, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Aur.list$Aur13, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(Aur.list$Aur10, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Aur.list$Aur10, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(Aur.list$Aur16, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Aur.list$Aur16, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(Aur.list$Aur19, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Aur.list$Aur19, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# low-quality cellのフィルタリング---------------------------------------------
# ミトコンドリアRNAが７以下、遺伝子数が200~2500の間の細胞のみ抽出


Aur.list$Aur1Endo <- subset(Aur.list$Aur1Endo, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 8)
Aur.list$Aur1Epi <- subset(Aur.list$Aur1Epi, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 8)
Aur.list$Aur2 <- subset(Aur.list$Aur2, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 8)
Aur.list$Aur3Endo <- subset(Aur.list$Aur3Endo, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 8)
Aur.list$Aur3Epi <- subset(Aur.list$Aur3Epi, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 8)
Aur.list$Aur4Endo <- subset(Aur.list$Aur4Endo, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 8)
Aur.list$Aur4Epi <- subset(Aur.list$Aur4Epi, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 8)
Aur.list$Aur3Blood <- subset(Aur.list$Aur3Blood, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 8)
Aur.list$Aur4Blood <- subset(Aur.list$Aur4Blood, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 8)
Aur.list$Aur12 <- subset(Aur.list$Aur12, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 8)
Aur.list$Aur13 <- subset(Aur.list$Aur13, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 8)

Aur.list$Aur10 <- subset(Aur.list$Aur10, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 8)
Aur.list$Aur16 <- subset(Aur.list$Aur16, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 8)
Aur.list$Aur19 <- subset(Aur.list$Aur19, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 8)


# 先程と比べlow-qualityデータがフィルターされていることを確認


plot1 <- FeatureScatter(Aur.list$Aur1Endo, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Aur.list$Aur1Epi, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(Aur.list$Aur2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Aur.list$Aur2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


plot1 <- FeatureScatter(Aur.list$Aur3Endo, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Aur.list$Aur3Epi, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(Aur.list$Aur4Endo, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Aur.list$Aur4Epi, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(Aur.list$Aur3Blood, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Aur.list$Aur3Blood, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(Aur.list$Aur4Blood, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Aur.list$Aur4Blood, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(Aur.list$Aur12, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Aur.list$Aur12, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(Aur.list$Aur13, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Aur.list$Aur13, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(Aur.list$Aur10, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Aur.list$Aur10, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(Aur.list$Aur16, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Aur.list$Aur16, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(Aur.list$Aur19, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Aur.list$Aur19, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


# データの正規化（normalization）-----------------------------------------------

Aur.list$Aur1Endo<- NormalizeData(Aur.list$Aur1Endo, normalization.method = "LogNormalize", scale.factor = 10000)
Aur.list$Aur1Epi <- NormalizeData(Aur.list$Aur1Epi, normalization.method = "LogNormalize", scale.factor = 10000)
Aur.list$Aur2<- NormalizeData(Aur.list$Aur2, normalization.method = "LogNormalize", scale.factor = 10000)
Aur.list$Aur3Endo<- NormalizeData(Aur.list$Aur3Endo, normalization.method = "LogNormalize", scale.factor = 10000)
Aur.list$Aur3Epi <- NormalizeData(Aur.list$Aur3Epi, normalization.method = "LogNormalize", scale.factor = 10000)
Aur.list$Aur4Endo<- NormalizeData(Aur.list$Aur4Endo, normalization.method = "LogNormalize", scale.factor = 10000)
Aur.list$Aur4Epi <- NormalizeData(Aur.list$Aur4Epi, normalization.method = "LogNormalize", scale.factor = 10000)
Aur.list$Aur3Blood <- NormalizeData(Aur.list$Aur3Blood, normalization.method = "LogNormalize", scale.factor = 10000)
Aur.list$Aur4Blood <- NormalizeData(Aur.list$Aur4Blood, normalization.method = "LogNormalize", scale.factor = 10000)
Aur.list$Aur12<- NormalizeData(Aur.list$Aur12, normalization.method = "LogNormalize", scale.factor = 10000)
Aur.list$Aur13<- NormalizeData(Aur.list$Aur13, normalization.method = "LogNormalize", scale.factor = 10000)

Aur.list$Aur10<- NormalizeData(Aur.list$Aur10, normalization.method = "LogNormalize", scale.factor = 10000)
Aur.list$Aur16<- NormalizeData(Aur.list$Aur16, normalization.method = "LogNormalize", scale.factor = 10000)
Aur.list$Aur19<- NormalizeData(Aur.list$Aur19, normalization.method = "LogNormalize", scale.factor = 10000)


all.genes <- rownames(Aur.list$Aur1Endo)
Aur.list$Aur1Endo <- ScaleData(Aur.list$Aur1Endo, features = all.genes)
all.genes <- rownames(Aur.list$Aur1Epi)
Aur.list$Aur1Epi <- ScaleData(Aur.list$Aur1Epi, features = all.genes)
all.genes <- rownames(Aur.list$Aur2)
Aur.list$Aur2 <- ScaleData(Aur.list$Aur2, features = all.genes)
all.genes <- rownames(Aur.list$Aur3Endo)
Aur.list$Aur3Endo <- ScaleData(Aur.list$Aur3Endo, features = all.genes)
all.genes <- rownames(Aur.list$Aur3Epi)
Aur.list$Aur3Epi <- ScaleData(Aur.list$Aur3Epi, features = all.genes)
all.genes <- rownames(Aur.list$Aur4Endo)
Aur.list$Aur4Endo <- ScaleData(Aur.list$Aur4Endo, features = all.genes)
all.genes <- rownames(Aur.list$Aur4Epi)
Aur.list$Aur4Epi <- ScaleData(Aur.list$Aur4Epi, features = all.genes)
all.genes <- rownames(Aur.list$Aur3Blood)
Aur.list$Aur3Blood <- ScaleData(Aur.list$Aur3Blood, features = all.genes)
all.genes <- rownames(Aur.list$Aur4Blood)
Aur.list$Aur4Blood <- ScaleData(Aur.list$Aur4Blood, features = all.genes)
all.genes <- rownames(Aur.list$Aur12)
Aur.list$Aur12 <- ScaleData(Aur.list$Aur12, features = all.genes)
all.genes <- rownames(Aur.list$Aur13)
Aur.list$Aur13 <- ScaleData(Aur.list$Aur13, features = all.genes)


all.genes <- rownames(Aur.list$Aur10)
Aur.list$Aur10 <- ScaleData(Aur.list$Aur10, features = all.genes)
all.genes <- rownames(Aur.list$Aur16)
Aur.list$Aur16 <- ScaleData(Aur.list$Aur16, features = all.genes)
all.genes <- rownames(Aur.list$Aur13)
Aur.list$Aur19 <- ScaleData(Aur.list$Aur19, features = all.genes)


# UMAPによるクラスタリング　2つ合わせて------------------------------------------------------------------------------

# normalize and identify variable features for each dataset independently
# data scaling変更231123

all.genes <- rownames(Aur.list["Aur1Endo"])+   rownames(Aur.list["Aur1Epi"])+   rownames(Aur.list["Aur2"])+   rownames(Aur.list["Aur3Endo"])+   rownames(Aur.list["Aur3Epi"])+   rownames(Aur.list["Aur4Endo"])+   rownames(Aur.list["Aur4Epi"])+   rownames(Aur.list["Aur3Blood"])+   rownames(Aur.list["Aur4Blood"])+   rownames(Aur.list["Aur12"])+   rownames(Aur.list["Aur13"])+   rownames(Aur.list["Aur10"])+   rownames(Aur.list["Aur16"])+   rownames(Aur.list["Aur19"])



Aur.list <- lapply(X = Aur.list, FUN = function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
  x <- ScaleData(x, features = all.genes)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- RunPCA(x, features = VariableFeatures(object = x), verbose = FALSE)
})


# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = Aur.list)


# Integrationの実行
Aur.anchors <- FindIntegrationAnchors(object.list = Aur.list, anchor.features = features, reduction = "rpca")

# this command creates an 'integrated' data assay
Aur.combined <- IntegrateData(anchorset =Aur.anchors)

# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(Aur.combined) <- "integrated"

Aur.combined <- ScaleData(Aur.combined, verbose = FALSE)
Aur.combined <- RunPCA(Aur.combined, npcs = 30, verbose = FALSE)
Aur.combined <- RunUMAP(Aur.combined, reduction = "pca", dims = 1:30)
Aur.combined <- FindNeighbors(Aur.combined, reduction = "pca", dims = 1:30)
Aur.combined <- FindClusters(Aur.combined, resolution = 0.10)

DimPlot(Aur.combined, reduction = "umap", group.by = "stim")
DimPlot(Aur.combined, reduction = "umap", split.by = "stim")
DimPlot(Aur.combined, reduction = "umap", label = TRUE, repel = TRUE)

DimPlot(Aur.combined, reduction = "umap", split.by = "AurBlood", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined, reduction = "umap", split.by = "sample1ex")
DimPlot(Aur.combined, reduction = "umap", split.by = "Auronly1ex")

###最初　first UMAP
saveRDS(Aur.combined, "~/Desktop/Human_LA_LAA/RDS/Aur_cluster_scAll_first_UMAP_231123.rds")
Aur.combined<- readRDS ("~/Desktop/Human_LA_LAA/RDS/Aur_cluster_scAll_first_UMAP_231123.rds")




#グループ分けAurBlood
# stim列を抽出する
cell.stim <- Aur.combined@meta.data[["stim"]]
# 新しいグループ（疾患ありなしの2群）とstim列との対応関係を定義する
newgroup <- c("Aur","Aur","Aur","Aur","Aur","Aur","Aur","Blood","Blood","Aur","Aur","Aur","Aur","Aur")
names(newgroup) <- c("Aur1Endo","Aur1Epi","Aur2","Aur3Endo","Aur3Epi","Aur4Endo","Aur4Epi","Aur3Blood","Aur4Blood","Aur12","Aur13","Aur10","Aur16","Aur19")
# 抽出したstim列を新しいグループに変換する
cell.newgroup <- newgroup[cell.stim]
# 新しいグループを新しいMetaDataとして追加する
Aur.combined <- AddMetaData(Aur.combined, cell.newgroup, col.name = "AurBlood")

### AFDuration ####
# stim列を抽出する
cell.stim <- Aur.combined@meta.data[["stim"]]
# 新しいグループ（疾患ありなしの2群）とstim列との対応関係を定義する
newgroup <- c("1PerAF","1PerAF","1PerAF","2LSPerAF","2LSPerAF","2LSPerAF","2LSPerAF","3Blood","3Blood","2LSPerAF","1PerAF","2LSPerAF","1PerAF","2LSPerAF")
names(newgroup) <- c("Aur1Endo","Aur1Epi","Aur2","Aur3Endo","Aur3Epi","Aur4Endo","Aur4Epi","Aur3Blood","Aur4Blood","Aur12","Aur13","Aur10","Aur16","Aur19")
# 抽出したstim列を新しいグループに変換する
cell.newgroup <- newgroup[cell.stim]
# 新しいグループを新しいMetaDataとして追加する
Aur.combined <- AddMetaData(Aur.combined, cell.newgroup, col.name = "AFDuration")


#グループ分けAuronly1ex
# stim列を抽出する
cell.stim <- Aur.combined@meta.data[["stim"]]
# 新しいグループ（疾患ありなしの2群）とstim列との対応関係を定義する
newgroup <- c("Other","Other","Aur","Aur","Aur","Aur","Aur","Other","Other","Aur","Aur","Aur","Aur","Aur")
names(newgroup) <- c("Aur1Endo","Aur1Epi","Aur2","Aur3Endo","Aur3Epi","Aur4Endo","Aur4Epi","Aur3Blood","Aur4Blood","Aur12","Aur13","Aur10","Aur16","Aur19")
# 抽出したstim列を新しいグループに変換する
cell.newgroup <- newgroup[cell.stim]
# 新しいグループを新しいMetaDataとして追加する
Aur.combined <- AddMetaData(Aur.combined, cell.newgroup, col.name = "Auronly1ex")


#グループ分けsample1ex
# stim列を抽出する
cell.stim <- Aur.combined@meta.data[["stim"]]
# 新しいグループ（疾患ありなしの2群）とstim列との対応関係を定義する
newgroup <- c("Other","Other","1.2","4.3","4.3","5.4","5.4","Other","Other","6.12","2.13","7.10","3.16","8.19")
names(newgroup) <- c("Aur1Endo","Aur1Epi","Aur2","Aur3Endo","Aur3Epi","Aur4Endo","Aur4Epi","Aur3Blood","Aur4Blood","Aur12","Aur13","Aur10","Aur16","Aur19")
# 抽出したstim列を新しいグループに変換する
cell.newgroup <- newgroup[cell.stim]
# 新しいグループを新しいMetaDataとして追加する
Aur.combined <- AddMetaData(Aur.combined, cell.newgroup, col.name = "sample1ex")


Aur.combined_trans <- SCTransform (Aur.combined)

cd_genes <- c("CD3E","NKG7","CD8A","CD68","CD14","S100A8","HLA-DQA1","CD79A","CD79B","GNLY","TYROBP","KIT", "MYH11","PECAM1","CDH5","ACTC1","TNNC1","MKI67")
DotPlot(Aur.combined_trans,features = cd_genes)+RotatedAxis()+coord_flip()



###最初　first UMAP
saveRDS(Aur.combined, "~/Desktop/Human_LA_LAA/RDS/Aur_cluster_scAll_first_UMAP_231123.rds")
Aur.combined<- readRDS ("~/Desktop/Human_LA_LAA/RDS/Aur_cluster_scAll_first_UMAP_231123.rds")

 
#CD45----
# 新しいグループ分けで可視化
DimPlot(Aur.combined, reduction = "umap", label = FALSE, repel = TRUE)
DimPlot(Aur.combined, reduction = "umap", label = TRUE,  repel = TRUE)

DimPlot(Aur.combined, reduction = "umap", group.by = "AurBlood")
DimPlot(Aur.combined, reduction = "umap", split.by = "AurBlood", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined, reduction = "umap", split.by = "AurBlood", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined, reduction = "umap", split.by = "stim",　label = TRUE,  repel = TRUE)


# Myeloid----
VlnPlot(Aur.combined, pt.size = 0, features = c("CD68","CD14","CSF1R","IL1B","S100A8", "S100A9") )
VlnPlot(Aur.combined, pt.size = 0, features = c("TREM2","CD9","APOE","CXCL3","TNF","CXCL8") )
VlnPlot(Aur.combined, pt.size = 0, features = c("HLA-DPA1","HLA-DPB1","HLA-DQA1","CLEC10A","FCER1A","CD1C") )
VlnPlot(Aur.combined, pt.size = 0, features = c("PTPRC","HLA-DPB1","HLA-DQA1","CLEC10A","FCER1A","CD1C") )

# Bcells-----
VlnPlot(Aur.combined, pt.size = 0, features = c("CD79A","CD79B","FCER2","CD22","MYH11","CD34") )
VlnPlot(Aur.combined, pt.size = 0, features = c("LUM","TAGLN","CD34","CDH5","PECAM1","TYROBP") )

# Tcells-----
VlnPlot(Aur.combined, pt.size = 0, features = c("CD3E","CD4","CD8A","NKG7","GNLY","GZMB") )
VlnPlot(Aur.combined, pt.size = 0, features = c("KLRB1","CCR6","CD1d1","MKI67","IL2RB","GZMB") )


VlnPlot(Aur.combined, features = c("CLEC10A","FCER1A","CD1D","FCGR3A") )

VlnPlot(Aur.combined, pt.size = 0 , features = c("ACTC1","TNNC1","FCER2","CD22","MYH11","CD34") )
VlnPlot(Aur.combined, pt.size = 0 , features = c("LUM","TAGLN","CD34","CDH5","PECAM1","TYROBP") )

VlnPlot(Aur.combined, pt.size = 0 , features = c("CD68","CD14","IL1B","S100A8","VCAN","CD3E") )
VlnPlot(Aur.combined, pt.size = 0 ,  features = c("CD4","CD8A","NKG7","GNLY","CD79A","CD79B") )
VlnPlot(Aur.combined, pt.size = 0 ,  features = c("FCER2","KIT","HDC","MKI67") )

cDC1
VlnPlot(Aur.combined, pt.size = 0, features = c("THBD","CLEC9A","XCR1","CADM1","IRF4","IRF8","BATF3") )
cDC2
VlnPlot(Aur.combined, pt.size = 0, features = c("THBD","CD1C","CLEC10A","FCER1A","IRF4","CD2") )
pDC
VlnPlot(Aur.combined, pt.size = 0, features = c("IL3RA","CLEC4C") )


FeaturePlot(Aur.combined, features=c('ACTC1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined, features=c('LUM'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined, features=c('CD34'), min.cutoff=0.5, max.cutoff='q90')



Aur.combined_trans <- SCTransform (Aur.combined)

cd_genes <- c("CD3E","NKG7","CD8A","CD68","CD14","S100A8","HLA-DQA1","CD79A","CD79B","GNLY","TYROBP","KIT", "MYH11","PECAM1","CDH5","ACTC1","TNNC1","MKI67")
DotPlot(Aur.combined_trans,features = cd_genes)+RotatedAxis()+coord_flip()



#CD45(Endothelial,SMC,myocardium除去)
Aur.combined_CD45 <- subset(Aur.combined, idents = c(0,1,2,3,4,5,6,8,9,10,11,12,14,16))
Aur.combined_CD45 <- ScaleData(Aur.combined_CD45, verbose = FALSE)
Aur.combined_CD45<- RunPCA(Aur.combined_CD45, npcs = 30, verbose = FALSE)
Aur.combined_CD45 <- RunUMAP(Aur.combined_CD45 , reduction = "pca", dims = 1:30)
Aur.combined_CD45 <- FindNeighbors(Aur.combined_CD45 , reduction = "pca", dims = 1:30)
Aur.combined_CD45 <- FindClusters(Aur.combined_CD45, resolution = 0.10)

DimPlot(Aur.combined_CD45, reduction = "umap", label = FALSE, repel = TRUE)
DimPlot(Aur.combined_CD45, reduction = "umap", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined_CD45, reduction = "umap", split.by = "AurBlood", label = FALSE,  repel = FALSE)
DimPlot(Aur.combined_CD45, reduction = "umap", split.by = "AurBlood", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined_CD45, reduction = "umap", split.by = "AFDuration")
DimPlot(Aur.combined_CD45, reduction = "umap", split.by = "sample1ex")
DimPlot(Aur.combined_CD45, reduction = "umap", split.by = "stim")
DimPlot(Aur.combined_CD45, reduction = "umap", split.by = "Auronly1ex")



#並べ替え
levels(Aur.combined_CD45) # [1] "0""1""2""3""4""5""6""7"8""9""10""11""12""13"
levels(Aur.combined_CD45) <- c("1", "4", "5","10","11","12", "0", "2", "6", "7","9","3","8","13")
new.cluster.ids <- c("0","1","2","3","4", "5", "6","7","8","9","10","11","12","13")
names(new.cluster.ids) <- levels(Aur.combined_CD45)
Aur.combined_CD45<- RenameIdents(Aur.combined_CD45, new.cluster.ids)




#CD45(Endothelial,SMC,myocardium除去)----並べ替え後final
saveRDS(Aur.combined_CD45, "~/Desktop/Human_LA_LAA/RDS/Aur_cluster_scCD45_final.rds")
Aur.combined_CD45<- readRDS ("~/Desktop/Human_LA_LAA/RDS/Aur_cluster_scCD45_final.rds")


# Olink----
VlnPlot(Aur.combined_CD45, pt.size = 0, features = c("IL17A","KCND2","RASGRF1","B4GALNT3","MTERF3", "FABP6") )
VlnPlot(Aur.combined_CD45, pt.size = 0, features = c("RNF4","ARNT","BIRC3","B4GALNT3","MTERF3", "FABP6") )
VlnPlot(Aur.combined_CD45, pt.size = 0, features = c("IL17A","IL17C","CXCL11","MMP12","MTERF3", "FABP6") )
# Myeloid----
VlnPlot(Aur.combined_CD45, pt.size = 0, features = c("CD68","CD14","CSF1R","IL1B","S100A8", "S100A9") )
VlnPlot(Aur.combined_CD45, pt.size = 0, features = c("TREM2","CD9","APOE","CXCL3","TNF","CXCL8") )
VlnPlot(Aur.combined_CD45, pt.size = 0, features = c("HLA-DPA1","HLA-DPB1","HLA-DQA1","CLEC10A","FCER1A","CD1C") )
VlnPlot(Aur.combined_CD45, pt.size = 0, features = c("PTPRC","HLA-DPB1","HLA-DQA1","CLEC10A","FCER1A","CD1C") )

# Bcells-----
VlnPlot(Aur.combined_CD45, pt.size = 0, features = c("CD79A","CD79B","FCER2","CD22","MYH11","CD34") )
VlnPlot(Aur.combined_CD45, pt.size = 0, features = c("LUM","TAGLN","CD34","CDH5","PECAM1","TYROBP") )

# Tcells-----
VlnPlot(Aur.combined_CD45, pt.size = 0, features = c("CD3E","CD4","CD8A","NKG7","GNLY","GZMB") )
VlnPlot(Aur.combined_CD45, pt.size = 0, features = c("KLRB1","CCR6","CD1d1","MKI67","IL2RB","GZMB") )


VlnPlot(Aur.combined_CD45, pt.size = 0, features = c("CLEC10A","FCER1A","CD1D","FCGR3A","CXCR2") )

VlnPlot(Aur.combined_CD45, pt.size = 0, features = c("ACTC1","TNNC1","FCER2","CD22","MYH11","CD34") )
VlnPlot(Aur.combined_CD45, pt.size = 0, features = c("LUM","TAGLN","CD34","CDH5","PECAM1","TYROBP") )

VlnPlot(Aur.combined_CD45, pt.size = 0 ,  features = c("FCER2","KIT","HDC","MKI67") )

cDC1
VlnPlot(Aur.combined_CD45, pt.size = 0 , features = c("THBD","CLEC9A","XCR1","CADM1","IRF4","IRF8","BATF3") )
cDC2
VlnPlot(Aur.combined_CD45, pt.size = 0 , features = c("THBD","CD1C","CLEC10A","FCER1A","IRF4","CD2") )
pDC
VlnPlot(Aur.combined_CD45, pt.size = 0 , features = c("IL3RA","CLEC4C") )

VlnPlot(Aur.combined_CD45,pt.size = 0, features = c("LYVE1","HLA-DOA","HLA-DQA1","HLA-DQA2","HLA-DQB1","TIMD4","TREM2") )

Manuscript
VlnPlot(Aur.combined_CD45,pt.size = 0, features = c("CD68","CD14","IL1B","S100A8","VCAN","HLA-DPA1") )
VlnPlot(Aur.combined_CD45,pt.size = 0, features = c("CD3E","NKG7","GNLY","CD79A","KIT","MKI67") )

VlnPlot(Aur.combined_CD45,pt.size = 0, features = c("CD68","CD14","CSF1R","CD163","LYVE1","CD9","TREM2","IL1B","S100A8","VCAN","LYZ","FCGR3A","CXCR2","FCGR3B","CSF3R","HLA-DPA1","CLEC9A","CD1C","IL3RA"))


Aur.combined_CD45_trans <- SCTransform (Aur.combined_CD45)

cd_genes <- c("CD68","CD14","VCAN","CXCR2","CSF1R","KIT","HLA-DQA2","CLEC9A","CLEC4C","CD3E","NKG7","GNLY","MKI67","CD79A","FCER2")
DotPlot(Aur.combined_CD45_trans,features = cd_genes)+RotatedAxis()+coord_flip()






Aur.combined_CD45_prop. <- prop.table(table(Idents(Aur.combined_CD45),Aur.combined_CD45$Auronly1ex))
write.table(Aur.combined_CD45_prop., "Aur.combined_CD45_prop.Auronly1ex.txt", quote=F, col.names=F, append=T)
write.csv(Aur.combined_CD45_prop., "Aur.combined_CD45_prop.Auronly1ex.csv")



Aur.combined_CD45_prop. <- prop.table(table(Idents(Aur.combined_CD45),Aur.combined_CD45$sample1ex))
write.table(Aur.combined_CD45_prop., "Aur.combined_CD45_prop.sample1ex.txt", quote=F, col.names=F, append=T)
write.csv(Aur.combined_CD45_prop., "Aur.combined_CD45_prop.sample1ex.csv")


# ヒートマップの作成------------------------------------------------------------------------------
Aur.combined_CD45.markers <- FindAllMarkers(Aur.combined_CD45, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Aur.combined_CD45.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- Aur.combined_CD45.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Aur.combined_CD45, features = top10$gene) + NoLegend()
write.table(top10, "table.txt", quote=F, col.names=F, append=T)


#myeloid-------
Aur.combined_1 <- subset(Aur.combined_CD45, idents = c(0,1,2,3,4,5))
Aur.combined_1 <- ScaleData(Aur.combined_1, verbose = FALSE)
Aur.combined_1<- RunPCA(Aur.combined_1, npcs = 30, verbose = FALSE)
Aur.combined_1 <- RunUMAP(Aur.combined_1 , reduction = "pca", dims = 1:30)
Aur.combined_1 <- FindNeighbors(Aur.combined_1 , reduction = "pca", dims = 1:30)
Aur.combined_1 <- FindClusters(Aur.combined_1, resolution = 0.25)


DimPlot(Aur.combined_1, reduction = "umap", label = FALSE, repel = TRUE)
DimPlot(Aur.combined_1, reduction = "umap", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined_1, reduction = "umap", split.by = "AurBlood", label = FALSE,  repel = FALSE)
DimPlot(Aur.combined_1, reduction = "umap", split.by = "AurBlood", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined_1, reduction = "umap", split.by = "AFDuration")
DimPlot(Aur.combined_1, reduction = "umap", split.by = "sample1ex")
DimPlot(Aur.combined_1, reduction = "umap", split.by = "stim")
DimPlot(Aur.combined_1, reduction = "umap", split.by = "Auronly1ex")

VlnPlot(Aur.combined_1, pt.size = 0 , features = c("FCGR3A","FCGR2A","ITGB2","CD16","CEACAM8","CD55","CD44","ELANE") )
VlnPlot(Aur.combined_1, pt.size = 0 , features = c("PTPRC","CD14","MPO","LYZ","ITGB2","CD16","FCGR3A","FCGR3B","CD55","CD44","ELANE","VCAN", "CD52","FCN1","S100A12","IL1B","CSF1R","S100A8","S100A9","CAMP","DEFB4A","CTSS","ITGAM","CD68","CYBB") )
VlnPlot(Aur.combined_1, pt.size = 0 , features = c("DEFB1","DEFA1","DEFB103B","S100A9","LYZ") )
VlnPlot(Aur.combined_1, pt.size = 0 , features = c("S100A6","FOS","S100A12","S100A9","HMGB2","KLF6"))



VlnPlot(Aur.combined_1, pt.size = 0 , features = c("HLA-DPA1","HLA-DPB1","HLA-DQA1","CLEC10A","FCER1A","CD1C") )
VlnPlot(Aur.combined_1, pt.size = 0 , features = c("TLR4","TLR2","TLR7","TLR3","TLR12", "TLR9") )
VlnPlot(Aur.combined_1, pt.size = 0 , features = c("HDC","KIT","TPSAB1","CCL5") )
VlnPlot(Aur.combined_1, pt.size = 0 , features = c("IL1B","TNF","CXCL8","NR1H3", "CEBPA","CEBPB") )
VlnPlot(Aur.combined_1, pt.size = 0 , features = c("CD68","CSF1R","ITGAM","S100A8","S100A9","S100A12") )
VlnPlot(Aur.combined_1, pt.size = 0 , features = c("FCGR3A","CD14","LYZ","VCAN","FCN1","CTSS") )
VlnPlot(Aur.combined_1, pt.size = 0 , features = c("CASP1","CASP4","DEFB4A","CYBB","IL1B","VEGFA") )
VlnPlot(Aur.combined_1, pt.size = 0 , features = c("CCL2","CCL20","CCL4","HMGB2","KLF6","VEGFA") )


VlnPlot(Aur.combined_1, pt.size = 0 , features = c("IL1B","TNF","CXCL8","CXCL3", "ABCG1","ABCA1") )
VlnPlot(Aur.combined_1, pt.size = 0 , features = c("CD9","TREM2","TGFB1","APOC1","APOE","C1QA", "C1QB") )
VlnPlot(Aur.combined_1, pt.size = 0 , features = c("C1QC","CD59","APOE","","NUPR1", "C1QB") )

VlnPlot(Aur.combined_1, pt.size = 0 , features = c("CD3E","CD4","CD8A","CCR2") )
VlnPlot(Aur.combined_1,pt.size = 0, features = c("CD3E","CD4","CD8A","IL7R","CCR7"))
FeaturePlot(Aur.combined_1, features=c('CD3E'), min.cutoff=0.5, max.cutoff='q90')


Aur.combined_1_trans <- SCTransform (Aur.combined_1)

cd_genes <- c("CD68","CD14","CSF1R","LYVE1","TREM2","CD9","IL1B","TNF","S100A8","VCAN","FCN1","FCGR3A","CXCR2","HLA-DPA1","CLEC9A","FCER1A")
DotPlot(Aur.combined_1_trans,features = cd_genes)+RotatedAxis()+coord_flip()

cd_genes <- c("CD68","CD14","CSF1R","CD163","LYVE1","FOLR2","C1QA","IL1B","TLR4","AREG","PDGFB","TREM2","TIMP1","CLEC5A","EREG","IL10","VEGFA","S100A8","VCAN","FCGR3A","CXCR2","FCGR3B","CSF3R","HLA-DPA1","HLA-DRA","CD1C","CLEC9A","IL3RA","KIT","CCR3","IL4")
DotPlot(Aur.combined_1_trans,features = cd_genes)+RotatedAxis()+coord_flip()


VlnPlot(Aur.combined_1, pt.size = 0 , features = c("CD3E","CD4","CD8A","CCR2") )
VlnPlot(Aur.combined_1, pt.size = 0 , features = c("CCL2","CCL20","CCL4","CCL5","KLF6","VEGFA") )
VlnPlot(Aur.combined_1,pt.size = 0, features = c("CD3E","NKG7","GNLY","CD79A","KIT","MKI67") )
VlnPlot(Aur.combined_1, pt.size = 0, features = c("ACTC1","TNNC1","FCER2","CD22","MYH11","CD34") )
VlnPlot(Aur.combined_1, pt.size = 0, features = c("LUM","TAGLN","CD34","CDH5","PECAM1","TYROBP") )
VlnPlot(Aur.combined_1, pt.size = 0, features = c("IFIT1","IFIT2","IFIT3","CD163","TREM2","LYVE1") )

# ヒートマップの作成------------------------------------------------------------------------------
Aur.combined_1.markers <- FindAllMarkers(Aur.combined_1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Aur.combined_1.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- Aur.combined_1.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Aur.combined_1, features = top10$gene) + NoLegend()
write.table(top10, "table.txt", quote=F, col.names=F, append=T)


Aur.combined_1_trans <- SCTransform (Aur.combined_1)

cd_genes <- c("CD68","CD14","CSF1R","CD163","LYVE1","CD9","TREM2","IL1B","S100A8","VCAN","LYZ","FCGR3A","CXCR2","FCGR3B","CSF3R","HLA-DPA1","CLEC9A","CD1C","IL3RA")
DotPlot(Aur.combined_1_trans,features = cd_genes)+RotatedAxis()+coord_flip()

# RDS仮
saveRDS(Aur.combined_1,  "~/Desktop/Human_LA_LAA/RDS/Aur_cluster_scMyeloid_first.rds")
Aur.combined_1 <- readRDS ("~/Desktop/Human_LA_LAA/RDS/Aur_cluster_scMyeloid_first.rds")



#(Myeloid cell Tcell doublet除去)--------
Aur.combined_1_1 <- subset(Aur.combined_1, idents = c(0,1,2,4,5,6,7,8,9,10,11,12,14))
Aur.combined_1_1 <- ScaleData(Aur.combined_1_1, verbose = FALSE)
Aur.combined_1_1<- RunPCA(Aur.combined_1_1, npcs = 30, verbose = FALSE)
Aur.combined_1_1 <- RunUMAP(Aur.combined_1_1 , reduction = "pca", dims = 1:30)
Aur.combined_1_1 <- FindNeighbors(Aur.combined_1_1 , reduction = "pca", dims = 1:30)
Aur.combined_1_1 <- FindClusters(Aur.combined_1_1, resolution = 0.5)





DimPlot(Aur.combined_1_1, reduction = "umap", label = FALSE, repel = TRUE)
DimPlot(Aur.combined_1_1, reduction = "umap", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined_1_1, reduction = "umap", split.by = "AurBlood", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined_1_1, reduction = "umap", split.by = "AFDuration")
DimPlot(Aur.combined_1_1, reduction = "umap", split.by = "stim")







# RDS Myeloid T cell doublet除去後
saveRDS(Aur.combined_1_1, "~/Desktop/Human_LA_LAA/RDS/Aur_cluster_scMyeloid_Tdoublet除去後.rds")
Aur.combined_1_1 <- readRDS ("~/Desktop/Human_LA_LAA/RDS/Aur_cluster_scMyeloid_Tdoublet除去後.rds")


VlnPlot(Aur.combined_1_1, pt.size = 0 , features = c("ACTC1","TNNC1","MYL7","NPPA","MYH11") )

FeaturePlot(Aur.combined_1_1, features=c('MYL7'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_1, features=c('CD3E'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_1, features=c('IL7R'), min.cutoff=0.1, max.cutoff='q90')
FeaturePlot(Aur.combined_1_1, features=c('CCR7'), min.cutoff=0.01, max.cutoff='q90')
FeaturePlot(Aur.combined_1_1, features=c('HLA-DRA'), min.cutoff=4, max.cutoff='q90')
FeaturePlot(Aur.combined_1_1, features=c('TNF'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_1, features=c('CCL3'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_1, features=c('CD3E'), min.cutoff=0.5, max.cutoff='q90')

VlnPlot(Aur.combined_1_1, pt.size = 0, features = c("CD3E","CD4","CD8A","NKG7","GNLY","GZMB") )
VlnPlot(Aur.combined_1_1, pt.size = 0, features = c("ACTC1","TNNC1","FCER2","CD22","MYH11","CD34") )
VlnPlot(Aur.combined_1_1, pt.size = 0, features = c("LUM","TAGLN","CD34","CDH5","PECAM1","TYROBP") )
VlnPlot(Aur.combined_1_1, pt.size = 0, features = c("CD3E","CD4","CD8A","NKG7","GNLY","MKI67") )
VlnPlot(Aur.combined_1_1, pt.size = 0, features = c("CD79A","CD79B","CD8A","NKG7","GNLY","MKI67") )
VlnPlot(Aur.combined_1_1, pt.size = 0, features = c("IL7R","CCR7","CD8A","NKG7","GNLY","MKI67") )
VlnPlot(Aur.combined_1_1, pt.size = 0, features = c("CCDR7","FSCN1","CD207","CD200","PDCD1LG2","CCL17","CCL22","TNFRSF4","CD274","CD40","CD80","CD86") )

# コンタミ探し
VlnPlot(Aur.combined_1_1, pt.size = 0, features = c("NPPA","MYL7","MYL4","TTN","MYBPC3","TNNT2"))
VlnPlot(Aur.combined_1_1, pt.size = 0, features = c("RYR2","PLN","SLC8A1","DCN","GSN","PDGFRA"))
VlnPlot(Aur.combined_1_1, pt.size = 0, features = c("VWF","PECAM1","CDH5","RGS5","ABCC9","KCNJ8"))
VlnPlot(Aur.combined_1_1, pt.size = 0, features = c("MYH11","TAGLN","ACTA2","PLP1","NRXN1","NRXN3"))
VlnPlot(Aur.combined_1_1, pt.size = 0, features = c("GPAM","FASN","LEP","MSLN","WT1","BNC1"))
VlnPlot(Aur.combined_1_1, pt.size = 0, features = c("PTPRC","CD163","CD79B","CD14","IL7R","CD3E","CD4","CD8A","NKG7","GNLY","GZMB"))


VlnPlot(Aur.combined_1_1, features = c("FCGR3A","FCGR2A","ITGB2","CD16","CEACAM8","CD55","CD44","ELANE") )
VlnPlot(Aur.combined_1_1, features = c("PTPRC","CD14","MPO","LYZ","ITGB2","CD16","FCGR3A","FCGR3B","CD55","CD44","ELANE","VCAN", "CD52","FCN1","S100A12","IL1B","CSF1R","S100A8","S100A9","CAMP","DEFB4A","CTSS","ITGAM","CD68","CYBB") )
VlnPlot(Aur.combined_1_1, features = c("DEFB1","DEFA1","DEFB103B") )
VlnPlot(Aur.combined_1_1, features = c("S100A6","FOS","S100A12","S100A9","VCAN","KLF6"))

VlnPlot(Aur.combined_1_1, pt.size = 0, features = c("CD14","FCGR3A") )

VlnPlot(Aur.combined_1_1, features = c("HLA-DPA1","HLA-DPB1","HLA-DQA1","CLEC10A","FCER1A","CD1C") )
VlnPlot(Aur.combined_1_1, features = c("TLR4","TLR2","TLR7","TLR3","TLR12", "TLR9") )
VlnPlot(Aur.combined_1_1, features = c("HDC","KIT","TPSAB1","CCL5") )
VlnPlot(Aur.combined_1_1, features = c("IL1B","TNF","CXCL8","NR1H3", "CEBPA","CEBPB") )
VlnPlot(Aur.combined_1_1, features = c("CD68","CSF1R","ITGAM","S100A8","S100A9","S100A12") )
VlnPlot(Aur.combined_1_1, features = c("FCGR3A","CD14","LYZ","VCAN","FCN1","CTSS") )
VlnPlot(Aur.combined_1_1, features = c("CASP1","CASP4","DEFB4A","CYBB","IL1B","VEGFA") )
VlnPlot(Aur.combined_1_1, features = c("CCL2","CCL20","CCL4","HMGB2","KLF6","VEGFA") )


VlnPlot(Aur.combined_1_1, features = c("IL1B","TNF","CXCL8","CXCL3", "ABCG1","ABCA1") )
VlnPlot(Aur.combined_1_1, features = c("CD9","TREM2","TGFB1","APOC1","APOE","C1QA", "C1QB") )
VlnPlot(Aur.combined_1_1, features = c("C1QC","CD59","APOE","","NUPR1", "C1QB") )

VlnPlot(Aur.combined_1_1,pt.size = 0, features = c("PTPRC","CD3E","CD4","CD8A","CCR2") )

LYVE1+ macrophage
VlnPlot(Aur.combined_1_1,pt.size = 0, features = c("LYVE1","HLA-DOA","HLA-DQA1","HLA-DQA2","HLA-DQB1","TIMD4","TREM2") )

Monocyte-derived macrophage
VlnPlot(Aur.combined_1_1,pt.size = 0, features = c("LYVE1","FOLR2","CEBPB","S100A8","CCL13","CCL18") )
VlnPlot(Aur.combined_1_1, features = c("LYVE1","FOLR2","CEBPB","S100A8","CCL13","CCL18") )
Antigen-presenting macrophage
VlnPlot(Aur.combined_1_1,pt.size = 0, features = c("FOLR2","LYVE1","MERTK","HLA-DRA","HLA-DMA","HLA-DMB","HLA-DPA1","TREM2") )

Dock4+ macrophage
VlnPlot(Aur.combined_1_1,pt.size = 0, features = c("DOCK4","IL4R","STAT3","ITGAM","C1QA","FOLR2") )

fibrotic
VlnPlot(Aur.combined_1_1,pt.size = 0, features = c("TREM2","CD9") )


VlnPlot(Aur.combined_1_1, features = c("CEACAM8","MPO","IL6","CD63","S100A8","S100A9") )
VlnPlot(Aur.combined_1_1, features = c("MKI67","ASP175","CSF1R","IL1B","S100A8","S100A9") )
VlnPlot(Aur.combined_1_1, features = c("CD80","CD86","CD40","CD209","IL10","CD83") )
VlnPlot(Aur.combined_1_1, features = c("MT2A","MALAT1","HNRNPU","CD3E","KLF4","CCL5") )
VlnPlot(Aur.combined_1_1, features = c("CD80","CD86","CD40","CD209","IL10","CD83") )
VlnPlot(Aur.combined_1_1, features = c("CD36","CCR2","CD14","CD61","CX3CR1","CCL5","PTPRC") )
VlnPlot(Aur.combined_1_1, features = c("HNRNPU","MALAT1","CD14","ICAM1","CX3CR1","CCL5","PTPRC") )
VlnPlot(Aur.combined_1_1, features = c("HDC","KIT","TPSAB1","GZMB","CCL5","TPSB2") )
VlnPlot(Aur.combined_1_1, features = c("ENPP3","MITF","CD63","ITGA2","CD33","ITGAM") )
VlnPlot(Aur.combined_1_1, features = c("TNF","TNFSF13","CXCL8","CXCL3","CCL3") )
VlnPlot(Aur.combined_1_1, features = c("CD9","MRC1","FCER2","CCL22","AREG","EREG") )
VlnPlot(Aur.combined_1_1, features = c("MKI67","TUBB","STMN1","TYMS","NFKB1","STAT4") )
VlnPlot(Aur.combined_1_1, features = c("C1QA","C1QB","C1QC","TYMS","NFKB1","STAT4") )
VlnPlot(Aur.combined_1_1, features = c("TET2","DNMT3A","ASXL1","PPM1D","NLRP3","STING1","CGAS" ))
VlnPlot(Aur.combined_1_1, features = c("OLR1","SCARB1","SCARB2","MSR1","MARCO","SRSF2") )
VlnPlot(Aur.combined_1_1, features = c("NR1H3","CTSD","CTSL","SPP1","MARCO","FABP4") )
VlnPlot(Aur.combined_1_1, features = c("TET2","TET1","TET3","SPP1","MARCO","SLC25A44") )
VlnPlot(Aur.combined_1_1, features = c("TNF","TNFSF13"),split.by = "stim",  cols = c("sa" ="deepskyblue", "Aur1Epi"="red"))
VlnPlot(Aur.combined_1_1, features = c("CD163","F13A1","MS4A4A","SPP1","MARCO","SLC25A44") )

manuscript
VlnPlot(Aur.combined_1_1,pt.size = 0, features = c("CD68","CD14","CSF1R","LYVE1","CD9","TREM2"))
VlnPlot(Aur.combined_1_1,pt.size = 0, features = c("IL1B","S100A8","VCAN","FCN1","LYZ","FCGR3A"))
VlnPlot(Aur.combined_1_1,pt.size = 0, features = c("FCGR3B","CSF3R","CLEC10A","FCER1A","HLA-DPA1"))

VlnPlot(Aur.combined_1_1,pt.size = 0, features = c("CCR2","CD163","CLEC10A","FCER1A","HLA-DPA1"))

VlnPlot(Aur.combined_1_1,pt.size = 0, features = c("CD3E","CD4","CD8A","IL7R","CCR7"))

VlnPlot(Aur.combined_1_1,pt.size = 0, features = c("HLA-DPA1","HLA-DPB1","HLA-DQA1","CLEC10A","FCER1A","CD1C") )
VlnPlot(Aur.combined_1_1,pt.size = 0, features = c("ITGAX","IFNG","CLEC9A","CADM1","CLEC10A","FCER1A") )

VlnPlot(Aur.combined_1_1,pt.size = 0, features = c("THBD","CLEC9A","XCR1","CADM1","CD1C","CLEC10A") )
VlnPlot(Aur.combined_1_1,pt.size = 0, features = c("FCER1A","CD2","IL3RA","CLEC4C") )



cDC1
VlnPlot(Aur.combined_1_1, pt.size = 0, features = c("THBD","CLEC9A","XCR1","CADM1","IRF4","IRF8","BATF3") )
cDC2
VlnPlot(Aur.combined_1_1, pt.size = 0, features = c("THBD","CD1C","CLEC10A","FCER1A","IRF4","CD2") )
pDC
VlnPlot(Aur.combined_1_1, pt.size = 0, features = c("IL3RA","CLEC4C") )

basophil
VlnPlot(Aur.combined_1_1,pt.size = 0, features = c("CD34","FCER1A","CD9","CLEC12A","ITGA2") )
VlnPlot(Aur.combined_1_1,pt.size = 0, features = c("KLF5","CSF2RB","MCPT8","") )

VlnPlot(Aur.combined_1_1, pt.size = 0, features = c("ACTC1","TNNC1","FCER2","CD22","MYH11","CD34") )


Aur.combined_1_1_trans <- SCTransform (Aur.combined_1_1)


cd_genes <- c("CD68","CD14","CSF1R","CD163","LYVE1","FOLR2","C1QA","IL1B","TLR4","AREG","PDGFB","TREM2","TIMP1","CLEC5A","EREG","IL10","VEGFA","S100A8","VCAN","FCGR3A","CXCR2","FCGR3B","CSF3R","ISG15","HLA-DPA1","HLA-DRA","CD1C","CLEC9A","IL3RA","KIT","CCR3","IL4","MYL7","CD79A","CCR7","IL7R")
DotPlot(Aur.combined_1_1_trans,features = cd_genes)+RotatedAxis()+coord_flip()




# ヒートマップの作成
Aur.combined_1_1.markers <- FindAllMarkers(Aur.combined_1_1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Aur.combined_1_1.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- Aur.combined_1_1.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Aur.combined_1_1, features = top10$gene) + NoLegend()
write.table(top10, "table.txt", quote=F, col.names=F, append=T)

#Myeloid(myocardium/Bcellコンタミ除去)--------
Aur.combined_1_2 <- subset(Aur.combined_1_1, idents = c(0,1,2,3,4,5,6,7,8,10,11,12,14,15))
Aur.combined_1_2 <- ScaleData(Aur.combined_1_2, verbose = FALSE)
Aur.combined_1_2<- RunPCA(Aur.combined_1_2, npcs = 30, verbose = FALSE)
Aur.combined_1_2 <- RunUMAP(Aur.combined_1_2 , reduction = "pca", dims = 1:30)
Aur.combined_1_2 <- FindNeighbors(Aur.combined_1_2 , reduction = "pca", dims = 1:30)
Aur.combined_1_2 <- FindClusters(Aur.combined_1_2, resolution = 0.6)

DimPlot(Aur.combined_1_2, reduction = "umap", label = FALSE, repel = TRUE)
DimPlot(Aur.combined_1_2, reduction = "umap", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined_1_2, reduction = "umap", split.by = "AurBlood", label = FALSE,  repel = FALSE)
DimPlot(Aur.combined_1_2, reduction = "umap", split.by = "AurBlood", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined_1_2, reduction = "umap", split.by = "AFDuration")
DimPlot(Aur.combined_1_2, reduction = "umap", split.by = "sample1ex")
DimPlot(Aur.combined_1_2, reduction = "umap", split.by = "stim")
DimPlot(Aur.combined_1_2, reduction = "umap", split.by = "Auronly1ex")


FeaturePlot(Aur.combined_1_2, features=c('TREM2'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_2, features=c('LYVE1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_2, features=c('IL1B'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_2, features=c('TIMP1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_2, features=c('AREG'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_2, features=c('CCR4'), min.cutoff=0.001, max.cutoff='q90')



# RDS
saveRDS(Aur.combined_1_2, "~/Desktop/Human_LA_LAA/RDS/Aur_cluster_scMyeloid_心筋B除去後.rds")
Aur.combined_1_2 <- readRDS ("~/Desktop/Human_LA_LAA/RDS/Aur_cluster_scMyeloid_心筋B除去後.rds")


# ヒートマップの作成
Aur.combined_1_2.markers <- FindAllMarkers(Aur.combined_1_2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Aur.combined_1_2.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- Aur.combined_1_2.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Aur.combined_1_2, features = top10$gene) + NoLegend()
write.table(top10, "table.txt", quote=F, col.names=F, append=T)




Manuscript候補
VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("CD68","CD14","CSF1R","CD163","IL1B","TLR4"))
VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("TIMP1","IL10","LYVE1","FOLR2","CD9","C1QA"))
VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("TREM2","S100A8","VCAN","LYZ","FCGR3A","ISG15"))
VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("CLEC9A","HLA-DRA","HLA-DQA1","CSF3R","HLA-DPA1","CD1C"))
VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("CXCR2","FCGR3B","CSF3R","HLA-DPA1","CD1C","IL3RA"))

VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("IL10","TGFB1","CLEC5A","TIMP1","AREG","SPP1"))
VlnPlot(Aur.combined_1_2, pt.size = 0.1, features = c("IL4","IL5","CD163","CLEC10A","AREG","MRC1"))

VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("IL4","IL13","STAT6","IRF4","SPP1"))

VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("NAMPT","GOS2","NEAT1"))

VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("TNF","TNFSF13","CXCL8","CXCL3","CCL3") )

VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("MMP9","MMP1","TIMP1","TIMP2","TIMP3") )
VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("MMP2","MMP12","TIMP1","TIMP2","TIMP3") )

VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("LYVE1","IL1B","TREM2","TIMP1","AREG","CLEC5A") )
VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("TLR4","IL1B","TREM2","TIMP1","IL10","CLEC5A") )

eosinophil
VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("CXCR4","IL3RA","PECAM1","PTGDR2","CD69","CCR3") )



FeaturePlot(Aur.combined_1_2, features=c('ISG15'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_2, features=c('CCR7'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_2, features=c('IL7R'), min.cutoff=0.2, max.cutoff='q90')
FeaturePlot(Aur.combined_1_2, features=c('CD79A'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_2, features=c('CD79B'), min.cutoff=0.5, max.cutoff='q90')


VlnPlot(Aur.combined_1_2,pt.size = 0,features = c("CCL2","CSF1","CSF2","CSF3","IFNG","IL4") )
VlnPlot(Aur.combined_1_2,pt.size = 0,features = c("CD44","IL1B","AREG","EGFR","CSF1R","IL34") )
VlnPlot(Aur.combined_1_2,pt.size = 0,features = c("CCR2","IL1B","AREG","EGFR","CSF1R","IL33") )
VlnPlot(Aur.combined_1_2,pt.size = 0,features = c("PDGFA","PDGRRA","PDGFB","PDGFRB","PDGFC","PDGFRC") )
VlnPlot(Aur.combined_1_2,pt.size = 0,features = c("PDGFRL","PDAP1","TGFB1","TGFB2","TGFBR1","TGFBR2") )
VlnPlot(Aur.combined_1_2,pt.size = 0,features = c("TGFB3","TGFB3R","IL17A","IL17RA","IL17B","IL17RB") )
VlnPlot(Aur.combined_1_2,pt.size = 0,features = c("IL17C","IL17RC","IL17D","IL17RD","IL17E","IL17RE") )
VlnPlot(Aur.combined_1_2,pt.size = 0,features = c("IL17REL","IL17F","IL17RF","MMP12","MMP13","SLCA12") )
VlnPlot(Aur.combined_1_2,pt.size = 0,features = c("HFE","CD63","ITGB1","CTRH1","IL6","IL11") )
VlnPlot(Aur.combined_1_2,pt.size = 0,features = c("BRD4","CX3CR1","MEOX1","CXCL10","CCL19","IL11") )
VlnPlot(Aur.combined_1_2,pt.size = 0,features = c("IL1R1","IL1R2","IL1RL1","IL1RAP","IL1RAP2","IL11") )
VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("TIMP1","TIMP2","TIMP3","MMP2","MMP9","IL10"))
VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("C1QA","C1QB","C1QC","IL10","CD9","TREM2"))


cDC1
VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("THBD","CLEC9A","XCR1","CADM1","IRF4","IRF8","BATF3") )
cDC2
VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("THBD","CD1C","CLEC10A","FCER1A","IRF4","CD2") )
pDC
VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("IL3RA","CLEC4C") )
moDC
VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("ITGAX","ITGAM","CD1A","CD1C","MRC1","CD209","SIRPA") )

VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("CD5"))
FeaturePlot(Aur.combined_1_2, features=c('CCR7'), min.cutoff=0.01, max.cutoff='q90')
FeaturePlot(Aur.combined_1_2, features=c('FSCN1'), min.cutoff=0.01, max.cutoff='q90')

#EGF
VlnPlot(Aur.combined_1_2,pt.size = 0,features = c("EGF","TGFA","AREG","BTC","HBEGF","EREG") )
VlnPlot(Aur.combined_1_2,pt.size = 0,features = c("EGFR","ERBB2","ERBB4","C1QA","C1QB","C1QC") )

#FN1
VlnPlot(Aur.combined_1_2,pt.size = 0,features = c("FN1","CCL2","CCR2","CSF1","CSF1R","IL34") )
VlnPlot(Aur.combined_1_2,pt.size = 0,features = c("ITGA2","ITGB1","ITGA4","ITGA5","ITGAV","ITGB7") )
VlnPlot(Aur.combined_1_2,pt.size = 0,features = c("ITGA2B","ITGB6","ITGBN8","CD44","SDC1","SDC4") )

#SPP1
VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("SPP1","NPPA","MYL7","MYL4","TTN","MYBPC3"))
VlnPlot(Aur.combined_1_2, pt.size = 0,features = c("CD44","ITGB1","ITGB3","ITGB5","ITGB6") )
FeaturePlot(Aur.combined_1_2, features=c('SPP1'), min.cutoff=0.5, max.cutoff='q90')


#PDGF
VlnPlot(Aur.combined_1_2,pt.size = 0,features = c("PDGFA","PDGFB","PDGFC","PDGFD","NPPA","MYL7") )
VlnPlot(Aur.combined_1_2,pt.size = 0,features = c("PDGFRA","PDGFRB","NPPA","MYL7","MYL4","TTN") )

VlnPlot(Aur.combined_1_2,pt.size = 0,features = c("STAT3","JAK2","NPPA","MYL7","MYL4","TTN") )


#COLLAGEN receptorのみ
VlnPlot(Aur.combined_1_2,pt.size = 0,features = c("ITGA1","ITGA3","ITGA9","ITGA10","ITGA11","ITGB1") )
VlnPlot(Aur.combined_1_2,pt.size = 0,features = c("ITGB8","ITGAV","CD44","SDC1","SDC4") )
#THBS
VlnPlot(Aur.combined_1_2,pt.size = 0,features = c("THBS1","THBS2","THBS3","THBS4","COMP") )
VlnPlot(Aur.combined_1_2,pt.size = 0,features = c("ITGA3","ITGAV","ITGB1","ITGB3","SDC1","SDC4") )
VlnPlot(Aur.combined_1_2,pt.size = 0,features = c("CD36","CD47","ITGB1","ITGB3","SDC1","SDC4") )

#VEGF
VlnPlot(Aur.combined_1_2,pt.size = 0,features = c("VEGFA","VEGFB","VEGFC","VEGFD","PGF") )
VlnPlot(Aur.combined_1_2,pt.size = 0,features = c("FLT1","KDR","FLT4","MYL7","MYL4","TTN") )


#BMP
VlnPlot(Aur.combined_1_2,pt.size = 0,features = c("BMP2","BMP4","GDF5","GDF6","GDF7","BMP15") )
VlnPlot(Aur.combined_1_2,pt.size = 0,features = c("BMP5","BMP6","BMP7","BMP8A","BMP8B","BMP10") )
VlnPlot(Aur.combined_1_2,pt.size = 0.1,features = c("BMPR1A","BMPR1B","ACVR2A","ACVR2B","BMPR2","ACVR1") )



#IGF
VlnPlot(Aur.combined_1_2,pt.size = 0,features = c("IGF1","IGF2","IGFL3","MYL7","MYL4","TTN") )
VlnPlot(Aur.combined_1_2,pt.size = 0,features = c("IGF1R","IGF2R","ITGAV","ITGB3","ITGA6","ITGB4") )
VlnPlot(Aur.combined_1_2,pt.size = 0,features = c("IGFLR1","IGF2R","ITGAV","ITGB3","ITGA6","ITGB4") )


VlnPlot(Aur.combined_1_2,pt.size = 0,features = c("CHI3L1","IGF2R","ITGAV","ITGB3","ITGA6","ITGB4") )

VlnPlot(Aur.combined_1_2,pt.size = 0,features = c("TNFAIP3","IGF2R","ITGAV","ITGB3","ITGA6","ITGB4","KIT","JAK2") )


FeaturePlot(Aur.combined_1_2, features=c('TIMP1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_2, features=c('IL1B'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_2, features=c('LYVE1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_2, features=c('TREM2'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_2, features=c('C1QC'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_2, features=c('AREG'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_2, features=c('IL10'), min.cutoff=0.5, max.cutoff='q90', split.by = "AFDuration")

FeaturePlot(Aur.combined_1_2, features=c('PDGFB'), min.cutoff=0.5, max.cutoff='q90')



Aur.combined_1_2_trans <- SCTransform (Aur.combined_1_2)

cd_genes <- c("CD68","CD14","CSF1R","CD163","IL1B","TLR4","TIMP1","IL10","LYVE1","FOLR2","CD9","C1QA","TREM2","S100A8","VCAN","LYZ","FCGR3A","ISG15","CXCR2","FCGR3B","CSF3R","HLA-DPA1","CD1C","IL3RA","CLEC9A","KIT","")
DotPlot(Aur.combined_1_2_trans,features = cd_genes)+RotatedAxis()+coord_flip()


Aur.combined_1_2_trans <- SCTransform (Aur.combined_1_2)

cd_genes <- c("CD68","CD14","CSF1R","CD163","IL1B","TLR4","TIMP1","LYVE1","FOLR2","C1QA","TREM2","S100A8","VCAN","FCGR3A","CXCR2","FCGR3B","CSF3R","HLA-DPA1","HLA-DRA","CD1C","CLEC9A","IL3RA","KIT","CCR3","IL4")
DotPlot(Aur.combined_1_2_trans,features = cd_genes)+RotatedAxis()+coord_flip()


New
Aur.combined_1_2_trans <- SCTransform (Aur.combined_1_2)

cd_genes <- c("CD68","CD14","CSF1R","CD163","LYVE1","FOLR2","C1QA","IL1B","TLR4","AREG","PDGFB","TREM2","TIMP1","CLEC5A","EREG","IL10","VEGFA","S100A8","VCAN","FCGR3A","CXCR2","FCGR3B","CSF3R","ISG15","HLA-DPA1","HLA-DRA","CD1C","CLEC9A","IL3RA","KIT","CCR3","IL4")
DotPlot(Aur.combined_1_2_trans,features = cd_genes)+RotatedAxis()+coord_flip()

cd_genes <- c("CCR7","FSCN1","CD5","CD68","CD14","CSF1R","CD163","LYVE1","FOLR2","C1QA","IL1B","TLR4","AREG","PDGFB","TREM2","TIMP1","CLEC5A","EREG","IL10","VEGFA","S100A8","VCAN","FCGR3A","CXCR2","FCGR3B","CSF3R","ISG15","HLA-DPA1","HLA-DRA","CD1C","CLEC9A","IL3RA","KIT","CCR3","IL4","KLF5","CD69")
DotPlot(Aur.combined_1_2_trans,features = cd_genes)+RotatedAxis()+coord_flip()



#Myeloid(baso除去)--------final
Aur.combined_1_3 <- subset(Aur.combined_1_2, idents = c(0,1,2,3,4,5,6,7,8,9,11,12,13,14,15))
Aur.combined_1_3 <- ScaleData(Aur.combined_1_3, verbose = FALSE)
Aur.combined_1_3<- RunPCA(Aur.combined_1_3, npcs = 30, verbose = FALSE)
Aur.combined_1_3 <- RunUMAP(Aur.combined_1_3 , reduction = "pca", dims = 1:30)
Aur.combined_1_3 <- FindNeighbors(Aur.combined_1_3 , reduction = "pca", dims = 1:30)
Aur.combined_1_3 <- FindClusters(Aur.combined_1_3, resolution = 0.6)

DimPlot(Aur.combined_1_3, reduction = "umap", label = FALSE, repel = TRUE)
DimPlot(Aur.combined_1_3, reduction = "umap", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined_1_3, reduction = "umap", split.by = "AurBlood", label = FALSE,  repel = FALSE)
DimPlot(Aur.combined_1_3, reduction = "umap", split.by = "AurBlood", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined_1_3, reduction = "umap", split.by = "sample1ex")
DimPlot(Aur.combined_1_3, reduction = "umap", split.by = "stim")
DimPlot(Aur.combined_1_3, reduction = "umap", split.by = "Auronly1ex")
DimPlot(Aur.combined_1_3, reduction = "umap", split.by = "EndoEpiBlood")
DimPlot(Aur.combined_1_3, reduction = "umap", split.by = "AFDuration1ex")



FeaturePlot(Aur.combined_1_3, features=c('TREM2'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_3, features=c('LYVE1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_3, features=c('IL1B'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_3, features=c('TIMP1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_3, features=c('AREG'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_3, features=c('CCL4'), min.cutoff=0.5, max.cutoff='q90')

# 並べ替える
levels(Aur.combined_1_3) # [1] "0" "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14"
levels(Aur.combined_1_3) <- c("0", "1",  "2", "3", "9", "5","6","13","7","8","11","12","14","4","10")

new.cluster.ids <- c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14")
names(new.cluster.ids) <- levels(Aur.combined_1_3)
Aur.combined_1_3<- RenameIdents(Aur.combined_1_3, new.cluster.ids)



# RDS
saveRDS(Aur.combined_1_3, "~/Desktop/Human_LA_LAA/RDS/Aur_cluster_scMyeloid_final.rds")
Aur.combined_1_3 <- readRDS ("~/Desktop/Human_LA_LAA/RDS/Aur_cluster_scMyeloid_final.rds")



#グループ分けEndoEpiBlood1ex
# stim列を抽出する
cell.stim <- Aur.combined_1_3@meta.data[["stim"]]
# 新しいグループ（疾患ありなしの2群）とstim列との対応関係を定義する
newgroup <- c("Other","Other","Other","1Endo","2Epi","1Endo","2Epi","3Blood","3Blood","Other","Other","Other","Other","Other")
names(newgroup) <- c("Aur1Endo","Aur1Epi","Aur2","Aur3Endo","Aur3Epi","Aur4Endo","Aur4Epi","Aur3Blood","Aur4Blood","Aur12","Aur13","Aur10","Aur16","Aur19")
# 抽出したstim列を新しいグループに変換する
cell.newgroup <- newgroup[cell.stim]
# 新しいグループを新しいMetaDataとして追加する
Aur.combined_1_3 <- AddMetaData(Aur.combined_1_3, cell.newgroup, col.name = "EndoEpiBlood")

### AFDuration1ex ####
# stim列を抽出する
cell.stim <- Aur.combined_1_3@meta.data[["stim"]]
# 新しいグループ（疾患ありなしの2群）とstim列との対応関係を定義する
newgroup <- c("3Other","3Other","1PerAF","2LSPerAF","2LSPerAF","2LSPerAF","2LSPerAF","3Other","3Other","2LSPerAF","1PerAF","2LSPerAF","1PerAF","2LSPerAF")
names(newgroup) <- c("Aur1Endo","Aur1Epi","Aur2","Aur3Endo","Aur3Epi","Aur4Endo","Aur4Epi","Aur3Blood","Aur4Blood","Aur12","Aur13","Aur10","Aur16","Aur19")
# 抽出したstim列を新しいグループに変換する
cell.newgroup <- newgroup[cell.stim]
# 新しいグループを新しいMetaDataとして追加する
Aur.combined_1_3 <- AddMetaData(Aur.combined_1_3, cell.newgroup, col.name = "AFDuration1ex")

# ヒートマップの作成
Aur.combined_1_3.markers <- FindAllMarkers(Aur.combined_1_3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Aur.combined_1_3.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- Aur.combined_1_3.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Aur.combined_1_3, features = top10$gene) + NoLegend()
write.table(top10, "table.txt", quote=F, col.names=F, append=T)


Aur.combined_1_3_prop. <- prop.table(table(Idents(Aur.combined_1_3),Aur.combined_1_3$AurBlood))
write.table(Aur.combined_1_3_prop., "Aur.combined_1_3_prop.AurBlood.txt", quote=F, col.names=F, append=T)
write.csv(Aur.combined_1_3_prop., "Aur.combined_1_3_prop.AurBlood.csv")

Aur.combined_1_3_prop. <- prop.table(table(Idents(Aur.combined_1_3),Aur.combined_1_3$AFDuration1ex))
write.table(Aur.combined_1_3_prop., "Aur.combined_1_3_prop.AFDuration1ex.txt", quote=F, col.names=F, append=T)
write.csv(Aur.combined_1_3_prop., "Aur.combined_1_3_prop.AFDuration1ex.csv")

Aur.combined_1_3_prop. <- prop.table(table(Idents(Aur.combined_1_3),Aur.combined_1_3$EndoEpiBlood))
write.table(Aur.combined_1_3_prop., "Aur.combined_1_3_prop.EndoEpiBlood.txt", quote=F, col.names=F, append=T)
write.csv(Aur.combined_1_3_prop., "Aur.combined_1_3_prop.EndoEpiBlood.csv")


Aur.combined_1_3_prop. <- prop.table(table(Idents(Aur.combined_1_3),Aur.combined_1_3$sample1ex))
write.table(Aur.combined_1_3_prop., "Aur.combined_1_3_prop.sample1ex.txt", quote=F, col.names=F, append=T)
write.csv(Aur.combined_1_3_prop., "Aur.combined_1_3_prop.sample1ex.csv")

Aur.combined_1_3_prop. <- prop.table(table(Idents(Aur.combined_1_3),Aur.combined_1_3$Auronly1ex))
write.table(Aur.combined_1_3_prop., "Aur.combined_1_3_prop.Auronly1ex.txt", quote=F, col.names=F, append=T)
write.csv(Aur.combined_1_3_prop., "Aur.combined_1_3_prop.Auronly1ex.csv")

FeaturePlot(Aur.combined_1_3, features=c('CD5'), min.cutoff=0.01, max.cutoff='q90')
FeaturePlot(Aur.combined_1_3, features=c('CD14'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_3, features=c('CD68'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_3, features=c('CD163'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_3, features=c('C1QA'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_3, features=c('HLA-DRA'), min.cutoff=3, max.cutoff='q90')
FeaturePlot(Aur.combined_1_3, features=c('CLEC10A'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_3, features=c('TNF'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_3, features=c('CCL3'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_3, features=c('CCL4'), min.cutoff=0.5, max.cutoff='q90')

FeaturePlot(Aur.combined_1_3, features=c('CD9'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_3, features=c('TREM2'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_3, features=c('LYVE1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_3, features=c('IL1B'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_3, features=c('AREG'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_3, features=c('EREG'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_3, features=c('TIMP1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_3, features=c('P2RY6'), min.cutoff=0.5, max.cutoff='q90')

# Olink----
VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("IL17A","KCND2","RASGRF1","B4GALNT3","MTERF3", "FABP6") )
VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("RNF4","ARNT","BIRC3","B4GALNT3","MTERF3", "FABP6") )
VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("IL17A","IL17C","CXCL11","MMP12","MTERF3", "FABP6") )
FeaturePlot(Aur.combined_1_3, features=c('MMP12'), min.cutoff=0.05, max.cutoff='q90')


#4 ISG Mo
VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("ISG15","ISG20","IFIT1","IFIT2","IFIT3","CXCL10","CCL2","CCL8","TNFSF10","IL1RN","STAT1","IRF7","TNFSF13B"))

#7
VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("RGS1","HLA-DPB1","HLA-DRB1","CD163"))

#2
VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("HES1","C1QA","C1QB","C1QC","IGF1","CCL3","CCL4","HBEGF","APOE","LIPA","CCL18"))

VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("TREM2","S100A8","VCAN","LYZ","FCGR3A","ISG15"))
VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("CLEC9A","HLA-DRA","HLA-DQA1","CSF3R","HLA-DPA1","CD1C"))
VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("CXCR2","FCGR3B","CSF3R","HLA-DPA1","CD1C","IL3RA"))
VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("KIT","FCGR3B","CSF3R","HLA-DPA1","CD1C","IL3RA"))

VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("IL10","TGFB1","CLEC5A","TIMP1","AREG","SPP1"))
VlnPlot(Aur.combined_1_3, pt.size = 0.1, features = c("IL4","IL5","CD163","CLEC10A","AREG","MRC1"))

VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("IL4","IL13","STAT6","IRF4","SPP1"))

VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("NAMPT","GOS2","NEAT1"))

VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("TNF","TNFSF13","CXCL8","CXCL3","CCL3") )

VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("MMP9","MMP1","TIMP1","TIMP2","TIMP3") )
VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("MMP2","MMP12","TIMP1","TIMP2","TIMP3") )

VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("LYVE1","IL1B","TREM2","TIMP1","AREG","CLEC5A") )
VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("TLR4","IL1B","TREM2","TIMP1","IL10","CLEC5A") )

eosinophil
VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("CXCR4","IL3RA","PECAM1","PTGDR2","CD69","CCR3") )



FeaturePlot(Aur.combined_1_3, features=c('ISG15'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_3, features=c('CCR7'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_3, features=c('IL7R'), min.cutoff=0.2, max.cutoff='q90')
FeaturePlot(Aur.combined_1_3, features=c('CD79A'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_3, features=c('CD79B'), min.cutoff=0.5, max.cutoff='q90')


VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("CCL2","CSF1","CSF2","CSF3","IFNG","IL4") )
VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("CD44","IL1B","AREG","EGFR","CSF1R","IL34") )
VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("CCR2","IL1B","AREG","EGFR","CSF1R","IL33") )
VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("PDGFA","PDGRRA","PDGFB","PDGFRB","PDGFC","PDGFRC") )
VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("PDGFRL","PDAP1","TGFB1","TGFB2","TGFBR1","TGFBR2") )
VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("TGFB3","TGFB3R","IL17A","IL17RA","IL17B","IL17RB") )
VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("IL17C","IL17RC","IL17D","IL17RD","IL17E","IL17RE") )
VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("IL17REL","IL17F","IL17RF","MMP12","MMP13","SLCA12") )
VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("HFE","CD63","ITGB1","CTRH1","IL6","IL11") )
VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("BRD4","CX3CR1","MEOX1","CXCL10","CCL19","IL11") )
VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("IL1R1","IL1R2","IL1RL1","IL1RAP","IL1RAP2","IL11") )
VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("TIMP1","TIMP2","TIMP3","MMP2","MMP9","IL10"))
VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("C1QA","C1QB","C1QC","IL10","CD9","TREM2"))


cDC1
VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("THBD","CLEC9A","XCR1","CADM1","IRF4","IRF8","BATF3") )
cDC2
VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("THBD","CD1C","CLEC10A","FCER1A","IRF4","CD2") )
pDC
VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("IL3RA","CLEC4C") )
moDC
VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("ITGAX","ITGAM","CD1A","CD1C","MRC1","CD209","SIRPA") )

mregDC
VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("CCR7","FSCN1","CCL17","CCL22","TNFRS4","CD200","CD274","CD40","CD80","CD86","PDCD11G2") )
regulatory
VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("CD200","CD274","PDCD1LG2") )
maturation
VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("CD40","CCR7","IL12B") )


VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("IL4","IL12A") )

VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("CD5"))
FeaturePlot(Aur.combined_1_3, features=c('CCR7'), min.cutoff=0.01, max.cutoff='q90')
FeaturePlot(Aur.combined_1_3, features=c('FSCN1'), min.cutoff=0.01, max.cutoff='q90')

#EGF
VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("EGF","TGFA","AREG","BTC","HBEGF","EREG") )
VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("EGFR","ERBB2","ERBB4","C1QA","C1QB","C1QC") )

#FN1
VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("FN1","CCL2","CCR2","CSF1","CSF1R","IL34") )
VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("ITGA2","ITGB1","ITGA4","ITGA5","ITGAV","ITGB7") )
VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("ITGA2B","ITGB6","ITGBN8","CD44","SDC1","SDC4") )

#SPP1
VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("SPP1","NPPA","MYL7","MYL4","TTN","MYBPC3"))
VlnPlot(Aur.combined_1_3, pt.size = 0,features = c("CD44","ITGB1","ITGB3","ITGB5","ITGB6") )
FeaturePlot(Aur.combined_1_3, features=c('SPP1'), min.cutoff=0.5, max.cutoff='q90')


#PDGF
VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("PDGFA","PDGFB","PDGFC","PDGFD","NPPA","MYL7") )
VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("PDGFRA","PDGFRB","NPPA","MYL7","MYL4","TTN") )

VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("STAT3","JAK2","NPPA","MYL7","MYL4","TTN") )


#COLLAGEN receptorのみ
VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("ITGA1","ITGA3","ITGA9","ITGA10","ITGA11","ITGB1") )
VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("ITGB8","ITGAV","CD44","SDC1","SDC4") )
#THBS
VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("THBS1","THBS2","THBS3","THBS4","COMP") )
VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("ITGA3","ITGAV","ITGB1","ITGB3","SDC1","SDC4") )
VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("CD36","CD47","ITGB1","ITGB3","SDC1","SDC4") )

#VEGF
VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("VEGFA","VEGFB","VEGFC","VEGFD","PGF") )
VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("FLT1","KDR","FLT4","MYL7","MYL4","TTN") )


#BMP
VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("BMP2","BMP4","GDF5","GDF6","GDF7","BMP15") )
VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("BMP5","BMP6","BMP7","BMP8A","BMP8B","BMP10") )
VlnPlot(Aur.combined_1_3,pt.size = 0.1,features = c("BMPR1A","BMPR1B","ACVR2A","ACVR2B","BMPR2","ACVR1") )



#IGF
VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("IGF1","IGF2","IGFL3","MYL7","MYL4","TTN") )
VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("IGF1R","IGF2R","ITGAV","ITGB3","ITGA6","ITGB4") )
VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("IGFLR1","IGF2R","ITGAV","ITGB3","ITGA6","ITGB4") )


VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("CHI3L1","IGF2R","ITGAV","ITGB3","ITGA6","ITGB4") )

VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("TNFAIP3","IGF2R","ITGAV","ITGB3","ITGA6","ITGB4") )
VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("P2RY6","GJA1") )


FeaturePlot(Aur.combined_1_3, features=c('TIMP1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_3, features=c('IL1B'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_3, features=c('LYVE1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_3, features=c('TREM2'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_3, features=c('C1QC'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_3, features=c('AREG'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_3, features=c('IL10'), min.cutoff=0.5, max.cutoff='q90', split.by = "AFDuration")

FeaturePlot(Aur.combined_1_3, features=c('PDGFB'), min.cutoff=0.5, max.cutoff='q90')



Aur.combined_1_3_trans <- SCTransform (Aur.combined_1_3)

cd_genes <- c("CD68","CD14","CSF1R","CD163","C1QA","IL1B","TNFAIP3","TLR4","AREG","EREG","TIMP1","CLEC5A","IL10","VEGFA","TREM2","CD9","LYVE1","FOLR2","S100A8","S100A12","LYZ","FCGR3A","ISG15","HLA-DPA1","HLA-DRA","TNF","CCL3","CCL4","CLEC10A","CLEC4A","CD1C","PDGFB","CLEC9A","IL3RA","CCR7","FSCN1","CCL19","CXCR2","FCGR3B","CSF3R","KIT","HDC")
DotPlot(Aur.combined_1_3_trans,features = cd_genes)+RotatedAxis()+coord_flip()


cd_genes <- c("CD68","CD14","CSF1R","CD163","CCL3","CCL4","TNF","IL1B","AREG","EREG","HBEGF","TREM2","CD9","PDGFB","IGF1","LYVE1","FOLR2","TIMP1","THBS1","CLEC5A","IL10","VEGFA","VEGFB","S100A8","S100A12","LYZ","FCGR3A","ISG15","IFIT3","HLA-DPA1","HLA-DRA","CLEC10A","CLEC4A","CD1C","CLEC9A","IL3RA","CCR7","FSCN1","CXCR2","FCGR3B","CSF3R","KIT","HDC")
DotPlot(Aur.combined_1_3_trans,features = cd_genes)+RotatedAxis()+coord_flip()


# ヒートマップの作成
Aur.combined_1_3.markers <- FindAllMarkers(Aur.combined_1_3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Aur.combined_1_3.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- Aur.combined_1_3.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Aur.combined_1_3, features = top10$gene) 
write.table(top10, "table.txt", quote=F, col.names=F, append=T)



#Macsのみ--------
Aur.combined_1_4 <- subset(Aur.combined_1_3, idents = c(0,1,2,3,4))
# RDS
saveRDS(Aur.combined_1_4, "~/Desktop/Human_LA_LAA/RDS/Aur_cluster_scMacs.rds")
Aur.combined_1_4 <- readRDS ("~/Desktop/Human_LA_LAA/RDS/Aur_cluster_scMacs.rds")

DimPlot(Aur.combined_1_4, reduction = "umap", split.by = "AurBlood")

DimPlot(Aur.combined_1_4, reduction = "umap", split.by = "AurBlood",cols = c("0" = "#F8766D","1" = "#E58700","2" = "#C99800","3" = "#A3A500","4" = "#6BB100"))


VlnPlot(Aur.combined_1_4, pt.size = 0, features = c("IL1B","TNFAIP3","AREG","THBS1","EREG","TIMP1"), cols = c("0" = "#F8766D","1" = "#E58700","2" = "#C99800","3" = "#A3A500","4" = "#6BB100"))
VlnPlot(Aur.combined_1_4, pt.size = 0, features = c("CLEC5A","IL10","VEGFA","TREM2","CD9","LYVE1"), cols = c("0" = "#F8766D","1" = "#E58700","2" = "#C99800","3" = "#A3A500","4" = "#6BB100"))
VlnPlot(Aur.combined_1_4, pt.size = 0, features = c("PDGFC","IGF1","LYVE1","C1QA","C1QB","C1QC"), cols = c("0" = "#F8766D","1" = "#E58700","2" = "#C99800","3" = "#A3A500","4" = "#6BB100"))

VlnPlot(Aur.combined_1_4, pt.size = 0, features = c("CCL3","CCL4","TNF","IL1B","AREG","EREG"), cols = c("0" = "#F8766D","1" = "#E58700","2" = "#C99800","3" = "#A3A500","4" = "#6BB100"))
VlnPlot(Aur.combined_1_4, pt.size = 0, features = c("HBEGF","TREM2","CD9","PDGFB","PDGFC","IGF1"), cols = c("0" = "#F8766D","1" = "#E58700","2" = "#C99800","3" = "#A3A500","4" = "#6BB100"))
VlnPlot(Aur.combined_1_4, pt.size = 0, features = c("LYVE1","FOLR2","TIMP1","THBS1","CLEC5A","IL10"), cols = c("0" = "#F8766D","1" = "#E58700","2" = "#C99800","3" = "#A3A500","4" = "#6BB100"))
VlnPlot(Aur.combined_1_4, pt.size = 0, features = c("VEGFA","VEGFB","CCL3","CCL4","TNF","IL1B"), cols = c("0" = "#F8766D","1" = "#E58700","2" = "#C99800","3" = "#A3A500","4" = "#6BB100"))


VlnPlot(Aur.combined_1_4, pt.size = 0,assay = "RNA", features = c("TNF","IL1B","TIMP1","TREM2","LYVE1","IL10"), cols = c("0" = "#F8766D","1" = "#E58700","2" = "#C99800","3" = "#A3A500","4" = "#6BB100"))
VlnPlot(Aur.combined_1_4, pt.size = 0,assay = "RNA", features = c("VEGFA","EREG","CLEC5A","AREG","CD9","LYVE1"), cols = c("0" = "#F8766D","1" = "#E58700","2" = "#C99800","3" = "#A3A500","4" = "#6BB100"))
VlnPlot(Aur.combined_1_4, pt.size = 0,assay = "RNA", features = c("CCL3","CCL4","FOLR2","PDGFB","PDGFC","IGF1"), cols = c("0" = "#F8766D","1" = "#E58700","2" = "#C99800","3" = "#A3A500","4" = "#6BB100"))
VlnPlot(Aur.combined_1_4, pt.size = 0,assay = "RNA", features = c("THBS1","CD9","TNF","VEGFB","HBEGF","PDGFA"), cols = c("0" = "#F8766D","1" = "#E58700","2" = "#C99800","3" = "#A3A500","4" = "#6BB100"))
VlnPlot(Aur.combined_1_4, pt.size = 0,assay = "RNA", features = c("MRC1","IGF1","LYVE1","FOLR2","HBEGF","PDGFA"), cols = c("0" = "#F8766D","1" = "#E58700","2" = "#C99800","3" = "#A3A500","4" = "#6BB100"))


色変更
library(scales)
show_col(hue_pal()(4))
show_col(hue_pal()(15))
cols = c("0" = "#F8766D","1" = "#E58700","2" = "#C99800","3" = "#A3A500","4" = "#6BB100")

cols = c("0" = "#F8766D","1" = "#E58700","2" = "#C99800","3" = "#A3A500","4" = "#6BB100","5" = "#00BA38","6" = "#00BF7D","7" = "#00BF7D")




# RDS
saveRDS(Aur.combined_1_4, "~/Desktop/Human_LA_LAA/RDS/Aur_cluster_scMacs.rds")
Aur.combined_1_4 <- readRDS ("~/Desktop/Human_LA_LAA/RDS/Aur_cluster_scMacs.rds")


#Macs+Monosのみ--------
Aur.combined_1_5 <- subset(Aur.combined_1_3, idents = c(0,1,2,3,4,5,6,7))
Aur.combined_1_5 <- ScaleData(Aur.combined_1_5, verbose = FALSE)
Aur.combined_1_5<- RunPCA(Aur.combined_1_5, npcs = 30, verbose = FALSE)
Aur.combined_1_5 <- RunUMAP(Aur.combined_1_5 , reduction = "pca", dims = 1:30)
Aur.combined_1_5 <- FindNeighbors(Aur.combined_1_5 , reduction = "pca", dims = 1:30)
Aur.combined_1_5 <- FindClusters(Aur.combined_1_5, resolution = 0.4)

# RDS
saveRDS(Aur.combined_1_5, "~/Desktop/Human_LA_LAA/RDS/Aur_cluster_scMacsMonos.rds")
Aur.combined_1_5 <- readRDS ("~/Desktop/Human_LA_LAA/RDS/Aur_cluster_scMacsMonos.rds")

DimPlot(Aur.combined_1_5, reduction = "umap", label = FALSE, repel = TRUE)
DimPlot(Aur.combined_1_5, reduction = "umap", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined_1_5, reduction = "umap", split.by = "AurBlood", label = FALSE,  repel = FALSE)
DimPlot(Aur.combined_1_5, reduction = "umap", split.by = "AurBlood", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined_1_5, reduction = "umap", split.by = "sample1ex")
DimPlot(Aur.combined_1_5, reduction = "umap", split.by = "stim")
DimPlot(Aur.combined_1_5, reduction = "umap", split.by = "Auronly1ex")
DimPlot(Aur.combined_1_5, reduction = "umap", split.by = "EndoEpiBlood")

DimPlot(Aur.combined_1_5, reduction = "umap", split.by = "Auronly1ex",cols = c("0" = "#F8766D","1" = "#E58700","2" = "#C99800","3" = "#A3A500","4" = "#6BB100","5" = "#00BA38","6" = "#00BF7D","7" = "#00C0AF"))
DimPlot(Aur.combined_1_5, reduction = "umap", split.by = "AurBlood",cols = c("0" = "#F8766D","1" = "#E58700","2" = "#C99800","3" = "#A3A500","4" = "#6BB100","5" = "#00BA38","6" = "#00BF7D","7" = "#00C0AF"))
# 並べ替える
levels(Aur.combined_1_5) # [1] "0" "1" "2" "3" "4" "5" "6" "7" 
levels(Aur.combined_1_5) <- c("0", "2",  "1", "3", "6", "4","5","7")

new.cluster.ids <- c("0","1","2","3","4","5","6","7")
names(new.cluster.ids) <- levels(Aur.combined_1_5)
Aur.combined_1_5<- RenameIdents(Aur.combined_1_5, new.cluster.ids)



Aur.combined_1_5_trans <- SCTransform (Aur.combined_1_5)

cd_genes <- c("CD68","CD14","CSF1R","CD163","C1QA","IL1B","TNFAIP3","TLR4","AREG","EREG","TIMP1","CLEC5A","IL10","VEGFA","TREM2","CD9","LYVE1","FOLR2","S100A8","S100A12","LYZ","FCGR3A","ISG15","HLA-DPA1","HLA-DRA","TNF","CCL3","CCL4","CLEC10A","CLEC4A","CD1C","PDGFB","CLEC9A","IL3RA","CCR7","FSCN1","CCL19","CXCR2","FCGR3B","CSF3R","KIT","HDC")
DotPlot(Aur.combined_1_5_trans,features = cd_genes)+RotatedAxis()+coord_flip()


cd_genes <- c("CD68","CD14","CSF1R","CD163","CCL3","CCL4","TNF","IL1B","AREG","EREG","HBEGF","TREM2","CD9","PDGFB","PDGFC","IGF1","LYVE1","FOLR2","TIMP1","THBS1","CLEC5A","IL10","VEGFA","VEGFB","S100A8","S100A12","LYZ","FCGR3A","ISG15","IFIT3","HLA-DPA1","HLA-DRA","CLEC10A","CLEC4A","CD1C","CLEC9A","IL3RA","CCR7","FSCN1","CXCR2","FCGR3B","CSF3R","KIT","HDC")
DotPlot(Aur.combined_1_5_trans,features = cd_genes)+RotatedAxis()+coord_flip()

Aur.combined_1_5_prop. <- prop.table(table(Idents(Aur.combined_1_5),Aur.combined_1_5$AurBlood))
write.table(Aur.combined_1_5_prop., "Aur.combined_1_5_prop.AurBlood.txt", quote=F, col.names=F, append=T)
write.csv(Aur.combined_1_5_prop., "Aur.combined_1_5_prop.AurBlood.csv")
Aur.combined_1_5_prop. <- prop.table(table(Idents(Aur.combined_1_5),Aur.combined_1_5$Auronly1ex))
write.table(Aur.combined_1_5_prop., "Aur.combined_1_5_prop.Auronly1ex.txt", quote=F, col.names=F, append=T)
write.csv(Aur.combined_1_5_prop., "Aur.combined_1_5_prop.Auronly1ex.csv")


#Inf Macsのみ--------
Aur.combined_1_6 <- subset(Aur.combined_1_3, idents = c(1,4))

VlnPlot(Aur.combined_1_6, pt.size = 0, features = c("AREG"), split.by = "sample1ex")
DimPlot(Aur.combined_1_6, reduction = "umap", split.by = "sample1ex")
DimPlot(Aur.combined_1_6, reduction = "umap", split.by = "Auronly1ex")

#Tcells--------
Aur.combined_2 <- subset(Aur.combined_CD45, idents = c(6,7,8,9,10))
Aur.combined_2 <- ScaleData(Aur.combined_2, verbose = FALSE)
Aur.combined_2<- RunPCA(Aur.combined_2, npcs = 30, verbose = FALSE)
Aur.combined_2 <- RunUMAP(Aur.combined_2 , reduction = "pca", dims = 1:30)
Aur.combined_2 <- FindNeighbors(Aur.combined_2, reduction = "pca", dims = 1:30)
Aur.combined_2 <- FindClusters(Aur.combined_2, resolution = 0.7)

DimPlot(Aur.combined_2, reduction = "umap", label = FALSE, repel = TRUE)
DimPlot(Aur.combined_2, reduction = "umap", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined_2, reduction = "umap", split.by = "AurBlood", label = FALSE,  repel = FALSE)
DimPlot(Aur.combined_2, reduction = "umap", split.by = "AurBlood", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined_2, reduction = "umap", split.by = "AFDuration")
DimPlot(Aur.combined_2, reduction = "umap", split.by = "sample1ex")
DimPlot(Aur.combined_2, reduction = "umap", split.by = "stim")
DimPlot(Aur.combined_2, reduction = "umap", split.by = "Auronly1ex")

# 並べ替える
levels(Aur.combined_2) # [1] "0" "1" "2" "3" "4" "5" "6" "7" "8" "9" "10"
levels(Aur.combined_2) <- c("2","4","1","3","6","8","9","0","5","7","10")

new.cluster.ids <- c("0","1","2","3","4","5","6","7","8","9","10")
names(new.cluster.ids) <- levels(Aur.combined_2)
Aur.combined_2<- RenameIdents(Aur.combined_2, new.cluster.ids)

# RDS
saveRDS(Aur.combined_2, "~/Desktop/Human_LA_LAA/RDS/Aur_cluster_scT.rds")
Aur.combined_2 <- readRDS ("~/Desktop/Human_LA_LAA/RDS/Aur_cluster_scT.rds")


# コンタミ探し
VlnPlot(Aur.combined_2, pt.size = 0, features = c("NPPA","MYL7","MYL4","TTN","MYBPC3","TNNT2"))
VlnPlot(Aur.combined_2, pt.size = 0, features = c("RYR2","PLN","SLC8A1","DCN","GSN","PDGFRA"))
VlnPlot(Aur.combined_2, pt.size = 0, features = c("VWF","PECAM1","CDH5","RGS5","ABCC9","KCNJ8"))
VlnPlot(Aur.combined_2, pt.size = 0, features = c("MYH11","TAGLN","ACTA2","PLP1","NRXN1","NRXN3"))
VlnPlot(Aur.combined_2, pt.size = 0, features = c("GPAM","FASN","LEP","MSLN","WT1","BNC1"))
VlnPlot(Aur.combined_2, pt.size = 0, features = c("PTPRC","CD163","CD79B","CD14","IL7R","CD3E","CD4","CD8A","NKG7","GNLY","GZMB"))



# Myeloid----
VlnPlot(Aur.combined_2, pt.size = 0, features = c("CD68","CD14","CSF1R","IL1B","S100A8", "S100A9") )
VlnPlot(Aur.combined_2, pt.size = 0, features = c("TREM2","CD9","APOE","CXCL3","TNF","CXCL8") )
VlnPlot(Aur.combined_2, pt.size = 0, features = c("HLA-DPA1","HLA-DPB1","HLA-DQA1","CLEC10A","FCER1A","CD1C") )
VlnPlot(Aur.combined_2, pt.size = 0, features = c("PTPRC","HLA-DPB1","HLA-DQA1","CLEC10A","FCER1A","CD1C") )

# Bcells-----
VlnPlot(Aur.combined_2, pt.size = 0, features = c("CD79A","CD79B","FCER2","CD22","MYH11","CD34") )
VlnPlot(Aur.combined_2, pt.size = 0, features = c("LUM","TAGLN","CD34","CDH5","PECAM1","TYROBP") )

# Tcells-----
VlnPlot(Aur.combined_2, pt.size = 0, features = c("CD3E","CD4","CD8A","NKG7","GNLY","GZMB") )
VlnPlot(Aur.combined_2, pt.size = 0, features = c("KLRB1","CCR6","CD1d1","MKI67","IL2RB","GZMB") )


VlnPlot(Aur.combined_2, features = c("CLEC10A","FCER1A","CD1D","FCGR3A") )

VlnPlot(Aur.combined_2, pt.size = 0 , features = c("ACTC1","TNNC1","FCER2","CD22","MYH11","CD34") )
VlnPlot(Aur.combined_2, pt.size = 0 , features = c("LUM","TAGLN","CD34","CDH5","PECAM1","TYROBP") )

VlnPlot(Aur.combined_2, pt.size = 0 , features = c("CD68","CD14","IL1B","S100A8","VCAN","CD3E") )
VlnPlot(Aur.combined_2, pt.size = 0 ,  features = c("CD4","CD8A","NKG7","GNLY","CD79A","CD79B") )
VlnPlot(Aur.combined_2, pt.size = 0 ,  features = c("FCER2","KIT","HDC","MKI67") )

FeaturePlot(Aur.combined_2, features=c('CD4'), min.cutoff=0.5, max.cutoff='q90')

FeaturePlot(Aur.combined_2, features=c('CD8A'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_2, features=c('AREG'), min.cutoff=0.5, max.cutoff='q90')

FeaturePlot(Aur.combined_2, features=c('NPPA'), min.cutoff=0.5, max.cutoff='q90')


# ヒートマップの作成
Aur.combined_2.markers <- FindAllMarkers(Aur.combined_2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Aur.combined_2.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- Aur.combined_2.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Aur.combined_2, features = top10$gene) + NoLegend()
write.table(top10, "table.txt", quote=F, col.names=F, append=T)




VlnPlot(Aur.combined_2,pt.size = 0, features = c("CD3E","CD4","CD8A","IL7R","GZMB","TYROBP") )
VlnPlot(Aur.combined_2,pt.size = 0, features = c("GZMA","GZMK","PRF1","CD28","LEF1","SELL") )
VlnPlot(Aur.combined_2,pt.size = 0, features = c("RORA","GATA3","CD40LG","IL10","IL4","IFNG") )
VlnPlot(Aur.combined_2,pt.size = 0, features = c("IL10","IL4","IL6","IFNG","IL17A") )
VlnPlot(Aur.combined_2,pt.size = 0, features = c("CCR7","GATA3","RORC","FOXP3") )
VlnPlot(Aur.combined_2,pt.size = 0, features = c("GZMB","TBX21","NKG7","GNLY","CX3CR1","CD69") )
VlnPlot(Aur.combined_2,pt.size = 0, features = c("CD68","S100A9","CSF1R","CD14","IL1B","EGR1") )

VlnPlot(Aur.combined_2,pt.size = 0, features = c("CD68","CSF1R","ITGAM","S100A8","S100A9","CXCR2") )
VlnPlot(Aur.combined_2,pt.size = 0, features = c("FCGR3A","CD14","LYZ","VCAN","FCN1","CTSS") )

VlnPlot(Aur.combined_2,pt.size = 0, features = c("RORA","GATA3","PD1","CD40LG","IL10","IL4","LILRB1R") ,split.by = "stim")

VlnPlot(Aur.combined_2,pt.size = 0, features = c("CD3E","CD4","CD8A","GZMB","NKG7","GNLY","TYROBP","FCGR3A","NCAM1"))


VlnPlot(Aur.combined_2,pt.size = 0, features = c("CD3E","CD4","CD8A","GZMB","NKG7","GNLY"))
VlnPlot(Aur.combined_2,pt.size = 0, features = c("GZMB","NKG7","IFNG","FCGR3A","NCAM1","MKI67"))

Aur.combined_2_trans <- SCTransform (Aur.combined_2)

cd_genes <- c("CD3E","CD4","CD8A","GZMB","NKG7","GNLY","TYROBP","IFNG","FCGR3A","NCAM1","MKI67")
DotPlot(Aur.combined_2_trans,features = cd_genes)+RotatedAxis()+coord_flip()




Aur.combined_2_prop. <- prop.table(table(Idents(Aur.combined_2),Aur.combined_2_0$AFControl))
write.table(Aur.combined_2_prop., "Aur.combined_2_prop..txt", quote=F, col.names=F, append=T)
write.csv(Aur.combined_2_prop., "Aur.combined_2_prop..csv")

Aur.combined_2_prop. <- prop.table(table(Idents(Aur.combined_2),Aur.combined_2_0$AFDuration))
write.table(Aur.combined_2_prop., "Aur.combined_2_prop..txt", quote=F, col.names=F, append=T)
write.csv(Aur.combined_2_prop., "Aur.combined_2_prop..csv")


#CD8----  未さきにCD4解析
Aur.combined_2_1 <- subset(Aur.combined_2, idents = c(2,3,4,5,6))
Aur.combined_2_1 <- ScaleData(Aur.combined_2_1, verbose = FALSE)
Aur.combined_2_1<- RunPCA(Aur.combined_2_1, npcs = 30, verbose = FALSE)
Aur.combined_2_1 <- RunUMAP(Aur.combined_2_1, reduction = "pca", dims = 1:30)
Aur.combined_2_1 <- FindNeighbors(Aur.combined_2_1, reduction = "pca", dims = 1:30)
Aur.combined_2_1 <- FindClusters(Aur.combined_2_1, resolution = 0.3)


# RDS
saveRDS(Aur.combined_2_1, "~/Desktop/Human_LA_LAA/RDS/Aur_2_1_cluster_sc5.rds")
Aur.combined_2_1 <- readRDS ("~/Desktop/Human_LA_LAA/RDS/Aur_2_1_cluster_sc5.rds")




# 新しいグループ分けで可視化
DimPlot(Aur.combined_2_1, reduction = "umap", label = FALSE, repel = TRUE)
DimPlot(Aur.combined_2_1, reduction = "umap", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined_2_1, reduction = "umap", split.by = "AurBlood")
DimPlot(Aur.combined_2_1, reduction = "umap", split.by = "AFDuration")
DimPlot(Aur.combined_2_1, reduction = "umap", split.by = "EndoEpiBlood")
DimPlot(Aur.combined_2_1, reduction = "umap", split.by = "stim")


VlnPlot(Aur.combined_2_1, pt.size = 0 , features = c("ACTC1","TNNC1","MYL7","NPPA","MYH11") )

FeaturePlot(Aur.combined_2_1, features=c('MYL7'), min.cutoff=0.5, max.cutoff='q90')

VlnPlot(Aur.combined_2_1, pt.size = 0 , features = c("ACTC1","TNNC1","FCER2","CD22","MYH11","CD34") )
VlnPlot(Aur.combined_2_1, pt.size = 0 , features = c("LUM","TAGLN","CD34","CDH5","PECAM1","TYROBP") )

VlnPlot(Aur.combined_2_1, pt.size = 0 , features = c("CD68","CD14","IL1B","S100A8","VCAN","CD3E") )
VlnPlot(Aur.combined_2_1, pt.size = 0 ,  features = c("CD4","CD8A","NKG7","GNLY","CD79A","CD79B") )
VlnPlot(Aur.combined_2_1, pt.size = 0 ,  features = c("FCER2","KIT","HDC","MKI67") )



VlnPlot(Aur.combined_2_1, pt.size = 0,    features = c("GZMA","GZMK","GZMB","TBX21","CX3CR1","NKG7") )
VlnPlot(Aur.combined_2_1, pt.size = 0,   features = c("CD69","IL7R","IFNG","CD103","SELL","CD27") )
VlnPlot(Aur.combined_2_1, pt.size = 0,   features = c("PRF1","IL2","IFNG","EOMES","SELL","CD27") )
VlnPlot(Aur.combined_2_1, pt.size = 0,  features = c("TNF","CD27","IFNG","MKI67","CD74","CXCR6") )
VlnPlot(Aur.combined_2_1, pt.size = 0, features = c("XCL1","CD7","LGALS3","ITGAE","CLEC10A","GNLY") )
VlnPlot(Aur.combined_2_1, pt.size = 0,features = c("GZMB","GZMH","PRF1","LEF1") )
VlnPlot(Aur.combined_2_1, pt.size = 0, features = c("KLF10","CCR7","IFNG","GZMH","CD40LG","CD4") )

VlnPlot(Aur.combined_2_1, pt.size = 0, features = c("KLF10","CCR7","IFNG","GZMH","CD40LG","CD4") )

VlnPlot(Aur.combined_2_1,pt.size = 0, features = c("CD8A","CD4","GNLY","TYROBP","TGFB1") )

VlnPlot(Aur.combined_2_1,pt.size = 0, features = c("LEF1","SELL","ACTN1","BACH2","CX3CR1","GZMK","SIPRG","IFNG","PDCD1","LAYN","SLC4A10") )
VlnPlot(Aur.combined_2_1,pt.size = 0, features = c("CD27","IL7R","GZMA","GZMK","CD69","CD74","GZMB","TBX21","NKG7","CX3CR1","GZMH","PRF1","TNF","IFNG","CD40LG","CXCR6") )

VlnPlot(Aur.combined_2_1,pt.size = 0, features = c("NPPA","MYL7") )



Naive
VlnPlot(Aur.combined_2_1, pt.size = 0, features = c("LEF1","LTB","CCR7") )
Effector
VlnPlot(Aur.combined_2_1, pt.size = 0, features = c("CX3CR1","FCGR3A","FGFBP2") )
MAIT
VlnPlot(Aur.combined_2_1, pt.size = 0, features = c("SLC4A10","ZBTB16","RORC") )
Tex
VlnPlot(Aur.combined_2_1, pt.size = 0, features = c("CTLA4","PDCD1","HAVCR2","LAYN") )


Manuscript
VlnPlot(Aur.combined_2_1,pt.size = 0, features = c("IL7R","LTB","GZMB","TBX21","NKG7","CX3CR1"))

VlnPlot(Aur.combined_2_1,pt.size = 0, features = c("GZMA","GZMK","CD69","CD74","IFNG","TNF","CXCR6"))
VlnPlot(Aur.combined_2_1,pt.size = 0, features = c("GZMK","CD69","CD74","IFNG","TNF","CXCR6"))


# ヒートマップの作成
Aur.combined_2_1.markers <- FindAllMarkers(Aur.combined_2_1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Aur.combined_2_1.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- Aur.combined_2_1.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Aur.combined_2_1, features = top10$gene) + NoLegend()
write.table(top10, "table.txt", quote=F, col.names=F, append=T)


Aur.combined_2_1_trans <- SCTransform (Aur.combined_2_1)

cd_genes <- c("CD27","IL7R","GZMA","GZMK","CD69","CD74","GZMB","TBX21","NKG7","CX3CR1","GZMH","PRF1","TNF","IFNG","CD40LG","CXCR6")
DotPlot(Aur.combined_2_1_trans,features = cd_genes)+RotatedAxis()+coord_flip()



New
Aur.combined_2_1_trans <- SCTransform (Aur.combined_2_1)

cd_genes <- c("IL7R","LTB","GZMB","TBX21","NKG7","CX3CR1","GZMA","GZMK","CD69","CD74","IFNG","TNF","CXCR6")
DotPlot(Aur.combined_2_1_trans,features = cd_genes)+RotatedAxis()+coord_flip()




#CD8(myocardiumコンタミ除去)----
Aur.combined_2_1_1 <- subset(Aur.combined_2_1, idents = c(0,1,2,3,4,5,6,7,8,10))
Aur.combined_2_1_1 <- ScaleData(Aur.combined_2_1_1, verbose = FALSE)
Aur.combined_2_1_1<- RunPCA(Aur.combined_2_1_1, npcs = 30, verbose = FALSE)
Aur.combined_2_1_1 <- RunUMAP(Aur.combined_2_1_1, reduction = "pca", dims = 1:30)
Aur.combined_2_1_1 <- FindNeighbors(Aur.combined_2_1_1, reduction = "pca", dims = 1:30)
Aur.combined_2_1_1 <- FindClusters(Aur.combined_2_1_1, resolution = 0.05)


# RDS
saveRDS(Aur.combined_2_1_1, "~/Desktop/Human_LA_LAA/RDS/Aur_2_1_1_cluster_sc5.rds")
Aur.combined_2_1_1 <- readRDS ("~/Desktop/Human_LA_LAA/RDS/Aur_2_1_1_cluster_sc5.rds")




# 新しいグループ分けで可視化
DimPlot(Aur.combined_2_1_1, reduction = "umap", label = FALSE, repel = TRUE)
DimPlot(Aur.combined_2_1_1, reduction = "umap", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined_2_1_1, reduction = "umap", split.by = "AurBlood")
DimPlot(Aur.combined_2_1_1, reduction = "umap", split.by = "AFDuration")
DimPlot(Aur.combined_2_1_1, reduction = "umap", split.by = "EndoEpiBlood")
DimPlot(Aur.combined_2_1_1, reduction = "umap", split.by = "stim")


VlnPlot(Aur.combined_2_1_1, pt.size = 0 , features = c("ACTC1","TNNC1","MYL7","NPPA","MYH11") )

FeaturePlot(Aur.combined_2_1_1, features=c('MYL7'), min.cutoff=0.5, max.cutoff='q90')

VlnPlot(Aur.combined_2_1_1, pt.size = 0 , features = c("ACTC1","TNNC1","FCER2","CD22","MYH11","CD34") )
VlnPlot(Aur.combined_2_1_1, pt.size = 0 , features = c("LUM","TAGLN","CD34","CDH5","PECAM1","TYROBP") )

VlnPlot(Aur.combined_2_1_1, pt.size = 0 , features = c("CD68","CD14","IL1B","S100A8","VCAN","CD3E","HLA-DR1") )
VlnPlot(Aur.combined_2_1_1, pt.size = 0 ,  features = c("CD4","CD8A","NKG7","GNLY","CD79A","CD79B") )
VlnPlot(Aur.combined_2_1_1, pt.size = 0 ,  features = c("FCER2","KIT","HDC","MKI67") )
VlnPlot(Aur.combined_2_1_1, pt.size = 0, features = c("CLEC9A","HLA-DRA","HLA-DQA1","CSF3R","HLA-DPA1","CD1C"))
VlnPlot(Aur.combined_2_1_1, pt.size = 0, features = c("CXCR2","FCGR3B","CSF3R","HLA-DPA1","CD1C","IL3RA"))



VlnPlot(Aur.combined_2_1_1, pt.size = 0,    features = c("GZMA","GZMK","GZMB","TBX21","CX3CR1","NKG7") )
VlnPlot(Aur.combined_2_1_1, pt.size = 0,   features = c("CD69","IL7R","IFNG","CD103","SELL","CD27") )
VlnPlot(Aur.combined_2_1_1, pt.size = 0,   features = c("PRF1","IL2","IFNG","EOMES","SELL","CD27") )
VlnPlot(Aur.combined_2_1_1, pt.size = 0,  features = c("TNF","CD27","IFNG","MKI67","CD74","CXCR6") )
VlnPlot(Aur.combined_2_1_1, pt.size = 0, features = c("XCL1","CD7","LGALS3","ITGAE","CLEC10A","GNLY") )
VlnPlot(Aur.combined_2_1_1, pt.size = 0,features = c("GZMB","GZMH","PRF1","LEF1") )
VlnPlot(Aur.combined_2_1_1, pt.size = 0, features = c("KLF10","CCR7","IFNG","GZMH","CD40LG","CD4") )

VlnPlot(Aur.combined_2_1_1, pt.size = 0, features = c("KLF10","CCR7","IFNG","GZMH","CD40LG","CD4") )

VlnPlot(Aur.combined_2_1_1,pt.size = 0, features = c("CD8A","CD4","GNLY","TYROBP","TGFB1") )

VlnPlot(Aur.combined_2_1_1,pt.size = 0, features = c("LEF1","SELL","ACTN1","BACH2","CX3CR1","GZMK","SIPRG","IFNG","PDCD1","LAYN","SLC4A10") )
VlnPlot(Aur.combined_2_1_1,pt.size = 0, features = c("CD27","IL7R","GZMA","GZMK","CD69","CD74","GZMB","TBX21","NKG7","CX3CR1","GZMH","PRF1","TNF","IFNG","CD40LG","CXCR6") )

VlnPlot(Aur.combined_2_1_1,pt.size = 0, features = c("NPPA","MYL7") )



Naive
VlnPlot(Aur.combined_2_1_1, pt.size = 0, features = c("LEF1","LTB","CCR7") )
Effector
VlnPlot(Aur.combined_2_1_1, pt.size = 0, features = c("CX3CR1","FCGR3A","FGFBP2") )
MAIT
VlnPlot(Aur.combined_2_1_1, pt.size = 0, features = c("SLC4A10","ZBTB16","RORC") )
Tex
VlnPlot(Aur.combined_2_1_1, pt.size = 0, features = c("CTLA4","PDCD1","HAVCR2","LAYN") )


Basic Research in Cardiology 2021
resident
VlnPlot(Aur.combined_2_1_1, pt.size = 0, features = c("CD69") )
naive
VlnPlot(Aur.combined_2_1_1, pt.size = 0, features = c("LEF1","LTB","CCR7") )
cytotoxic
VlnPlot(Aur.combined_2_1_1, pt.size = 0, features = c("NKG7","PRF1","GZMB","GNLY") )
Tex
VlnPlot(Aur.combined_2_1_1, pt.size = 0, features = c("PDCD1","HAVCR2","LAG1") )


FeaturePlot (Aur.combined_2_1_1, features = c("LTB"), min.cutoff=0.00005, max.cutoff='q90')
FeaturePlot (Aur.combined_2_1_1, features = c("HAVCR2"), min.cutoff=0.00005, max.cutoff='q90')
FeaturePlot (Aur.combined_2_1_1, features = c("PRF1"), min.cutoff=0.00005, max.cutoff='q90')
FeaturePlot (Aur.combined_2_1_1, features = c("SLC4A10"), min.cutoff=0.00005, max.cutoff='q90')
FeaturePlot (Aur.combined_2_1_1, features = c("CX3CR1"), min.cutoff=0.00005, max.cutoff='q90')
FeaturePlot (Aur.combined_2_1_1, features = c("CXCR6"), min.cutoff=0.00005, max.cutoff='q90')



Manuscript
VlnPlot(Aur.combined_2_1_1,pt.size = 0, features = c("IL7R","LTB","GZMB","TBX21","NKG7","CX3CR1"))

VlnPlot(Aur.combined_2_1_1,pt.size = 0, features = c("GZMA","GZMK","CD69","CD74","IFNG","TNF","CXCR6"))
VlnPlot(Aur.combined_2_1_1,pt.size = 0, features = c("GZMK","CD69","CD74","IFNG","TNF","CXCR6"))


# ヒートマップの作成
Aur.combined_2_1_1.markers <- FindAllMarkers(Aur.combined_2_1_1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Aur.combined_2_1_1.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- Aur.combined_2_1_1.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Aur.combined_2_1_1, features = top10$gene) + NoLegend()
write.table(top10, "table.txt", quote=F, col.names=F, append=T)





New
Aur.combined_2_1_1_trans <- SCTransform (Aur.combined_2_1_1)

cd_genes <- c("IL7R","LTB","GZMB","TBX21","NKG7","CX3CR1","GZMA","GZMK","CD69","CD74","IFNG","TNF","CXCR6")
DotPlot(Aur.combined_2_1_1_trans,features = cd_genes)+RotatedAxis()+coord_flip()





Aur.combined_2_1_1_prop. <- prop.table(table(Idents(Aur.combined_2_1_1),Aur.combined_2_1_1$AFControl))
write.table(Aur.combined_2_1_1_prop., "Aur.combined_2_1_1_prop..txt", quote=F, col.names=F, append=T)
write.csv(Aur.combined_2_1_1_prop., "Aur.combined_2_1_1_prop..csv")

Aur.combined_2_1_prop. <- prop.table(table(Idents(Aur.combined_2_1),Aur.combined_2_1$AFDuration))
write.table(Aur.combined_2_1_prop., "Aur.combined_2_1_1_prop..txt", quote=F, col.names=F, append=T)
write.csv(Aur.combined_2_1_prop., "Aur.combined_2_1_1_prop..csv")















#CD4----
Aur.combined_2_2 <- subset(Aur.combined_2, idents = c(1,6,8,13,16))
Aur.combined_2_2 <- ScaleData(Aur.combined_2_2, verbose = FALSE)
Aur.combined_2_2<- RunPCA(Aur.combined_2_2, npcs = 30, verbose = FALSE)
Aur.combined_2_2 <- RunUMAP(Aur.combined_2_2, reduction = "pca", dims = 1:30)
Aur.combined_2_2 <- FindNeighbors(Aur.combined_2_2, reduction = "pca", dims = 1:30)
Aur.combined_2_2 <- FindClusters(Aur.combined_2_2, resolution = 0.5)



DimPlot(Aur.combined_2_2, reduction = "umap", label = FALSE, repel = TRUE)
DimPlot(Aur.combined_2_2, reduction = "umap", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined_2_2, reduction = "umap", split.by = "AurBlood", label = FALSE,  repel = FALSE)
DimPlot(Aur.combined_2_2, reduction = "umap", split.by = "AurBlood", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined_2_2, reduction = "umap", split.by = "AFDuration")
DimPlot(Aur.combined_2_2, reduction = "umap", split.by = "sample1ex")
DimPlot(Aur.combined_2_2, reduction = "umap", split.by = "stim")
DimPlot(Aur.combined_2_2, reduction = "umap", split.by = "Auronly1ex")

Aur.combined_2_2_trans <- SCTransform (Aur.combined_2_2)

cd_genes <- c("CD3E","CD4","CD8A","GZMB","NKG7","GNLY","FOXP3","CTLA4","SELL","CCR7","LEF1","CD28","PRF1","GZMA","GZMK","CD69","IL7R","CCL4","CCL5","TNFAIP3","TBX21","GATA3","RORC")
DotPlot(Aur.combined_2_2_trans,features = cd_genes)+RotatedAxis()+coord_flip()


# RDS
saveRDS(Aur.combined_2_2, "~/Desktop/Human_LA_LAA/RDS/Aur_cluster_scAll_Aur2_2.rds")
Aur.combined_2_2 <- readRDS ("~/Desktop/Human_LA_LAA/RDS/Aur_cluster_scAll_Aur2_2.rds")




# コンタミ探し
VlnPlot(Aur.combined_2_2, pt.size = 0, features = c("NPPA","MYL7","MYL4","TTN","MYBPC3","TNNT2"))
VlnPlot(Aur.combined_2_2, pt.size = 0, features = c("RYR2","PLN","SLC8A1","DCN","GSN","PDGFRA"))
VlnPlot(Aur.combined_2_2, pt.size = 0, features = c("VWF","PECAM1","CDH5","RGS5","ABCC9","KCNJ8"))
VlnPlot(Aur.combined_2_2, pt.size = 0, features = c("MYH11","TAGLN","ACTA2","PLP1","NRXN1","NRXN3"))
VlnPlot(Aur.combined_2_2, pt.size = 0, features = c("GPAM","FASN","LEP","MSLN","WT1","BNC1"))
VlnPlot(Aur.combined_2_2, pt.size = 0, features = c("PTPRC","CD163","CD79B","CD14","IL7R","CD3E","CD4","CD8A","NKG7","GNLY","GZMB"))



VlnPlot(Aur.combined_2_2, pt.size = 0 , features = c("ACTC1","TNNC1","MYL7","NPPA","MYH11") )

FeaturePlot(Aur.combined_2_2, features=c('NPPA'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_2_2, features=c('CD8A'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_2_2, features=c('GNLY'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_2_2, features=c('CD4'), min.cutoff=0.5, max.cutoff='q90')


VlnPlot(Aur.combined_2_2, pt.size = 0 , features = c("ACTC1","TNNC1","FCER2","CD22","MYH11","CD34") )
VlnPlot(Aur.combined_2_2, pt.size = 0 , features = c("LUM","TAGLN","CD34","CDH5","PECAM1","TYROBP") )

VlnPlot(Aur.combined_2_2, pt.size = 0 , features = c("CD68","CD14","IL1B","S100A8","VCAN","CD3E") )
VlnPlot(Aur.combined_2_2, pt.size = 0 ,  features = c("CD4","CD8A","NKG7","GNLY","CD79A","CD79B") )
VlnPlot(Aur.combined_2_2, pt.size = 0 ,  features = c("FCER2","KIT","HDC","MKI67") )

VlnPlot(Aur.combined_2_2,pt.size = 0, features = c("HLA-DPA1","HLA-DPB1","HLA-DQA1","CLEC10A","FCER1A","CD1C") )
VlnPlot(Aur.combined_2_2,pt.size = 0, features = c("THBD","CLEC9A","XCR1","CADM1","CD1C","CLEC10A") )
VlnPlot(Aur.combined_2_2,pt.size = 0, features = c("FCER1A","CD2","IL3RA","CLEC4C") )



cDC1
VlnPlot(Aur.combined_2_2, pt.size = 0, features = c("THBD","CLEC9A","XCR1","CADM1","IRF4","IRF8","BATF3") )
cDC2
VlnPlot(Aur.combined_2_2, pt.size = 0, features = c("THBD","CD1C","CLEC10A","FCER1A","IRF4","CD2") )
pDC
VlnPlot(Aur.combined_2_2, pt.size = 0, features = c("IL3RA","CLEC4C") )


# ヒートマップの作成
Aur.combined_2_2.markers <- FindAllMarkers(Aur.combined_2_2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Aur.combined_2_2.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- Aur.combined_2_2.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Aur.combined_2_2, features = top10$gene) + NoLegend()
write.table(top10, "table.txt", quote=F, col.names=F, append=T)




#CD4(最終、CD8除去)--------
Aur.combined_2_3 <- subset(Aur.combined_2_2, idents = c(0,1,2,3,4,6,7,8,9,10,11,12))
Aur.combined_2_3 <- ScaleData(Aur.combined_2_3, verbose = FALSE)
Aur.combined_2_3<- RunPCA(Aur.combined_2_3, npcs = 30, verbose = FALSE)
Aur.combined_2_3 <- RunUMAP(Aur.combined_2_3 , reduction = "pca", dims = 1:30)
Aur.combined_2_3 <- FindNeighbors(Aur.combined_2_3, reduction = "pca", dims = 1:30)
Aur.combined_2_3 <- FindClusters(Aur.combined_2_3, resolution = 0.5)


FeaturePlot(Aur.combined_2_3, features=c('FOXP3'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_2_3, features=c('PRF1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_2_3, features=c('GATA3'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_2_3, features=c('GNLY'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_2_3, features=c('AREG'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_2_3, features=c('CD8A'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_2_3, features=c('CD4'), min.cutoff=0.5, max.cutoff='q90')

# 新しいグループ分けで可視化
DimPlot(Aur.combined_2_3, reduction = "umap", label = FALSE, repel = TRUE)
DimPlot(Aur.combined_2_3, reduction = "umap", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined_2_3, reduction = "umap", split.by = "AurBlood", label = FALSE,  repel = FALSE)
DimPlot(Aur.combined_2_3, reduction = "umap", split.by = "AurBlood", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined_2_3, reduction = "umap", split.by = "AFDuration")
DimPlot(Aur.combined_2_3, reduction = "umap", split.by = "sample1ex")
DimPlot(Aur.combined_2_3, reduction = "umap", split.by = "stim")
DimPlot(Aur.combined_2_3, reduction = "umap", split.by = "Auronly1ex")

# RDS
saveRDS(Aur.combined_2_3, "~/Desktop/Human_LA_LAA/RDS/Aur_cluster_scAll_Aur2_3.rds")
Aur.combined_2_3 <- readRDS ("~/Desktop/Human_LA_LAA/RDS/Aur_cluster_scAll_Aur2_3.rds")



VlnPlot(Aur.combined_2_3,pt.size = 0, features = c("GZMB","GZMK","PRF1","LEF1","IL7R","SELL") )
VlnPlot(Aur.combined_2_3,pt.size = 0, features = c("TBX21","GATA3","CCR7","RORC","FOXP3","CTLA4") )
VlnPlot(Aur.combined_2_3,pt.size = 0, features = c("IL17A","IL23","IL2","IL10","IL4","IFNG") )
VlnPlot(Aur.combined_2_3,pt.size = 0, features = c("GZMB","CD28","PRF1","CD69","IL7R","SELL") )
VlnPlot(Aur.combined_2_3,pt.size = 0, features = c("CXCR3","CCR4","CCR6","TER1","TNF")) 

VlnPlot(Aur.combined_2_3,pt.size = 0, features = c("CCL4","CCL5","TNFAIP3","RORC","CD28","GZMB") )
VlnPlot(Aur.combined_2_3,pt.size = 0, features = c("LTB","NKG7","IFNG","CD69","CD28","GZMB") )
VlnPlot(Aur.combined_2_3,pt.size = 0, features = c("CD44","NKG7","IFNG","CD69","CD28","KLF10") )

VlnPlot(Aur.combined_2_3,pt.size = 0, features = c("CD3E","CD4","CD8A","CD68","CD14","IL1B","S100A9") )
VlnPlot(Aur.combined_2_3,pt.size = 0, features = c("CD3E","CD4","CD8A","GNLY","GZMB","TYROBP") )
VlnPlot(Aur.combined_2_3, pt.size = 0, features = c("AREG","CXCL8","TIMP1","IL4","IL10","IL33") )


VlnPlot(Aur.combined_2_3,pt.size = 0,features = c("CCL2","CSF1","CSF2","CSF3","IFNG","IL4") )
VlnPlot(Aur.combined_2_3,pt.size = 0,features = c("CD44","IL1B","AREG","EGFR","CSF1R","IL34") )
VlnPlot(Aur.combined_2_3,pt.size = 0,features = c("CCR2","IL1B","AREG","EGFR","CSF1R","IL34") )
VlnPlot(Aur.combined_2_3,pt.size = 0,features = c("PDGFA","PDGRRA","PDGFB","PDGFRB","PDGFC","PDGFRC") )
VlnPlot(Aur.combined_2_3,pt.size = 0,features = c("PDGFRL","PDAP1","TGFB1","TGFB2","TGFBR1","TGFBR2") )
VlnPlot(Aur.combined_2_3,pt.size = 0,features = c("TGFB3","TGFB3R","IL17A","IL17RA","IL17B","IL17RB") )
VlnPlot(Aur.combined_2_3,pt.size = 0,features = c("IL17C","IL17RC","IL17D","IL17RD","IL17E","IL17RE") )
VlnPlot(Aur.combined_2_3,pt.size = 0,features = c("IL17REL","IL17F","IL17RF","MMP12","MMP13","SLCA12") )
VlnPlot(Aur.combined_2_3,pt.size = 0,features = c("HFE","CD63","ITGB1","CTRH1","IL6","IL11") )
VlnPlot(Aur.combined_2_3,pt.size = 0,features = c("BRD4","CX3CR1","MEOX1","CXCL10","CCL19","IL11") )
VlnPlot(Aur.combined_2_3,pt.size = 0,features = c("IL1R1","IL1R2","IL1RL1","IL1RAP","IL1RAP2","IL11") )
VlnPlot(Aur.combined_2_3, pt.size = 0, features = c("TIMP1","TIMP2","TIMP3","MMP2","MMP9","IL10"))
VlnPlot(Aur.combined_2_3, pt.size = 0, features = c("C1QA","C1QB","C1QC","IL10","CD9","TREM2"))



#EGF
VlnPlot(Aur.combined_2_3,pt.size = 0.1,features = c("EGF","TGFA","AREG","BTC","HBEGF","EREG") )
VlnPlot(Aur.combined_2_3,pt.size = 1,features = c("EGFR","ERBB2","ERBB4","C1QA","C1QB","C1QC") )



FeaturePlot(Aur.combined_2_3,pt.size = 0, features=c('CD69'), min.cutoff=2.0, max.cutoff='q90',split.by = "stim" )
FeaturePlot(Aur.combined_2_3,pt.size = 0, features=c('CD69'), min.cutoff=2.0, max.cutoff='q90')

pt.size = 0, 

VlnPlot(Aur.combined_2_3,pt.size = 0, features = c("CD28","CD40LG") ,  split.by = "AFDuration" , cols = c("Control" ="yellow", "LSPerAF"="blue","PerAF" ="red" ),pt.size = 0.1)

VlnPlot(Aur.combined_2_3,pt.size = 0, features = c("CD28","CD40LG") ,  split.by = "AFDuration" , cols = c("Control" ="yellow", "LSPerAF"="blue","PerAF" ="red" ),pt.size = 0)



Manuscript
VlnPlot(Aur.combined_2_3,pt.size = 0, features = c("CD3E","GZMA","GZMK","CCL4","CCL5","TNFAIP3"))
VlnPlot(Aur.combined_2_3,pt.size = 0, features = c("IL7R","SELL","CCR7","CD28","PRF1","GATA3"))
VlnPlot(Aur.combined_2_3,pt.size = 0, features = c("CD69","CCR7","LEF1","FOXP3","CTLA4"))



Aur.combined_2_3_trans <- SCTransform (Aur.combined_2_3)

cd_genes <- c("CD3E","CD4","CD8A","GZMB","NKG7","GNLY","FOXP3","CTLA4","SELL","CCR7","LEF1","CD28","PRF1","GZMA","GZMK","CD69","IL7R","CCL4","CCL5","TNFAIP3","TBX21","GATA3","RORC","MKI67","AREG")
DotPlot(Aur.combined_2_3_trans,features = cd_genes)+RotatedAxis()+coord_flip()




New
Aur.combined_2_3_trans <- SCTransform (Aur.combined_2_3)

cd_genes <- c("CD3E","CD28","PRF1","GZMA","GZMK","CCL4","CCL5","TNFAIP3","IL7R","SELL","CCR7","LEF1","AREG","FOXP3","CTLA4")
DotPlot(Aur.combined_2_3_trans,features = cd_genes)+RotatedAxis()+coord_flip()



Aur.combined_2_3_trans <- SCTransform (Aur.combined_2_3)

cd_genes <- c("CD3E","GZMA","GZMK","CCL4","CCL5","TNFAIP3","IFNG","CD28","PRF1","IL7R","SELL","CCR7","LEF1","AREG","FOXP3","CTLA4")
DotPlot(Aur.combined_2_3_trans,features = cd_genes)+RotatedAxis()+coord_flip()



# ヒートマップの作成
Aur.combined_2_3.markers <- FindAllMarkers(Aur.combined_2_3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Aur.combined_2_3.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- Aur.combined_2_3.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Aur.combined_2_3, features = top10$gene) + NoLegend()
write.table(top10, "table.txt", quote=F, col.names=F, append=T)





Aur.combined_2_3_prop. <- prop.table(table(Idents(Aur.combined_2_3),Aur.combined_2_3$AFControl))
write.table(Aur.combined_2_3_prop., "Aur.combined_2_3_prop..txt", quote=F, col.names=F, append=T)
write.csv(Aur.combined_2_3_prop., "Aur.combined_2_3_prop..csv")

Aur.combined_2_3_prop. <- prop.table(table(Idents(Aur.combined_2_3),Aur.combined_2_3$AFDuration))
write.table(Aur.combined_2_3_prop., "Aur.combined_2_3_prop..txt", quote=F, col.names=F, append=T)
write.csv(Aur.combined_2_3_prop., "Aur.combined_2_3_prop..csv")



#NK cell----
Aur.combined_2_4 <- subset(Aur.combined_2_0, idents = c(2,3))
Aur.combined_2_4 <- ScaleData(Aur.combined_2_4, verbose = FALSE)
Aur.combined_2_4<- RunPCA(Aur.combined_2_4, npcs = 30, verbose = FALSE)
Aur.combined_2_4 <- RunUMAP(Aur.combined_2_4, reduction = "pca", dims = 1:30)
Aur.combined_2_4 <- FindNeighbors(Aur.combined_2_4, reduction = "pca", dims = 1:30)
Aur.combined_2_4 <- FindClusters(Aur.combined_2_4, resolution = 0.2)
p1 <- DimPlot(Aur.combined_2_4, reduction = "umap", split.by = "stim")
p1
p2 <- DimPlot(Aur.combined_2_4, reduction = "umap", label = TRUE, repel = TRUE)
p2


### AFDuration ####
# stim列を抽出する
cell.stim <- Aur.combined_2_4@meta.data[["stim"]]
# 新しいグループ（疾患ありなしの2群）とstim列との対応関係を定義する
newgroup <- c("PerAF","PerAF","PerAF","LSPerAF","LSPerAF","LSPerAF","LSPerAF","Control","Control","Control","Control","Control","Control")
names(newgroup) <- c("Aur1Endo","Aur1Epi","Aur2","Aur3Endo","Aur3Epi","Aur4Endo","Aur4Epi","Control1","Control2","Control3","Control4","Control5","Control6")
# 抽出したstim列を新しいグループに変換する
cell.newgroup <- newgroup[cell.stim]
# 新しいグループを新しいMetaDataとして追加する
Aur.combined_2_4 <- AddMetaData(Aur.combined_2_4, cell.newgroup, col.name = "AFDuration")



# 新しいグループ分けで可視化
DimPlot(Aur.combined_2_4, reduction = "umap")
DimPlot(Aur.combined_2_4, reduction = "umap", label = TRUE,  repel = TRUE)

DimPlot(Aur.combined_2_4, reduction = "umap", group.by = "AFControl")
DimPlot(Aur.combined_2_4, reduction = "umap",split.by = "AFControl")

DimPlot(Aur.combined_2_4, reduction = "umap",split.by = "AFDuration")

DimPlot(Aur.combined_2_4, reduction = "umap", split.by = "stim")

# RDS
saveRDS(Aur.combined_2_4, "~/Desktop/Human_LA_LAA/RDS/Aur_2_4_cluster_withcontorl.rds")
Aur.combined_2_4 <- readRDS ("~/Desktop/Human_LA_LAA/RDS/Aur_2_4_cluster_withcontorl.rds")

FCGR3A=CD16 NCAM1=CD56 CSF2=GM-CSF CCL3=MIP-1α CCL4=MIP-1β
VlnPlot(Aur.combined_2_4,pt.size = 0, features = c(c"IL10"))

NKが産生 Il-5と13はなし
VlnPlot(Aur.combined_2_4,pt.size = 0, features = c("IFNG","CSF2","IL10","TNF","CCL3","CCL4"))

# ヒートマップの作成
Aur.combined_2_4.markers <- FindAllMarkers(Aur.combined_2_4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Aur.combined_2_4.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- Aur.combined_2_4.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Aur.combined_2_4, features = top10$gene) + NoLegend()
write.table(top10, "table.txt", quote=F, col.names=F, append=T)





Aur.combined_2_4_prop. <- prop.table(table(Idents(Aur.combined_2_4),Aur.combined_2_4$AFControl))
write.table(Aur.combined_2_4_prop., "Aur.combined_2_4_prop..txt", quote=F, col.names=F, append=T)
write.csv(Aur.combined_2_4_prop., "Aur.combined_2_4_prop..csv")

Aur.combined_2_4_prop. <- prop.table(table(Idents(Aur.combined_2_4),Aur.combined_2_4$AFDuration))
write.table(Aur.combined_2_4_prop., "Aur.combined_2_4_prop..txt", quote=F, col.names=F, append=T)
write.csv(Aur.combined_2_4_prop., "Aur.combined_2_4_prop..csv")













#Bcells---------
Aur.combined_3 <- subset(Aur.combined_CD45_1, idents = c(9))
Aur.combined_3 <- ScaleData(Aur.combined_3, verbose = FALSE)
Aur.combined_3 <- RunPCA(Aur.combined_3, npcs = 30, verbose = FALSE)
Aur.combined_3 <- RunUMAP(Aur.combined_3 , reduction = "pca", dims = 1:30)
Aur.combined_3 <- FindNeighbors(Aur.combined_3 , reduction = "pca", dims = 1:30)
Aur.combined_3 <- FindClusters(Aur.combined_3, resolution = 0.5)
p1 <- DimPlot(Aur.combined_3, reduction = "umap", split.by = "stim")
p1
p2 <- DimPlot(Aur.combined_3, reduction = "umap", label = FALSE, repel = TRUE)
p2 



### AFDuration ####
# stim列を抽出する
cell.stim <- Aur.combined_3@meta.data[["stim"]]
# 新しいグループ（疾患ありなしの2群）とstim列との対応関係を定義する
newgroup <- c("PerAF","PerAF","PerAF","LSPerAF","LSPerAF","LSPerAF","LSPerAF","Control","Control","Control","Control","Control","Control")
names(newgroup) <- c("Aur1Endo","Aur1Epi","Aur2","Aur3Endo","Aur3Epi","Aur4Endo","Aur4Epi","Control1","Control2","Control3","Control4","Control5","Control6")
# 抽出したstim列を新しいグループに変換する
cell.newgroup <- newgroup[cell.stim]
# 新しいグループを新しいMetaDataとして追加する
Aur.combined_3 <- AddMetaData(Aur.combined_3, cell.newgroup, col.name = "AFDuration")




# 新しいグループ分けで可視化
DimPlot(Aur.combined_3, reduction = "umap")
DimPlot(Aur.combined_3, reduction = "umap", label = TRUE,  repel = TRUE)

DimPlot(Aur.combined_3, reduction = "umap", group.by = "AFControl")
DimPlot(Aur.combined_3, reduction = "umap",split.by = "AFControl")

DimPlot(Aur.combined_3, reduction = "umap",split.by = "AFDuration")

DimPlot(Aur.combined_3, reduction = "umap", split.by = "stim")




VlnPlot(Aur.combined_3, features = c("CD79A","CD79B","CD27","FCER2","CD22") )
VlnPlot(Aur.combined_3, features = c("CD38","PTPRC","CD40","CD24","CD79B","FCER2") )
VlnPlot(Aur.combined_3, features = c("IGHG1","IGHG3","IGHG","IGHM","IGHD","IGHE"))
VlnPlot(Aur.combined_3,pt.size = 0, features = c("IGHG1","IGHG2","IGHG3","IGHG4","IGHD","IGHM","IGHE","IGHA1"))
VlnPlot(Aur.combined_3, features = c("IGHG1","IGHG3","IGHG","IGHM","IGHD","IGHE"),split.by = "stim",  cols = c("sa" ="deepskyblue", "Aur1Epi"="red"))
VlnPlot(Aur.combined_3,pt.size = 0, features = c("CD40","IL4R","IL21R","PAX5","BACH2","ICOSL","CD80","CD86","HLA-DQB1","HLA-DQB1","HLA-DPA1"))
VlnPlot(Aur.combined_3, features = c("CD69","IL7R","LTB","SELL","CCL4","CCL3","GZMB","GZLY","NKG7","CD22","CD33","CD59","LGHG4","CD19","MS4A1","HLA-DRA","HLA-DRB1","IGHA1","IGHM","IGHD","IGHG1","IGHG2","IGHG3","FCER2","CXCR4","CD24","CD86","CD27","CD38","TNFRSF17"))


VlnPlot(Aur.combined_3,pt.size = 0, features = c("IL-4","IL-5","IL-6","TGFB","IFNG"))

VlnPlot(Aur.combined_3,pt.size = 0, features = c("CD19","ITGAM","MS4A1","ITGAX","IGHM","IGHD","CD27","CD38","CD14"))
VlnPlot(Aur.combined_3,pt.size = 0, features = c("CCR7","LEF1","LTB"))


Breg
VlnPlot(Aur.combined_3,pt.size = 0, features = c("IL10","GZMB","PDCD1"))

memory marker
VlnPlot(Aur.combined_3,pt.size = 0, features = c("MS4A1"))
plasma marker
VlnPlot(Aur.combined_3,pt.size = 0, features = c("PRDM1"))

CD20=MS4A1 
△CD11b=ITGAM CD11c=ITGAX

Tcell
VlnPlot(Aur.combined_3, features = c("CD3E","CD4","CD8A","GNLY","GZMB","TYROBP") )


# ヒートマップの作成------------------------------------------------------------------------------
Aur.combined_3.markers <- FindAllMarkers(Aur.combined_3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Aur.combined_3.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- Aur.combined_3.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Aur.combined_3, features = top10$gene) + NoLegend()
write.table(top10, "table.txt", quote=F, col.names=F, append=T)


Aur.combined_3_trans <- SCTransform (Aur.combined_3)

cd_genes <- c("CD27","CD69","IL7R","LTB","SELL","CCL4","CCL3","GZMA","GNLY","NKG7","CD22","CD33","CD59","CD19","MS4A1","ITGAM","HLA-DRA","HLA-DRB1","IGHA1","IGHM","IGHD","IGHG1","IGHG2","IGHG3","IGHG4","FCER2","CXCR4","CD24","CD86","CD38","TNFRSF17")
DotPlot(Aur.combined_3_trans,features = cd_genes)+RotatedAxis()+coord_flip()




# RDS Bcells(T cell doublet除去前)
saveRDS(Aur.combined_3, "~/Desktop/Human_LA_LAA/RDS/Aur_3_cluster_withcontorl.rds")
Aur.combined_3 <- readRDS ("~/Desktop/Human_LA_LAA/RDS/Aur_3_cluster_withcontorl.rds")



#Bcells(T cell doublet除去後)---------
Aur.combined_3_1 <- subset(Aur.combined_3, idents = c(0,1,2,4,5))
Aur.combined_3_1 <- ScaleData(Aur.combined_3_1, verbose = FALSE)
Aur.combined_3_1 <- RunPCA(Aur.combined_3_1, npcs = 30, verbose = FALSE)
Aur.combined_3_1 <- RunUMAP(Aur.combined_3_1 , reduction = "pca", dims = 1:30)
Aur.combined_3_1 <- FindNeighbors(Aur.combined_3_1 , reduction = "pca", dims = 1:30)
Aur.combined_3_1 <- FindClusters(Aur.combined_3_1, resolution = 0.4)
p1 <- DimPlot(Aur.combined_3_1, reduction = "umap", split.by = "stim")
p1
p2 <- DimPlot(Aur.combined_3_1, reduction = "umap", label = FALSE, repel = TRUE)
p2 

# 並べ替え
levels(Aur.combined_3_1) # [1] "0" "1" "2" "3" 
levels(Aur.combined_3_1) <- c("0", "2", "1", "3")

new.cluster.ids <- c("0","1","2","3")
names(new.cluster.ids) <- levels(Aur.combined_3_1)
Aur.combined_3_1<- RenameIdents(Aur.combined_3_1, new.cluster.ids)

# RDS Bcells(T cell doublet除去後)
saveRDS(Aur.combined_3_1, "~/Desktop/Human_LA_LAA/RDS/Aur_3_1_cluster_withcontorl.rds")
Aur.combined_3_1 <- readRDS ("~/Desktop/Human_LA_LAA/RDS/Aur_3_1_cluster_withcontorl.rds")



### AFDuration ####
# stim列を抽出する
cell.stim <- Aur.combined_3_1@meta.data[["stim"]]
# 新しいグループ（疾患ありなしの2群）とstim列との対応関係を定義する
newgroup <- c("PerAF","PerAF","PerAF","LSPerAF","LSPerAF","LSPerAF","LSPerAF","Control","Control","Control","Control","Control","Control")
names(newgroup) <- c("Aur1Endo","Aur1Epi","Aur2","Aur3Endo","Aur3Epi","Aur4Endo","Aur4Epi","Control1","Control2","Control3","Control4","Control5","Control6")
# 抽出したstim列を新しいグループに変換する
cell.newgroup <- newgroup[cell.stim]
# 新しいグループを新しいMetaDataとして追加する
Aur.combined_3_1 <- AddMetaData(Aur.combined_3_1, cell.newgroup, col.name = "AFDuration")




# 新しいグループ分けで可視化
DimPlot(Aur.combined_3_1, reduction = "umap")
DimPlot(Aur.combined_3_1, reduction = "umap", label = TRUE,  repel = TRUE)

DimPlot(Aur.combined_3_1, reduction = "umap", group.by = "AFControl")
DimPlot(Aur.combined_3_1, reduction = "umap",split.by = "AFControl")

DimPlot(Aur.combined_3_1, reduction = "umap",split.by = "AFDuration")

DimPlot(Aur.combined_3_1, reduction = "umap", split.by = "stim")




VlnPlot(Aur.combined_3_1, features = c("CD79A","CD79B","CD27","FCER2","CD22") )
VlnPlot(Aur.combined_3_1, features = c("CD38","PTPRC","CD40","CD24","CD79B","FCER2") )
VlnPlot(Aur.combined_3_1, features = c("IGHG1","IGHG3","IGHG","IGHM","IGHD","IGHE"))
VlnPlot(Aur.combined_3_1, features = c("IGHG1","IGHG2","IGHG3","IGHG4","IGHD","IGHM","IGHE","IGHA1"))
VlnPlot(Aur.combined_3_1, features = c("IGHG1","IGHG3","IGHG","IGHM","IGHD","IGHE"),split.by = "stim",  cols = c("sa" ="deepskyblue", "Aur1Epi"="red"))
VlnPlot(Aur.combined_3_1,pt.size = 0, features = c("CD40","IL4R","IL21R","PAX5","BACH2","ICOSL","CD80","CD86","HLA-DQB1","HLA-DQB1","HLA-DPA1"))
VlnPlot(Aur.combined_3_1, features = c("CD69","IL7R","LTB","SELL","CCL4","CCL3","GZMB","GNLY","NKG7","CD22","CD33","CD59","IGHG4","CD19","MS4A1","HLA-DRA","HLA-DRB1","IGHA1","IGHM","IGHD","IGHG1","IGHG2","IGHG3","FCER2","CXCR4","CD24","CD86","CD27","CD38","TNFRSF17"))


VlnPlot(Aur.combined_3_1,pt.size = 0, features = c("IL-4","IL-5","IL-6","TGFB","IFNG"))

VlnPlot(Aur.combined_3_1,pt.size = 0, features = c("CD19","ITGAM","MS4A1","ITGAX","IGHM","IGHD","CD27","CD38","CD14"))
VlnPlot(Aur.combined_3_1,pt.size = 0, features = c("CCR7","LEF1","LTB"))


VlnPlot(Aur.combined_3_1,pt.size = 0, features = c("CXCR3","TNF","KLF10","ZEB2","TEFC","ZBTB32","YBX3"))




Breg
VlnPlot(Aur.combined_3_1,pt.size = 0, features = c("IL10","GZMB","PDCD1"))

memory marker
VlnPlot(Aur.combined_3_1,pt.size = 0, features = c("MS4A1"))
plasma marker
VlnPlot(Aur.combined_3_1,pt.size = 0, features = c("PRDM1"))

naive
VlnPlot(Aur.combined_3_1,pt.size = 0, features = c("TCL1A","IL4R","CCR7","LEF1","LTB"))
VlnPlot(Aur.combined_3_1, features = c("IL4R","CCR7","LEF1","LTB"))

CD20=MS4A1 
△CD11b=ITGAM CD11c=ITGAX

Tcell
VlnPlot(Aur.combined_3_1, features = c("CD79B","CD79A","CD20","CD19") )


Manuscript
VlnPlot(Aur.combined_3_1,pt.size = 0, features = c("CD79A","CD79B","FCER2","CD22","MS4A1","AIM2"))
VlnPlot(Aur.combined_3_1,pt.size = 0, features = c("CD27","PRDM1","IGHM","IGHD","IGHG1","IGHG2","IGHG3","IGHG4","IGHA1"))
VlnPlot(Aur.combined_3_1,pt.size = 0, features = c("IGHD","IGHG1","IGHG2","IGHG3","IGHG4","IGHA1"))



# ヒートマップの作成------------------------------------------------------------------------------
Aur.combined_3_1.markers <- FindAllMarkers(Aur.combined_3_1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Aur.combined_3_1.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- Aur.combined_3_1.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Aur.combined_3_1, lines.width = NULL, features = top10$gene) 
write.table(top10, "table.txt", quote=F, col.names=F, append=T)


Aur.combined_3_1_trans <- SCTransform (Aur.combined_3_1)

cd_genes <- c("CD27","CD69","IL7R","LTB","SELL","CCL4","CCL3","GZMA","GNLY","NKG7","CD22","CD33","CD59","CD19","MS4A1","ITGAM","HLA-DRA","HLA-DRB1","IGHA1","IGHM","IGHD","IGHG1","IGHG2","IGHG3","IGHG4","FCER2","CXCR4","CD24","CD86","CD38","TNFRSF17")
DotPlot(Aur.combined_3_1_trans,features = cd_genes)+RotatedAxis()+coord_flip()


Aur.combined_3_1_trans <- SCTransform (Aur.combined_3_1)

cd_genes <- c("CD79A","CD79B","FCER2","CD22","MS4A1","CD27","HLA-DRA","HLA-DRB1","IGHM","IGHD","IGHG1","IGHG2","IGHG3","IGHG4","IGHA1")
DotPlot(Aur.combined_3_1_trans,features = cd_genes)+RotatedAxis()+coord_flip()



OK
Aur.combined_3_1_trans <- SCTransform (Aur.combined_3_1)

cd_genes <- c("CD79A","CD79B","FCER2","CD22","MS4A1","AIM2","CD27","PRDM1","IGHM","IGHD","IGHG1","IGHG2","IGHG3","IGHG4","IGHA1")
DotPlot(Aur.combined_3_1_trans,features = cd_genes)+RotatedAxis()+coord_flip()



Aur.combined_3_1_prop. <- prop.table(table(Idents(Aur.combined_3_1),Aur.combined_3_1$AFControl))
write.table(Aur.combined_3_1_prop., "Aur.combined_3_1_prop..txt", quote=F, col.names=F, append=T)
write.csv(Aur.combined_3_1_prop., "Aur.combined_3_1_prop..csv")

Aur.combined_3_1_prop. <- prop.table(table(Idents(Aur.combined_3_1),Aur.combined_3_1$AFDuration))
write.table(Aur.combined_3_1_prop., "Aur.combined_3_1_prop..txt", quote=F, col.names=F, append=T)
write.csv(Aur.combined_3_1_prop., "Aur.combined_3_1_prop..csv")
