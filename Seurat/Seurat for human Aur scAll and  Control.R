# ライブラリの読み込み----------------------------------------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

# split the dataset into a list of two seurat objects (stim and CTRL)----
Aur.list<-list()
data.Aur1Endo <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur1/outs/per_sample_outs/Aur1_Endo/count/sample_feature_bc_matrix/")
data.Aur1Epi <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur1/outs/per_sample_outs/Aur1_Epi/count/sample_feature_bc_matrix/")
data.Aur2 <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur2_without_CellPlex/outs/filtered_feature_bc_matrix/")
data.Aur3Endo <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur3/outs/per_sample_outs/Aur3_Endo/count/sample_feature_bc_matrix/")
data.Aur3Epi <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur3/outs/per_sample_outs/Aur3_Epi/count/sample_feature_bc_matrix/")
data.Aur4Endo <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur4/outs/per_sample_outs/Aur4_Endo/count/sample_feature_bc_matrix/")
data.Aur4Epi <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur4/outs/per_sample_outs/Aur4_Epi/count/sample_feature_bc_matrix/")
data.Aur12 <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur12/outs/filtered_feature_bc_matrix")
data.Aur13 <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur13/outs/filtered_feature_bc_matrix")

data.Aur10 <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur10/outs/filtered_feature_bc_matrix")
data.Aur16 <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur16/outs/filtered_feature_bc_matrix")
data.Aur19 <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur19/outs/filtered_feature_bc_matrix")


data.Control1 <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/Nature2020control/Control1_ERR6449788_cellranger/ERR6449788/outs/filtered_feature_bc_matrix/")
data.Control2 <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/Nature2020control/Control2_ERR6449807_cellranger/outs/filtered_feature_bc_matrix/")
data.Control3 <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/Nature2020control/Control3_ERR6449822_cellranger/outs/filtered_feature_bc_matrix/")
data.Control4 <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/Nature2020control/Control4_ERR7423275_cellranger/outs/filtered_feature_bc_matrix/")
data.Control5 <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/Nature2020control/Control5_ERR7423294_cellranger/outs/filtered_feature_bc_matrix/")
data.Control6 <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/Nature2020control/Control6_ERR7423309_cellranger/outs/filtered_feature_bc_matrix/")


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



Aur.list$Control1 <- CreateSeuratObject(counts = data.Control1, project = "A",min.cells = 3, min.features = 200)
Aur.list$Control1 <- AddMetaData(Aur.list$Control1 , "Control1", col.name = "stim")

Aur.list$Control2 <- CreateSeuratObject(counts = data.Control2, project = "A",min.cells = 3, min.features = 200)
Aur.list$Control2 <- AddMetaData(Aur.list$Control2 , "Control2", col.name = "stim")

Aur.list$Control3 <- CreateSeuratObject(counts = data.Control3, project = "A",min.cells = 3, min.features = 200)
Aur.list$Control3 <- AddMetaData(Aur.list$Control3 , "Control3", col.name = "stim")

Aur.list$Control4 <- CreateSeuratObject(counts = data.Control4, project = "A",min.cells = 3, min.features = 200)
Aur.list$Control4 <- AddMetaData(Aur.list$Control4 , "Control4", col.name = "stim")

Aur.list$Control5 <- CreateSeuratObject(counts = data.Control5, project = "A",min.cells = 3, min.features = 200)
Aur.list$Control5 <- AddMetaData(Aur.list$Control5 , "Control5", col.name = "stim")


Aur.list$Control6 <- CreateSeuratObject(counts = data.Control6, project = "A",min.cells = 3, min.features = 200)
Aur.list$Control6 <- AddMetaData(Aur.list$Control6 , "Control6", col.name = "stim")

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
Aur.list$Aur12[["percent.mt"]] <- PercentageFeatureSet(Aur.list$Aur12, pattern = "^MT-")
Aur.list$Aur13[["percent.mt"]] <- PercentageFeatureSet(Aur.list$Aur13, pattern = "^MT-")

Aur.list$Aur10[["percent.mt"]] <- PercentageFeatureSet(Aur.list$Aur10, pattern = "^MT-")
Aur.list$Aur16[["percent.mt"]] <- PercentageFeatureSet(Aur.list$Aur16, pattern = "^MT-")
Aur.list$Aur19[["percent.mt"]] <- PercentageFeatureSet(Aur.list$Aur19, pattern = "^MT-")


Aur.list$Control1[["percent.mt"]] <- PercentageFeatureSet(Aur.list$Control1, pattern = "^MT-")
Aur.list$Control2[["percent.mt"]] <- PercentageFeatureSet(Aur.list$Control2, pattern = "^MT-")
Aur.list$Control3[["percent.mt"]] <- PercentageFeatureSet(Aur.list$Control3, pattern = "^MT-")
Aur.list$Control4[["percent.mt"]] <- PercentageFeatureSet(Aur.list$Control4, pattern = "^MT-")
Aur.list$Control5[["percent.mt"]] <- PercentageFeatureSet(Aur.list$Control5, pattern = "^MT-")
Aur.list$Control6[["percent.mt"]] <- PercentageFeatureSet(Aur.list$Control6, pattern = "^MT-")


VlnPlot(Aur.list$Aur1Endo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Aur.list$Aur1Epi, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Aur.list$Aur2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Aur.list$Aur3Endo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Aur.list$Aur3Epi, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Aur.list$Aur4Endo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Aur.list$Aur4Epi, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Aur.list$Aur12, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Aur.list$Aur13, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

VlnPlot(Aur.list$Aur10, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Aur.list$Aur16, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Aur.list$Aur19, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

VlnPlot(Aur.list$Control1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Aur.list$Control2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Aur.list$Control3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Aur.list$Control4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Aur.list$Control5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Aur.list$Control6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


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


plot1 <- FeatureScatter(Aur.list$Control1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Aur.list$Control1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2



plot1 <- FeatureScatter(Aur.list$Control2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Aur.list$Control2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(Aur.list$Control3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Aur.list$Control3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(Aur.list$Control4, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Aur.list$Control4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(Aur.list$Control5, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Aur.list$Control5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(Aur.list$Control6, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Aur.list$Control6, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
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
Aur.list$Aur12 <- subset(Aur.list$Aur12, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 8)
Aur.list$Aur13 <- subset(Aur.list$Aur13, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 8)

Aur.list$Aur10 <- subset(Aur.list$Aur10, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 8)
Aur.list$Aur16 <- subset(Aur.list$Aur16, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 8)
Aur.list$Aur19 <- subset(Aur.list$Aur19, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 8)


Aur.list$Control1 <- subset(Aur.list$Control1, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 8)
Aur.list$Control2 <- subset(Aur.list$Control2, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 8)
Aur.list$Control3 <- subset(Aur.list$Control3, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 8)
Aur.list$Control4 <- subset(Aur.list$Control4, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 8)
Aur.list$Control5 <- subset(Aur.list$Control5, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 8)
Aur.list$Control6 <- subset(Aur.list$Control6, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 8)





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



plot1 <- FeatureScatter(Aur.list$Control1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Aur.list$Control1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(Aur.list$Control2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Aur.list$Control2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(Aur.list$Control3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Aur.list$Control3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(Aur.list$Control4, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Aur.list$Control4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(Aur.list$Control5, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Aur.list$Control5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(Aur.list$Control6, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Aur.list$Control6, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# データの正規化（normalization）-----------------------------------------------
Aur.list$Aur1Endo<- NormalizeData(Aur.list$Aur1Endo, normalization.method = "LogNormalize", scale.factor = 10000)
Aur.list$Aur1Epi <- NormalizeData(Aur.list$Aur1Epi, normalization.method = "LogNormalize", scale.factor = 10000)
Aur.list$Aur2<- NormalizeData(Aur.list$Aur2, normalization.method = "LogNormalize", scale.factor = 10000)
Aur.list$Aur3Endo<- NormalizeData(Aur.list$Aur3Endo, normalization.method = "LogNormalize", scale.factor = 10000)
Aur.list$Aur3Epi <- NormalizeData(Aur.list$Aur3Epi, normalization.method = "LogNormalize", scale.factor = 10000)
Aur.list$Aur4Endo<- NormalizeData(Aur.list$Aur4Endo, normalization.method = "LogNormalize", scale.factor = 10000)
Aur.list$Aur4Epi <- NormalizeData(Aur.list$Aur4Epi, normalization.method = "LogNormalize", scale.factor = 10000)
Aur.list$Aur12<- NormalizeData(Aur.list$Aur12, normalization.method = "LogNormalize", scale.factor = 10000)
Aur.list$Aur13<- NormalizeData(Aur.list$Aur13, normalization.method = "LogNormalize", scale.factor = 10000)

Aur.list$Aur10<- NormalizeData(Aur.list$Aur10, normalization.method = "LogNormalize", scale.factor = 10000)
Aur.list$Aur16<- NormalizeData(Aur.list$Aur16, normalization.method = "LogNormalize", scale.factor = 10000)
Aur.list$Aur19<- NormalizeData(Aur.list$Aur19, normalization.method = "LogNormalize", scale.factor = 10000)

Aur.list$Control1 <- NormalizeData(Aur.list$Control1, normalization.method = "LogNormalize", scale.factor = 10000)
Aur.list$Control2 <- NormalizeData(Aur.list$Control2, normalization.method = "LogNormalize", scale.factor = 10000)
Aur.list$Control3 <- NormalizeData(Aur.list$Control3, normalization.method = "LogNormalize", scale.factor = 10000)
Aur.list$Control4 <- NormalizeData(Aur.list$Control4, normalization.method = "LogNormalize", scale.factor = 10000)
Aur.list$Control5 <- NormalizeData(Aur.list$Control5, normalization.method = "LogNormalize", scale.factor = 10000)
Aur.list$Control6 <- NormalizeData(Aur.list$Control6, normalization.method = "LogNormalize", scale.factor = 10000)

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


all.genes <- rownames(Aur.list$Control1)
Aur.list$Control1 <- ScaleData(Aur.list$Control1, features = all.genes)
all.genes <- rownames(Aur.list$Control2)
Aur.list$Control2 <- ScaleData(Aur.list$Control2, features = all.genes)
all.genes <- rownames(Aur.list$Control3)
Aur.list$Control3 <- ScaleData(Aur.list$Control3, features = all.genes)
all.genes <- rownames(Aur.list$Control4)
Aur.list$Control4 <- ScaleData(Aur.list$Control4, features = all.genes)
all.genes <- rownames(Aur.list$Control5)
Aur.list$Control5 <- ScaleData(Aur.list$Control5, features = all.genes)
all.genes <- rownames(Aur.list$Control6)
Aur.list$Control6 <- ScaleData(Aur.list$Control6, features = all.genes)


Aur.list$Aur1Endo
Aur.list$Aur1Epi 
Aur.list$Aur2 
Aur.list$Aur3Endo
Aur.list$Aur3Epi 
Aur.list$Aur4Endo
Aur.list$Aur4Epi 
Aur.list$Aur12
Aur.list$Aur13

Aur.list$Aur10
Aur.list$Aur16
Aur.list$Aur19


Aur.list$Control1
Aur.list$Control2
Aur.list$Control3
Aur.list$Control4
Aur.list$Control5
Aur.list$Control6



# UMAPによるクラスタリング　2つ合わせて------------------------------------------------------------------------------

# normalize and identify variable features for each dataset independently




all.genes <- rownames(Aur.list["Aur1Endo"]) +  rownames(Aur.list["Aur1Epi"]) +   rownames(Aur.list["Aur2"]) +rownames(Aur.list["Aur3Endo"])+   rownames(Aur.list["Aur3Epi"]) +   rownames(Aur.list["Aur4Endo"])+   rownames(Aur.list["Aur4Epi"])+   rownames(Aur.list["Aur12"])+   rownames(Aur.list["Aur13"])+   rownames(Aur.list["Aur10"])+   rownames(Aur.list["Aur16"])+   rownames(Aur.list["Aur19"])+   rownames(Aur.list["Control1"])+   rownames(Aur.list["Control2"])+   rownames(Aur.list["Control3"])+   rownames(Aur.list["Control4"])+   rownames(Aur.list["Control5"])+   rownames(Aur.list["Control6"])
                                                                                                                                                                                              
                                                                                                                                                                                                                                                                                                                                                                  
Aur.list <- lapply(X = Aur.list, FUN = function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
  x <- ScaleData(x, features = all.genes)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = Aur.list)
                                                        

# Integrationの実行
Aur.anchors <- FindIntegrationAnchors(object.list = Aur.list, anchor.features = features)

# this command creates an 'integrated' data assay
Aur.combined <- IntegrateData(anchorset =Aur.anchors)

# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(Aur.combined) <- "integrated"

Aur.combined <- ScaleData(Aur.combined, verbose = FALSE)
Aur.combined <- RunPCA(Aur.combined, npcs = 30, verbose = FALSE)
Aur.combined <- RunUMAP(Aur.combined, reduction = "pca", dims = 1:30)
Aur.combined <- FindNeighbors(Aur.combined, reduction = "pca", dims = 1:30)
Aur.combined <- FindClusters(Aur.combined, resolution = 0.15)

DimPlot(Aur.combined, reduction = "umap", group.by = "stim")
DimPlot(Aur.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 <- DimPlot(Aur.combined, reduction = "umap", split.by = "stim")
p1
p2 <- DimPlot(Aur.combined, reduction = "umap", label = TRUE,  repel = TRUE)

# RDS_first UMAP
saveRDS(Aur.combined, "~/Desktop/Human_LA_LAA/RDS/Aur_cluster_withcontorl_AurscAll_最初.rds")
Aur.combined <- readRDS ("~/Desktop/Human_LA_LAA/RDS/Aur_cluster_withcontorl_AurscAll_最初.rds")



### AFControl1ex ##
# stim列を抽出する
cell.stim <- Aur.combined@meta.data[["stim"]]
# 新しいグループ（疾患ありなしの2群）とstim列との対応関係を定義する
newgroup <- c("3.Other","3.Other","2.AF","2.AF","2.AF","2.AF","2.AF","2.AF","2.AF","2.AF","2.AF","2.AF","1.Control","1.Control","1.Control","1.Control","1.Control","1.Control")
names(newgroup) <- c("Aur1Endo","Aur1Epi","Aur2","Aur3Endo","Aur3Epi","Aur4Endo","Aur4Epi","Aur12","Aur13","Aur10","Aur16","Aur19","Control1","Control2","Control3","Control4","Control5","Control6")
# 抽出したstim列を新しいグループに変換する
cell.newgroup <- newgroup[cell.stim]
# 新しいグループを新しいMetaDataとして追加する
Aur.combined <- AddMetaData(Aur.combined, cell.newgroup, col.name = "AFControl1ex")

### AFDuration1ex ####
# stim列を抽出する
cell.stim <- Aur.combined@meta.data[["stim"]]
# 新しいグループ（疾患ありなしの2群）とstim列との対応関係を定義する
newgroup <- c("4Other","4Other","2PerAF","3LSPerAF","3LSPerAF","3LSPerAF","3LSPerAF","3LSPerAF","2PerAF","3LSPerAF","2PerAF","3LSPerAF","1Control","1Control","1Control","1Control","1Control","1Control")
names(newgroup) <- c("Aur1Endo","Aur1Epi","Aur2","Aur3Endo","Aur3Epi","Aur4Endo","Aur4Epi","Aur12","Aur13","Aur10","Aur16","Aur19","Control1","Control2","Control3","Control4","Control5","Control6")
# 抽出したstim列を新しいグループに変換する
cell.newgroup <- newgroup[cell.stim]
# 新しいグループを新しいMetaDataとして追加する
Aur.combined <- AddMetaData(Aur.combined, cell.newgroup, col.name = "AFDuration1ex")

#グループ分けsamp1ex
# stim列を抽出する
cell.stim <- Aur.combined@meta.data[["stim"]]
# 新しいグループ（疾患ありなしの2群）とstim列との対応関係を定義する
newgroup <- c("Other","Other","1.2","4.3","4.3","5.4","5.4","Other","Other","6.12","2.13","7.10","3.16","8.19","9.Control1","9.Control2","9.Control3","9.Control1","9.Control2","9.Control3")
names(newgroup) <- c("Aur1Endo","Aur1Epi","Aur2","Aur3Endo","Aur3Epi","Aur4Endo","Aur4Epi","Aur3Blood","Aur4Blood","Aur12","Aur13","Aur10","Aur16","Aur19","Control1","Control2","Control3","Control4","Control5","Control6")
# 抽出したstim列を新しいグループに変換する
cell.newgroup <- newgroup[cell.stim]
# 新しいグループを新しいMetaDataとして追加する
Aur.combined <- AddMetaData(Aur.combined, cell.newgroup, col.name = "sample1ex")

#グループ分けfinal for reveise 
# stim列を抽出する
cell.stim <- Aur.combined_CD45_1@meta.data[["stim"]]
# 新しいグループ（疾患ありなしの2群）とstim列との対応関係を定義する
newgroup <- c("Other","Other","1.2","4.3","4.3","5.4","5.4","Other","Other","6.12","2.13","7.10","3.16","8.19","9.Control1.1","9.Control2.1","Other","9.Control1.2","9.Control2.2","Other")
names(newgroup) <- c("Aur1Endo","Aur1Epi","Aur2","Aur3Endo","Aur3Epi","Aur4Endo","Aur4Epi","Aur3Blood","Aur4Blood","Aur12","Aur13","Aur10","Aur16","Aur19","Control1","Control2","Control3","Control4","Control5","Control6")
# 抽出したstim列を新しいグループに変換する
cell.newgroup <- newgroup[cell.stim]
# 新しいグループを新しいMetaDataとして追加する
Aur.combined_CD45_1 <- AddMetaData(Aur.combined_CD45_1, cell.newgroup, col.name = "sample1ex1cotrlex")
# Visualization 別々------

DimPlot(Aur.combined, reduction = "umap", label = FALSE, repel = TRUE)
DimPlot(Aur.combined, reduction = "umap", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined, reduction = "umap", split.by = "AFControl1ex", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined, reduction = "umap", split.by = "AFControl1ex", label = FALSE,  repel = FALSE)
DimPlot(Aur.combined, reduction = "umap", split.by = "AFDuration1ex")
DimPlot(Aur.combined, reduction = "umap", split.by = "sample1ex")
DimPlot(Aur.combined, reduction = "umap", split.by = "stim")

# Myeloid----
VlnPlot(Aur.combined, pt.size = 0, features = c("CD68","CD14","CSF1R","IL1B","S100A8", "S100A9") )
VlnPlot(Aur.combined, pt.size = 0, features = c("TREM2","CD9","APOE","CXCL3","C1QA","TNF","CXCL8") )
VlnPlot(Aur.combined, pt.size = 0, features = c("HLA-DPA1","HLA-DPB1","HLA-DQA1","CLEC10A","FCER1A","CD1C") )
VlnPlot(Aur.combined, pt.size = 0, features = c("PTPRC","HLA-DPB1","HLA-DQA1","CLEC10A","FCER1A","CD1C") )

# Bcells-----
VlnPlot(Aur.combined, pt.size = 0, features = c("CD79A","CD79B","FCER2","CD22","MYH11","CD34") )
VlnPlot(Aur.combined, pt.size = 0, features = c("LUM","TAGLN","CD34","CDH5","PECAM1","TYROBP") )

# Tcells-----
VlnPlot(Aur.combined, pt.size = 0, features = c("CD3E","CD4","CD8A","NKG7","GNLY","GZMB") )
VlnPlot(Aur.combined, pt.size = 0, features = c("KLRB1","CCR6","CD1d1","MKI67","IL2RB","GZMB") )

VlnPlot(Aur.combined, pt.size = 0, features = c("CD3E","NKG7","CD8A","CD68","CD14","S100A8") )
VlnPlot(Aur.combined, pt.size = 0, features = c("HLA-DQA1","CD79A","CD79B","GNLY","TYROBP","KIT") )
VlnPlot(Aur.combined, pt.size = 0, features = c("MYH11","PECAM1","CDH5","ACTC1","TNNC1","MKI67") )

VlnPlot(Aur.combined, pt.size = 0 , features = c("CD68","CD14","IL1B","S100A8","VCAN","CD3E") )
VlnPlot(Aur.combined, pt.size = 0 ,  features = c("CD4","CD8A","NKG7","GNLY","CD79A","CD79B") )
VlnPlot(Aur.combined, pt.size = 0 ,  features = c("FCER2","KIT","HDC","MKI67") )

# コンタミ探し-----
VlnPlot(Aur.combined, pt.size = 0, features = c("NPPA","MYL7","MYL4","TTN","MYBPC3","TNNT2"))
VlnPlot(Aur.combined, pt.size = 0, features = c("RYR2","PLN","SLC8A1","DCN","GSN","PDGFRA"))
VlnPlot(Aur.combined, pt.size = 0, features = c("VWF","PECAM1","CDH5","RGS5","ABCC9","KCNJ8"))
VlnPlot(Aur.combined, pt.size = 0, features = c("MYH11","TAGLN","ACTA2","PLP1","NRXN1","NRXN3"))
VlnPlot(Aur.combined, pt.size = 0, features = c("GPAM","FASN","LEP","MSLN","WT1","BNC1"))
VlnPlot(Aur.combined, pt.size = 0, features = c("PTPRC","CD163","CD79B","CD14","IL7R","CD3E","CD4","CD8A","NKG7","GNLY","GZMB"))

FeaturePlot(Aur.combined, features=c('ACTC1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined, features=c('LUM'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined, features=c('CD34'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined, features=c('CD3E'), min.cutoff=0.5, max.cutoff='q90')

Aur.combined_trans <- SCTransform (Aur.combined)

cd_genes <- c("CD3E","NKG7","CD8A","CD68","CD14","S100A8","HLA-DQA1","CD79A","CD79B","GNLY","TYROBP","KIT", "MYH11","PECAM1","CDH5","ACTC1","TNNC1","MKI67")
DotPlot(Aur.combined_trans,features = cd_genes)+RotatedAxis()+coord_flip()

DimPlot(Aur.combined, reduction = "umap", label = TRUE,  repel = TRUE)

#CD45 MyocariumおよびEndothelial cellおよびSMC除去----
Aur.combined_CD45 <- subset(Aur.combined, idents = c(0,1,2,4,5,7,9,11,13,14,15,17,18,19))
Aur.combined_CD45 <- ScaleData(Aur.combined_CD45, verbose = FALSE)
Aur.combined_CD45<- RunPCA(Aur.combined_CD45, npcs = 30, verbose = FALSE)
Aur.combined_CD45 <- RunUMAP(Aur.combined_CD45 , reduction = "pca", dims = 1:30)
Aur.combined_CD45 <- FindNeighbors(Aur.combined_CD45 , reduction = "pca", dims = 1:30)
Aur.combined_CD45 <- FindClusters(Aur.combined_CD45, resolution = 0.1)




# Visualization 別々------
DimPlot(Aur.combined_CD45, reduction = "umap", label = FALSE, repel = TRUE)
DimPlot(Aur.combined_CD45, reduction = "umap", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined_CD45, reduction = "umap", split.by = "AFControl1ex", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined_CD45, reduction = "umap", split.by = "AFControl1ex", label = FALSE,  repel = FALSE)
DimPlot(Aur.combined_CD45, reduction = "umap", split.by = "AFDuration1ex")
DimPlot(Aur.combined_CD45, reduction = "umap", split.by = "sample1ex")
DimPlot(Aur.combined_CD45, reduction = "umap", split.by = "stim")

# RDS_CD45
saveRDS(Aur.combined_CD45, "~/Desktop/Human_LA_LAA/RDS/Aur_cluster_withcontorl_AurscCD45_最初.rds")
Aur.combined_CD45 <- readRDS ("~/Desktop/Human_LA_LAA/RDS/Aur_cluster_withcontorl_AurscCD45_最初.rds")

Aur.combined_CD45_trans <- SCTransform (Aur.combined_CD45)

cd_genes <- c("CD68","CD14","VCAN","CXCR2","CSF1R","KIT","HLA-DQA2","CLEC9A","CLEC4C","CD3E","NKG7","GNLY","MKI67","CD79A","FCER2")
DotPlot(Aur.combined_CD45_trans,features = cd_genes)+RotatedAxis()+coord_flip()




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

# コンタミ探し
VlnPlot(Aur.combined_CD45, pt.size = 0, features = c("NPPA","MYL7","MYL4","TTN","MYBPC3","TNNT2"))
VlnPlot(Aur.combined_CD45, pt.size = 0, features = c("RYR2","PLN","SLC8A1","DCN","GSN","PDGFRA"))
VlnPlot(Aur.combined_CD45, pt.size = 0, features = c("VWF","PECAM1","CDH5","RGS5","ABCC9","KCNJ8"))
VlnPlot(Aur.combined_CD45, pt.size = 0, features = c("MYH11","TAGLN","ACTA2","PLP1","NRXN1","NRXN3"))
VlnPlot(Aur.combined_CD45, pt.size = 0, features = c("GPAM","FASN","LEP","MSLN","WT1","BNC1"))
VlnPlot(Aur.combined_CD45, pt.size = 0, features = c("PTPRC","CD163","CD79B","CD14","IL7R","CD3E","CD4","CD8A","NKG7","GNLY","GZMB"))



Aur.combined_CD45_trans <- SCTransform (Aur.combined_CD45)

cd_genes <- c("CD3E","NKG7","CD8A","CD68","CD14","S100A8","HLA-DQA1","CD79A","CD79B","GNLY","TYROBP","KIT", "MYH11","PECAM1","CDH5","ACTC1","TNNC1","MKI67","PLP1","NRXN1","NRXN3")
DotPlot(Aur.combined_CD45_trans,features = cd_genes)+RotatedAxis()+coord_flip()




# ヒートマップの作成------------------------------------------------------------------------------
Aur.combined_CD45.markers <- FindAllMarkers(Aur.combined_CD45, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Aur.combined_CD45.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- Aur.combined_CD45.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Aur.combined_CD45, features = top10$gene) + NoLegend()
write.table(top10, "table.txt", quote=F, col.names=F, append=T)

#CD45 final(Neuronal除去)----
Aur.combined_CD45_1 <- subset(Aur.combined_CD45, idents = c(0,1,2,3,4,5,6,7,8,9,10,11,12,14))
Aur.combined_CD45_1 <- ScaleData(Aur.combined_CD45_1, verbose = FALSE)
Aur.combined_CD45_1<- RunPCA(Aur.combined_CD45_1, npcs = 30, verbose = FALSE)
Aur.combined_CD45_1 <- RunUMAP(Aur.combined_CD45_1 , reduction = "pca", dims = 1:30)
Aur.combined_CD45_1 <- FindNeighbors(Aur.combined_CD45_1 , reduction = "pca", dims = 1:30)
Aur.combined_CD45_1 <- FindClusters(Aur.combined_CD45_1, resolution = 0.1)

### AFControl1ex ##
# stim列を抽出する
cell.stim <- Aur.combined_CD45_1@meta.data[["stim"]]
# 新しいグループ（疾患ありなしの2群）とstim列との対応関係を定義する
newgroup <- c("3.Other","3.Other","2.AF","2.AF","2.AF","2.AF","2.AF","2.AF","2.AF","2.AF","2.AF","2.AF","1.Control","1.Control","1.Control","1.Control","1.Control","1.Control")
names(newgroup) <- c("Aur1Endo","Aur1Epi","Aur2","Aur3Endo","Aur3Epi","Aur4Endo","Aur4Epi","Aur12","Aur13","Aur10","Aur16","Aur19","Control1","Control2","Control3","Control4","Control5","Control6")
# 抽出したstim列を新しいグループに変換する
cell.newgroup <- newgroup[cell.stim]
# 新しいグループを新しいMetaDataとして追加する
Aur.combined_CD45_1 <- AddMetaData(Aur.combined_CD45_1, cell.newgroup, col.name = "AFControl1ex")
### Aur1ex ##
# stim列を抽出する
cell.stim <- Aur.combined_CD45_1@meta.data[["stim"]]
# 新しいグループ（疾患ありなしの2群）とstim列との対応関係を定義する
newgroup <- c("2.Other","2.Other","1.Aur","1.Aur","1.Aur","1.Aur","1.Aur","1.Aur","1.Aur","1.Aur","1.Aur","1.Aur","1.Aur","1.Aur","1.Aur","1.Aur","1.Aur","1.Aur")
names(newgroup) <- c("Aur1Endo","Aur1Epi","Aur2","Aur3Endo","Aur3Epi","Aur4Endo","Aur4Epi","Aur12","Aur13","Aur10","Aur16","Aur19","Control1","Control2","Control3","Control4","Control5","Control6")
# 抽出したstim列を新しいグループに変換する
cell.newgroup <- newgroup[cell.stim]
# 新しいグループを新しいMetaDataとして追加する
Aur.combined_CD45_1 <- AddMetaData(Aur.combined_CD45_1, cell.newgroup, col.name = "Aur1ex")



# Visualization 別々------
DimPlot(Aur.combined_CD45_1, reduction = "umap", label = FALSE, repel = TRUE)
DimPlot(Aur.combined_CD45_1, reduction = "umap", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined_CD45_1, reduction = "umap", split.by = "Aur1ex")
DimPlot(Aur.combined_CD45_1, reduction = "umap", split.by = "AFControl1ex", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined_CD45_1, reduction = "umap", split.by = "AFControl1ex", label = FALSE,  repel = FALSE)
DimPlot(Aur.combined_CD45_1, reduction = "umap", split.by = "AFDuration1ex")
DimPlot(Aur.combined_CD45_1, reduction = "umap", split.by = "sample1ex1cotrlex")
DimPlot(Aur.combined_CD45_1, reduction = "umap", split.by = "stim")


test <- subset(Aur.combined_CD45_1, subset = sample1ex1cotrlex == "9.Control1.1")
test

# RDS_CD45
saveRDS(Aur.combined_CD45_1, "~/Desktop/Human_LA_LAA/RDS/Aur_cluster_withcontorl_AurscCD45_1_revision20241018.rds")
DimPlot(Aur.combined_CD45_1, reduction = "umap", split.by = "AFDuration1ex")


Aur.combined_CD45_1 <- readRDS ("~/Desktop/Human_LA_LAA/RDS/Aur_cluster_withcontorl_AurscCD45_1.rds")



Aur.combined_CD45_1_trans <- SCTransform (Aur.combined_CD45_1)

cd_genes <- c("CD68","CD14","VCAN","CXCR2","CSF1R","KIT","HLA-DQA2","CLEC9A","CLEC4C","CD3E","NKG7","GNLY","MKI67","CD79A","FCER2")
DotPlot(Aur.combined_CD45_1_trans,features = cd_genes)+RotatedAxis()+coord_flip()

#並べ替え
levels(Aur.combined_CD45_1) # [1] "0""1""2""3""4""5""6""7"8""9""10""11""12"
levels(Aur.combined_CD45_1) <- c("2", "3", "5","10","11","12", "0", "1", "6", "7","9","4","8")
new.cluster.ids <- c("0","1","2","3","4", "5", "6","7","8","9","10","11","12")
names(new.cluster.ids) <- levels(Aur.combined_CD45_1)
Aur.combined_CD45_1<- RenameIdents(Aur.combined_CD45_1, new.cluster.ids)




# Myeloid----
VlnPlot(Aur.combined_CD45_1, pt.size = 0, features = c("CD68","CD14","CSF1R","IL1B","S100A8", "S100A9") )
VlnPlot(Aur.combined_CD45_1, pt.size = 0, features = c("TREM2","CD9","APOE","CXCL3","TNF","CXCL8") )
VlnPlot(Aur.combined_CD45_1, pt.size = 0, features = c("HLA-DPA1","HLA-DPB1","HLA-DQA1","CLEC10A","FCER1A","CD1C") )
VlnPlot(Aur.combined_CD45_1, pt.size = 0, features = c("PTPRC","HLA-DPB1","HLA-DQA1","CLEC10A","FCER1A","CD1C") )

# Bcells-----
VlnPlot(Aur.combined_CD45_1, pt.size = 0, features = c("CD79A","CD79B","FCER2","CD22","MYH11","CD34") )
VlnPlot(Aur.combined_CD45_1, pt.size = 0, features = c("LUM","TAGLN","CD34","CDH5","PECAM1","TYROBP") )

# Tcells-----
VlnPlot(Aur.combined_CD45_1, pt.size = 0, features = c("CD3E","CD4","CD8A","NKG7","GNLY","GZMB") )
VlnPlot(Aur.combined_CD45_1, pt.size = 0, features = c("KLRB1","CCR6","CD1d1","MKI67","IL2RB","GZMB") )


VlnPlot(Aur.combined_CD45_1, pt.size = 0, features = c("CLEC10A","FCER1A","CD1D","FCGR3A","CXCR2") )

VlnPlot(Aur.combined_CD45_1, pt.size = 0, features = c("ACTC1","TNNC1","FCER2","CD22","MYH11","CD34") )
VlnPlot(Aur.combined_CD45_1, pt.size = 0, features = c("LUM","TAGLN","CD34","CDH5","PECAM1","TYROBP") )

VlnPlot(Aur.combined_CD45_1, pt.size = 0 ,  features = c("FCER2","KIT","HDC","MKI67") )

cDC1
VlnPlot(Aur.combined_CD45_1, pt.size = 0 , features = c("THBD","CLEC9A","XCR1","CADM1","IRF4","IRF8","BATF3") )
cDC2
VlnPlot(Aur.combined_CD45_1, pt.size = 0 , features = c("THBD","CD1C","CLEC10A","FCER1A","IRF4","CD2") )
pDC
VlnPlot(Aur.combined_CD45_1, pt.size = 0 , features = c("IL3RA","CLEC4C") )

VlnPlot(Aur.combined_CD45_1,pt.size = 0, features = c("LYVE1","HLA-DOA","HLA-DQA1","HLA-DQA2","HLA-DQB1","TIMD4","TREM2") )

Manuscript
VlnPlot(Aur.combined_CD45_1,pt.size = 0, features = c("CD68","CD14","IL1B","S100A8","VCAN","HLA-DPA1") )
VlnPlot(Aur.combined_CD45_1,pt.size = 0, features = c("CD3E","NKG7","GNLY","CD79A","KIT","MKI67") )

VlnPlot(Aur.combined_CD45_1,pt.size = 0, features = c("CD68","CD14","CSF1R","CD163","LYVE1","CD9","TREM2","IL1B","S100A8","VCAN","LYZ","FCGR3A","CXCR2","FCGR3B","CSF3R","HLA-DPA1","CLEC9A","CD1C","IL3RA"))

# コンタミ探し
VlnPlot(Aur.combined_CD45_1, pt.size = 0, features = c("NPPA","MYL7","MYL4","TTN","MYBPC3","TNNT2"))
VlnPlot(Aur.combined_CD45_1, pt.size = 0, features = c("RYR2","PLN","SLC8A1","DCN","GSN","PDGFRA"))
VlnPlot(Aur.combined_CD45_1, pt.size = 0, features = c("VWF","PECAM1","CDH5","RGS5","ABCC9","KCNJ8"))
VlnPlot(Aur.combined_CD45_1, pt.size = 0, features = c("MYH11","TAGLN","ACTA2","PLP1","NRXN1","NRXN3"))
VlnPlot(Aur.combined_CD45_1, pt.size = 0, features = c("GPAM","FASN","LEP","MSLN","WT1","BNC1"))
VlnPlot(Aur.combined_CD45_1, pt.size = 0, features = c("PTPRC","CD163","CD79B","CD14","IL7R","CD3E","CD4","CD8A","NKG7","GNLY","GZMB"))

Aur.combined_CD45_1_trans <- SCTransform (Aur.combined_CD45_1)


cd_genes <- c("CD68","CD14","CSF1R","VCAN","FCGR3B","CXCR2","KIT","HDC","HLA-DPA1","CLEC9A","CLEC4C","CD3E","NKG7","GNLY","MKI67","CD79A","FCER2")
DotPlot(Aur.combined_CD45_1_trans,features = cd_genes)+RotatedAxis()+coord_flip()


Aur.combined_CD45_1_prop. <- prop.table(table(Idents(Aur.combined_CD45_1),Aur.combined_CD45_1$AFControl1ex))
write.table(Aur.combined_CD45_1_prop., "Aur.combined_CD45_1_prop.txt", quote=F, col.names=F, append=T)
write.csv(Aur.combined_CD45_1_prop., "Aur.combined_1_1_prop..csv")




# ヒートマップの作成------------------------------------------------------------------------------
Aur.combined_CD45_1.markers <- FindAllMarkers(Aur.combined_CD45_1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Aur.combined_CD45_1.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- Aur.combined_CD45_1.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Aur.combined_CD45_1, features = top10$gene) + NoLegend()
write.table(top10, "table.txt", quote=F, col.names=F, append=T)


#myeloid-------
Aur.combined_1 <- subset(Aur.combined_CD45_1, idents = c(0,1,2,3,4,5))
Aur.combined_1 <- ScaleData(Aur.combined_1, verbose = FALSE)
Aur.combined_1<- RunPCA(Aur.combined_1, npcs = 30, verbose = FALSE)
Aur.combined_1 <- RunUMAP(Aur.combined_1 , reduction = "pca", dims = 1:30)
Aur.combined_1 <- FindNeighbors(Aur.combined_1 , reduction = "pca", dims = 1:30)
Aur.combined_1 <- FindClusters(Aur.combined_1, resolution = 0.5)

# RDS_Myeloid Aur.combined_1
saveRDS(Aur.combined_1, "~/Desktop/Human_LA_LAA/RDS/Aur_cluster_withcontorl_AurscMy_1.rds")
Aur.combined_1 <- readRDS ("~/Desktop/Human_LA_LAA/RDS/Aur_cluster_withcontorl_AurscMy_1.rds")

# Visualization 別々------
DimPlot(Aur.combined_1, reduction = "umap", label = FALSE, repel = TRUE)
DimPlot(Aur.combined_1, reduction = "umap", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined_1, reduction = "umap", split.by = "Aur1ex")
DimPlot(Aur.combined_1, reduction = "umap", split.by = "AFControl1ex", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined_1, reduction = "umap", split.by = "AFControl1ex", label = FALSE,  repel = FALSE)
DimPlot(Aur.combined_1, reduction = "umap", split.by = "AFDuration1ex")
DimPlot(Aur.combined_1, reduction = "umap", split.by = "sample1ex")
DimPlot(Aur.combined_1, reduction = "umap", split.by = "stim")

Aur.combined_1_trans <- SCTransform (Aur.combined_1)

cd_genes <- c("CD68","CD14","CSF1R","LYVE1","TREM2","CD9","IL1B","TNF","S100A8","VCAN","FCN1","FCGR3A","CXCR2","HLA-DPA1","CLEC9A","FCER1A")
DotPlot(Aur.combined_1_trans,features = cd_genes)+RotatedAxis()+coord_flip()

cd_genes <- c("CD68","CD14","CSF1R","CD163","C1QA","IL1B","AREG","EREG","TIMP1","THBS1","CLEC5A","IL10","VEGFA","VEGFB","TREM2","CD9","PDGFB","PDGFC","IGF1","LYVE1","FOLR2","S100A8","S100A12","LYZ","FCGR3A","ISG15","HLA-DPA1","HLA-DRA","TNF","CCL3","CCL4","CLEC10A","CLEC4A","CD1C","CLEC9A","IL3RA","CCR7","FSCN1","CCL19","CXCR2","FCGR3B","CSF3R","KIT","HDC","MYL7","MYL4","TAGLN","ACTA2","CD79A","CD79B")
DotPlot(Aur.combined_1_trans,features = cd_genes)+RotatedAxis()+coord_flip()

# ヒートマップの作成------------------------------------------------------------------------------
Aur.combined_1.markers <- FindAllMarkers(Aur.combined_1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Aur.combined_1.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- Aur.combined_1.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Aur.combined_1, features = top10$gene) + NoLegend()
write.table(top10, "table.txt", quote=F, col.names=F, append=T)

# コンタミ探し
VlnPlot(Aur.combined_1, pt.size = 0, features = c("NPPA","MYL7","MYL4","TTN","MYBPC3","TNNT2"))
VlnPlot(Aur.combined_1, pt.size = 0, features = c("RYR2","PLN","SLC8A1","DCN","GSN","PDGFRA"))
VlnPlot(Aur.combined_1, pt.size = 0, features = c("VWF","PECAM1","CDH5","RGS5","ABCC9","KCNJ8"))
VlnPlot(Aur.combined_1, pt.size = 0, features = c("MYH11","TAGLN","ACTA2","PLP1","NRXN1","NRXN3"))
VlnPlot(Aur.combined_1, pt.size = 0, features = c("GPAM","FASN","LEP","MSLN","WT1","BNC1"))
VlnPlot(Aur.combined_1, pt.size = 0, features = c("PTPRC","CD163","CD79A","CD79B","FCER2","CD14","IL7R","CD3E","CD4","CD8A","NKG7","GNLY","GZMB"))

FeaturePlot(Aur.combined_1, features=c('MYL7'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1, features=c('SLC8A1'), min.cutoff=0.2, max.cutoff='q90')

FeaturePlot(Aur.combined_1, features=c('RGS5'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1, features=c('KCNJ8'), min.cutoff=0.2, max.cutoff='q90')
FeaturePlot(Aur.combined_1, features=c('TAGLN'), min.cutoff=0.2, max.cutoff='q90')



VlnPlot(Aur.combined_1, pt.size = 0, features = c("CD3E","NKG7","CD8A","CD68","CD14","S100A8") )
VlnPlot(Aur.combined_1, pt.size = 0, features = c("HLA-DQA1","CD79A","CD79B","GNLY","TYROBP","KIT") )
VlnPlot(Aur.combined_1, pt.size = 0, features = c("MYH11","PECAM1","CDH5","ACTC1","TNNC1","MKI67") )

VlnPlot(Aur.combined_1, pt.size = 0, features = c("CD79A","CD79B","FCER2","CD22","MYH11","CD34") )
VlnPlot(Aur.combined_1, pt.size = 0, features = c("LUM","TAGLN","CD34","CDH5","PECAM1","TYROBP") )


VlnPlot(Aur.combined_1, pt.size = 0, features = c("CD68","S100A8","CD14","CSF1R","CD163","LYVE1"))
VlnPlot(Aur.combined_1, pt.size = 0, features = c("CD9","C1QA","C1QB","TREM2","IL1B","TNF"))
VlnPlot(Aur.combined_1, pt.size = 0, features = c("CXCL3","CXCL8","TLR2","TLR4","VCAN","FCN1"))
VlnPlot(Aur.combined_1, pt.size = 0, features = c("LYZ","FCGR3A","CSF3R","FCGR3B","CLEC10A","FCER1A"))
VlnPlot(Aur.combined_1, pt.size = 0, features = c("CD1D","HLA-DRA","HLA-DQA1","HLA-DPA1","CD1C","CD2"))
VlnPlot(Aur.combined_1, pt.size = 0, features = c("CLEC9A","XCR1","CADM1","CLEC4C","IL3RA","SPP1"))





#Myeloid(MyocariumおよびEC/PCおよびT・Bcell doublet除去)--------
Aur.combined_1_1 <- subset(Aur.combined_1, idents = c(0,1,2,3,4,5,6,7,8,13,14,16,17,18))
Aur.combined_1_1 <- ScaleData(Aur.combined_1_1, verbose = FALSE)
Aur.combined_1_1<- RunPCA(Aur.combined_1_1, npcs = 30, verbose = FALSE)
Aur.combined_1_1 <- RunUMAP(Aur.combined_1_1 , reduction = "pca", dims = 1:30)
Aur.combined_1_1 <- FindNeighbors(Aur.combined_1_1 , reduction = "pca", dims = 1:30)
Aur.combined_1_1 <- FindClusters(Aur.combined_1_1, resolution = 0.50)


# RDS_Myeloid Aur.combined_1_1
saveRDS(Aur.combined_1_1, "~/Desktop/Human_LA_LAA/RDS/Aur_cluster_withcontorl_AurscMy_1_1.rds")
Aur.combined_1_1 <- readRDS ("~/Desktop/Human_LA_LAA/RDS/Aur_cluster_withcontorl_AurscMy_1_1.rds")


# Visualization 別々------
DimPlot(Aur.combined_1_1, reduction = "umap", label = FALSE, repel = TRUE)
DimPlot(Aur.combined_1_1, reduction = "umap", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined_1_1, reduction = "umap", split.by = "Aur1ex")
DimPlot(Aur.combined_1_1, reduction = "umap", split.by = "AFControl1ex", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined_1_1, reduction = "umap", split.by = "AFControl1ex", label = FALSE,  repel = FALSE)
DimPlot(Aur.combined_1_1, reduction = "umap", split.by = "AFDuration1ex")
DimPlot(Aur.combined_1_1, reduction = "umap", split.by = "sample1ex")
DimPlot(Aur.combined_1_1, reduction = "umap", split.by = "stim")


# 並べ替える未
levels(Aur.combined_1_1) # [1] "0" "1" "2" "3" "4" "5" "6" "7" "8" 
levels(Aur.combined_1_1) <- c("0", "1", "4", "3", "5", "2", "6", "7", "8" )

new.cluster.ids <- c("0","1","2","3","4", "5","6","7","8")
names(new.cluster.ids) <- levels(Aur.combined_1_1)
Aur.combined_1_1<- RenameIdents(Aur.combined_1_1, new.cluster.ids)


Aur.combined_1_1_trans <- SCTransform (Aur.combined_1_1)
cd_genes <- c("CD68","CD14","CSF1R","CD163","C1QA","IL1B","AREG","EREG","TIMP1","THBS1","CLEC5A","IL10","VEGFA","VEGFB","TREM2","CD9","PDGFB","PDGFC","IGF1","LYVE1","FOLR2","S100A8","S100A12","LYZ","FCGR3A","ISG15","HLA-DPA1","HLA-DRA","TNF","CCL3","CCL4","CLEC10A","CLEC4A","CD1C","CLEC9A","IL3RA","CCR7","FSCN1","CCL19","CXCR2","FCGR3B","CSF3R","KIT","HDC","CD3E","CD79A","DCN","GSN","PDGFRA","NPPA","MYL7","MYL4","VWF","PECAM1","CDH5","RGS5","ABCC9","KCNJ8","MYH11","TAGLN","ACTA2","MKI67")
DotPlot(Aur.combined_1_1_trans,features = cd_genes)+RotatedAxis()+coord_flip()


"DCN","GSN","PDGFRA","NPPA","MYL7","MYL4","VWF","PECAM1","CDH5","RGS5","ABCC9","KCNJ8","MYH11","TAGLN","ACTA2","CD3E","CD4","CD8A","NKG7","GNLY","GZMB","IL7R","CD40LG","CD14","C1QA","CD68","HLA-DPA1","HLA-DPB1","HLA-DQA1","CLEC10A","FCER1A","CD1C","MKI67","CD79A","CD79B","KIT","HDC"


# 新しいグループ分けで可視化
DimPlot(Aur.combined_1_1, reduction = "umap", label = FALSE, repel = TRUE)
DimPlot(Aur.combined_1_1, reduction = "umap", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined_1_1, reduction = "umap", split.by = "AFControl")
DimPlot(Aur.combined_1_1, reduction = "umap", split.by = "stim")
DimPlot(Aur.combined_1_1, reduction = "umap", split.by = "AFDuration")





# ヒートマップの作成------------------------------------------------------------------------------
Aur.combined_1_1.markers <- FindAllMarkers(Aur.combined_1_1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Aur.combined_1_1.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- Aur.combined_1_1.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Aur.combined_1_1, features = top10$gene) + NoLegend()
write.table(top10, "table.txt", quote=F, col.names=F, append=T)



VlnPlot(Aur.combined_1_1, pt.size = 0, features = c("CD3E","NKG7","CD8A","CD68","CD14","S100A8") )
VlnPlot(Aur.combined_1_1, pt.size = 0, features = c("HLA-DQA1","CD79A","CD79B","GNLY","TYROBP","KIT") )
VlnPlot(Aur.combined_1_1, pt.size = 0, features = c("MYH11","PECAM1","CDH5","ACTC1","TNNC1","MKI67") )

VlnPlot(Aur.combined_1_1, pt.size = 0, features = c("CD79A","CD79B","FCER2","CD22","MYH11","CD34") )
VlnPlot(Aur.combined_1_1, pt.size = 0, features = c("LUM","TAGLN","CD34","CDH5","PECAM1","TYROBP") )


VlnPlot(Aur.combined_1_1, pt.size = 0, features = c("CD68","S100A8","CD14","CSF1R","CD163","LYVE1"))
VlnPlot(Aur.combined_1_1, pt.size = 0, features = c("CD9","C1QA","C1QB","TREM2","IL1B","TNF"))
VlnPlot(Aur.combined_1_1, pt.size = 0, features = c("CXCL3","CXCL8","TLR2","TLR4","VCAN","FCN1"))
VlnPlot(Aur.combined_1_1, pt.size = 0, features = c("LYZ","FCGR3A","CSF3R","FCGR3B","CLEC10A","FCER1A"))
VlnPlot(Aur.combined_1_1, pt.size = 0, features = c("CD1D","HLA-DRA","HLA-DQA1","HLA-DPA1","CD1C","CD2"))
VlnPlot(Aur.combined_1_1, pt.size = 0, features = c("CLEC9A","XCR1","CADM1","CLEC4C","IL3RA","SPP1"))




VlnPlot(Aur.combined_1_1, features = c("FCGR3A","FCGR2A","ITGB2","CD16","CEACAM8","CD55","CD44","ELANE") )
VlnPlot(Aur.combined_1_1, features = c("PTPRC","CD14","MPO","LYZ","ITGB2","CD16","FCGR3A","FCGR3B","CD55","CD44","ELANE","VCAN", "CD52","FCN1","S100A12","IL1B","CSF1R","S100A8","S100A9","CAMP","DEFB4A","CTSS","ITGAM","CD68","CYBB") )
VlnPlot(Aur.combined_1_1, features = c("DEFB1","DEFA1","DEFB103B","S100A9","LYZ") )
VlnPlot(Aur.combined_1_1, features = c("S100A6","FOS","S100A12","S100A9","HMGB2","KLF6"))



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

VlnPlot(Aur.combined_1_1, features = c("CD3E","CD4","CD8A","CCR2") )

LYVE1+ macrophage
VlnPlot(Aur.combined_1_1,pt.size = 0, features = c("LYVE1","HLA-DOA","HLA-DQA1","HLA-DQA2","HLA-DQB1","TIMD4") )

Monocyte-derived macrophage
VlnPlot(Aur.combined_1_1,pt.size = 0, features = c("LYVE1","FOLR2","CEBPB","S100A8","CCL13","CCL18") )

Antigen-presenting macrophage
VlnPlot(Aur.combined_1_1,pt.size = 0, features = c("FOLR2","LYVE1","MERTK","HLA-DRA","HLA-DMA","HLA-DMB","HLA-DPA1","TREM2") )

Dock4+ macrophage
VlnPlot(Aur.combined_1_1,pt.size = 0, features = c("DOCK4","IL4R","STAT3","ITGAM","C1QA","FOLR2") )

fibrotic
VlnPlot(Aur.combined_1_1,pt.size = 0, features = c("TREM2","CD9","C1QA") )


VlnPlot(Aur.combined_1_1, features = c("CD163","MKI67","TRPG1","FGF13","CCL7") )
VlnPlot(Aur.combined_1_1, features = c("CD163","MKI67","TRPG1","FGF13","CCL7") )



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
VlnPlot(Aur.combined_1_1,pt.size = 0,  features = c("TET2","DNMT3A","ASXL1","PPM1D","NLRP3","STING1","CGAS" ))
VlnPlot(Aur.combined_1_1, features = c("OLR1","SCARB1","SCARB2","MSR1","MARCO","SRSF2") )
VlnPlot(Aur.combined_1_1, features = c("NR1H3","CTSD","CTSL","SPP1","MARCO","FABP4") )
VlnPlot(Aur.combined_1_1, features = c("TET2","TET1","TET3","SPP1","MARCO","SLC25A44") )
VlnPlot(Aur.combined_1_1, features = c("TNF","TNFSF13"),split.by = "stim",  cols = c("sa" ="deepskyblue", "Aur1Epi"="red"))
VlnPlot(Aur.combined_1_1, features = c("CD163","F13A1","MS4A4A","SPP1","MARCO","SLC25A44") )

VlnPlot(Aur.combined_1_1, features = c("MYL7","ACTC1","TNNT2","MYL4","ACTC1","TNNC1") )

DC

VlnPlot(Aur.combined_1_1, features = c("HLA-DPA1","HLA-DPB1","HLA-DQA1","CLEC10A","FCER1A","CD1C") )
VlnPlot(Aur.combined_1_1, features = c("ITGAX","IFNG","CLEC9A","CADM1","CLEC10A","FCER1A") )

VlnPlot(Aur.combined_1_1, features = c("THBD","CLEC9A","XCR1","CADM1","CD1C","CLEC10A") )
VlnPlot(Aur.combined_1_1, features = c("FCER1A","CD2","IL3RA","CLEC4C") )



cDC1→My.8
VlnPlot(Aur.combined_1_1, features = c("THBD","CLEC9A","XCR1","CADM1","IRF4","IRF8","BATF3") )
cDC2→My.7
VlnPlot(Aur.combined_1_1, features = c("THBD","CD1C","CLEC10A","FCER1A","IRF4","CD2") )
pDC→My.9
VlnPlot(Aur.combined_1_1, features = c("IL3RA","CLEC4C",) )

AXL..7

VlnPlot(Aur.combined_1_1, features = c("CD1C","CD141","CD163","CD14","CD5") )
VlnPlot(Aur.combined_1_1, features = c("CCR7","LAMP3","THBD","AXL","SIGLEC6") )

☓
VlnPlot(Aur.combined_1_1, features = c("CD141","CD123","CD303","CD2","CD172A") )

Neutrophil
VlnPlot(Aur.combined_1_1, features = c("CSF3R","FCGR3B","NAMPT","CXCR2","CCR2") )
VlnPlot(Aur.combined_1_1, features = c("CSF3R","FCGR3B","NAMPT","CXCR2","CCR2") )
VlnPlot(Aur.combined_1_1,pt.size = 0,  features = c("NAMPT","SRGN","TNFRSF10C","CD14","CLEC4E","MXD1","CGAS" ))


Aur.combined_1_1_trans <- SCTransform (Aur.combined_1_1)

cd_genes <- c("CD68","CD14","CSF1R","LYVE1","TREM2","CD9","IL1B","TNF","S100A8","VCAN","FCN1","FCGR3A","FCGR3B","CXCR2","HLA-DPA1","CLEC9A","CD1C","IL3RA")
DotPlot(Aur.combined_1_1_trans,features = cd_genes)+RotatedAxis()+coord_flip()

New
Aur.combined_1_1_trans <- SCTransform (Aur.combined_1_1)

cd_genes <- c("CD68","CD14","CSF1R","LYVE1","TREM2","CD9","IL1B","TNF","S100A8","VCAN","FCN1","LYZ","FCGR3A","FCGR3B","CSF3R","HLA-DPA1","CLEC9A","CD1C","IL3RA")
DotPlot(Aur.combined_1_1_trans,features = cd_genes)+RotatedAxis()+coord_flip()



Aur.combined_1_1_prop. <- prop.table(table(Idents(Aur.combined_1_1),Aur.combined_1_1$AFControl))
write.table(Aur.combined_1_1_prop., "Aur.combined_1_1_prop..txt", quote=F, col.names=F, append=T)
write.csv(Aur.combined_1_1_prop., "Aur.combined_1_1_prop..csv")

Aur.combined_1_1_prop. <- prop.table(table(Idents(Aur.combined_1_1),Aur.combined_1_1$AFDuration))
write.table(Aur.combined_1_1_prop., "Aur.combined_1_1_prop..txt", quote=F, col.names=F, append=T)
write.csv(Aur.combined_1_1_prop., "Aur.combined_1_1_prop..csv")







#macs+mono--------
Aur.combined_1_2 <- subset(Aur.combined_1_1, idents = c(0,1,3,5,7,8,9,15,16,18))
Aur.combined_1_2 <- ScaleData(Aur.combined_1_2, verbose = FALSE)
Aur.combined_1_2<- RunPCA(Aur.combined_1_2, npcs = 30, verbose = FALSE)
Aur.combined_1_2 <- RunUMAP(Aur.combined_1_2 , reduction = "pca", dims = 1:30)
Aur.combined_1_2 <- FindNeighbors(Aur.combined_1_2 , reduction = "pca", dims = 1:30)
Aur.combined_1_2 <- FindClusters(Aur.combined_1_2, resolution = 0.4)

# RDS_Myeloid Aur.combined_1_2
saveRDS(Aur.combined_1_2, "~/Desktop/Human_LA_LAA/RDS/Aur_cluster_withcontorl_AurscMy_1_2.rds")
Aur.combined_1_2 <- readRDS ("~/Desktop/Human_LA_LAA/RDS/Aur_cluster_withcontorl_AurscMy_1_2.rds")


# Visualization 別々------
DimPlot(Aur.combined_1_2, reduction = "umap", label = FALSE, repel = TRUE)
DimPlot(Aur.combined_1_2, reduction = "umap", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined_1_2, reduction = "umap", split.by = "Aur1ex")
DimPlot(Aur.combined_1_2, reduction = "umap", split.by = "AFControl1ex", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined_1_2, reduction = "umap", split.by = "AFControl1ex", label = FALSE,  repel = FALSE)
DimPlot(Aur.combined_1_2, reduction = "umap", split.by = "AFDuration1ex")
DimPlot(Aur.combined_1_2, reduction = "umap", split.by = "sample1ex")
DimPlot(Aur.combined_1_2, reduction = "umap", split.by = "stim")


# コンタミ探し
VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("NPPA","MYL7","MYL4","TTN","MYBPC3","TNNT2"))
VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("RYR2","PLN","SLC8A1","DCN","GSN","PDGFRA"))
VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("VWF","PECAM1","CDH5","RGS5","ABCC9","KCNJ8"))
VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("MYH11","TAGLN","ACTA2","PLP1","NRXN1","NRXN3"))
VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("GPAM","FASN","LEP","MSLN","WT1","BNC1"))
VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("PTPRC","CD163","CD79A","CD79B","FCER2","CD14","IL7R","CD3E","CD4","CD8A","NKG7","GNLY","GZMB"))



# ヒートマップの作成------------------------------------------------------------------------------
Aur.combined_1_2.markers <- FindAllMarkers(Aur.combined_1_2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Aur.combined_1_2.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- Aur.combined_1_2.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Aur.combined_1_2, features = top10$gene) + NoLegend()
write.table(top10, "table.txt", quote=F, col.names=F, append=T)


Aur.combined_1_2_trans <- SCTransform (Aur.combined_1_2)

cd_genes <- c("CD68","CD14","CSF1R","CD163","LYVE1","FOLR2","IL1B","TLR4","TIMP1","IL10","CD9","C1QA","TREM2","S100A8","VCAN","LYZ","FCGR3A")
DotPlot(Aur.combined_1_2_trans,features = cd_genes)+RotatedAxis()+coord_flip()

cd_genes <- c("CD68","CD14","CSF1R","CD163","C1QA","IL1B","AREG","EREG","TIMP1","THBS1","CLEC5A","IL10","VEGFA","VEGFB","TREM2","CD9","PDGFB","PDGFC","IGF1","LYVE1","FOLR2","S100A8","S100A12","LYZ","FCGR3A","ISG15","IFIT3","IFIT5","HLA-DPA1","HLA-DRA","TNF","CCL3","CCL4","CLEC10A","CLEC4A","CD1C","CLEC9A","IL3RA","CCR7","FSCN1","CCL19","CXCR2","FCGR3B","CSF3R","KIT","HDC","CD3E","CD79A","DCN","GSN","PDGFRA","NPPA","MYL7","MYL4","VWF","PECAM1","CDH5","RGS5","ABCC9","KCNJ8","MYH11","TAGLN","ACTA2","MKI67")
DotPlot(Aur.combined_1_2_trans,features = cd_genes)+RotatedAxis()+coord_flip()



DimPlot(AF, reduction = "umap", label = FALSE, repel = TRUE)
DimPlot(Aur.combined_1_2, reduction = "umap", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined_1_2, reduction = "umap", split.by = "AFControl",cols = c("0" = "#F8766D","1" = "#DE8C00","2" = "#B79F00","3" = "#00BA3B","4" = "#00C08B"))
DimPlot(Aur.combined_1_2, reduction = "umap", split.by = "stim")
DimPlot(Aur.combined_1_2, reduction = "umap", split.by = "AFDuration")



VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("CD3E","NKG7","CD8A","CD68","CD14","S100A8") )
VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("HLA-DQA1","CD79A","CD79B","GNLY","TYROBP","KIT") )
VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("MYH11","PECAM1","CDH5","ACTC1","TNNC1","MKI67") )

VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("CD79A","CD79B","FCER2","CD22","MYH11","CD34") )
VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("LUM","TAGLN","CD34","CDH5","PECAM1","TYROBP") )


VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("TIMP1","IL10","AREG","CLEC5A","PECAM1","TYROBP") )

Manuscript
VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("CD68","S100A8","CD14","CSF1R","CD163","LYVE1"))
VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("CD9","C1QA","C1QB","TREM2","IL1B","TNF"))
VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("CXCL3","CXCL8","TLR2","TLR4","VCAN","FCN1"))
VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("LYZ","FCGR3A","CSF3R","FCGR3B","CLEC10A","FCER1A"))
VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("CD1D","HLA-DRA","HLA-DQA1","HLA-DPA1","CD1C","CD2"))
VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("CLEC9A","XCR1","CADM1","CLEC4C","IL3RA","SPP1"))
VlnPlot(Aur.combined_1_2, pt.size = 0, features = c("TIMP1","IL10","ISG15","PDGFB"))

FeaturePlot(Aur.combined_1_2, features=c('LYVE1'), min.cutoff=0.4, max.cutoff='q90')
FeaturePlot(Aur.combined_1_2, features=c('LYVE1'), min.cutoff=0.4, max.cutoff='q90',split.by = "AFControl")
FeaturePlot(Aur.combined_1_2, features=c('CD163'), min.cutoff=0.4, max.cutoff='q90',split.by = "AFControl")
FeaturePlot(Aur.combined_1_2, features=c('AREG'), min.cutoff=0.4, max.cutoff='q90',split.by = "AFControl")
FeaturePlot(Aur.combined_1_2, features=c('AREG'), min.cutoff=0.4, max.cutoff='q90')
FeaturePlot(Aur.combined_1_2, features=c('EREG'), min.cutoff=0.4, max.cutoff='q90',split.by = "AFControl")
FeaturePlot(Aur.combined_1_2, features=c('EREG'), min.cutoff=0.4, max.cutoff='q90')
FeaturePlot(Aur.combined_1_2, features=c('TIMP1'), min.cutoff=2.0, max.cutoff='q90')

FeaturePlot(Aur.combined_1_2, features=c('ACTC1'), min.cutoff=0.2, max.cutoff='q90')
FeaturePlot(Aur.combined_1_2, features=c('PDGFB'), min.cutoff=0.4, max.cutoff='q90',split.by = "AFControl")
FeaturePlot(Aur.combined_1_2, features=c("LYVE1","FOLR2","C1QA","IL1B","TLR4","AREG","PDGFB","TREM2","TIMP1","CLEC5A","EREG","IL10","VEGFA"), min.cutoff=0.4, max.cutoff='q90',split.by = "AFControl")


Aur.combined_1_2_trans <- SCTransform (Aur.combined_1_2)

cd_genes <- c("CD68","CD14","CSF1R","CD163","LYVE1","FOLR2","IL1B","TLR4","AREG","EREG","CD9","C1QA","TREM2","PDGFB","S100A8","VCAN","LYZ","FCGR3A")
DotPlot(Aur.combined_1_2_trans,features = cd_genes)+RotatedAxis()+coord_flip()

Aur.combined_1_2_prop. <- prop.table(table(Idents(Aur.combined_1_2),Aur.combined_1_2$AFControl))
write.table(Aur.combined_1_2_prop., "Aur.combined_1_2_prop..txt", quote=F, col.names=F, append=T)
write.csv(Aur.combined_1_2_prop., "Aur.combined_1_2_prop..csv")




#macs+mono_final(心筋\rericyteコンタミ除去)--------
Aur.combined_1_3 <- subset(Aur.combined_1_2, idents = c(0,1,2,3,4,5,6,7))
Aur.combined_1_3 <- ScaleData(Aur.combined_1_3, verbose = FALSE)
Aur.combined_1_3<- RunPCA(Aur.combined_1_3, npcs = 30, verbose = FALSE)
Aur.combined_1_3 <- RunUMAP(Aur.combined_1_3 , reduction = "pca", dims = 1:30)
Aur.combined_1_3 <- FindNeighbors(Aur.combined_1_3 , reduction = "pca", dims = 1:30)
Aur.combined_1_3 <- FindClusters(Aur.combined_1_3, resolution = 0.4)

# RDS_Myeloid Aur.combined_1_3(最終並び替えなど前)
)
saveRDS(Aur.combined_1_3, "~/Desktop/Human_LA_LAA/RDS/Aur_cluster_withcontorl_AurscMy_1_3.rds")
Aur.combined_1_3 <- readRDS ("~/Desktop/Human_LA_LAA/RDS/Aur_cluster_withcontorl_AurscMy_1_3.rds")



# Visualization 別々------
DimPlot(Aur.combined_1_3, reduction = "umap", label = FALSE, repel = TRUE)
DimPlot(Aur.combined_1_3, reduction = "umap", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined_1_3, reduction = "umap", split.by = "Aur1ex")
DimPlot(Aur.combined_1_3, reduction = "umap", split.by = "AFControl1ex", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined_1_3, reduction = "umap", split.by = "AFControlAur1Control3ex", label = FALSE,  repel = FALSE)
DimPlot(Aur.combined_1_3, reduction = "umap", split.by = "AFDurationAur1Control3ex")
DimPlot(Aur.combined_1_3, reduction = "umap", split.by = "sample1ex")
DimPlot(Aur.combined_1_3, reduction = "umap", split.by = "stim")




VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("CD68","S100A8","CD14","CSF1R","CD163","LYVE1"))
VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("CD9","C1QA","C1QB","TREM2","IL1B","TNF"))
VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("CXCL3","CXCL8","TLR2","TLR4","VCAN","FCN1"))
VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("LYZ","FCGR3A","CSF3R","FCGR3B","CLEC10A","FCER1A"))
VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("CD1D","HLA-DRA","HLA-DQA1","HLA-DPA1","CD1C","CD2"))
VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("CLEC9A","XCR1","CADM1","CLEC4C","IL3RA","SPP1"))
VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("TIMP1","IL10","ISG15","PDGFB"))
VlnPlot(Aur.combined_1_3, pt.size = 0, features = c("CCL3","CCL4","ISG15","PDGFB"))


#New cluster resident--------

levels(Aur.combined_1_3) # [1] "0" "1" "2" "3" "4" "5" "6" "7" "8"
levels(Aur.combined_1_3) <- c("0", "1", "2", "3", "4", "5", "6", "7","8" )

new.cluster.ids <- c("2", "3", "0", "4", "5", "1", "6", "0","7")
names(new.cluster.ids) <- levels(Aur.combined_1_3)
Aur.combined_1_3<- RenameIdents(Aur.combined_1_3, new.cluster.ids)
DimPlot(Aur.combined_1_3, reduction = "umap", label = TRUE,  repel = TRUE)
DimPlot(Aur.combined_1_3, reduction = "umap", split.by = "stim")

new.cluster.ids <- c("")

VlnPlot(Aur.combined_1_3, features = c("AREG"), pt.size = 0.1 , split.by = "AFControl1ex" )


# RDS_Myeloid final
saveRDS(Aur.combined_1_3, "~/Desktop/Human_LA_LAA/RDS/Aur_cluster_withcontorl_AurscMy_final.rds")
Aur.combined_1_3 <- readRDS ("~/Desktop/Human_LA_LAA/RDS/Aur_cluster_withcontorl_AurscMy_final.rds")

DimPlot(Aur.combined_1_3, reduction = "umap", split.by = "AFControl1ex",cols = c("0" = "#F8766D","1" = "#E58700","2" = "#C99800","3" = "#A3A500","4" = "#6BB100","5" = "#00BA38","6" = "#00BF7D","7" = "#00C0AF"))
DimPlot(Aur.combined_1_3, reduction = "umap", split.by = "Aur1ex",cols = c("0" = "#F8766D","1" = "#E58700","2" = "#C99800","3" = "#A3A500","4" = "#6BB100","5" = "#00BA38","6" = "#00BF7D","7" = "#00C0AF"))
cols = c("0" = "#F8766D","1" = "#E58700","2" = "#C99800","3" = "#A3A500","4" = "#6BB100","5" = "#00BA38","6" = "#00BF7D","7" = "#00BF7D")




### AFControlAur1Control3ex ##
# stim列を抽出する
cell.stim <- Aur.combined_1_3@meta.data[["stim"]]
# 新しいグループ（疾患ありなしの2群）とstim列との対応関係を定義する
newgroup <- c("3.Other","3.Other","2.AF","2.AF","2.AF","2.AF","2.AF","2.AF","2.AF","2.AF","2.AF","2.AF","1.Control","1.Control","3.Other","1.Control","1.Control","3.Other")
names(newgroup) <- c("Aur1Endo","Aur1Epi","Aur2","Aur3Endo","Aur3Epi","Aur4Endo","Aur4Epi","Aur12","Aur13","Aur10","Aur16","Aur19","Control1","Control2","Control3","Control4","Control5","Control6")
# 抽出したstim列を新しいグループに変換する
cell.newgroup <- newgroup[cell.stim]
# 新しいグループを新しいMetaDataとして追加する
Aur.combined_1_3 <- AddMetaData(Aur.combined_1_3, cell.newgroup, col.name = "AFControlAur1Control3ex")




### AFDurationAur1Control3ex ####
# stim列を抽出する
cell.stim <- Aur.combined_1_3@meta.data[["stim"]]
# 新しいグループ（疾患ありなしの2群）とstim列との対応関係を定義する
newgroup <- c("4Other","4Other","2PerAF","3LSPerAF","3LSPerAF","3LSPerAF","3LSPerAF","3LSPerAF","2PerAF","3LSPerAF","2PerAF","3LSPerAF","1Control","1Control","4Other","1Control","1Control","4Other")
names(newgroup) <- c("Aur1Endo","Aur1Epi","Aur2","Aur3Endo","Aur3Epi","Aur4Endo","Aur4Epi","Aur12","Aur13","Aur10","Aur16","Aur19","Control1","Control2","Control3","Control4","Control5","Control6")
# 抽出したstim列を新しいグループに変換する
cell.newgroup <- newgroup[cell.stim]
# 新しいグループを新しいMetaDataとして追加する
Aur.combined_1_3 <- AddMetaData(Aur.combined_1_3, cell.newgroup, col.name = "AFDurationAur1Control3ex")

VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("AREG"), split.by = "AFDurationAur1Control3ex" )
VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("EREG"), split.by = "AFDurationAur1Control3ex" )
VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("HBEGF"), split.by = "AFDurationAur1Control3ex" )
VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("IGF1"), split.by = "AFDurationAur1Control3ex" )
VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("PDGFB"), split.by = "AFDurationAur1Control3ex" )

VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("AREG","EREG","HBEGF","IGF1","PDGFB"), split.by = "AFDurationAur1Control3ex" )

Aur.combined_1_3_trans <- SCTransform (Aur.combined_1_3)

cd_genes <- c("CD68","CD14","CSF1R","CD163","CCL3","CCL4","TNF","IL1B","AREG","EREG","HBEGF","TREM2","CD9","PDGFB","PDGFC","IGF1","LYVE1","FOLR2","TIMP1","THBS1","CLEC5A","IL10","VEGFA","VEGFB","S100A8","S100A12","LYZ","FCGR3A","ISG15","IFIT3")
DotPlot(Aur.combined_1_3_trans,features = cd_genes)+RotatedAxis()+coord_flip()

Aur.combined_1_3_prop. <- prop.table(table(Idents(Aur.combined_1_3),Aur.combined_1_3$Aur1ex))
write.table(Aur.combined_1_3_prop., "Aur.combined_1_3_prop.Aur1ex.txt", quote=F, col.names=F, append=T)
write.csv(Aur.combined_1_3_prop., "Aur.combined_1_3_prop.Aur1ex.csv")
Aur.combined_1_3_prop. <- prop.table(table(Idents(Aur.combined_1_3),Aur.combined_1_3$AFControl1ex))
write.table(Aur.combined_1_3_prop., "Aur.combined_1_3_prop.AFControl1ex.txt", quote=F, col.names=F, append=T)
write.csv(Aur.combined_1_3_prop., "Aur.combined_1_3_prop.AFControl1ex.csv")

Aur.combined_1_3_prop. <- prop.table(table(Idents(Aur.combined_1_3),Aur.combined_1_3$sample1ex))
write.table(Aur.combined_1_3_prop., "Aur.combined_1_3_prop.sample1ex.txt", quote=F, col.names=F, append=T)
write.csv(Aur.combined_1_3_prop., "Aur.combined_1_3_prop.sample1ex.csv")

# ヒートマップの作成
Aur.combined_1_3.markers <- FindAllMarkers(Aur.combined_1_3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Aur.combined_1_3.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- Aur.combined_1_3.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Aur.combined_1_3, features = top10$gene) + NoLegend()
write.table(top10, "table.txt", quote=F, col.names=F, append=T)


FeaturePlot(Aur.combined_1_3, features=c('AREG'), min.cutoff=1, max.cutoff='q90',split.by = "AFControlAur1Control3ex")
FeaturePlot(Aur.combined_1_3, features=c('AREG'), min.cutoff=1, max.cutoff='q90')
FeaturePlot(Aur.combined_1_3, features=c('EREG'), min.cutoff=0.4, max.cutoff='q90',split.by = "AFControlAur1Control3ex")
FeaturePlot(Aur.combined_1_3, features=c('EREG'), min.cutoff=0.4, max.cutoff='q90')
FeaturePlot(Aur.combined_1_3, features=c('TIMP1'), min.cutoff=1.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_3, features=c('TIMP1'), min.cutoff=1.5, max.cutoff='q90',split.by = "AFControlAur1Control3ex")

FeaturePlot(Aur.combined_1_3, features=c('HBEGF'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(Aur.combined_1_3, features=c('HBEGF'), min.cutoff=0.5, max.cutoff='q90',split.by = "AFControlAur1Control3ex")
VlnPlot(Aur.combined_1_3, features=c('HBEGF'), split.by = "AFControlAur1Control3ex")
VlnPlot(Aur.combined_1_3,pt.size = 0, features=c('IGF1'), split.by = "AFControlAur1Control3ex")

VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("PDGFB"), split.by = "AFControlAur1Control3ex" )
VlnPlot(Aur.combined_1_3,pt.size = 0.1, features = c("HBEGF") , assay = "RNA", split.by = "AFControlAur1Control3ex" , cols = c("Control" ="deepskyblue", "AF"="red"))

png("vlnplot_NK.png", width=800, height=400) 

#EGF
VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("EGF","TGFA","AREG","BTC","HBEGF","EREG") )
VlnPlot(Aur.combined_1_3,pt.size = 0,features = c("EGFR","ERBB2","ERBB4") )


#Tcells--------
Aur.combined_2 <- subset(Aur.combined_CD45_1, idents = c(6,7,8,9))
Aur.combined_2 <- ScaleData(Aur.combined_2, verbose = FALSE)
Aur.combined_2<- RunPCA(Aur.combined_2, npcs = 30, verbose = FALSE)
Aur.combined_2 <- RunUMAP(Aur.combined_2 , reduction = "pca", dims = 1:30)
Aur.combined_2 <- FindNeighbors(Aur.combined_2, reduction = "pca", dims = 1:30)
Aur.combined_2 <- FindClusters(Aur.combined_2, resolution = 0.1)
p1 <- DimPlot(Aur.combined_2, reduction = "umap", split.by = "stim")
p1
p2 <- DimPlot(Aur.combined_2, reduction = "umap", label = TRUE, repel = TRUE)
p2
DimPlot(Aur.combined_2, reduction = "umap", label = TRUE, repel = TRUE)

FeaturePlot(Aur.combined_2, features=c('TRGV9'), min.cutoff=0.5, max.cutoff='q90', split.by = "AFControl")

# RDS
saveRDS(Aur.combined_2, "~/Desktop/Human_LA_LAA/RDS/Aur_2_cluster_withcontorl.rds")
Aur.combined_2<- readRDS ("~/Desktop/Human_LA_LAA/RDS/Aur_2_cluster_withcontorl.rds")

Aur.combined_2_trans <- SCTransform (Aur.combined_2)

cd_genes <- c("CD3E","CD4","CD8A","GZMB","NKG7","GNLY","TYROBP","FCGR3A","NCAM1","MKI67","CD79A")
DotPlot(Aur.combined_2_trans,features = cd_genes)+RotatedAxis()+coord_flip()

cd_genes <- c("CD3E","CD4","CD8A","IL7R","GATA3","XCL1","CD7","XCL2","CLIC3","KLRC2","TRDV2","TRGV9")
DotPlot(Aur.combined_2_trans,features = cd_genes)+RotatedAxis()+coord_flip()

ここから

#Tcells(最終、コンタミ除去)--------
Aur.combined_2_0 <- subset(Aur.combined_2, idents = c(0,1,2,3))
Aur.combined_2_0 <- ScaleData(Aur.combined_2_0, verbose = FALSE)
Aur.combined_2_0<- RunPCA(Aur.combined_2_0, npcs = 30, verbose = FALSE)
Aur.combined_2_0 <- RunUMAP(Aur.combined_2_0 , reduction = "pca", dims = 1:30)
Aur.combined_2_0 <- FindNeighbors(Aur.combined_2_0, reduction = "pca", dims = 1:30)
Aur.combined_2_0 <- FindClusters(Aur.combined_2_0, resolution = 0.1)

FeaturePlot(Aur.combined_2_0, features=c('TRGV9'), min.cutoff=0.05, max.cutoff='q90', split.by = "AFControl")

Aur.combined_2_0_trans <- SCTransform (Aur.combined_2_0)

cd_genes <- c("CD3E","CD4","CD8A","GZMB","NKG7","GNLY","TYROBP","FCGR3A","NCAM1","MKI67","CD79A")
DotPlot(Aur.combined_2_0_trans,features = cd_genes)+RotatedAxis()+coord_flip()

cd_genes <- c("CD3E","CD4","CD8A","IL7R","GATA3","XCL1","CD7","XCL2","CLIC3","KLRC2","TRDV2","TRGV9")
DotPlot(Aur.combined_2_0_trans,features = cd_genes)+RotatedAxis()+coord_flip()

### JRAS追記ここから ####
# stim列を抽出する
cell.stim <- Aur.combined_2_0@meta.data[["stim"]]
# 新しいグループ（疾患ありなしの2群）とstim列との対応関係を定義する
newgroup <- c("AF","AF","AF","AF","AF","AF","AF","Control","Control","Control","Control","Control","Control")
names(newgroup) <- c("Aur1Endo","Aur1Epi","Aur2","Aur3Endo","Aur3Epi","Aur4Endo","Aur4Epi","Control1","Control2","Control3","Control4","Control5","Control6")
# 抽出したstim列を新しいグループに変換する
cell.newgroup <- newgroup[cell.stim]
# 新しいグループを新しいMetaDataとして追加する
Aur.combined_2_0 <- AddMetaData(Aur.combined_2_0, cell.newgroup, col.name = "AFControl")

# RDS
saveRDS(Aur.combined_2_0, "~/Desktop/Human_LA_LAA/RDS/Aur_2_0_cluster_withcontorl.rds")
Aur.combined_2_0 <- readRDS ("~/Desktop/Human_LA_LAA/RDS/Aur_2_0_cluster_withcontorl.rds")

### AFDuration ####
# stim列を抽出する
cell.stim <- Aur.combined_2_0@meta.data[["stim"]]
# 新しいグループ（疾患ありなしの2群）とstim列との対応関係を定義する
newgroup <- c("PerAF","PerAF","PerAF","LSPerAF","LSPerAF","LSPerAF","LSPerAF","Control","Control","Control","Control","Control","Control")
names(newgroup) <- c("Aur1Endo","Aur1Epi","Aur2","Aur3Endo","Aur3Epi","Aur4Endo","Aur4Epi","Control1","Control2","Control3","Control4","Control5","Control6")
# 抽出したstim列を新しいグループに変換する
cell.newgroup <- newgroup[cell.stim]
# 新しいグループを新しいMetaDataとして追加する
Aur.combined_2_0 <- AddMetaData(Aur.combined_2_0, cell.newgroup, col.name = "AFDuration")


# 新しいグループ分けで可視化
DimPlot(Aur.combined_2_0, reduction = "umap", label = FALSE, repel = TRUE)
DimPlot(Aur.combined_2_0, reduction = "umap", label = TRUE,  repel = TRUE)

DimPlot(Aur.combined_2_0, reduction = "umap", group.by = "AFControl")
DimPlot(Aur.combined_2_0, reduction = "umap", split.by = "AFControl")
DimPlot(Aur.combined_2_0, reduction = "umap", split.by = "AFControl", label = TRUE,  repel = TRUE)

DimPlot(Aur.combined_2_0, reduction = "umap", split.by = "AFDuration")

DimPlot(Aur.combined_2_0, reduction = "umap", split.by = "stim")



VlnPlot(Aur.combined_2_0, pt.size = 0,features = c("CD3E","CD4","CD8A","IL7R","GZMB","TYROBP") )
VlnPlot(Aur.combined_2_0, features = c("CD3E","CD4","CD8A","IL7R","GZMB","TYROBP") )
VlnPlot(Aur.combined_2_0, features = c("GZMA","GZMK","PRF1","CD28","LEF1","SELL") )
VlnPlot(Aur.combined_2_0, features = c("RORA","GATA3","CD40LG","IL10","IL4","IFNG") )
VlnPlot(Aur.combined_2_0, features = c("IL10","IL4","IL6","IFNG","IL17A") )
VlnPlot(Aur.combined_2_0, features = c("CCR7","GATA3","RORC","FOXP3") )
VlnPlot(Aur.combined_2_0, features = c("GZMB","TBX21","NKG7","GNLY","CX3CR1","CD69") )
VlnPlot(Aur.combined_2_0, features = c("CD68","S100A9","CSF1R","CD14","IL1B","EGR1") )

VlnPlot(Aur.combined_2_0,pt.size = 0, features = c("CD68","CSF1R","ITGAM","S100A8","S100A9","CXCR2") )
VlnPlot(Aur.combined_2_0, features = c("FCGR3A","CD14","LYZ","VCAN","FCN1","CTSS") )

VlnPlot(Aur.combined_2_0, features = c("RORA","GATA3","PD1","CD40LG","IL10","IL4","LILRB1R") ,split.by = "stim")

VlnPlot(Aur.combined_2_0,pt.size = 0, features = c("CD3E","CD4","CD8A","GZMB","NKG7","GNLY","TYROBP","FCGB3A","NCAM1"))


VlnPlot(Aur.combined_2_0,pt.size = 0, features = c("CD3E","CD4","CD8A","GZMB","NKG7","GNLY"))
VlnPlot(Aur.combined_2_0,pt.size = 0, features = c("GZMB","NKG7","GNLY","TYROBP","FCGR3A","NCAM1"))


# ヒートマップの作成
Aur.combined_2_0.markers <- FindAllMarkers(Aur.combined_2_0, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Aur.combined_2_0.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- Aur.combined_2_0.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Aur.combined_2_0, features = top10$gene) + NoLegend()
write.table(top10, "table.txt", quote=F, col.names=F, append=T)





Aur.combined_2_0_trans <- SCTransform (Aur.combined_2_0)

cd_genes <- c("CD3E","CD4","CD8A","GZMB","NKG7","GNLY","TYROBP","FCGR3A","NCAM1")
DotPlot(Aur.combined_2_0_trans,features = cd_genes)+RotatedAxis()+coord_flip()


Aur.combined_2_0_prop. <- prop.table(table(Idents(Aur.combined_2_0),Aur.combined_2_0$AFControl))
write.table(Aur.combined_2_0_prop., "Aur.combined_2_0_prop..txt", quote=F, col.names=F, append=T)
write.csv(Aur.combined_2_0_prop., "Aur.combined_2_0_prop..csv")

Aur.combined_2_0_prop. <- prop.table(table(Idents(Aur.combined_2_0),Aur.combined_2_0$AFDuration))
write.table(Aur.combined_2_0_prop., "Aur.combined_2_0_prop..txt", quote=F, col.names=F, append=T)
write.csv(Aur.combined_2_0_prop., "Aur.combined_2_0_prop..csv")



#CD8----
Aur.combined_2_1 <- subset(Aur.combined_2_0, idents = c(0))
Aur.combined_2_1 <- ScaleData(Aur.combined_2_1, verbose = FALSE)
Aur.combined_2_1<- RunPCA(Aur.combined_2_1, npcs = 30, verbose = FALSE)
Aur.combined_2_1 <- RunUMAP(Aur.combined_2_1, reduction = "pca", dims = 1:30)
Aur.combined_2_1 <- FindNeighbors(Aur.combined_2_1, reduction = "pca", dims = 1:30)
Aur.combined_2_1 <- FindClusters(Aur.combined_2_1, resolution = 0.2)
p1 <- DimPlot(Aur.combined_2_1, reduction = "umap", split.by = "stim")
p1
p2 <- DimPlot(Aur.combined_2_1, reduction = "umap", label = TRUE, repel = TRUE)
p2

### JRAS追記ここから ####
# stim列を抽出する
cell.stim <- Aur.combined_2_1@meta.data[["stim"]]
# 新しいグループ（疾患ありなしの2群）とstim列との対応関係を定義する
newgroup <- c("Endo","Endo","Endo","Epi","Epi","Epi")
names(newgroup) <- c("Aur1Endo","Aur3Endo","Aur4Endo","Aur1Epi","Aur3Epi","Aur4Epi")
# 抽出したstim列を新しいグループに変換する
cell.newgroup <- newgroup[cell.stim]
# 新しいグループを新しいMetaDataとして追加する
Aur.combined_2_1 <- AddMetaData(Aur.combined_2_1, cell.newgroup, col.name = "EndoEpi")

# RDS
saveRDS(Aur.combined_2_1, "~/Desktop/Human_LA_LAA/RDS/Aur_2_1_cluster_withcontorl.rds")
Aur.combined_2_1 <- readRDS ("~/Desktop/Human_LA_LAA/RDS/Aur_2_1_cluster_withcontorl.rds")

### AFDuration ####
# stim列を抽出する
cell.stim <- Aur.combined_2_1@meta.data[["stim"]]
# 新しいグループ（疾患ありなしの2群）とstim列との対応関係を定義する
newgroup <- c("PerAF","PerAF","PerAF","LSPerAF","LSPerAF","LSPerAF","LSPerAF","Control","Control","Control","Control","Control","Control")
names(newgroup) <- c("Aur1Endo","Aur1Epi","Aur2","Aur3Endo","Aur3Epi","Aur4Endo","Aur4Epi","Control1","Control2","Control3","Control4","Control5","Control6")
# 抽出したstim列を新しいグループに変換する
cell.newgroup <- newgroup[cell.stim]
# 新しいグループを新しいMetaDataとして追加する
Aur.combined_2_1 <- AddMetaData(Aur.combined_2_1, cell.newgroup, col.name = "AFDuration")


# 新しいグループ分けで可視化
DimPlot(Aur.combined_2_1, reduction = "umap", label = FALSE, repel = TRUE)
DimPlot(Aur.combined_2_1, reduction = "umap", label = TRUE,  repel = TRUE)

DimPlot(Aur.combined_2_1, reduction = "umap", group.by = "AFControl")
DimPlot(Aur.combined_2_1, reduction = "umap", split.by = "AFControl")

DimPlot(Aur.combined_2_1, reduction = "umap", split.by = "AFDuration")

DimPlot(Aur.combined_2_1, reduction = "umap", split.by = "stim")


VlnPlot(Aur.combined_2_1,    features = c("GZMA","GZMK","GZMB","TBX21","CX3CR1","NKG7") )
VlnPlot(Aur.combined_2_1,   features = c("CD69","IL7R","IFNG","CD103","SELL","CD27") )
VlnPlot(Aur.combined_2_1,   features = c("PRF1","IL2","IFNG","EOMES","SELL","CD27") )
VlnPlot(Aur.combined_2_1,  features = c("TNF","CD27","IFNG","MKI67","CD74","CXCR6") )
VlnPlot(Aur.combined_2_1, features = c("XCL1","CD7","LGALS3","ITGAE","CLEC10A","GNLY") )
VlnPlot(Aur.combined_2_1,features = c("GZMB","GZMH","PRF1","LEF1") )
VlnPlot(Aur.combined_2_1, features = c("KLF10","CCR7","IFNG","GZMH","CD40LG","CD4") )

VlnPlot(Aur.combined_2_1, features = c("KLF10","CCR7","IFNG","GZMH","CD40LG","CD4") )

VlnPlot(Aur.combined_2_1,pt.size = 0, features = c("CD8A","CD4","GNLY","TYROBP","TGFB1") )

VlnPlot(Aur.combined_2_1,pt.size = 0, features = c("LEF1","SELL","ACTN1","BACH2","CX3CR1","GZMK","SIPRG","IFNG","PDCD1","LAYN","SLC4A10") )
VlnPlot(Aur.combined_2_1,pt.size = 0, features = c("CD27","IL7R","GZMA","GZMK","CD69","CD74","GZMB","TBX21","NKG7","CX3CR1","GZMH","PRF1","TNF","IFNG","CD40LG","CXCR6") )

Naive
VlnPlot(Aur.combined_2_1, features = c("LEF1","LTB","CCR7") )
Effector
VlnPlot(Aur.combined_2_1, features = c("CX3CR1","FCGR3A","FGFBP2") )
MAIT
VlnPlot(Aur.combined_2_1, features = c("SLC4A10","ZBTB16","RORC") )
Tex
VlnPlot(Aur.combined_2_1, features = c("CTLA4","PDCD1","HAVCR2","LAYN") )

VlnPlot(Aur.combined_2_1, features = c("CD40LG","CD28") , split.by = "AFDuration" , cols = c("Control" ="yellow", "LSPerAF"="blue","PerAF" ="red"),pt.size = 0)


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




Aur.combined_2_1_prop. <- prop.table(table(Idents(Aur.combined_2_1),Aur.combined_2_1$AFControl))
write.table(Aur.combined_2_1_prop., "Aur.combined_2_1_prop..txt", quote=F, col.names=F, append=T)
write.csv(Aur.combined_2_1_prop., "Aur.combined_2_1_prop..csv")

Aur.combined_2_1_prop. <- prop.table(table(Idents(Aur.combined_2_1),Aur.combined_2_1$AFDuration))
write.table(Aur.combined_2_1_prop., "Aur.combined_2_1_prop..txt", quote=F, col.names=F, append=T)
write.csv(Aur.combined_2_1_prop., "Aur.combined_2_1_prop..csv")

#CD4----
Aur.combined_2_2 <- subset(Aur.combined_2_0, idents = c(1))
Aur.combined_2_2 <- ScaleData(Aur.combined_2_2, verbose = FALSE)
Aur.combined_2_2<- RunPCA(Aur.combined_2_2, npcs = 30, verbose = FALSE)
Aur.combined_2_2 <- RunUMAP(Aur.combined_2_2, reduction = "pca", dims = 1:30)
Aur.combined_2_2 <- FindNeighbors(Aur.combined_2_2, reduction = "pca", dims = 1:30)
Aur.combined_2_2 <- FindClusters(Aur.combined_2_2, resolution = 0.2)

p1 <- DimPlot(Aur.combined_2_2, reduction = "umap", split.by = "stim")
p1
p2 <- DimPlot(Aur.combined_2_2, reduction = "umap", label = TRUE, repel = TRUE)
p2

#CD4(最終、コンタミ除去)--------
Aur.combined_2_3 <- subset(Aur.combined_2_2, idents = c(0,1,2,3,4,5))
Aur.combined_2_3 <- ScaleData(Aur.combined_2_3, verbose = FALSE)
Aur.combined_2_3<- RunPCA(Aur.combined_2_3, npcs = 30, verbose = FALSE)
Aur.combined_2_3 <- RunUMAP(Aur.combined_2_3 , reduction = "pca", dims = 1:30)
Aur.combined_2_3 <- FindNeighbors(Aur.combined_2_3, reduction = "pca", dims = 1:30)
Aur.combined_2_3 <- FindClusters(Aur.combined_2_3, resolution = 0.3)


### JRAS追記ここから ####
# stim列を抽出する
cell.stim <- Aur.combined_2_3@meta.data[["stim"]]
# 新しいグループ（疾患ありなしの2群）とstim列との対応関係を定義する
newgroup <- c("Endo","Endo","Endo","Epi","Epi","Epi")
names(newgroup) <- c("Aur1Endo","Aur3Endo","Aur4Endo","Aur1Epi","Aur3Epi","Aur4Epi")
# 抽出したstim列を新しいグループに変換する
cell.newgroup <- newgroup[cell.stim]
# 新しいグループを新しいMetaDataとして追加する
Aur.combined_2_3 <- AddMetaData(Aur.combined_2_3, cell.newgroup, col.name = "EndoEpi")


# 並べ替え
levels(Aur.combined_2_3) # [1] "0" "1" "2" "3" "4" "5" "6" 
levels(Aur.combined_2_3) <- c("0", "1", "3", "2", "4","5","6")

new.cluster.ids <- c("0","1","2","3","4","5","6")
names(new.cluster.ids) <- levels(Aur.combined_2_3)
Aur.combined_2_3<- RenameIdents(Aur.combined_2_3, new.cluster.ids)

# RDS
saveRDS(Aur.combined_2_3, "~/Desktop/Human_LA_LAA/RDS/Aur_2_3_cluster_withcontorl.rds")
Aur.combined_2_3 <- readRDS ("~/Desktop/Human_LA_LAA/RDS/Aur_2_3_cluster_withcontorl.rds")

### AFDuration ####
# stim列を抽出する
cell.stim <- Aur.combined_2_3@meta.data[["stim"]]
# 新しいグループ（疾患ありなしの2群）とstim列との対応関係を定義する
newgroup <- c("PerAF","PerAF","PerAF","LSPerAF","LSPerAF","LSPerAF","LSPerAF","Control","Control","Control","Control","Control","Control")
names(newgroup) <- c("Aur1Endo","Aur1Epi","Aur2","Aur3Endo","Aur3Epi","Aur4Endo","Aur4Epi","Control1","Control2","Control3","Control4","Control5","Control6")
# 抽出したstim列を新しいグループに変換する
cell.newgroup <- newgroup[cell.stim]
# 新しいグループを新しいMetaDataとして追加する
Aur.combined_2_3 <- AddMetaData(Aur.combined_2_3, cell.newgroup, col.name = "AFDuration")



# 新しいグループ分けで可視化
DimPlot(Aur.combined_2_3, reduction = "umap")
DimPlot(Aur.combined_2_3, reduction = "umap", label = TRUE,  repel = TRUE)

DimPlot(Aur.combined_2_3, reduction = "umap", group.by = "AFControl")
DimPlot(Aur.combined_2_3, reduction = "umap",split.by = "AFControl")

DimPlot(Aur.combined_2_3, reduction = "umap",split.by = "AFDuration")

DimPlot(Aur.combined_2_3, reduction = "umap", split.by = "stim")



VlnPlot(Aur.combined_2_3, features = c("GZMB","GZMK","PRF1","LEF1","IL7R","SELL") )
VlnPlot(Aur.combined_2_3, features = c("TBX21","GATA3","CCR7","RORC","FOXP3","CTLA4") )
VlnPlot(Aur.combined_2_3, features = c("IL17A","IL23","IL2","IL10","IL4","IFNG") )
VlnPlot(Aur.combined_2_3, features = c("GZMB","CD28","PRF1","CD69","IL7R","SELL") )
VlnPlot(Aur.combined_2_3, features = c("CXCR3","CCR4","CCR6","TER1"), split.by ="stim",cols = c("sa" ="deepskyblue", "Aur1Epi"="red") )

VlnPlot(Aur.combined_2_3, features = c("CCL4","CCL5","TNFAIP3","RORC","CD28","GZMB") )
VlnPlot(Aur.combined_2_3, features = c("LTB","NKG7","IFNG","CD69","CD28","GZMB") )
VlnPlot(Aur.combined_2_3, features = c("CD44","NKG7","IFNG","CD69","CD28","KLF10") )

VlnPlot(Aur.combined_2_3, features = c("CD3E","CD4","CD8A","CD68","CD14","IL1B","S100A9") )
VlnPlot(Aur.combined_2_3, features = c("CD3E","CD4","CD8A","GNLY","GZMB","TYROBP") )

FeaturePlot(Aur.combined_2_3, features=c('CD69'), min.cutoff=2.0, max.cutoff='q90',split.by = "stim" )
FeaturePlot(Aur.combined_2_3, features=c('CD69'), min.cutoff=2.0, max.cutoff='q90')

pt.size = 0, 

VlnPlot(Aur.combined_2_3, features = c("CD28","CD40LG") ,  split.by = "AFDuration" , cols = c("Control" ="yellow", "LSPerAF"="blue","PerAF" ="red" ),pt.size = 0.1)

VlnPlot(Aur.combined_2_3, features = c("CD28","CD40LG") ,  split.by = "AFDuration" , cols = c("Control" ="yellow", "LSPerAF"="blue","PerAF" ="red" ),pt.size = 0)



Manuscript
VlnPlot(Aur.combined_2_3,pt.size = 0, features = c("CD3E","GZMA","GZMK","CCL4","CCL5","TNFAIP3"))
VlnPlot(Aur.combined_2_3,pt.size = 0, features = c("IL7R","SELL","CCR7","LEF1","CD28","PRF1","GATA3","CD69"))
VlnPlot(Aur.combined_2_3,pt.size = 0, features = c("CCR7","LEF1","CD28","PRF1","GATA3","CD69"))



Aur.combined_2_3_trans <- SCTransform (Aur.combined_2_3)

cd_genes <- c("CD3E","CD4","CD8A","GZMB","NKG7","GNLY","FOXP3","CTLA4","SELL","CCR7","LEF1","CD28","PRF1","GZMA","GZMK","CD69","IL7R","CCL4","CCL5","TNFAIP3","TBX21","GATA3","RORC")
DotPlot(Aur.combined_2_3_trans,features = cd_genes)+RotatedAxis()+coord_flip()




New
Aur.combined_2_3_trans <- SCTransform (Aur.combined_2_3)

cd_genes <- c("CD3E","GZMA","GZMK","CCL4","CCL5","TNFAIP3","IL7R","SELL","CCR7","LEF1","CD28","PRF1","GATA3","CD69")
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




# VInplot----
VlnPlot(Aur.combined_1, features = c("CXCL10","S100A8","HMGB2","HLA-DPA1","HLA-DQA1","HLA-DPB1") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "Aur1Epi"="red"))
VlnPlot(Aur.combined_1, features = c("APOE") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "Aur1Epi"="red"))
VlnPlot(Aur.combined_1, features = c("CXCL10") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "Aur1Epi"="red"))

Aur.combined_1_DC <- subset(Aur.combined_1, idents = c("5_sa", "5_Aur1Epi"))
VlnPlot(Aur.combined_1_DC, features = c("CD40","CD80","CD86","HLA-DPA1","HLA-DQA1","HLA-DPB1") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "Aur1Epi"="red"))
VlnPlot(Aur.combined_1_DC, features = c("CLEC10A","FCER1A","CD1C","IL3RA","RELB") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "Aur1Epi"="red"))
VlnPlot(Aur.combined_1_DC, features = c("NECTIN2","PVR","TNFSF4","IL3RA","RELB") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "Aur1Epi"="red"))


Aur.combined_2_2_cmCD4 <- subset(Aur.combined_2_2, idents = c("3_sa", "3_Aur1Epi"))
VlnPlot(Aur.combined_2_2_cmCD4, features = c("CD28","CD40LG","CTLA4","TIGHT","FOXP3","TNFRSF4") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "Aur1Epi"="red"))
VlnPlot(Aur.combined_2_2_cmCD4, features = c("TIGIT","CD226" ,"IL17A","IL4","IL10") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "Aur1Epi"="red"))

Aur.combined_2_2_eCD4 <- subset(Aur.combined_2_2, idents = c("2_sa", "2_Aur1Epi"))
VlnPlot(Aur.combined_2_eCD4, features = c("CD28","CD40LG","PDCD1","LAG3","CTLA4","TNFRSF4") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "Aur1Epi"="red"))
VlnPlot(Aur.combined_2_eCD4, features = c("TIGIT","CD226","HAVCR2","PDCD1") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "Aur1Epi"="red"))

Aur.combined_2_2_nCD4 <- subset(Aur.combined_2_2, idents = c("0_sa","3_sa","1_sa", "0_Aur1Epi", "1_Aur1Epi", "3_Aur1Epi"))
VlnPlot(Aur.combined_2_2_nCD4, features = c("CD28","CD40LG","CTLA4","LAG3","CTLA4","TNFRSF4") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "Aur1Epi"="red"))
VlnPlot(Aur.combined_2_2_nCD4, features = c("TIGIT","CD226","HAVCR2","PDCD1") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "Aur1Epi"="red"))

Aur.combined_2_2_e1CD4 <- subset(Aur.combined_2_2, idents = c("0_sa", "0_Aur1Epi"))
VlnPlot(Aur.combined_2_2_e1CD4, features = c("CD28","CD40LG","PDCD1","LAG3","CTLA4","TNFRSF4") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "Aur1Epi"="red"))
VlnPlot(Aur.combined_2_2_e1CD4, features = c("TIGIT","CD226","HAVCR2","PDCD1") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "Aur1Epi"="red"))

Aur.combined_2_2_regCD4 <- subset(Aur.combined_2_2, idents = c("4_sa", "4_Aur1Epi"))
VlnPlot(Aur.combined_2_2_regCD4, features = c("CD28","CD40LG","PDCD1","LAG3","CTLA4","TNFRSF4") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "Aur1Epi"="red"))
VlnPlot(Aur.combined_2_2_regCD4, features = c("TIGIT","CD226","HAVCR2","PDCD1") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "Aur1Epi"="red"))


Aur.combined_2.markers <- FindAllMarkers(Aur.combined_2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Aur.combined_2.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- Aur.combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Aur.combined_2, features = top10$gene) + NoLegend()
write.table(top10, "table.txt", quote=F, col.names=F, append=T)


Aur.combined_3.markers <- FindAllMarkers(Aur.combined_3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Aur.combined_3.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- Aur.combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(Aur.combined_3, features = top10$gene) + NoLegend()
write.table(top10, "table.txt", quote=F, col.names=F, append=T)

# volcano all myeloid------

Aur.combined$celltype <- Idents(Aur.combined)
Aur.combined$celltype.stim <- paste(Idents(Aur.combined), Aur.combined$stim, sep="_")
Idents(Aur.combined) <- "celltype.stim" # Idents()関数で変更
levels(Aur.combined) # 変更できているか確認_1


Aur.combined_Myeloid_sa<- subset(Aur.combined, idents = c("5_sa","6_sa","7_sa"))
Aur.combined_Myeloid_Aur1Epi<- subset(Aur.combined, idents = c("5_Aur1Epi","6_Aur1Epi","7_Aur1Epi"))

Aur.combined_Tcells_sa<- subset(Aur.combined, idents = c("0_sa","1_sa","2_sa","3_sa","4_sa"))
Aur.combined_Tcells_Aur1Epi<- subset(Aur.combined, idents = c("0_Aur1Epi","1_Aur1Epi","2_Aur1Epi","3_Aur1Epi","4_Aur1Epi"))

Aur.combined_Bcells_sa<- subset(Aur.combined, idents = c("8_sa","9_sa"))
Aur.combined_Bcells_Aur1Epi<- subset(Aur.combined, idents = c("8_Aur1Epi","9_Aur1Epi"))



Aur.table_myeloid <- FindMarkers(Aur.combined, ident.1 = c("5_sa","6_sa","7_sa"), ident.2 =c("5_Aur1Epi","6_Aur1Epi","7_Aur1Epi"), verbose = FALSE, logfc.threshold = 0)

Aur.table_myeloid$logp <- -log10(Aur.table_myeloid$p_val)

Aur.table_myeloid_filtered_left = subset(Aur.table_myeloid, logp>=5 & avg_log2FC <= -1.0)
Aur.table_myeloid_filtered_right = subset(Aur.table_myeloid, logp>=5 & avg_log2FC >= 1.0)

genes.to.label.left <- rownames(Aur.table_myeloid_filtered_left)
genes.to.label.right <- rownames(Aur.table_myeloid_filtered_right)

p1 <- ggplot(Aur.table_myeloid, aes(avg_log2FC, logp, label)) + geom_point() 
p1 <- LabelPoints(plot = p1, points = genes.to.label.right,color="blue", repel = TRUE, xnudge=0)
p1 <- LabelPoints(plot = p1, points = genes.to.label.left,color="red", repel = TRUE, xnudge=0)
p1







# VOLCANO Trem2 macs----



# volcano all CD4------

Aur.combined_2$celltype <- Idents(Aur.combined_2)
Aur.combined_2$celltype.stim <- paste(Idents(Aur.combined_2), Aur.combined_2$stim, sep="_")
Idents(Aur.combined_2) <- "celltype.stim" # Idents()関数で変更
levels(Aur.combined_2) # 変更できているか確認_1

Aur.table_CD4 <- FindMarkers(Aur.combined_2, ident.1 = c("2_sa"), ident.2 =c("2_Aur1Epi"), verbose = FALSE, logfc.threshold = 0)

Aur.table_CD4$logp <- -log10(Aur.table_CD4$p_val)

Aur.table_CD4_filtered_left = subset(Aur.table_CD4, logp>=2 & avg_log2FC <= -0.5)
Aur.table_CD4_filtered_right = subset(Aur.table_CD4, logp>=2 & avg_log2FC >= 0.5)

genes.to.label.left <- rownames(Aur.table_CD4_filtered_left)
genes.to.label.right <- rownames(Aur.table_CD4_filtered_right)

p1 <- ggplot(Aur.table_CD4, aes(avg_log2FC, logp, label)) + geom_point() 
p1 <- LabelPoints(plot = p1, points = genes.to.label.right,color="blue", repel = TRUE, xnudge=0)
p1 <- LabelPoints(plot = p1, points = genes.to.label.left,color="red", repel = TRUE, xnudge=0)
p1

write.table(Aur.table_CD4, "effecCD4.txt", quote=F, col.names=F, append=T)

# central memory CD4

Aur.combined_2_2$celltype <- Idents(Aur.combined_2_2)
Aur.combined_2_2$celltype.stim <- paste(Idents(Aur.combined_2_2), Aur.combined_2_2$stim, sep="_")
Idents(Aur.combined_2_2) <- "celltype.stim" # Idents()関数で変更
levels(Aur.combined_2_2) # 変更できているか確認

Aur.table_cCD4 <- FindMarkers(Aur.combined_2_2, ident.1 = c("3_sa"), ident.2 =c("3_Aur1Epi"), verbose = FALSE, logfc.threshold = 0)

Aur.table_cCD4$logp <- -log10(Aur.table_cCD4$p_val)

Aur.table_cCD4_filtered_left = subset(Aur.table_cCD4, logp>=2 & avg_log2FC <= -0.4)
Aur.table_cCD4_filtered_right = subset(Aur.table_cCD4, logp>=2 & avg_log2FC >= 0.4)

genes.to.label.left <- rownames(Aur.table_cCD4_filtered_left)
genes.to.label.right <- rownames(Aur.table_cCD4_filtered_right)

p1 <- ggplot(Aur.table_cCD4, aes(avg_log2FC, logp, label)) + geom_point() 
p1 <- LabelPoints(plot = p1, points = genes.to.label.right,color="blue", repel = TRUE, xnudge=-0.1)
p1 <- LabelPoints(plot = p1, points = genes.to.label.left,color="red", repel = TRUE, xnudge=-0.1)
p1

write.table(Aur.table1, "table.txt", quote=F, col.names=F, append=T)

# volcano myeloid subset
Aur.combined_1$celltype <- Idents(Aur.combined_1)
Aur.combined_1$celltype.stim <- paste(Idents(Aur.combined_1), Aur.combined_1$stim, sep="_")
Idents(Aur.combined_1) <- "celltype.stim" # Idents()関数で変更
levels(Aur.combined_1)

Aur.table1 <- FindMarkers(Aur.combined_1, ident.1 = "5_Aur1Epi", ident.2 = "5_sa", verbose = FALSE, logfc.threshold = 0)
Aur.table2 <- FindMarkers(Aur.combined_1, ident.1 = "0_Aur1Epi", ident.2 = "0_sa", verbose = FALSE, logfc.threshold = 0)
Aur.table3 <- FindMarkers(Aur.combined_1, ident.1 = "1_Aur1Epi", ident.2 = "1_sa", verbose = FALSE, logfc.threshold = 0)

Aur.table1$logp <- -log10(Aur.table1$p_val)
Aur.table2$logp <- -log10(Aur.table2$p_val)
Aur.table3$logp <- -log10(Aur.table3$p_val)
Aur.table4$logp <- -log10(Aur.table4$p_val)

# cluster1----
Aur.table1_filtered_left = subset(Aur.table1, logp>=1.30103& avg_log2FC <= -0.8)
Aur.table1_filtered_right = subset(Aur.table1, logp>=1.30103 & avg_log2FC >= 0.8)

genes.to.label.left <- rownames(Aur.table1_filtered_left)
genes.to.label.right <- rownames(Aur.table1_filtered_right)

p1 <- ggplot(Aur.table1, aes(avg_log2FC, logp, label)) + geom_point() 
p1 <- LabelPoints(plot = p1, points = genes.to.label.right,color="red", repel = TRUE, xnudge=0)
p1 <- LabelPoints(plot = p1, points = genes.to.label.left,color="blue", repel = TRUE, xnudge=0)
p1

write.table(Aur.table1, "table.txt", quote=F, col.names=F, append=T)

VlnPlot(Aur.combined_1, ident.1 = "5_sa", ident.2 = "5_Aur1Epi", features = c("CXCL10") ,  split.by = "stim" , cols = c("sa" ="deepskyblue", "Aur1Epi"="red"))

# cluster2----
Aur.table2_filtered_left = subset(Aur.table2, logp>=2 & avg_log2FC <= -0.8)
Aur.table2_filtered_right = subset(Aur.table2, logp>=2 & avg_log2FC >= 0.8)

genes.to.label.left <- rownames(Aur.table2_filtered_left)
genes.to.label.right <- rownames(Aur.table2_filtered_right)

p1 <- ggplot(Aur.table2, aes(avg_log2FC, logp, label)) + geom_point() 
p1 <- LabelPoints(plot = p1, points = genes.to.label.right,color="red", repel = TRUE, xnudge=0)
p1 <- LabelPoints(plot = p1, points = genes.to.label.left,color="blue", repel = TRUE, xnudge=0)
p1

write.table(Aur.table2, "table.txt", quote=F, col.names=F, append=T)

# cluster3----

Aur.table3_filtered_left = subset(Aur.table3, logp>=2 & avg_log2FC <= -0.8)
Aur.table3_filtered_right = subset(Aur.table3, logp>=2 & avg_log2FC >= 0.8)

genes.to.label.left <- rownames(Aur.table3_filtered_left)
genes.to.label.right <- rownames(Aur.table3_filtered_right)

p1 <- ggplot(Aur.table3, aes(avg_log2FC, logp, label)) + geom_point() 
p1 <- LabelPoints(plot = p1, points = genes.to.label.right,color="red", repel = TRUE, xnudge=0)
p1 <- LabelPoints(plot = p1, points = genes.to.label.left,color="blue", repel = TRUE, xnudge=0)
p1

write.table(Aur.table2, "table.txt", quote=F, col.names=F, append=T)


plot1 <- FeatureScatter(Aur.combined_3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Aur.combined_3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

Aur.combined_3



