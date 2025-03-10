# ライブラリの読み込み----------------------------------------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

# split the dataset into a list of two seurat objects (stim and CTRL)----
AurN.list<-list()
data.Aur4N <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur4N/outs/filtered_feature_bc_matrix")
data.Aur9_12N <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur9_12N/outs/filtered_feature_bc_matrix")
data.Aur13_16N <- Read10X(data.dir = "~/Desktop/Human_LA_LAA/CellRanger/Aur13_16N/outs/filtered_feature_bc_matrix")

AurN.list$Aur4N<- CreateSeuratObject(counts = data.Aur4N, project = "A", min.cells = 3, min.features = 200)
AurN.list$Aur4N <- AddMetaData(AurN.list$Aur4N , "Aur4N", col.name = "stim")

Aur9_12N<- CreateSeuratObject(counts = data.Aur9_12N, project = "A", min.cells = 3, min.features = 200)
a <- read.table("~/Desktop/Human_LA_LAA/CellRanger/Aur9_12N/CellSNP0911_vireo/donor_ids.tsv", header = TRUE,sep = "\t")
for(i in 1:length(Aur9_12N@assays[["RNA"]]@data@Dimnames[[2]])){
  while(Aur9_12N@assays[["RNA"]]@data@Dimnames[[2]][i] != a[i,1]){
    a <- a[-i,]
  }
}
Aur9_12N <- AddMetaData(Aur9_12N , a[["donor_id"]], col.name = "donor_id")
AurN.list$Aur9N <- subset(Aur9_12N, subset = donor_id == c("donor0"))
AurN.list$Aur9N <- AddMetaData(AurN.list$Aur9N , "Aur9N", col.name = "stim")
AurN.list$Aur12N <- subset(Aur9_12N, subset = donor_id == c("donor1"))
AurN.list$Aur12N <- AddMetaData(AurN.list$Aur12N , "Aur12N", col.name = "stim")

Aur13_16N<- CreateSeuratObject(counts = data.Aur13_16N, project = "A", min.cells = 3, min.features = 200)
a <- read.table("~/Desktop/Human_LA_LAA/CellRanger/Aur13_16N/CellSNP0911_vireo/donor_ids.tsv", header = TRUE,sep = "\t")
for(i in 1:length(Aur13_16N@assays[["RNA"]]@data@Dimnames[[2]])){
  while(Aur13_16N@assays[["RNA"]]@data@Dimnames[[2]][i] != a[i,1]){
    a <- a[-i,]
  }
}
Aur13_16N <- AddMetaData(Aur13_16N , a[["donor_id"]], col.name = "donor_id")
AurN.list$Aur13N <- subset(Aur13_16N, subset = donor_id == c("donor0"))
AurN.list$Aur13N <- AddMetaData(AurN.list$Aur13N , "Aur13N", col.name = "stim")
AurN.list$Aur16N <- subset(Aur13_16N, subset = donor_id == c("donor1"))
AurN.list$Aur16N <- AddMetaData(AurN.list$Aur16N , "Aur16N", col.name = "stim")



rm(data.Aur13_16N)
rm(data.Aur4N)
rm(data.Aur9_12N)
rm(a)
rm(Aur13_16N)
rm(Aur9_12N)



# low-quality cellの確認①------------------------------------------------------
# MT-から始まるミトコンドリアRNAを"percent.mt"列としてデータに追加

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
AurN.list$Aur4N[["percent.mt"]] <- PercentageFeatureSet(AurN.list$Aur4N, pattern = "^MT-")
AurN.list$Aur9N[["percent.mt"]] <- PercentageFeatureSet(AurN.list$Aur9N, pattern = "^MT-")
AurN.list$Aur12N[["percent.mt"]] <- PercentageFeatureSet(AurN.list$Aur12N, pattern = "^MT-")
AurN.list$Aur13N[["percent.mt"]] <- PercentageFeatureSet(AurN.list$Aur13N, pattern = "^MT-")
AurN.list$Aur16N[["percent.mt"]] <- PercentageFeatureSet(AurN.list$Aur16N, pattern = "^MT-")


VlnPlot(AurN.list$Aur4N, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(AurN.list$Aur9N, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(AurN.list$Aur12N, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(AurN.list$Aur13N, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(AurN.list$Aur16N, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

# percent.mtが極端に大きい（10%以上）、またはnFeature_RNAが小さすぎる（200以下）・大きすぎる(2500以上)細胞は、
# 死細胞やdoubletの可能性が高い
plot1 <- FeatureScatter(AurN.list$Aur4N, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(AurN.list$Aur4N, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(AurN.list$Aur9N, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(AurN.list$Aur9N, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(AurN.list$Aur12N, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(AurN.list$Aur12N, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(AurN.list$Aur13N, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(AurN.list$Aur13N, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(AurN.list$Aur16N, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(AurN.list$Aur16N, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# low-quality cellのフィルタリング---------------------------------------------
# ミトコンドリアRNAが７以下、遺伝子数が200~2500の間の細胞のみ抽出
AurN.list$Aur4N <- subset(AurN.list$Aur4N, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mt < 30)
AurN.list$Aur9N <- subset(AurN.list$Aur9N, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mt < 30)
AurN.list$Aur12N <- subset(AurN.list$Aur12N, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mt < 30)
AurN.list$Aur13N <- subset(AurN.list$Aur13N, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mt < 30)
AurN.list$Aur16N <- subset(AurN.list$Aur16N, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mt < 30)


# 先程と比べlow-qualityデータがフィルターされていることを確認
plot1 <- FeatureScatter(AurN.list$Aur4N, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(AurN.list$Aur4N, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(AurN.list$Aur9N, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(AurN.list$Aur9N, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(AurN.list$Aur12N, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(AurN.list$Aur12N, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(AurN.list$Aur13N, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(AurN.list$Aur13N, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot1 <- FeatureScatter(AurN.list$Aur16N, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(AurN.list$Aur16N, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# データの正規化（normalization）-----------------------------------------------
AurN.list$Aur4N<- NormalizeData(AurN.list$Aur4N, normalization.method = "LogNormalize", scale.factor = 10000)
AurN.list$Aur9N<- NormalizeData(AurN.list$Aur9N, normalization.method = "LogNormalize", scale.factor = 10000)
AurN.list$Aur12N<- NormalizeData(AurN.list$Aur12N, normalization.method = "LogNormalize", scale.factor = 10000)
AurN.list$Aur13N<- NormalizeData(AurN.list$Aur13N, normalization.method = "LogNormalize", scale.factor = 10000)
AurN.list$Aur16N<- NormalizeData(AurN.list$Aur16N, normalization.method = "LogNormalize", scale.factor = 10000)


all.genes <- rownames(AurN.list$Aur4N)
AurN.list$Aur4N <- ScaleData(AurN.list$Aur4N, features = all.genes)

all.genes <- rownames(AurN.list$Aur9N)
AurN.list$Aur9N <- ScaleData(AurN.list$Aur9N, features = all.genes)

all.genes <- rownames(AurN.list$Aur12N)
AurN.list$Aur12N <- ScaleData(AurN.list$Aur12N, features = all.genes)

all.genes <- rownames(AurN.list$Aur13N)
AurN.list$Aur13N <- ScaleData(AurN.list$Aur13N, features = all.genes)

all.genes <- rownames(AurN.list$Aur16N)
AurN.list$Aur16N <- ScaleData(AurN.list$Aur16N, features = all.genes)


# UMAPによるクラスタリング　2つ合わせて------------------------------------------------------------------------------

# normalize and identify variable features for each dataset independently
all.genes <- rownames(AurN.list["Aur4N"])+   rownames(AurN.list["Aur9N"])+   rownames(AurN.list["Aur12N"])+   rownames(AurN.list["Aur13N"])+   rownames(AurN.list["Aur16N"])



AurN.list <- lapply(X = AurN.list, FUN = function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
  x <- ScaleData(x, features = all.genes)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = AurN.list)

# Integrationの実行
AurN.anchors <- FindIntegrationAnchors(object.list = AurN.list, anchor.features = features)

# this command creates an 'integrated' data assay
AurN.combined <- IntegrateData(anchorset =AurN.anchors)

# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(AurN.combined) <- "integrated"

AurN.combined <- ScaleData(AurN.combined, verbose = FALSE)
AurN.combined <- RunPCA(AurN.combined, npcs = 30, verbose = FALSE)
AurN.combined <- RunUMAP(AurN.combined, reduction = "pca", dims = 1:30)
AurN.combined <- FindNeighbors(AurN.combined, reduction = "pca", dims = 1:30)
AurN.combined <- FindClusters(AurN.combined, resolution = 0.35)






# 並べ替える
levels(AurN.combined) # [1] "0" "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" 
levels(AurN.combined) <- c("0", "1",  "7", "2", "5", "10","3","6","4","8","9")

new.cluster.ids <- c("0","1","2","3","4","5","6","7","8","9","10")
names(new.cluster.ids) <- levels(AurN.combined)
AurN.combined<- RenameIdents(AurN.combined, new.cluster.ids)



# ヒートマップの作成------------------------------------------------------------------------------

AurN.combined.markers <- FindAllMarkers(AurN.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
AurN.combined.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- AurN.combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(AurN.combined, features = top10$gene) 
write.table(top10, "table.txt", quote=F, col.names=F, append=T)







# stim列を抽出する
cell.stim <- AurN.combined_1@meta.data[["stim"]]
# 新しいグループ（疾患ありなしの2群）とstim列との対応関係を定義する
newgroup <- c("2PermanentAF","1EarlyAF","1EarlyAF","2PermanentAF","2PermanentAF")
names(newgroup) <- c("Aur12N","Aur13N","Aur16N","Aur4N","Aur9N")
# 抽出したstim列を新しいグループに変換する
cell.newgroup <- newgroup[cell.stim]
# 新しいグループを新しいMetaDataとして追加する
AurN.combined_1 <- AddMetaData(AurN.combined_1, cell.newgroup, col.name = "AFDuration")



DimPlot(AurN.combined, reduction = "umap", label = TRUE, repel = FALSE)
DimPlot(AurN.combined, reduction = "umap", split.by = "stim")
DimPlot(AurN.combined, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(AurN.combined, reduction = "umap", split.by = "AFDuration")



saveRDS(AurN.combined, "~/Desktop/Human_LA_LAA/RDS/AurN_cluster2_sn_LAA_All.rds")
AurN.combined<- readRDS ("~/Desktop/Human_LA_LAA/RDS/AurN_cluster2_sn_LAA_All.rds")





# Cells of the adult human heart 参照-----

"Ventricular Carrdiomyocyte"
VlnPlot(AurN.combined,pt.size = 0,features = c("MYH7","MYL2","FHL2") )

"Atrial Carrdiomyocyte"
VlnPlot(AurN.combined,pt.size = 0, features = c("NPPA","MYL7","MYL4") )
FeaturePlot (AurN.combined, features = c("ACTA2"), min.cutoff=0.00005, max.cutoff='q90')
FeaturePlot (AurN.combined, features = c("TCF21"), min.cutoff=0.00005, max.cutoff='q90')
VlnPlot(AurN.combined, features = c("FAP","TCF21") )
FeaturePlot (AurN.combined, features = c("FAP"), min.cutoff=0.00005, max.cutoff='q90')

VlnPlot(AurN.combined,pt.size = 0,features = c("TTN","MYBPC3","TNNT2","RYR2","PLN","SLC8A1") )

"Fibroblast"
VlnPlot(AurN.combined,pt.size = 0, features = c("DCN","GSN","PDGFRA"))

VlnPlot(AurN.combined,pt.size = 0, features = c("FAP","TCF21") ,split.by = "AFControl")
VlnPlot(AurN.combined, pt.size = 0, features = c("POSTN","NOX4","COL1A1","COL1A2") )

"Endothelial"
VlnPlot(AurN.combined,pt.size = 0,features = c("VWF","PECAM1","CDH5") )

"Pericyte"
VlnPlot(AurN.combined,pt.size = 0, features = c("RGS5","ABCC9","KCNJ8") )

"Smooth muscle"
VlnPlot(AurN.combined,pt.size = 0, features = c("MYH11","TAGLN","ACTA2") )

"Neuronal"
VlnPlot(AurN.combined,pt.size = 0, features = c("PLP1","NRXN1","NRXN3") )

"Lymphoid"
VlnPlot(AurN.combined,pt.size = 0, features = c("CD8A","IL7R","CD40LG") )
VlnPlot(AurN.combined,pt.size = 0, features = c("CD3E","CD4","CD8A","NKG7","GNLY","GZMB"))
VlnPlot(AurN.combined,pt.size = 0, features = c("MKI67","CD79A","CD79B"))

"Myeloid"
VlnPlot(AurN.combined,pt.size = 0, features = c("CD14","C1QA","CD68") )
VlnPlot(AurN.combined,pt.size = 0, features = c("CD68","CD14","CSF1R","IL1B","S100A8", "S100A9") )
VlnPlot(AurN.combined,pt.size = 0, features = c("TREM2","CD9","APOE","CXCL3","C1QA","TNF","CXCL8") )
VlnPlot(AurN.combined,pt.size = 0, features = c("HLA-DPA1","HLA-DPB1","HLA-DQA1","CLEC10A","FCER1A","CD1C") )
VlnPlot(AurN.combined,pt.size = 0, features = c("PTPRC","HLA-DPB1","HLA-DQA1","CLEC10A","FCER1A","CD1C") )
VlnPlot(AurN.combined,pt.size = 0, features = c("KIT") )

"Adipocytes"
VlnPlot(AurN.combined,pt.size = 0, features = c("GPAM","FASN","LEP") )

"Mesothelial"
VlnPlot(AurN.combined,pt.size = 0, features = c("MSLN","WT1","BNC1") )

VlnPlot(AurN.combined,pt.size = 0, features = c("BMP10") )
VlnPlot(AurN.combined,pt.size = 0, features = c("ACLY") ,group.by = "AFControl")
VlnPlot(AurN.combined,pt.size = 0, features = c("TIMP1") ,split.by = "AFControl")


AurN.combined_trans <- SCTransform (AurN.combined)

cd_genes <- c("DCN","GSN","PDGFRA","NPPA","MYL7","MYL4","VWF","PECAM1","CDH5","RGS5","ABCC9","KCNJ8","MYH11","TAGLN","ACTA2","CD3E","CD4","CD8A","NKG7","GNLY","GZMB","IL7R","CD40LG","CD14","C1QA","CD68","HLA-DPA1","HLA-DPB1","HLA-DQA1","CLEC10A","FCER1A","CD1C","MKI67","CD79A","CD79B","KIT","HDC")
DotPlot(AurN.combined_trans,features = cd_genes)+RotatedAxis()+coord_flip()


VlnPlot(AurN.combined,pt.size = 0, features = c("FKBP5","HIF1A","SLC8A1") )



AurN.combined_prop. <- prop.table(table(Idents(AurN.combined),AurN.combined))
write.table(AurN.combined_prop., "AurN.combined_prop.txt", quote=F, col.names=F, append=T)
write.csv(AurN.combined_prop., "AurN.combined_prop.csv")

AurN.combined_prop. <- prop.table(table(Idents(AurN.combined),AurN.combined$AFControl))
write.table(AurN.combined_prop., "AurN.combined_prop.AFControl.txt", quote=F, col.names=F, append=T)
write.csv(AurN.combined_prop., "AurN.combined_prop.Control.csv")

AurN.combined_prop. <- prop.table(table(Idents(AurN.combined),AurN.combined$AFDuration))
write.table(AurN.combined_prop., "AurN.combined_prop.AFDuration.txt", quote=F, col.names=F, append=T)
write.csv(AurN.combined_prop., "AurN.combined_prop.AFDuration.csv")



#All(コンタミ除去後)----
AurN.combined_1 <- subset(AurN.combined, idents = c(0,1,3,4,5,6,7,8,9,10,11,12,13))
AurN.combined_1 <- ScaleData(AurN.combined_1, verbose = FALSE)
AurN.combined_1 <- RunPCA(AurN.combined_1, npcs = 30, verbose = FALSE)
AurN.combined_1 <- RunUMAP(AurN.combined_1 , reduction = "pca", dims = 1:30)
AurN.combined_1 <- FindNeighbors(AurN.combined_1 , reduction = "pca", dims = 1:30)
AurN.combined_1 <- FindClusters(AurN.combined_1, resolution = 0.25)

DimPlot(AurN.combined_1, reduction = "umap", label = FALSE, repel = FALSE)
DimPlot(AurN.combined_1, reduction = "umap", split.by = "stim")
DimPlot(AurN.combined_1, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(AurN.combined_1, reduction = "umap", split.by = "AFDuration")


# Olink----
VlnPlot(AurN.combined_1, pt.size = 0, features = c("IL17A","KCND2","RASGRF1","B4GALNT3","MTERF3", "FABP6") )
VlnPlot(AurN.combined_1, pt.size = 0, features = c("RNF4","ARNT","BIRC3","B4GALNT3","MTERF3", "FABP6") )
VlnPlot(AurN.combined_1, pt.size = 0, features = c("IL17A","IL17C","CXCL11","MMP12","MTERF3", "FABP6") )



# 並べ替える
levels(AurN.combined_1) # [1] "0" "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11"
levels(AurN.combined_1) <- c("0", "2" ,"4", "1", "5","10","3","7","6","8","9","11")

new.cluster.ids <- c("0","1","2","3","4","5","6","7","8","9","10","11")
names(new.cluster.ids) <- levels(AurN.combined_1)
AurN.combined_1<- RenameIdents(AurN.combined_1, new.cluster.ids)


saveRDS(AurN.combined_1, "~/Desktop/Human_LA_LAA/RDS/AurN_cluster_sn_LAA_All最終.rds")
AurN.combined_1<- readRDS ("~/Desktop/Human_LA_LAA/RDS/AurN_cluster_sn_LAA_All最終.rds")


VlnPlot(AurN.combined_1,pt.size = 0,features = c("PDGFA","PDGRRA","PDGFB","PDGFRB","PDGFC","PDGFRC") )
VlnPlot(AurN.combined_1,pt.size = 0,features = c("PDGFRL","PDAP1","TGFB1","TGFB2","TGFB1R","TGFB2R") )
VlnPlot(AurN.combined_1,pt.size = 0,features = c("TGFB3","TGFB3R","IL17A","IL17RA","IL17B","IL17RB") )
VlnPlot(AurN.combined_1,pt.size = 0,features = c("IL17C","IL17RC","IL17D","IL17RD","IL17E","IL17RE") )
VlnPlot(AurN.combined_1,pt.size = 0,features = c("IL17REL","IL17F","IL17RF","MMP12","MMP13","SLCA12") )
VlnPlot(AurN.combined_1,pt.size = 0,features = c("HFE","CD63","ITGB1","CTRH1","IL6","IL11") )
VlnPlot(AurN.combined_1,pt.size = 0,features = c("BRD4","CX3CR1","MEOX1","CXCL10","CCL19","IL11") )


VlnPlot(AurN.combined_1,pt.size = 0,features = c("COL6A3","IL1R1","MEOX1","CXCL10","CCL19","IL11") )
VlnPlot(AurN.combined_1,pt.size = 0,features = c("MEG3","BICC1","ABI3BP","LAM2","FGF14","IL11") )

FeaturePlot(AurN.combined_1, features=c('MEG3'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_1, features=c('TAGLN'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_1, features=c('TTN'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_1, features=c('MYH11'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_1, features=c('ACTA2'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_1, features=c('RGS5'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_1, features=c('ABCC9'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_1, features=c('KCNJ8'), min.cutoff=0.5, max.cutoff='q90')

FeaturePlot(AurN.combined_1, features=c('SPP1'), min.cutoff=0.5, max.cutoff='q90')


# Cells of the adult human heart 参照-----

"Atrial Carrdiomyocyte"
VlnPlot(AurN.combined_1,pt.size = 0, features = c("NPPA","MYL7","MYL4") )
VlnPlot(AurN.combined_1,pt.size = 0,features = c("TTN","MYBPC3","TNNT2","RYR2","PLN","SLC8A1") )

"Fibroblast"
VlnPlot(AurN.combined_1,pt.size = 0, features = c("DCN","GSN","PDGFRA"))

VlnPlot(AurN.combined_1,pt.size = 0, features = c("FAP","TCF21") ,split.by = "AFControl")
VlnPlot(AurN.combined_1, pt.size = 0, features = c("POSTN","NOX4","COL1A1","COL1A2") )

"Endothelial"
VlnPlot(AurN.combined_1,pt.size = 0,features = c("VWF","PECAM1","CDH5") )

"Pericyte"
VlnPlot(AurN.combined_1,pt.size = 0, features = c("RGS5","ABCC9","KCNJ8") )

"Smooth muscle"
VlnPlot(AurN.combined_1,pt.size = 0, features = c("MYH11","TAGLN","ACTA2") )
VlnPlot(AurN.combined_1,pt.size = 0, features = c("NTRK3","RGS5") )

"Neuronal"
VlnPlot(AurN.combined_1,pt.size = 0, features = c("PLP1","NRXN1","NRXN3") )

"Lymphoid"
VlnPlot(AurN.combined_1,pt.size = 0, features = c("CD8A","IL7R","CD40LG") )
VlnPlot(AurN.combined_1,pt.size = 0, features = c("CD3E","CD4","CD8A","NKG7","GNLY","GZMB"))
VlnPlot(AurN.combined_1,pt.size = 0, features = c("MKI67","CD79A","CD79B"))
VlnPlot(AurN.combined_1,pt.size = 0, features = c("CCR7","TCF7","PTPRC","SYTL3","CD69","SKAP1"))

"Myeloid"
VlnPlot(AurN.combined_1,pt.size = 0, features = c("CD14","C1QA","CD68") )
VlnPlot(AurN.combined_1,pt.size = 0, features = c("CD68","CD14","CSF1R","IL1B","S100A8", "S100A9") )
VlnPlot(AurN.combined_1,pt.size = 0, features = c("TREM2","CD9","APOE","CXCL3","C1QA","TNF","CXCL8") )
VlnPlot(AurN.combined_1,pt.size = 0, features = c("HLA-DPA1","HLA-DPB1","HLA-DQA1","CLEC10A","FCER1A","CD1C") )
VlnPlot(AurN.combined_1,pt.size = 0, features = c("PTPRC","HLA-DPB1","HLA-DQA1","CLEC10A","FCER1A","CD1C") )
VlnPlot(AurN.combined_1,pt.size = 0, features = c("KIT","CD163") )

"Adipocytes"
VlnPlot(AurN.combined_1,pt.size = 0, features = c("GPAM","FASN","LEP") )

"Mesothelial"
VlnPlot(AurN.combined_1,pt.size = 0, features = c("MSLN","WT1","BNC1") )

FeaturePlot(AurN.combined_1, features=c('GJA1'), min.cutoff=0.5, max.cutoff='q90')


FeaturePlot(AurN.combined_1, features=c('POSTN'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_1, features=c('COL3A1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_1, features=c('DCN'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_1, features=c('NPPA'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_1, features=c('VWF'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_1, features=c('SORBS2'), min.cutoff=0.5, max.cutoff='q90', split.by = "AFDuration")


FeaturePlot(AurN.combined_1, features=c('MYBPC3'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_1, features=c('CD163'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_1, features=c('IL7R'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_1, features=c('MKI67'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_1, features=c('VWF'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_1, features=c('RGS5'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_1, features=c('MYH11'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_1, features=c('GPAM'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_1, features=c('NRXN1'), min.cutoff=0.5, max.cutoff='q90')

FeaturePlot(AurN.combined_1, features=c('COL1A1','MYBPC3','CD163','IL7R','MKI67','VWF','RGS5','MYH11','GPAM','NRXN1'), min.cutoff=0.5, max.cutoff='q90')

AurN.combined_1.markers <- FindAllMarkers(AurN.combined_1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
AurN.combined_1.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- AurN.combined_1.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(AurN.combined_1, features = top10$gene) 
write.table(top10, "table.txt", quote=F, col.names=F, append=T)




AurN.combined_1_trans <- SCTransform (AurN.combined_1)

cd_genes <- c("DCN","COL1A1","PDGFRA","MYL4","MYBPC3","RYR2","CD163","CD14","PTPRC","IL7R","MKI67","VWF","PECAM1","CDH5","RGS5","ABCC9","KCNJ8","MYH11","ACTA2","GPAM","FASN","LEP","PLP1","NRXN1","NRXN3")
DotPlot(AurN.combined_1_trans,features = cd_genes)+RotatedAxis()+coord_flip()

AurN.combined_1_trans <- SCTransform (AurN.combined_1)

cd_genes <- c("COL1A1","PDGFRA","MEG3","MYL4","MYBPC3","RYR2","CD163","CD14","PTPRC","IL7R","MKI67","VWF","PECAM1","CDH5","RGS5","ABCC9","KCNJ8","MYH11","ACTA2","GPAM","FASN","LEP","PLP1","NRXN1","NRXN3")
DotPlot(AurN.combined_1_trans,features = cd_genes)+RotatedAxis()+coord_flip()



AurN.combined_1_prop. <- prop.table(table(Idents(AurN.combined_1)))
write.table(AurN.combined_1_prop., "AurN.combined_1_prop.txt", quote=F, col.names=F, append=T)
write.csv(AurN.combined_1_prop., "AurN.combined_1_prop.csv")


AurN.combined_1_prop. <- prop.table(table(Idents(AurN.combined_1),AurN.combined_1$AFDuration))
write.table(AurN.combined_1_prop., "AurN.combined_1_prop.AFDuration.txt", quote=F, col.names=F, append=T)
write.csv(AurN.combined_1_prop., "AurN.combined_1_prop.AFDuration.csv")


AurN.combined_1_prop. <- prop.table(table(Idents(AurN.combined_1),AurN.combined_1$stim))
write.table(AurN.combined_1_prop., "AurN.combined_1_prop.sep.txt", quote=F, col.names=F, append=T)
write.csv(AurN.combined_1_prop., "AurN.combined_1_prop.sep.csv")


#非免疫細胞----
AurN.combined_1_1 <- subset(AurN.combined_1, idents = c(0,1,2,6,7,8,9,10,11))

DimPlot(AurN.combined_1_1, reduction = "umap", label = FALSE, repel = TRUE)

#非免疫細胞 Cell chat----
saveRDS(AurN.combined_1_1, "~/Desktop/Human_LA_LAA/RDS/AurN_cluster_sn_LAA_非免疫細胞.rds")
AurN.combined_1_1<- readRDS ("~/Desktop/Human_LA_LAA/RDS/AurN_cluster_sn_LAA_非免疫細胞.rds")



#fibroblast----
AurN.combined_fibroblast <- subset(AurN.combined_1, idents = c(0))
AurN.combined_fibroblast <- ScaleData(AurN.combined_fibroblast, verbose = FALSE)
AurN.combined_fibroblast<- RunPCA(AurN.combined_fibroblast, npcs = 30, verbose = FALSE)
AurN.combined_fibroblast <- RunUMAP(AurN.combined_fibroblast , reduction = "pca", dims = 1:30)
AurN.combined_fibroblast <- FindNeighbors(AurN.combined_fibroblast , reduction = "pca", dims = 1:30)
AurN.combined_fibroblast <- FindClusters(AurN.combined_fibroblast, resolution = 0.35)



DimPlot(AurN.combined_fibroblast, reduction = "umap", label = FALSE, repel = TRUE)
DimPlot(AurN.combined_fibroblast, reduction = "umap", label = TRUE,  repel = TRUE)
DimPlot(AurN.combined_fibroblast, reduction = "umap", split.by = "AFDuration")
DimPlot(AurN.combined_fibroblast, reduction = "umap", split.by = "stim")


saveRDS(AurN.combined_fibroblast, "~/Desktop/Human_LA_LAA/RDS/AurN_cluster_sn_LAA_fibroblast.rds")
AurN.combined_fibroblast<- readRDS ("~/Desktop/Human_LA_LAA/RDS/AurN_cluster_sn_LAA_fibroblast.rds")


VlnPlot(AurN.combined_fibroblast,pt.size = 0, features = c("DCN","GSN","PDGFRA","POSTN","TNC","PDGFRA") )



# コンタミ探し
VlnPlot(AurN.combined_fibroblast, pt.size = 0, features = c("NPPA","MYL7","MYL4","TTN","MYBPC3","TNNT2"))
VlnPlot(AurN.combined_fibroblast, pt.size = 0, features = c("RYR2","PLN","SLC8A1","DCN","GSN","PDGFRA"))
VlnPlot(AurN.combined_fibroblast, pt.size = 0, features = c("VWF","PECAM1","CDH5","RGS5","ABCC9","KCNJ8"))
VlnPlot(AurN.combined_fibroblast, pt.size = 0, features = c("MYH11","TAGLN","ACTA2","PLP1","NRXN1","NRXN3"))
VlnPlot(AurN.combined_fibroblast, pt.size = 0, features = c("GPAM","FASN","LEP","MSLN","WT1","BNC1"))
VlnPlot(AurN.combined_fibroblast, pt.size = 0, features = c("PTPRC","CD163","CD79B","CD14","IL7R","CD3E","CD4","CD8A","NKG7","GNLY","GZMB"))

FeaturePlot(AurN.combined_fibroblast, features=c('MYL7'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblast, features=c('NPPA'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblast, features=c('TTN'), min.cutoff=0.5, max.cutoff='q90')

# ヒートマップの作成------------------------------------------------------------------------------

AurN.combined_fibroblast.markers <- FindAllMarkers(AurN.combined_fibroblast, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
AurN.combined_fibroblast.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- AurN.combined_fibroblast.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(AurN.combined_fibroblast, features = top10$gene) 
write.table(top10, "table.txt", quote=F, col.names=F, append=T)


VlnPlot(AurN.combined_fibroblast, pt.size = 0, features = c("TIMP1","TIMP2","TIMP3","MMP2","MMP9","TGFB1"))
VlnPlot(AurN.combined_fibroblast, pt.size = 0, features = c("C1QA","C1QB","C1QC","IL10","CD9","TREM2"))


# aFB Nature 2020 cells of the adult human heart
VlnPlot(AurN.combined_fibroblast, pt.size = 0, features = c("GSN","PCOLCE2","FBLN1","ABCA10","SCN7A","ELL2"))
VlnPlot(AurN.combined_fibroblast, pt.size = 0, features = c("GSN","PCOLCE2","NAMPT","ATP13A3","VWF","PTPRB"))

FB4　TGFB
VlnPlot(AurN.combined_fibroblast, pt.size = 0, features = c("TGFB","POSTN","TNC"))


FB6　ECM
VlnPlot(AurN.combined_fibroblast, pt.size = 0, features = c("OSMR","ILST6"))


# Clin Trans Med 2023
VlnPlot(AurN.combined_fibroblast, pt.size = 0, features = c("PDRFRL","MFAP5","CLIP","CXCL1","ACKR3","TNFAIP6"))
VlnPlot(AurN.combined_fibroblast, pt.size = 0, features = c("PTX3","SERINA3","ITLN1","LGBP3","ID3","RNASE1"))

VlnPlot(AurN.combined_fibroblast, pt.size = 0, features = c("TNF","TLR4","IL1B","CXCL3","CCL2","CCL5") )

VlnPlot(AurN.combined_fibroblast, pt.size = 0, features = c("ACTA2","FAP","POSTN","LUM","COL3A1") )

# Nature 2022 Single-nucleus profiling of human dilated and hypertrophic cardiomyopathy LV
Activated fibroblast
VlnPlot(AurN.combined_fibroblast, pt.size = 0, features = c("POSTN","NOX4","COL1A1","COL1A2") )
FB-TLL2
VlnPlot(AurN.combined_fibroblast, pt.size = 0, features = c("TLL2","ERI2","ACSM1","THUMPD1") )
Act
VlnPlot(AurN.combined_fibroblast, pt.size = 0, features = c("POSTN","TNC","THBS4","MEOX01") )
FB-ZBTB7C
VlnPlot(AurN.combined_fibroblast, pt.size = 0, features = c("ZBTB7C","PCOLCE2","SEMA3C","KAZN") )

FB-X1
VlnPlot(AurN.combined_fibroblast, pt.size = 0, features = c("TTN","MT-CO2","MT-ATP6","KMT-CO1"))

FB-PTCHD4
VlnPlot(AurN.combined_fibroblast, pt.size = 0, features = c("PTCHD4","ZMAT3","DOB2","PAPPA") )



#fibroblast(心筋doublet除去後)----
AurN.combined_fibroblast1 <- subset(AurN.combined_fibroblast, idents = c(0,2,3,4))
AurN.combined_fibroblast1 <- ScaleData(AurN.combined_fibroblast1, verbose = FALSE)
AurN.combined_fibroblast1<- RunPCA(AurN.combined_fibroblast1, npcs = 30, verbose = FALSE)
AurN.combined_fibroblast1 <- RunUMAP(AurN.combined_fibroblast1 , reduction = "pca", dims = 1:30)
AurN.combined_fibroblast1 <- FindNeighbors(AurN.combined_fibroblast1 , reduction = "pca", dims = 1:30)
AurN.combined_fibroblast1 <- FindClusters(AurN.combined_fibroblast1, resolution = 0.35)



DimPlot(AurN.combined_fibroblast1, reduction = "umap", label = FALSE, repel = TRUE)
DimPlot(AurN.combined_fibroblast1, reduction = "umap", label = TRUE,  repel = TRUE)
DimPlot(AurN.combined_fibroblast1, reduction = "umap", split.by = "AFDuration")
DimPlot(AurN.combined_fibroblast1, reduction = "umap", split.by = "stim")


saveRDS(AurN.combined_fibroblast1, "~/Desktop/Human_LA_LAA/RDS/AurN_cluster_sn_LAA_fibroblast1.rds")
AurN.combined_fibroblast1<- readRDS ("~/Desktop/Human_LA_LAA/RDS/AurN_cluster_sn_LAA_fibroblast1.rds")




VlnPlot(AurN.combined_fibroblast1,pt.size = 0, features = c("DCN","GSN","PDGFRA","POSTN","TNC","PDGFRA") )



# コンタミ探し
VlnPlot(AurN.combined_fibroblast1, pt.size = 0, features = c("NPPA","MYL7","MYL4","TTN","MYBPC3","TNNT2"))
VlnPlot(AurN.combined_fibroblast1, pt.size = 0, features = c("RYR2","PLN","SLC8A1","DCN","GSN","PDGFRA"))
VlnPlot(AurN.combined_fibroblast1, pt.size = 0, features = c("VWF","PECAM1","CDH5","RGS5","ABCC9","KCNJ8"))
VlnPlot(AurN.combined_fibroblast1, pt.size = 0, features = c("MYH11","TAGLN","ACTA2","PLP1","NRXN1","NRXN3"))
VlnPlot(AurN.combined_fibroblast1, pt.size = 0, features = c("GPAM","FASN","LEP","MSLN","WT1","BNC1"))
VlnPlot(AurN.combined_fibroblast1, pt.size = 0, features = c("PTPRC","CD163","CD79B","CD14","IL7R","CD3E","CD4","CD8A","NKG7","GNLY","GZMB"))


# ヒートマップの作成------------------------------------------------------------------------------

AurN.combined_fibroblast1.markers <- FindAllMarkers(AurN.combined_fibroblast1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
AurN.combined_fibroblast1.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- AurN.combined_fibroblast1.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(AurN.combined_fibroblast1, features = top10$gene) 
write.table(top10, "table.txt", quote=F, col.names=F, append=T)


VlnPlot(AurN.combined_fibroblast1, pt.size = 0, features = c("TIMP1","TIMP2","TIMP3","MMP2","MMP9","TGFB1"))
VlnPlot(AurN.combined_fibroblast1, pt.size = 0, features = c("C1QA","C1QB","C1QC","IL10","CD9","TREM2"))


# aFB Nature 2020 cells of the adult human heart
VlnPlot(AurN.combined_fibroblast1, pt.size = 0, features = c("GSN","PCOLCE2","FBLN1","ABCA10","SCN7A","ELL2"))
VlnPlot(AurN.combined_fibroblast1, pt.size = 0, features = c("GSN","PCOLCE2","NAMPT","ATP13A3","VWF","PTPRB"))

FB4　TGFB
VlnPlot(AurN.combined_fibroblast1, pt.size = 0, features = c("TGFB","POSTN","TNC"))


FB6　ECM
VlnPlot(AurN.combined_fibroblast1, pt.size = 0, features = c("OSMR","ILST6"))


# Clin Trans Med 2023
VlnPlot(AurN.combined_fibroblast1, pt.size = 0, features = c("PDRFRL","MFAP5","CLIP","CXCL1","ACKR3","TNFAIP6"))
VlnPlot(AurN.combined_fibroblast1, pt.size = 0, features = c("PTX3","SERINA3","ITLN1","LGBP3","ID3","RNASE1"))

VlnPlot(AurN.combined_fibroblast1, pt.size = 0, features = c("TNF","TLR4","IL1B","CXCL3","CCL2","CCL5") )

VlnPlot(AurN.combined_fibroblast1, pt.size = 0, features = c("ACTA2","FAP","POSTN","LUM","COL3A1") )

# Nature 2022 Single-nucleus profiling of human dilated and hypertrophic cardiomyopathy LV
Activated fibroblast
VlnPlot(AurN.combined_fibroblast1, pt.size = 0, features = c("POSTN","NOX4","COL1A1","COL1A2") )
FB-TLL2
VlnPlot(AurN.combined_fibroblast1, pt.size = 0, features = c("TLL2","ERI2","ACSM1","THUMPD1") )
Act
VlnPlot(AurN.combined_fibroblast1, pt.size = 0, features = c("POSTN","TNC","THBS4","MEOX01") )
FB-ZBTB7C
VlnPlot(AurN.combined_fibroblast1, pt.size = 0, features = c("ZBTB7C","PCOLCE2","SEMA3C","KAZN") )

FB-X1
VlnPlot(AurN.combined_fibroblast1, pt.size = 0, features = c("TTN","MT-CO2","MT-ATP6","KMT-CO1"))

FB-PTCHD4
VlnPlot(AurN.combined_fibroblast1, pt.size = 0, features = c("PTCHD4","ZMAT3","DOB2","PAPPA") )




#最終fibroblast(再度心筋doublet除去後)----
AurN.combined_fibroblast2 <- subset(AurN.combined_fibroblast1, idents = c(0,1,2,3))
AurN.combined_fibroblast2 <- ScaleData(AurN.combined_fibroblast2, verbose = FALSE)
AurN.combined_fibroblast2 <- RunPCA(AurN.combined_fibroblast2, npcs = 30, verbose = FALSE)
AurN.combined_fibroblast2 <- RunUMAP(AurN.combined_fibroblast2 , reduction = "pca", dims = 1:30)
AurN.combined_fibroblast2 <- FindNeighbors(AurN.combined_fibroblast2 , reduction = "pca", dims = 1:30)
AurN.combined_fibroblast2 <- FindClusters(AurN.combined_fibroblast2, resolution = 0.35)



DimPlot(AurN.combined_fibroblast2, reduction = "umap", label = FALSE, repel = TRUE)
DimPlot(AurN.combined_fibroblast2, reduction = "umap", label = TRUE,  repel = TRUE)
DimPlot(AurN.combined_fibroblast2, reduction = "umap", split.by = "AFDuration")
DimPlot(AurN.combined_fibroblast2, reduction = "umap", split.by = "stim")


saveRDS(AurN.combined_fibroblast2, "~/Desktop/Human_LA_LAA/RDS/AurN_cluster_sn_LAA_fibroblast2.rds")
AurN.combined_fibroblast2<- readRDS ("~/Desktop/Human_LA_LAA/RDS/AurN_cluster_sn_LAA_fibroblast2.rds")




VlnPlot(AurN.combined_fibroblast2,pt.size = 0, features = c("DCN","GSN","PDGFRA","POSTN","TNC","PDGFRA") )



# コンタミ探し
VlnPlot(AurN.combined_fibroblast2, pt.size = 0, features = c("NPPA","MYL7","MYL4","TTN","MYBPC3","TNNT2"))
VlnPlot(AurN.combined_fibroblast2, pt.size = 0, features = c("RYR2","PLN","SLC8A1","DCN","GSN","PDGFRA"))
VlnPlot(AurN.combined_fibroblast2, pt.size = 0, features = c("VWF","PECAM1","CDH5","RGS5","ABCC9","KCNJ8"))
VlnPlot(AurN.combined_fibroblast2, pt.size = 0, features = c("MYH11","TAGLN","ACTA2","PLP1","NRXN1","NRXN3"))
VlnPlot(AurN.combined_fibroblast2, pt.size = 0, features = c("GPAM","FASN","LEP","MSLN","WT1","BNC1"))
VlnPlot(AurN.combined_fibroblast2, pt.size = 0, features = c("PTPRC","CD163","CD79B","CD14","IL7R","CD3E","CD4","CD8A","NKG7","GNLY","GZMB"))







# ヒートマップの作成------------------------------------------------------------------------------

AurN.combined_fibroblast2.markers <- FindAllMarkers(AurN.combined_fibroblast2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
AurN.combined_fibroblast2.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- AurN.combined_fibroblast2.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(AurN.combined_fibroblast2, features = top10$gene) 
write.table(top10, "table.txt", quote=F, col.names=F, append=T)


VlnPlot(AurN.combined_fibroblast2, pt.size = 0, features = c("TIMP1","TIMP2","TIMP3","MMP2","MMP9","IL10"))
VlnPlot(AurN.combined_fibroblast2, pt.size = 0, features = c("C1QA","C1QB","C1QC","IL10","CD9","TREM2"))
VlnPlot(AurN.combined_fibroblast2, pt.size = 0, features = c("CNN1"))


# aFB Nature 2020 cells of the adult human heart
VlnPlot(AurN.combined_fibroblast2, pt.size = 0, features = c("GSN","PCOLCE2","FBLN1","ABCA10","SCN7A","ELL2"))
VlnPlot(AurN.combined_fibroblast2, pt.size = 0, features = c("GSN","PCOLCE2","NAMPT","ATP13A3","VWF","PTPRB"))

FB4　TGFB
VlnPlot(AurN.combined_fibroblast2, pt.size = 0, features = c("TGFB","POSTN","TNC"))


FB6　ECM
VlnPlot(AurN.combined_fibroblast2, pt.size = 0, features = c("OSMR","ILST6"))


# Clin Trans Med 2023
VlnPlot(AurN.combined_fibroblast2, pt.size = 0, features = c("PDRFRL","MFAP5","CLIP","CXCL1","ACKR3","TNFAIP6"))
VlnPlot(AurN.combined_fibroblast2, pt.size = 0, features = c("PTX3","SERINA3","ITLN1","LGBP3","ID3","RNASE1"))
VlnPlot(AurN.combined_fibroblast2, pt.size = 0, features = c("TNF","TLR4","IL1B","CXCL3","CCL2","CCL5") )
VlnPlot(AurN.combined_fibroblast2, pt.size = 0.1, features = c("TNF","TLR4","IL1B","CXCL3","CCL2","CCL5") )

VlnPlot(AurN.combined_fibroblast2, pt.size = 0, features = c("ACTA2","FAP","POSTN","LUM","COL3A1") )

# Nature 2022 Single-nucleus profiling of human dilated and hypertrophic cardiomyopathy LV
Activated fibroblast
VlnPlot(AurN.combined_fibroblast2, pt.size = 0, features = c("POSTN","NOX4","COL1A1","COL1A2") )
FB-TLL2
VlnPlot(AurN.combined_fibroblast2, pt.size = 0, features = c("TLL2","ERI2","ACSM1","THUMPD1") )
Act
VlnPlot(AurN.combined_fibroblast2, pt.size = 0, features = c("POSTN","TNC","THBS4","MEOX01") )
FB-ZBTB7C
VlnPlot(AurN.combined_fibroblast2, pt.size = 0, features = c("ZBTB7C","PCOLCE2","SEMA3C","KAZN") )

FB-X1
VlnPlot(AurN.combined_fibroblast2, pt.size = 0, features = c("TTN","MT-CO2","MT-ATP6","KMT-CO1"))

FB-PTCHD4
VlnPlot(AurN.combined_fibroblast2, pt.size = 0, features = c("PTCHD4","ZMAT3","DOB2","PAPPA") )


VlnPlot(AurN.combined_fibroblast2,pt.size = 0, features = c("ABCA9","ABCA6","ABCA10","TIMP1") )

VlnPlot(AurN.combined_fibroblast2, pt.size = 0, features = c("PDGF","PDRFRA","SPP1","IL1B","IL1R1","AREG"))
VlnPlot(AurN.combined_fibroblast2, pt.size = 0, features = c("EGFR","CCL2","CCR2","CSF1","CSF1R","IL34"))


# J Adv Res 2022 China AF snRNAseq
FB1 myofibroblast
VlnPlot(AurN.combined_fibroblast2, pt.size = 0, features = c("ACTG2","MYH11","ITGA8"))
FB2 
VlnPlot(AurN.combined_fibroblast2, pt.size = 0, features = c("POSTN","ELN"))
FB3 lipid trafficking
VlnPlot(AurN.combined_fibroblast2, pt.size = 0, features = c("NEGR1","ABCA8","ABCA6") )
FB4 cell adhesion and proliferaion
VlnPlot(AurN.combined_fibroblast2, pt.size = 0, features = c("KCND2","KCTD16","FGF7","FLRT2") )
FB5 macropharge
VlnPlot(AurN.combined_fibroblast2, pt.size = 0, features = c("CD163","MRC1") )
FB5 Tcell-like
VlnPlot(AurN.combined_fibroblast2, pt.size = 0, features = c("SKAPI","ITK") )


VlnPlot(AurN.combined_fibroblast2, pt.size = 0, features = c("FAP","THBS4","TSHZ2","FAM155A","COL1A2","MIR100HG","NOX4","COL22A1","COL1A1") )
VlnPlot(AurN.combined_fibroblast2, pt.size = 0, features = c("AEBP1","CHAD","COMP") )
VlnPlot(AurN.combined_fibroblast2, pt.size = 0.1, features = c("COL3A1","LUM","COMP") )
VlnPlot(AurN.combined_fibroblast2, pt.size = 0.1, features = c("ACTA2","SMA") )
VlnPlot(CTEPH.combined_SMC_EC_1, pt.size = 0, features = c("PDGF","PDRFRA","SPP1","IL1B","IL1R1","AREG"))
VlnPlot(CTEPH.combined_SMC_EC_1, pt.size = 0, features = c("EGFR","CCL2","CCR2","CSF1","CSF1R","IL34"))


AurN.combined_fibroblast2_trans <- SCTransform (AurN.combined_fibroblast2)

cd_genes <- c("ABCA6","ABCA8","CD9","FBLN1","FBLN2","TIMP1","MMP2","COL3A1","LUM","POSTN","CCL2","EGR1","ELN","FBN1","PRG4")
DotPlot(AurN.combined_fibroblast2_trans,features = cd_genes)+RotatedAxis()+coord_flip()

FeaturePlot(AurN.combined_fibroblast2, features=c('MYH11'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblast2, features=c('LUM'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblast2, features=c('COL3A1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblast2, features=c('POSTN'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblast2, features=c('FBLN1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblast2, features=c('FBN1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblast2, features=c('CCL2'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblast2, features=c('ABCA6'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblast2, features=c('EGR1'), min.cutoff=0.5, max.cutoff='q90')

FeaturePlot(AurN.combined_fibroblast2, features=c('PRG4'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblast2, features=c('MMP2'), min.cutoff=0.5, max.cutoff='q90')

VlnPlot(AurN.combined_fibroblast2, pt.size = 0, features = c("IL1R1","CXCR2") )


AurN.combined_fibroblast2_prop. <- prop.table(table(Idents(AurN.combined_fibroblast2)))
write.table(AurN.combined_fibroblast2_prop., "AurN.combined_fibroblast2_prop.txt", quote=F, col.names=F, append=T)
write.csv(AurN.combined_fibroblast2_prop., "AurN.combined_fibroblast2_prop.csv")


AurN.combined_fibroblast2_prop. <- prop.table(table(Idents(AurN.combined_fibroblast2),AurN.combined_fibroblast2$AFDuration))
write.table(AurN.combined_fibroblast2_prop., "AurN.combined_fibroblast2_prop.AFDuration.txt", quote=F, col.names=F, append=T)
write.csv(AurN.combined_fibroblast2_prop., "AurN.combined_fibroblast2_prop.AFDuration.csv")













#fibroblast+myofibroblast(cardiomyocytweコンタミ除去)----
AurN.combined_fibroblastmyofibroblast1 <- subset(AurN.combined_fibroblastmyofibroblast, idents = c(1,2,3,4))
AurN.combined_fibroblastmyofibroblast1 <- ScaleData(AurN.combined_fibroblastmyofibroblast, verbose = FALSE)
AurN.combined_fibroblastmyofibroblast1 <- RunPCA(AurN.combined_fibroblastmyofibroblast, npcs = 30, verbose = FALSE)
AurN.combined_fibroblastmyofibroblast1 <- RunUMAP(AurN.combined_fibroblastmyofibroblast , reduction = "pca", dims = 1:30)
AurN.combined_fibroblastmyofibroblast1 <- FindNeighbors(AurN.combined_fibroblastmyofibroblast , reduction = "pca", dims = 1:30)
AurN.combined_fibroblastmyofibroblast1 <- FindClusters(AurN.combined_fibroblastmyofibroblast, resolution = 0.35)



DimPlot(AurN.combined_fibroblastmyofibroblast, reduction = "umap", label = FALSE, repel = TRUE)
DimPlot(AurN.combined_fibroblastmyofibroblast, reduction = "umap", label = TRUE,  repel = TRUE)
DimPlot(AurN.combined_fibroblastmyofibroblast, reduction = "umap", split.by = "AFDuration")
DimPlot(AurN.combined_fibroblastmyofibroblast, reduction = "umap", split.by = "stim")


saveRDS(AurN.combined_fibroblastmyofibroblast, "~/Desktop/Human_LA_LAA/RDS/AurN_cluster_sn_LAA_fibroblastmyofibroblast.rds")
AurN.combined_fibroblastmyofibroblast<- readRDS ("~/Desktop/Human_LA_LAA/RDS/AurN_cluster_sn_LAA_fibroblastmyofibroblast.rds")




VlnPlot(AurN.combined_fibroblastmyofibroblast,pt.size = 0, features = c("DCN","GSN","PDGFRA","POSTN","TNC","PDGFRA") )



# コンタミ探し
VlnPlot(AurN.combined_fibroblastmyofibroblast, pt.size = 0, features = c("NPPA","MYL7","MYL4","TTN","MYBPC3","TNNT2"))
VlnPlot(AurN.combined_fibroblastmyofibroblast, pt.size = 0, features = c("RYR2","PLN","SLC8A1","DCN","GSN","PDGFRA"))
VlnPlot(AurN.combined_fibroblastmyofibroblast, pt.size = 0, features = c("VWF","PECAM1","CDH5","RGS5","ABCC9","KCNJ8"))
VlnPlot(AurN.combined_fibroblastmyofibroblast, pt.size = 0, features = c("MYH11","TAGLN","ACTA2","PLP1","NRXN1","NRXN3"))
VlnPlot(AurN.combined_fibroblastmyofibroblast, pt.size = 0, features = c("GPAM","FASN","LEP","MSLN","WT1","BNC1"))
VlnPlot(AurN.combined_fibroblastmyofibroblast, pt.size = 0, features = c("PTPRC","CD163","CD79B","CD14","IL7R","CD3E","CD4","CD8A","NKG7","GNLY","GZMB"))







# ヒートマップの作成------------------------------------------------------------------------------

AurN.combined_fibroblastmyofibroblast.markers <- FindAllMarkers(AurN.combined_fibroblastmyofibroblast, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
AurN.combined_fibroblastmyofibroblast.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- AurN.combined_fibroblastmyofibroblast.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(AurN.combined_fibroblastmyofibroblast, features = top10$gene) 
write.table(top10, "table.txt", quote=F, col.names=F, append=T)


VlnPlot(AurN.combined_fibroblastmyofibroblast, pt.size = 0, features = c("TIMP1","TIMP2","TIMP3","MMP2","MMP9","IL10"))
VlnPlot(AurN.combined_fibroblastmyofibroblast, pt.size = 0, features = c("C1QA","C1QB","C1QC","IL10","CD9","TREM2"))
VlnPlot(AurN.combined_fibroblastmyofibroblast, pt.size = 0, features = c("CNN1"))


# aFB Nature 2020 cells of the adult human heart
VlnPlot(AurN.combined_fibroblastmyofibroblast, pt.size = 0, features = c("GSN","PCOLCE2","FBLN1","ABCA10","SCN7A","ELL2"))
VlnPlot(AurN.combined_fibroblastmyofibroblast, pt.size = 0, features = c("GSN","PCOLCE2","NAMPT","ATP13A3","VWF","PTPRB"))

FB4　TGFB
VlnPlot(AurN.combined_fibroblastmyofibroblast, pt.size = 0, features = c("TGFB","POSTN","TNC"))


FB6　ECM
VlnPlot(AurN.combined_fibroblastmyofibroblast, pt.size = 0, features = c("OSMR","ILST6"))


# Clin Trans Med 2023
VlnPlot(AurN.combined_fibroblastmyofibroblast, pt.size = 0, features = c("PDRFRL","MFAP5","CLIP","CXCL1","ACKR3","TNFAIP6"))
VlnPlot(AurN.combined_fibroblastmyofibroblast, pt.size = 0, features = c("PTX3","SERINA3","ITLN1","LGBP3","ID3","RNASE1"))
VlnPlot(AurN.combined_fibroblastmyofibroblast, pt.size = 0, features = c("TNF","TLR4","IL1B","CXCL3","CCL2","CCL5") )
VlnPlot(AurN.combined_fibroblastmyofibroblast, pt.size = 0.1, features = c("TNF","TLR4","IL1B","CXCL3","CCL2","CCL5") )

VlnPlot(AurN.combined_fibroblast+myofibroblast, pt.size = 0, features = c("ACTA2","FAP","POSTN","LUM","COL3A1") )

# Nature 2022 Single-nucleus profiling of human dilated and hypertrophic cardiomyopathy LV
Activated fibroblast
VlnPlot(AurN.combined_fibroblastmyofibroblast, pt.size = 0, features = c("POSTN","NOX4","COL1A1","COL1A2") )
FB-TLL2
VlnPlot(AurN.combined_fibroblastmyofibroblast, pt.size = 0, features = c("TLL2","ERI2","ACSM1","THUMPD1") )
Act
VlnPlot(AurN.combined_fibroblastmyofibroblast, pt.size = 0, features = c("POSTN","TNC","THBS4","MEOX01") )
FB-ZBTB7C
VlnPlot(AurN.combined_fibroblastmyofibroblast, pt.size = 0, features = c("ZBTB7C","PCOLCE2","SEMA3C","KAZN") )

FB-X1
VlnPlot(AurN.combined_fibroblastmyofibroblast, pt.size = 0, features = c("TTN","MT-CO2","MT-ATP6","KMT-CO1"))

FB-PTCHD4
VlnPlot(AurN.combined_fibroblastmyofibroblast, pt.size = 0, features = c("PTCHD4","ZMAT3","DOB2","PAPPA") )


VlnPlot(AurN.combined_fibroblastmyofibroblast,pt.size = 0, features = c("ABCA9","ABCA6","ABCA10","TIMP1") )

VlnPlot(AurN.combined_fibroblastmyofibroblast, pt.size = 0, features = c("PDGF","PDRFRA","SPP1","IL1B","IL1R1","AREG"))
VlnPlot(AurN.combined_fibroblastmyofibroblast, pt.size = 0, features = c("EGFR","CCL2","CCR2","CSF1","CSF1R","IL34"))


# J Adv Res 2022 China AF snRNAseq
FB1 myofibroblast
VlnPlot(AurN.combined_fibroblastmyofibroblast, pt.size = 0, features = c("ACTG2","MYH11","ITGA8"))
FB2 
VlnPlot(AurN.combined_fibroblastmyofibroblast, pt.size = 0, features = c("POSTN","ELN"))
FB3 lipid trafficking
VlnPlot(AurN.combined_fibroblastmyofibroblast, pt.size = 0, features = c("NEGR1","ABCA8","ABCA6") )
FB4 cell adhesion and proliferaion
VlnPlot(AurN.combined_fibroblastmyofibroblast, pt.size = 0, features = c("KCND2","KCTD16","FGF7","FLRT2") )
FB5 macropharge
VlnPlot(AurN.combined_fibroblastmyofibroblast, pt.size = 0, features = c("CD163","MRC1") )
FB5 Tcell-like
VlnPlot(AurN.combined_fibroblastmyofibroblast, pt.size = 0, features = c("SKAPI","ITK") )


VlnPlot(AurN.combined_fibroblastmyofibroblast, pt.size = 0, features = c("FAP","THBS4","TSHZ2","FAM155A","COL1A2","MIR100HG","NOX4","COL22A1","COL1A1") )
VlnPlot(AurN.combined_fibroblastmyofibroblast, pt.size = 0, features = c("AEBP1","CHAD","COMP") )
VlnPlot(AurN.combined_fibroblastmyofibroblast, pt.size = 0.1, features = c("COL3A1","LUM","COMP") )
VlnPlot(AurN.combined_fibroblastmyofibroblast, pt.size = 0.1, features = c("ACTA2","SMA") )
VlnPlot(AurN.combined_fibroblastmyofibroblast, pt.size = 0, features = c("PDGF","PDRFRA","SPP1","IL1B","IL1R1","AREG"))
VlnPlot(AurN.combined_fibroblastmyofibroblast, pt.size = 0, features = c("EGFR","CCL2","CCR2","CSF1","CSF1R","IL34"))


AurN.combined_fibroblastmyofibroblast_trans <- SCTransform (AurN.combined_fibroblastmyofibroblast)

cd_genes <- c("ABCA6","ABCA8","CD9","FBLN1","FBLN2","TIMP1","MMP2","COL3A1","LUM","POSTN","CCL2","EGR1","ELN","FBN1","PRG4")
DotPlot(AurN.combined_fibroblastmyofibroblast_trans,features = cd_genes)+RotatedAxis()+coord_flip()


FeaturePlot(AurN.combined_fibroblastmyofibroblast, features=c('TTN'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast, features=c('NPPA'), min.cutoff=0.5, max.cutoff='q90')

FeaturePlot(AurN.combined_fibroblastmyofibroblast, features=c('MYH11'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast, features=c('ACTA2'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast, features=c('GSN'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast, features=c('LUM'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast, features=c('COL3A1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast, features=c('POSTN'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast, features=c('FBLN1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast, features=c('FBN1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast, features=c('CCL2'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast, features=c('ABCA6'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast, features=c('EGR1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast, features=c('MYH11'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast, features=c('PRG4'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast, features=c('MMP2'), min.cutoff=0.5, max.cutoff='q90')

VlnPlot(AurN.combined_fibroblastmyofibroblast, pt.size = 0, features = c("IL1R1","CXCR2") )


AurN.combined_fibroblast+myofibroblast_prop. <- prop.table(table(Idents(AurN.combined_fibroblast+myofibroblast)))
write.table(AurN.combined_fibroblast+myofibroblast_prop., "AurN.combined_fibroblast+myofibroblast_prop.txt", quote=F, col.names=F, append=T)
write.csv(AurN.combined_fibroblast+myofibroblast_prop., "AurN.combined_fibroblast+myofibroblast_prop.csv")


AurN.combined_fibroblast+myofibroblast_prop. <- prop.table(table(Idents(AurN.combined_fibroblast+myofibroblast),AurN.combined_fibroblast2$AFDuration))
write.table(AurN.combined_fibroblast+myofibroblast_prop., "AurN.combined_fibroblast+myofibroblast_prop.AFDuration.txt", quote=F, col.names=F, append=T)
write.csv(AurN.combined_fibroblast+myofibroblast_prop., "AurN.combined_fibroblast+myofibroblast_prop.AFDuration.csv")










#fibroblast+myofibroblast(myocardium doublet除去Aur4Nのみ認める)----
AurN.combined_fibroblastmyofibroblast1 <- subset(AurN.combined_fibroblastmyofibroblast, idents = c(1,2,3,4))
AurN.combined_fibroblastmyofibroblast1 <- ScaleData(AurN.combined_fibroblastmyofibroblast1, verbose = FALSE)
AurN.combined_fibroblastmyofibroblast1 <- RunPCA(AurN.combined_fibroblastmyofibroblast1, npcs = 30, verbose = FALSE)
AurN.combined_fibroblastmyofibroblast1 <- RunUMAP(AurN.combined_fibroblastmyofibroblast1 , reduction = "pca", dims = 1:30)
AurN.combined_fibroblastmyofibroblast1 <- FindNeighbors(AurN.combined_fibroblastmyofibroblast1 , reduction = "pca", dims = 1:30)
AurN.combined_fibroblastmyofibroblast1 <- FindClusters(AurN.combined_fibroblastmyofibroblast1, resolution = 0.35)



DimPlot(AurN.combined_fibroblastmyofibroblast1, reduction = "umap", label = FALSE, repel = TRUE)
DimPlot(AurN.combined_fibroblastmyofibroblast1, reduction = "umap", label = TRUE,  repel = TRUE)
DimPlot(AurN.combined_fibroblastmyofibroblast1, reduction = "umap", split.by = "AFDuration")
DimPlot(AurN.combined_fibroblastmyofibroblast1, reduction = "umap", split.by = "stim")


saveRDS(AurN.combined_fibroblastmyofibroblast1, "~/Desktop/Human_LA_LAA/RDS/AurN_cluster_sn_LAA_fibroblastmyofibroblast1.rds")
AurN.combined_fibroblastmyofibroblast1<- readRDS ("~/Desktop/Human_LA_LAA/RDS/AurN_cluster_sn_LAA_fibroblastmyofibroblast1.rds")




VlnPlot(AurN.combined_fibroblastmyofibroblast1,pt.size = 0, features = c("DCN","GSN","PDGFRA","POSTN","TNC","PDGFRA") )



# コンタミ探し
VlnPlot(AurN.combined_fibroblastmyofibroblast1, pt.size = 0, features = c("NPPA","MYL7","MYL4","TTN","MYBPC3","TNNT2"))
VlnPlot(AurN.combined_fibroblastmyofibroblast1, pt.size = 0, features = c("VWF","PECAM1","CDH5","RGS5","ABCC9","KCNJ8"))
VlnPlot(AurN.combined_fibroblastmyofibroblast1, pt.size = 0, features = c("MYH11","TAGLN","ACTA2","PLP1","NRXN1","NRXN3"))
VlnPlot(AurN.combined_fibroblastmyofibroblast1, pt.size = 0, features = c("GPAM","FASN","LEP","MSLN","WT1","BNC1"))
VlnPlot(AurN.combined_fibroblastmyofibroblast1, pt.size = 0, features = c("PTPRC","CD163","CD79B","CD14","IL7R","CD3E","CD4","CD8A","NKG7","GNLY","GZMB"))







# ヒートマップの作成------------------------------------------------------------------------------

AurN.combined_fibroblastmyofibroblast1.markers <- FindAllMarkers(AurN.combined_fibroblastmyofibroblast1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
AurN.combined_fibroblastmyofibroblast1.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- AurN.combined_fibroblastmyofibroblast1.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(AurN.combined_fibroblastmyofibroblast1, features = top10$gene) 
write.table(top10, "table.txt", quote=F, col.names=F, append=T)


VlnPlot(AurN.combined_fibroblastmyofibroblast1, pt.size = 0, features = c("TIMP1","TIMP2","TIMP3","MMP2","MMP9","IL10"))
VlnPlot(AurN.combined_fibroblastmyofibroblast1, pt.size = 0, features = c("C1QA","C1QB","C1QC","IL10","CD9","TREM2"))
VlnPlot(AurN.combined_fibroblastmyofibroblast1, pt.size = 0, features = c("CNN1"))


# aFB Nature 2020 cells of the adult human heart
VlnPlot(AurN.combined_fibroblastmyofibroblast1, pt.size = 0, features = c("GSN","PCOLCE2","FBLN1","ABCA10","SCN7A","ELL2"))
VlnPlot(AurN.combined_fibroblastmyofibroblast1, pt.size = 0, features = c("GSN","PCOLCE2","NAMPT","ATP13A3","VWF","PTPRB"))

FB4　TGFB
VlnPlot(AurN.combined_fibroblastmyofibroblast1, pt.size = 0, features = c("TGFB","POSTN","TNC"))


FB6　ECM
VlnPlot(AurN.combined_fibroblastmyofibroblast1, pt.size = 0, features = c("OSMR","ILST6"))


# Clin Trans Med 2023
VlnPlot(AurN.combined_fibroblastmyofibroblast1, pt.size = 0, features = c("PDRFRL","MFAP5","CLIP","CXCL1","ACKR3","TNFAIP6"))
VlnPlot(AurN.combined_fibroblastmyofibroblast1, pt.size = 0, features = c("PTX3","SERINA3","ITLN1","LGBP3","ID3","RNASE1"))
VlnPlot(AurN.combined_fibroblastmyofibroblast1, pt.size = 0, features = c("TNF","TLR4","IL1B","CXCL3","CCL2","CCL5") )
VlnPlot(AurN.combined_fibroblastmyofibroblast1, pt.size = 0.1, features = c("TNF","TLR4","IL1B","CXCL3","CCL2","CCL5") )

VlnPlot(AurN.combined_fibroblast+myofibroblast1, pt.size = 0, features = c("ACTA2","FAP","POSTN","LUM","COL3A1") )

# Nature 2022 Single-nucleus profiling of human dilated and hypertrophic cardiomyopathy LV
Activated fibroblast
VlnPlot(AurN.combined_fibroblastmyofibroblast1, pt.size = 0, features = c("POSTN","NOX4","COL1A1","COL1A2") )
FB-TLL2
VlnPlot(AurN.combined_fibroblastmyofibroblast1, pt.size = 0, features = c("TLL2","ERI2","ACSM1","THUMPD1") )
Act
VlnPlot(AurN.combined_fibroblastmyofibroblast1, pt.size = 0, features = c("POSTN","TNC","THBS4","MEOX01") )
FB-ZBTB7C
VlnPlot(AurN.combined_fibroblastmyofibroblast1, pt.size = 0, features = c("ZBTB7C","PCOLCE2","SEMA3C","KAZN") )

FB-X1
VlnPlot(AurN.combined_fibroblastmyofibroblast1, pt.size = 0, features = c("TTN","MT-CO2","MT-ATP6","KMT-CO1"))

FB-PTCHD4
VlnPlot(AurN.combined_fibroblastmyofibroblast1, pt.size = 0, features = c("PTCHD4","ZMAT3","DOB2","PAPPA") )


VlnPlot(AurN.combined_fibroblastmyofibroblast1,pt.size = 0, features = c("ABCA9","ABCA6","ABCA10","TIMP1") )

VlnPlot(AurN.combined_fibroblastmyofibroblast1, pt.size = 0, features = c("PDGF","PDRFRA","SPP1","IL1B","IL1R1","AREG"))
VlnPlot(AurN.combined_fibroblastmyofibroblast1, pt.size = 0, features = c("EGFR","CCL2","CCR2","CSF1","CSF1R","IL34"))


# J Adv Res 2022 China AF snRNAseq
FB1 myofibroblast
VlnPlot(AurN.combined_fibroblastmyofibroblast1, pt.size = 0, features = c("ACTG2","MYH11","ITGA8"))
FB2 
VlnPlot(AurN.combined_fibroblastmyofibroblast1, pt.size = 0, features = c("POSTN","ELN"))
FB3 lipid trafficking
VlnPlot(AurN.combined_fibroblastmyofibroblast1, pt.size = 0, features = c("NEGR1","ABCA8","ABCA6") )
FB4 cell adhesion and proliferaion
VlnPlot(AurN.combined_fibroblastmyofibroblast1, pt.size = 0, features = c("KCND2","KCTD16","FGF7","FLRT2") )
FB5 macropharge
VlnPlot(AurN.combined_fibroblastmyofibroblast1, pt.size = 0, features = c("CD163","MRC1") )
FB5 Tcell-like
VlnPlot(AurN.combined_fibroblastmyofibroblast1, pt.size = 0, features = c("SKAPI","ITK") )


VlnPlot(AurN.combined_fibroblastmyofibroblast1, pt.size = 0, features = c("FAP","THBS4","TSHZ2","FAM155A","COL1A2","MIR100HG","NOX4","COL22A1","COL1A1") )
VlnPlot(AurN.combined_fibroblastmyofibroblast1, pt.size = 0, features = c("AEBP1","CHAD","COMP") )
VlnPlot(AurN.combined_fibroblastmyofibroblast1, pt.size = 0.1, features = c("COL3A1","LUM","COMP") )
VlnPlot(AurN.combined_fibroblastmyofibroblast1, pt.size = 0.1, features = c("ACTA2","SMA") )
VlnPlot(AurN.combined_fibroblastmyofibroblast1, pt.size = 0, features = c("PDGF","PDRFRA","SPP1","IL1B","IL1R1","AREG"))
VlnPlot(AurN.combined_fibroblastmyofibroblast1, pt.size = 0, features = c("EGFR","CCL2","CCR2","CSF1","CSF1R","IL34"))


AurN.combined_fibroblastmyofibroblast1_trans <- SCTransform (AurN.combined_fibroblastmyofibroblast1)

cd_genes <- c("ABCA6","ABCA8","CD9","FBLN1","FBLN2","TIMP1","MMP2","COL3A1","LUM","POSTN","CCL2","EGR1","ELN","FBN1","PRG4")
DotPlot(AurN.combined_fibroblastmyofibroblast1_trans,features = cd_genes)+RotatedAxis()+coord_flip()


FeaturePlot(AurN.combined_fibroblastmyofibroblast1, features=c('TTN'), min.cutoff=0.5, max.cutoff='q90')

FeaturePlot(AurN.combined_fibroblastmyofibroblast1, features=c('MYH11'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast1, features=c('ACTA2'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast1, features=c('GSN'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast1, features=c('LUM'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast1, features=c('COL3A1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast1, features=c('POSTN'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast1, features=c('FBLN1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast1, features=c('FBN1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast1, features=c('CCL2'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast1, features=c('ABCA6'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast1, features=c('EGR1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast1, features=c('MYH11'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast1, features=c('PRG4'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast1, features=c('MMP2'), min.cutoff=0.5, max.cutoff='q90')

VlnPlot(AurN.combined_fibroblastmyofibroblast1, pt.size = 0, features = c("IL1R1","CXCR2") )


AurN.combined_fibroblastmyofibroblast1_prop. <- prop.table(table(Idents(AurN.combined_fibroblastmyofibroblast1)))
write.table(AurN.combined_fibroblastmyofibroblast1_prop., "AurN.combined_fibroblastmyofibroblast1_prop.txt", quote=F, col.names=F, append=T)
write.csv(AurN.combined_fibroblastmyofibroblast1_prop., "AurN.combined_fibroblastmyofibroblast1_prop.csv")


AurN.combined_fibroblastmyofibroblast1_prop. <- prop.table(table(Idents(AurN.combined_fibroblastmyofibroblast1),AurN.combined_fibroblast2$AFDuration))
write.table(AurN.combined_fibroblastmyofibroblast1_prop., "AurN.combined_fibroblastmyofibroblast1_prop.AFDuration.txt", quote=F, col.names=F, append=T)
write.csv(AurN.combined_fibroblastmyofibroblast1_prop., "AurN.combined_fibroblastmyofibroblast1_prop.AFDuration.csv")










#最終fibroblast+myofibroblast(myocardium doublet除去)----
AurN.combined_fibroblastmyofibroblast2 <- subset(AurN.combined_fibroblastmyofibroblast1, idents = c(0,1,2,3,4))
AurN.combined_fibroblastmyofibroblast2 <- ScaleData(AurN.combined_fibroblastmyofibroblast2, verbose = FALSE)
AurN.combined_fibroblastmyofibroblast2 <- RunPCA(AurN.combined_fibroblastmyofibroblast2, npcs = 30, verbose = FALSE)
AurN.combined_fibroblastmyofibroblast2 <- RunUMAP(AurN.combined_fibroblastmyofibroblast2 , reduction = "pca", dims = 1:30)
AurN.combined_fibroblastmyofibroblast2 <- FindNeighbors(AurN.combined_fibroblastmyofibroblast2 , reduction = "pca", dims = 1:30)
AurN.combined_fibroblastmyofibroblast2 <- FindClusters(AurN.combined_fibroblastmyofibroblast2, resolution = 0.35)



DimPlot(AurN.combined_fibroblastmyofibroblast2, reduction = "umap", label = FALSE, repel = TRUE)
DimPlot(AurN.combined_fibroblastmyofibroblast2, reduction = "umap", label = TRUE,  repel = TRUE)
DimPlot(AurN.combined_fibroblastmyofibroblast2, reduction = "umap", split.by = "AFDuration")
DimPlot(AurN.combined_fibroblastmyofibroblast2, reduction = "umap", split.by = "stim")


saveRDS(AurN.combined_fibroblastmyofibroblast2, "~/Desktop/Human_LA_LAA/RDS/AurN_cluster_sn_LAA_fibroblastmyofibroblast2.rds")
AurN.combined_fibroblastmyofibroblast2<- readRDS ("~/Desktop/Human_LA_LAA/RDS/AurN_cluster_sn_LAA_fibroblastmyofibroblast2.rds")

FeaturePlot(AurN.combined_fibroblastmyofibroblast2, features=c("COL6A1","COL6A2","COL6A3","ULBP2"), min.cutoff=0.5, max.cutoff='q90')



# Olink----
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("IL17A","KCND2","RASGRF1","B4GALNT3","MTERF3", "FABP6") )
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("RNF4","ARNT","BIRC3","B4GALNT3","MTERF3", "FABP6") )
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("ULBP2") )


MHC2
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DMA","HLA-DMB","HLA-DQA2"))
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("HLA-DOA","HLA-DOB","HLA-DQB1","HLA-DRA","HLA-DRB1","HLA-DRB3"))
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("HLA-DOA","HLA-DOB","HLA-DQB1","HLA-DRA","HLA-DRB4","HLA-DRB5"))
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("COL6A1","COL6A2","COL6A3","HLA-DRA","HLA-DRB4","HLA-DRB5"))

VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0.1, features = c("SFRP1","SFRP2","FRZB","SFRP4","SFRP5"))
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("SFRP1","SFRP2","FRZB","SFRP4","SFRP5"))
FeaturePlot(AurN.combined_fibroblastmyofibroblast2, features=c("SFRP1","SFRP2","FRZB","SFRP4","SFRP5"), min.cutoff=0.2, max.cutoff='q90')


COL6A1 COL6A2 COL6A3

VlnPlot(AurN.combined_fibroblastmyofibroblast2,pt.size = 0, features = c("DCN","GSN","PDGFRA","POSTN","TNC","PDGFRA") )



# コンタミ探し
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("NPPA","MYL7","MYL4","TTN","MYBPC3","TNNT2"))
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("RYR2","PLN","SLC8A1","DCN","GSN","PDGFRA"))
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("VWF","PECAM1","CDH5","RGS5","ABCC9","KCNJ8"))
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("MYH11","TAGLN","ACTA2","PLP1","NRXN1","NRXN3"))
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("GPAM","FASN","LEP","MSLN","WT1","BNC1"))
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("PTPRC","CD163","CD79B","CD14","IL7R","CD3E","CD4","CD8A","NKG7","GNLY","GZMB"))






# ヒートマップの作成------------------------------------------------------------------------------

AurN.combined_fibroblastmyofibroblast2.markers <- FindAllMarkers(AurN.combined_fibroblastmyofibroblast2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
AurN.combined_fibroblastmyofibroblast2.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- AurN.combined_fibroblastmyofibroblast2.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(AurN.combined_fibroblastmyofibroblast2, features = top10$gene) 
write.table(top10, "table.txt", quote=F, col.names=F, append=T)


VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("TIMP1","TIMP2","TIMP3","MMP2","MMP9","IL10"))
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("C1QA","C1QB","C1QC","IL10","CD9","TREM2"))
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("CNN1"))


# aFB Nature 2020 cells of the adult human heart
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("GSN","PCOLCE2","FBLN1","ABCA10","SCN7A","ELL2"))
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("GSN","PCOLCE2","NAMPT","ATP13A3","VWF","PTPRB"))

FB4　TGFB
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("TGFB","POSTN","TNC"))


FB6　ECM
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("OSMR","ILST6"))


# Clin Trans Med 2023
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("PDRFRL","MFAP5","CLIP","CXCL1","ACKR3","TNFAIP6"))
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("PTX3","SERINA3","ITLN1","LGBP3","ID3","RNASE1"))
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("TNF","TLR4","IL1B","CXCL3","CCL2","CCL5") )
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0.1, features = c("TNF","TLR4","IL1B","CXCL3","CCL2","CCL5") )

VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("ACTA2","FAP","POSTN","LUM","COL3A1") )

# Nature 2022 Single-nucleus profiling of human dilated and hypertrophic cardiomyopathy LV
Activated fibroblast
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("POSTN","NOX4","COL1A1","COL1A2") )
FB-TLL2
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("TLL2","ERI2","ACSM1","THUMPD1") )
Act
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("POSTN","TNC","THBS4","MEOX01") )
FB-ZBTB7C
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("ZBTB7C","PCOLCE2","SEMA3C","KAZN") )

FB-X1
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("TTN","MT-CO2","MT-ATP6","KMT-CO1"))

FB-PTCHD4
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("PTCHD4","ZMAT3","DOB2","PAPPA") )


VlnPlot(AurN.combined_fibroblastmyofibroblast2,pt.size = 0, features = c("ABCA9","ABCA6","ABCA10","TIMP1") )

VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("PDGF","PDRFRA","SPP1","IL1B","IL1R1","AREG"))
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("EGFR","CCL2","CCR2","CSF1","CSF1R","IL34"))


# J Adv Res 2022 China AF snRNAseq
FB1 myofibroblast
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("ACTG2","MYH11","ITGA8"))
FB2 
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("POSTN","ELN"))
FB3 lipid trafficking
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("NEGR1","ABCA8","ABCA6") )
FB4 cell adhesion and proliferaion
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("KCND2","KCTD16","FGF7","FLRT2") )
FB5 macropharge
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("CD163","MRC1") )
FB5 Tcell-like
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("SKAPI","ITK") )


VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("FAP","THBS4","TSHZ2","FAM155A","COL1A2","MIR100HG","NOX4","COL22A1","COL1A1") )
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("AEBP1","CHAD","COMP") )
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0.1, features = c("COL3A1","LUM","COMP") )
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0.1, features = c("ACTA2","SMA") )
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("PDGF","PDRFRA","SPP1","IL1B","IL1R1","AREG"))
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("EGFR","CCL2","CCR2","CSF1","CSF1R","IL34"))


VlnPlot(AurN.combined_fibroblastmyofibroblast2,pt.size = 0,features = c("CCL2","CSF1","CSF2","CSF3","IFNG","IL4") )
VlnPlot(AurN.combined_fibroblastmyofibroblast2,pt.size = 0,features = c("CD44","IL1B","AREG","EGFR","CSF1R","IL34") )
VlnPlot(AurN.combined_fibroblastmyofibroblast2,pt.size = 0,features = c("CCR2","IL1B","AREG","EGFR","CSF1R","IL34") )
VlnPlot(AurN.combined_fibroblastmyofibroblast2,pt.size = 0,features = c("PDGFA","PDGRRA","PDGFB","PDGFRB","PDGFC","PDGFRC") )
VlnPlot(AurN.combined_fibroblastmyofibroblast2,pt.size = 0,features = c("PDGFRL","PDAP1","TGFB1","TGFB2","TGFBR1","TGFBR2") )
VlnPlot(AurN.combined_fibroblastmyofibroblast2,pt.size = 0,features = c("TGFB3","TGFB3R","IL17A","IL17RA","IL17B","IL17RB") )
VlnPlot(AurN.combined_fibroblastmyofibroblast2,pt.size = 0,features = c("IL17C","IL17RC","IL17D","IL17RD","IL17E","IL17RE") )
VlnPlot(AurN.combined_fibroblastmyofibroblast2,pt.size = 0,features = c("IL17REL","IL17F","IL17RF","MMP12","MMP13","SLCA12") )
VlnPlot(AurN.combined_fibroblastmyofibroblast2,pt.size = 0,features = c("HFE","CD63","ITGB1","CTRH1","IL6","IL11") )
VlnPlot(AurN.combined_fibroblastmyofibroblast2,pt.size = 0,features = c("BRD4","CX3CR1","MEOX1","CXCL10","CCL19","IL11") )
VlnPlot(AurN.combined_fibroblastmyofibroblast2,pt.size = 0,features = c("IL1R1","IL1R2","IL1RL1","IL1RAP","IL1RAP2","IL11") )
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("TIMP1","TIMP2","TIMP3","MMP2","MMP9","IL10"))
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("C1QA","C1QB","C1QC","IL10","CD9","TREM2"))


#cellchat_VlnPlot
#EGF
VlnPlot(AurN.combined_fibroblastmyofibroblast2,pt.size = 0.1,features = c("EGF","TGFA","AREG","BTC","HBEGF","EREG") )
VlnPlot(AurN.combined_fibroblastmyofibroblast2,pt.size = 1,features = c("EGFR","ERBB2","ERBB4","C1QA","C1QB","C1QC") )


FeaturePlot(AurN.combined_fibroblastmyofibroblast2, features=c("EGFR"), min.cutoff=1, max.cutoff='q90')

#FN1
VlnPlot(AurN.combined_fibroblastmyofibroblast2,pt.size = 0,features = c("FN1","CCL2","CCR2","CSF1","CSF1R","IL34") )
VlnPlot(AurN.combined_fibroblastmyofibroblast2,pt.size = 0,features = c("ITGA2","ITGB1","ITGA4","ITGA5","ITGAV","ITGB7") )
VlnPlot(AurN.combined_fibroblastmyofibroblast2,pt.size = 0,features = c("ITGA2B","ITGB6","ITGBN8","CD44","SDC1","SDC4") )

#SPP1
VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("SPP1","NPPA","MYL7","MYL4","TTN","MYBPC3"))
VlnPlot(AurN.combined_fibroblastmyofibroblast2,pt.size = 0,features = c("CD44","ITGB1","ITGB3","ITGB5","ITGB6") )

FeaturePlot(AurN.combined_fibroblastmyofibroblast2, features=c('SPP1'), min.cutoff=0.01, max.cutoff='q90')

#PDGF
VlnPlot(AurN.combined_fibroblastmyofibroblast2,pt.size = 0,features = c("PDGFA","PDGFB","PDGFC","PDGFD","NPPA","MYL7") )
VlnPlot(AurN.combined_fibroblastmyofibroblast2,pt.size = 0,features = c("PDGFRA","PDGFRB","NPPA","MYL7","MYL4","TTN") )

VlnPlot(AurN.combined_fibroblastmyofibroblast2,pt.size = 0,features = c("STAT3","JAK2","NPPA","MYL7","MYL4","TTN") )


#COLLAGEN receptorのみ
VlnPlot(AurN.combined_fibroblastmyofibroblast2,pt.size = 0,features = c("ITGA1","ITGA3","ITGA9","ITGA10","ITGA11","ITGB1") )
VlnPlot(AurN.combined_fibroblastmyofibroblast2,pt.size = 0,features = c("ITGB8","ITGAV","CD44","SDC1","SDC4") )

#THBS
VlnPlot(AurN.combined_fibroblastmyofibroblast2,pt.size = 0,features = c("THBS1","THBS2","THBS3","THBS4","COMP") )
VlnPlot(AurN.combined_fibroblastmyofibroblast2,pt.size = 0,features = c("ITGA3","ITGAV","ITGB1","ITGB3","SDC1","SDC4") )
VlnPlot(AurN.combined_fibroblastmyofibroblast2,pt.size = 0,features = c("CD36","CD47","ITGB1","ITGB3","SDC1","SDC4") )

#VEGF
VlnPlot(AurN.combined_fibroblastmyofibroblast2,pt.size = 0,features = c("VEGFA","VEGFB","VEGFC","VEGFD","PGF") )
VlnPlot(AurN.combined_fibroblastmyofibroblast2,pt.size = 0,features = c("FLT1","KDR","FLT4","MYL7","MYL4","TTN") )


#IGF
VlnPlot(AurN.combined_fibroblastmyofibroblast2,pt.size = 0,features = c("IGF1","IGF2","IGFL3","MYL7","MYL4","TTN") )
VlnPlot(AurN.combined_fibroblastmyofibroblast2,pt.size = 0,features = c("IGF1R","IGF2R","ITGAV","ITGB3","ITGA6","ITGB4") )
VlnPlot(AurN.combined_fibroblastmyofibroblast2,pt.size = 0,features = c("IGFLR1","IGF2R","ITGAV","ITGB3","ITGA6","ITGB4") )

VlnPlot(AurN.combined_fibroblastmyofibroblast2,pt.size = 0,assay = "RNA",features = c("IGF1R","IGF2R","ITGAV","ITGB3","ITGA6","ITGB4") )


FB
VlnPlot(AurN.combined_fibroblastmyofibroblast2,pt.size = 0,features = c("NKD2","CCL13","CCL19","LOX","LOXL1","LOXL2") )
VlnPlot(AurN.combined_fibroblastmyofibroblast2,pt.size = 0,features = c("DLL1","DLL3","DLL4","JAG1","JAG2","INHBA") )
VlnPlot(AurN.combined_fibroblastmyofibroblast2,pt.size = 0,features = c("CCR2","CCR5","DLL4","JAG1","JAG2","INHBA") )

VlnPlot(AurN.combined_fibroblastmyofibroblast2,pt.size = 0,features = c("ADAMTS6","MFAP5","LINC01133","GFPT6","FSTL1","HAS2","PRG4","CD55") )


FeaturePlot(AurN.combined_fibroblastmyofibroblast2, features=c('INHBA'), min.cutoff=0.5, max.cutoff='q90')


AurN.combined_fibroblastmyofibroblast2_trans <- SCTransform (AurN.combined_fibroblastmyofibroblast2)

cd_genes <- c("ABCA6","ABCA8","CD9","FBLN1","FBLN2","TIMP1","MMP2","COL3A1","LUM","POSTN","CCL2","EGR1","ELN","FBN1","HAS2","PRG4","MFAP5","ADAMTS6","ACTA2","MYH11")
DotPlot(AurN.combined_fibroblastmyofibroblast2_trans,features = cd_genes)+RotatedAxis()+coord_flip()

FeaturePlot(AurN.combined_fibroblastmyofibroblast2, features=c('TTN'), min.cutoff=0.5, max.cutoff='q90')

FeaturePlot(AurN.combined_fibroblastmyofibroblast2, features=c('MYH11'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast2, features=c('ACTA2'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast2, features=c('GSN'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast2, features=c('LUM'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast2, features=c('COL3A1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast2, features=c('POSTN'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast2, features=c('FBLN1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast2, features=c('FBN1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast2, features=c('CCL2'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast2, features=c('ABCA6'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast2, features=c('EGR1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast2, features=c('MYH11'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast2, features=c('PRG4'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_fibroblastmyofibroblast2, features=c('MMP2'), min.cutoff=0.5, max.cutoff='q90')

VlnPlot(AurN.combined_fibroblastmyofibroblast2, pt.size = 0, features = c("IL1R1","CXCR2") )


AurN.combined_fibroblastmyofibroblast2_prop. <- prop.table(table(Idents(AurN.combined_fibroblastmyofibroblast2)))
write.table(AurN.combined_fibroblastmyofibroblast2_prop., "AurN.combined_fibroblastmyofibroblast2_prop.txt", quote=F, col.names=F, append=T)
write.csv(AurN.combined_fibroblastmyofibroblast2_prop., "AurN.combined_fibroblastmyofibroblast2_prop.csv")


AurN.combined_fibroblastmyofibroblast2_prop. <- prop.table(table(Idents(AurN.combined_fibroblastmyofibroblast2),AurN.combined_fibroblastmyofibroblast2$AFDuration))
write.table(AurN.combined_fibroblastmyofibroblast2_prop., "AurN.combined_fibroblastmyofibroblast2_prop.AFDuration.txt", quote=F, col.names=F, append=T)
write.csv(AurN.combined_fibroblastmyofibroblast2_prop., "AurN.combined_fibroblastmyofibroblast2_prop.AFDuration.csv")















#Endothelialcellのみ抽出--------
AurN.combined_endothelialcells <- subset(AurN.combined_1, idents = c(6,7))
AurN.combined_endothelialcells <- ScaleData(AurN.combined_endothelialcells, verbose = FALSE)
AurN.combined_endothelialcells <- RunPCA(AurN.combined_endothelialcells, npcs = 30, verbose = FALSE)
AurN.combined_endothelialcells <- RunUMAP(AurN.combined_endothelialcells , reduction = "pca", dims = 1:30)
AurN.combined_endothelialcells <- FindNeighbors(AurN.combined_endothelialcells , reduction = "pca", dims = 1:30)
AurN.combined_endothelialcells <- FindClusters(AurN.combined_endothelialcells, resolution = 0.2)

DimPlot(AurN.combined_endothelialcells, reduction = "umap", label = FALSE, repel = TRUE)
DimPlot(AurN.combined_endothelialcells, reduction = "umap", label = TRUE,  repel = TRUE)
DimPlot(AurN.combined_endothelialcells, reduction = "umap", split.by = "AFDuration")
DimPlot(AurN.combined_endothelialcells, reduction = "umap", split.by = "stim")
# コンタミ探し
VlnPlot(AurN.combined_endothelialcells, pt.size = 0, features = c("NPPA","MYL7","MYL4","TTN","MYBPC3","TNNT2"))
VlnPlot(AurN.combined_endothelialcells, pt.size = 0, features = c("RYR2","PLN","SLC8A1","DCN","GSN","PDGFRA"))
VlnPlot(AurN.combined_endothelialcells, pt.size = 0, features = c("VWF","PECAM1","CDH5","RGS5","ABCC9","KCNJ8"))
VlnPlot(AurN.combined_endothelialcells, pt.size = 0, features = c("MYH11","TAGLN","ACTA2","PLP1","NRXN1","NRXN3"))
VlnPlot(AurN.combined_endothelialcells, pt.size = 0, features = c("GPAM","FASN","LEP","MSLN","WT1","BNC1"))
VlnPlot(AurN.combined_endothelialcells, pt.size = 0, features = c("PTPRC","CD163","CD79B","CD14","IL7R","CD3E","CD4","CD8A","NKG7","GNLY","GZMB"))

VlnPlot(AurN.combined_endothelialcells, pt.size = 0, features = c("PTPRC","CD163","CD79B","CD14","CD68"))
VlnPlot(AurN.combined_endothelialcells, pt.size = 0, features = c("PTPRC","CD163","CD79B","CD14","CD68"))

# Nature 2020 cells of the adult human heart
pan-EC marker
VlnPlot(AurN.combined_endothelialcells, pt.size = 0, features = c("PECAM1","CDH5","VWF"))
capillary EC
VlnPlot(AurN.combined_endothelialcells, pt.size = 0, features = c("RGCC","CA4"))
immune EC
VlnPlot(AurN.combined_endothelialcells, pt.size = 0, features = c("CXCL1","CCL2","IL6","ICAM1"))
arterial EC
VlnPlot(AurN.combined_endothelialcells, pt.size = 0, features = c("SEMA3G","EFNB2","DLL4","ICAM1"))
venous EC
VlnPlot(AurN.combined_endothelialcells, pt.size = 0, features = c("NR2F2","ACKR1"))
atrial EC
VlnPlot(AurN.combined_endothelialcells, pt.size = 0, features = c("SMOC1","NPR3"))
lymphatic EC
VlnPlot(AurN.combined_endothelialcells, pt.size = 0, features = c("PROX1","TBX1","PDPN"))

# △Nature 2022 Single-nucleus profiling of human dilated and hypertrophic cardiomyopathy LV
EC-PKD1L1
VlnPlot(AurN.combined_endothelialcells, pt.size = 0, features = c("PKD1L1","RGCC","BTNL9") )
EC-NEBL
VlnPlot(AurN.combined_endothelialcells, pt.size = 0, features = c("NEBL","PCSK5","PRDM16") )
EC-LHX6
VlnPlot(AurN.combined_endothelialcells, pt.size = 0, features = c("LHX6","KCNIP4","ZBTB7C") )
EC-TMEM163
VlnPlot(AurN.combined_endothelialcells, pt.size = 0, features = c("TMEM163","KIT","LOXHD1"))
L-EC
VlnPlot(AurN.combined_endothelialcells, pt.size = 0, features = c("MMRN1","RELN","CCL21"))
EC-DKK2
VlnPlot(AurN.combined_endothelialcells, pt.size = 0, features = c("DKK2","NKAIN2","TOX"))
EC-X1
VlnPlot(AurN.combined_endothelialcells, pt.size = 0, features = c("DMD2","PRKG1","NEGR1"))
EC-X2
VlnPlot(AurN.combined_endothelialcells, pt.size = 0, features = c("MYH7","MT-ND3","MT-ND4"))
EC-ADAMTS9
VlnPlot(AurN.combined_endothelialcells, pt.size = 0, features = c("ADAMTS9","PITPNC1"))


# ヒートマップの作成------------------------------------------------------------------------------

AurN.combined_endothelialcells.markers <- FindAllMarkers(AurN.combined_endothelialcells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
AurN.combined_endothelialcells.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- AurN.combined_endothelialcells.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(AurN.combined_endothelialcells, features = top10$gene) 
write.table(top10, "table.txt", quote=F, col.names=F, append=T)




saveRDS(AurN.combined_endothelialcells, "~/Desktop/Human_LA_LAA/RDS/AurN_cluster_sn_LAA_endothelialcells.rds")
AurN.combined_endothelialcells<- readRDS ("~/Desktop/Human_LA_LAA/RDS/AurN_cluster_sn_LAA_endothelialcells.rds")





#Endothelialcellのみ抽出(cardiomyocyte？doublet除去)--------

AurN.combined_endothelialcells1 <- subset(AurN.combined_endothelialcells, idents = c(0,1,2,3))
AurN.combined_endothelialcells1 <- ScaleData(AurN.combined_endothelialcells1, verbose = FALSE)
AurN.combined_endothelialcells1 <- RunPCA(AurN.combined_endothelialcells1, npcs = 30, verbose = FALSE)
AurN.combined_endothelialcells1 <- RunUMAP(AurN.combined_endothelialcells1 , reduction = "pca", dims = 1:30)
AurN.combined_endothelialcells1 <- FindNeighbors(AurN.combined_endothelialcells1 , reduction = "pca", dims = 1:30)
AurN.combined_endothelialcells1 <- FindClusters(AurN.combined_endothelialcells1, resolution = 0.2)

saveRDS(AurN.combined_endothelialcells1, "~/Desktop/Human_LA_LAA/RDS/AurN_cluster2_sn_LAA_endothelialcells1.rds")
AurN.combined_endothelialcells1<- readRDS ("~/Desktop/Human_LA_LAA/RDS/AurN_cluster2_sn_LAA_endothelialcells1.rds")




DimPlot(AurN.combined_endothelialcells1, reduction = "umap", label = FALSE, repel = TRUE)

DimPlot(AurN.combined_endothelialcells1, reduction = "umap", split.by = "AFDuration")
DimPlot(AurN.combined_endothelialcells1, reduction = "umap", split.by = "stim")


# Olink----
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("IL17A","KCND2","RASGRF1","B4GALNT3","MTERF3", "FABP6") )
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("RNF4","ARNT","BIRC3","B4GALNT3","MTERF3", "FABP6") )


# コンタミ探し
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("NPPA","MYL7","MYL4","TTN","MYBPC3","TNNT2"))
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("RYR2","PLN","SLC8A1","DCN","GSN","PDGFRA"))
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("VWF","PECAM1","CDH5","RGS5","ABCC9","KCNJ8"))
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("MYH11","TAGLN","ACTA2","PLP1","NRXN1","NRXN3"))
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("GPAM","FASN","LEP","MSLN","WT1","BNC1"))
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("PTPRC","CD163","CD79B","CD14","IL7R","CD3E","CD4","CD8A","NKG7","GNLY","GZMB"))


VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("ADAMTS1","NNMT","SLCO2A1","TFPI","TFPI2","A2M") )

VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("EGR1","MEG3","CARMN","PCSK7","FOS","ACTA2"))
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("RLP28","RPL41","COL3A1","CXCL2","EGFL7"))
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("TNFRSF1A","TNFRSF1B","MKI67","COL3A1","CXCL2","EGFL7"))
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("SEMA3G","CLIC3","DEPP1","ACKR1","VCAM1","SELE"))
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("ABCC9","AGT","RGS5","PROX1","PDPN","LYVE1"))
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("DCN","CFD","LUM","POSTN","NPR3","CDH11"))
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("RGCC","CA4","B3GNT5","COL4A1","C1orf54","NEURL1B"))
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("TFPI","TFPI2","ADAMTS1","COL4A1","C1orf54","NEURL1B"))
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("THBD","PTGIS","NOS1","PLAT","VWF","NOS2"))
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("IL1RN","MYD88","IL1R1","CXCR2","TLR4","TLR9"))


VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("MKI67"))
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("CASQ2","SYNE2"))



# Clin Trans Med 2023
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("VWF","SEMA3G","CLIC3","DEPP1","ACKR1","VCAM1"))
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("SELE","PDPN","LYVE1","PROX1","RGCC","CA4"))
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("B3GNT5","COL4A1","C1orf54","NEURL1B","POSTN","NPR3") )
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("CDH11","DCN","CFD","LUM","ABCC9","AGT") )
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("CDH11","DCN","CFD","LUM","ABCC9","RGS5") )


# Nature 2020 cells of the adult human heart
pan-EC marker
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("PECAM1","CDH5","VWF"))
capillary EC
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("RGCC","CA4"))
immune EC
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("CXCL1","CCL2","IL6","ICAM1"))
arterial EC
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("SEMA3G","EFNB2","DLL4","ICAM1"))
venous EC
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("NR2F2","ACKR1"))
atrial EC
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("SMOC1","NPR3"))
lymphatic EC
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("PROX1","TBX1","PDPN"))

# △Nature 2022 Single-nucleus profiling of human dilated and hypertrophic cardiomyopathy LV
EC-PKD1L1
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("PKD1L1","RGCC","BTNL9") )
EC-NEBL
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("NEBL","PCSK5","PRDM16") )
EC-LHX6
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("LHX6","KCNIP4","ZBTB7C") )
EC-TMEM163
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("TMEM163","KIT","LOXHD1"))
L-EC
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("MMRN1","RELN","CCL21"))
EC-DKK2
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("DKK2","NKAIN2","TOX"))
EC-X1
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("DMD2","PRKG1","NEGR1"))
EC-X2
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("MYH7","MT-ND3","MT-ND4"))
EC-ADAMTS9
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("ADAMTS9","PITPNC1"))

# J Adv Res 2022 China AF snRNAseq

EC1
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("ACKR1","VCAM1","PECAM1","MPLVAP"))
EC2
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0.1, features = c("CA4","PROM1","KDR"))
EC3
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("GJA5","MET","PLCG2"))
EC4
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("CA4"))


VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("PLVAP","ITGA6"))




VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("TNF","TLR4","IL1B","CXCL3","CCL2","CCL5") )

VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("IL1R1","CXCR2") )

VlnPlot(AurN.combined_endothelialcells1, pt.size = 0.1, features = c("ACTA2","MYH11") )

VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("PDGF","PDRFRA","SPP1","IL1B","IL1R1","AREG"))
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0.1, features = c("EGFR","CCL2","CCR2","CSF1","CSF1R","IL34"))
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0.1, features = c("EGFR","CCL2","CCR2","CSF1","CSF1R","IL34"))


VlnPlot(AurN.combined_endothelialcells1,pt.size = 0,features = c("CCL2","CSF1","CSF2","CSF3","IFNG","IL4") )
VlnPlot(AurN.combined_endothelialcells1,pt.size = 0,features = c("CD44","IL1B","AREG","EGFR","CSF1R","IL34") )
VlnPlot(AurN.combined_endothelialcells1,pt.size = 0,features = c("CCR2","IL1B","AREG","EGFR","CSF1R","IL34") )
VlnPlot(AurN.combined_endothelialcells1,pt.size = 0,features = c("PDGFA","PDGRRA","PDGFB","PDGFRB","PDGFC","PDGFRC") )
VlnPlot(AurN.combined_endothelialcells1,pt.size = 0,features = c("PDGFRL","PDAP1","TGFB1","TGFB2","TGFBR1","TGFBR2") )
VlnPlot(AurN.combined_endothelialcells1,pt.size = 0,features = c("TGFB3","TGFB3R","IL17A","IL17RA","IL17B","IL17RB") )
VlnPlot(AurN.combined_endothelialcells1,pt.size = 0,features = c("IL17C","IL17RC","IL17D","IL17RD","IL17E","IL17RE") )
VlnPlot(AurN.combined_endothelialcells1,pt.size = 0,features = c("IL17REL","IL17F","IL17RF","MMP12","MMP13","SLCA12") )
VlnPlot(AurN.combined_endothelialcells1,pt.size = 0,features = c("HFE","CD63","ITGB1","CTRH1","IL6","IL11") )
VlnPlot(AurN.combined_endothelialcells1,pt.size = 0,features = c("BRD4","CX3CR1","MEOX1","CXCL10","CCL19","IL11") )
VlnPlot(AurN.combined_endothelialcells1,pt.size = 0,features = c("IL1R1","IL1R2","IL1RL1","IL1RAP","IL1RAP2","IL11") )
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("TIMP1","TIMP2","TIMP3","MMP2","MMP9","IL10"))
VlnPlot(AurN.combined_endothelialcells1, pt.size = 0, features = c("C1QA","C1QB","C1QC","IL10","CD9","TREM2"))


#EGF
VlnPlot(AurN.combined_endothelialcells1,pt.size = 0,features = c("EGF","TGFA","AREG","BTC","HBEGF","EREG") )
VlnPlot(AurN.combined_endothelialcells1,pt.size = 0,features = c("EGFR","ERBB2","ERBB4","MYL7","MYL4","TTN") )

#FN1
VlnPlot(AurN.combined_endothelialcells1,pt.size = 0,features = c("FN1","CCL2","CCR2","CSF1","CSF1R","IL34") )
VlnPlot(AurN.combined_endothelialcells1,pt.size = 0,features = c("ITGA2","ITGB1","ITGA4","ITGA5","ITGAV","ITGB7") )
VlnPlot(AurN.combined_endothelialcells1,pt.size = 0,features = c("ITGA2B","ITGB6","ITGBN8","CD44","SDC1","SDC4") )


#BMP
VlnPlot(AurN.combined_endothelialcells1,pt.size = 0,features = c("BMP2","BMP4","GDF5","GDF6","GDF7","BMP15") )
VlnPlot(AurN.combined_endothelialcells1,pt.size = 0,features = c("BMP5","BMP6","BMP7","BMP8A","BMP8B","BMP10") )
VlnPlot(AurN.combined_endothelialcells1,pt.size = 0,features = c("BMPR1A","BMPR1B","ACVR2A","ACVR2B","BMPR2","ACVR1") )

#THBS
VlnPlot(AurN.combined_endothelialcells1,pt.size = 0,features = c("THBS1","THBS2","THBS3","THBS4","COMP") )
VlnPlot(AurN.combined_endothelialcells1,pt.size = 0,features = c("ITGA3","ITGAV","ITGB1","ITGB3","SDC1","SDC4") )
VlnPlot(AurN.combined_endothelialcells1,pt.size = 0,features = c("CD36","CD47","ITGB1","ITGB3","SDC1","SDC4") )

#VEGF
VlnPlot(AurN.combined_endothelialcells1,pt.size = 0,features = c("VEGFA","VEGFB","VEGFC","VEGFD","PGF") )
VlnPlot(AurN.combined_endothelialcells1,pt.size = 0,features = c("FLT1","KDR","FLT4","MYL7","MYL4","TTN") )
#IGF
VlnPlot(AurN.combined_endothelialcells1,pt.size = 0,features = c("IGF1","IGF2","IGFL3","MYL7","MYL4","TTN") )
VlnPlot(AurN.combined_endothelialcells1,pt.size = 0,features = c("IGF1R","IGF2R","ITGAV","ITGB3","ITGA6","ITGB4") )
VlnPlot(AurN.combined_endothelialcells1,pt.size = 0,features = c("IGFLR1","IGF2R","ITGAV","ITGB3","ITGA6","ITGB4") )

VlnPlot(AurN.combined_endothelialcells1,pt.size = 0,assay = "RNA",features = c("IGF1R","IGF2R","ITGAV","ITGB3","ITGA6","ITGB4") )


VlnPlot(AurN.combined_endothelialcells1,pt.size = 0,features = c("NKD2","CCL13","CCL19","LOX","LOXL1","LOXL2") )
VlnPlot(AurN.combined_endothelialcells1,pt.size = 0,features = c("DLL1","DLL3","DLL4","JAG1","JAG2","INHBA") )
VlnPlot(AurN.combined_endothelialcells1,pt.size = 0,features = c("CCR2","CCR5","DLL4","JAG1","JAG2","INHBA") )

VlnPlot(AurN.combined_endothelialcells1,pt.size = 0,features = c("IL6R","IL1B","TNFRSF1A","TNFRSF1B","JAG2","INHBA") )


# ヒートマップの作成------------------------------------------------------------------------------

AurN.combined_endothelialcells1.markers <- FindAllMarkers(AurN.combined_endothelialcells1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
AurN.combined_endothelialcells1.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- AurN.combined_endothelialcells1.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(AurN.combined_endothelialcells1, features = top10$gene) 
write.table(top10, "table.txt", quote=F, col.names=F, append=T)


AurN.combined_endothelialcells1_trans <- SCTransform (AurN.combined_endothelialcells1)

cd_genes <- c("RGCC","CA4","NOTCH4","EGFL7","ITGA6","KDR","MKI67","NR2F2","ACKR1","VCAM1","CXCL2","TLR4","EGR1","SMOC1","NPR3","POSTN","COL3A1","SEMA3G","EFNB2","MET","ADAMTS1","THBD")
DotPlot(AurN.combined_endothelialcells1_trans,features = cd_genes)+RotatedAxis()+coord_flip()

AurN.combined_endothelialcells1_trans <- SCTransform (AurN.combined_endothelialcells1)

cd_genes <- c("RGCC","CA4","NOTCH4","EGFL7","ITGA6","NR2F2","ACKR1","VCAM1","CXCL2","TLR4","EGR1","SMOC1","NPR3","POSTN","COL3A1","IL1R1","SEMA3G","EFNB2","ADAMTS1","THBD","JAG1")
DotPlot(AurN.combined_endothelialcells1_trans,features = cd_genes)+RotatedAxis()+coord_flip()


AurN.combined_endothelialcells1_prop. <- prop.table(table(Idents(AurN.combined_endothelialcells1)))
write.table(AurN.combined_endothelialcells1_prop., "AurN.combined_endothelialcells1_prop.txt", quote=F, col.names=F, append=T)
write.csv(AurN.combined_endothelialcells1_prop., "AurN.combined_endothelialcells1_prop.csv")


AurN.combined_endothelialcells1_prop. <- prop.table(table(Idents(AurN.combined_endothelialcells1),AurN.combined_endothelialcells1$AFDuration))
write.table(AurN.combined_endothelialcells1_prop., "AurN.combined_endothelialcells1_prop.AFDuration.txt", quote=F, col.names=F, append=T)
write.csv(AurN.combined_endothelialcells1_prop., "AurN.combined_endothelialcells1_prop.AFDuration.csv")







#fibroblast(including myofibroblast)+endothelial----
AurN.combined_FBEC <- subset(AurN.combined_1, idents = c(0,6,7,9))
AurN.combined_FBEC <- ScaleData(AurN.combined_FBEC, verbose = FALSE)
AurN.combined_FBEC <- RunPCA(AurN.combined_FBEC, npcs = 30, verbose = FALSE)
AurN.combined_FBEC <- RunUMAP(AurN.combined_FBEC , reduction = "pca", dims = 1:30)
AurN.combined_FBEC <- FindNeighbors(AurN.combined_FBEC , reduction = "pca", dims = 1:30)
AurN.combined_FBEC <- FindClusters(AurN.combined_FBEC, resolution = 0.35)



DimPlot(AurN.combined_FBEC, reduction = "umap", label = FALSE, repel = TRUE)
DimPlot(AurN.combined_FBEC, reduction = "umap", label = TRUE,  repel = TRUE)
DimPlot(AurN.combined_FBEC, reduction = "umap", split.by = "AFDuration")
DimPlot(AurN.combined_FBEC, reduction = "umap", split.by = "stim")


# 並べ替える
levels(AurN.combined_FBEC) # [1] "0" "1" "2" "3" "4" "5" "6" 
levels(AurN.combined_FBEC) <- c("0", "1" ,"4", "5", "2","3","6")

new.cluster.ids <- c("0","1","2","3","4","5","6")
names(new.cluster.ids) <- levels(AurN.combined_FBEC)
AurN.combined_FBEC<- RenameIdents(AurN.combined_FBEC, new.cluster.ids)


saveRDS(AurN.combined_FBEC, "~/Desktop/Human_LA_LAA/RDS/AurN_cluster_sn_LAA_FBEC.rds")
AurN.combined_FBEC<- readRDS ("~/Desktop/Human_LA_LAA/RDS/AurN_cluster_sn_LAA_FBEC.rds")




VlnPlot(AurN.combined_FBEC,pt.size = 0, features = c("DCN","GSN","PDGFRA","POSTN","TNC","PDGFRA") )



# コンタミ探し
VlnPlot(AurN.combined_FBEC, pt.size = 0, features = c("NPPA","MYL7","MYL4","TTN","MYBPC3","TNNT2"))
VlnPlot(AurN.combined_FBEC, pt.size = 0, features = c("RYR2","PLN","SLC8A1","DCN","GSN","PDGFRA"))
VlnPlot(AurN.combined_FBEC, pt.size = 0, features = c("VWF","PECAM1","CDH5","RGS5","ABCC9","KCNJ8"))
VlnPlot(AurN.combined_FBEC, pt.size = 0, features = c("MYH11","TAGLN","ACTA2","PLP1","NRXN1","NRXN3"))
VlnPlot(AurN.combined_FBEC, pt.size = 0, features = c("GPAM","FASN","LEP","MSLN","WT1","BNC1"))
VlnPlot(AurN.combined_FBEC, pt.size = 0, features = c("PTPRC","CD163","CD79B","CD14","IL7R","CD3E","CD4","CD8A","NKG7","GNLY","GZMB"))


FeaturePlot(AurN.combined_FBEC, features=c('NPPA'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_FBEC, features=c('TTN'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_FBEC, features=c('POSTN'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_FBEC, features=c('COL3A1'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_FBEC, features=c('LUM'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_FBEC, features=c('ACTA2'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_FBEC, features=c('MYH11'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_FBEC, features=c('TAGLIN'), min.cutoff=0.5, max.cutoff='q90')
FeaturePlot(AurN.combined_FBEC, features=c('DCN'), min.cutoff=0.5, max.cutoff='q90')



# ヒートマップの作成------------------------------------------------------------------------------

AurN.combined_FBEC.markers <- FindAllMarkers(AurN.combined_FBEC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
AurN.combined_FBEC.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- AurN.combined_FBEC.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(AurN.combined_FBEC, features = top10$gene) 
write.table(top10, "table.txt", quote=F, col.names=F, append=T)


VlnPlot(AurN.combined_FBEC, pt.size = 0, features = c("TIMP1","TIMP2","TIMP3","MMP2","MMP9","IL10"))
VlnPlot(AurN.combined_FBEC, pt.size = 0, features = c("C1QA","C1QB","C1QC","IL10","CD9","TREM2"))
VlnPlot(AurN.combined_FBEC, pt.size = 0, features = c("CNN1"))


# aFB Nature 2020 cells of the adult human heart
VlnPlot(AurN.combined_FBEC, pt.size = 0, features = c("GSN","PCOLCE2","FBLN1","ABCA10","SCN7A","ELL2"))
VlnPlot(AurN.combined_FBEC, pt.size = 0, features = c("GSN","PCOLCE2","NAMPT","ATP13A3","VWF","PTPRB"))

FB4　TGFB
VlnPlot(AurN.combined_FBEC, pt.size = 0, features = c("TGFB","POSTN","TNC"))


FB6　ECM
VlnPlot(AurN.combined_FBEC, pt.size = 0, features = c("OSMR","ILST6"))


# Clin Trans Med 2023
VlnPlot(AurN.combined_FBEC, pt.size = 0, features = c("PDRFRL","MFAP5","CLIP","CXCL1","ACKR3","TNFAIP6"))
VlnPlot(AurN.combined_FBEC, pt.size = 0, features = c("PTX3","SERINA3","ITLN1","LGBP3","ID3","RNASE1"))
VlnPlot(AurN.combined_FBEC, pt.size = 0, features = c("TNF","TLR4","IL1B","CXCL3","CCL2","CCL5") )
VlnPlot(AurN.combined_FBEC, pt.size = 0.1, features = c("TNF","TLR4","IL1B","CXCL3","CCL2","CCL5") )

VlnPlot(AurN.combined_FBEC, pt.size = 0, features = c("ACTA2","FAP","POSTN","LUM","COL3A1") )

# Nature 2022 Single-nucleus profiling of human dilated and hypertrophic cardiomyopathy LV
Activated fibroblast
VlnPlot(AurN.combined_FBEC, pt.size = 0, features = c("POSTN","NOX4","COL1A1","COL1A2") )
FB-TLL2
VlnPlot(AurN.combined_FBEC, pt.size = 0, features = c("TLL2","ERI2","ACSM1","THUMPD1") )
Act
VlnPlot(AurN.combined_FBEC, pt.size = 0, features = c("POSTN","TNC","THBS4","MEOX01") )
FB-ZBTB7C
VlnPlot(AurN.combined_FBEC, pt.size = 0, features = c("ZBTB7C","PCOLCE2","SEMA3C","KAZN") )

FB-X1
VlnPlot(AurN.combined_FBEC, pt.size = 0, features = c("TTN","MT-CO2","MT-ATP6","KMT-CO1"))

FB-PTCHD4
VlnPlot(AurN.combined_FBEC, pt.size = 0, features = c("PTCHD4","ZMAT3","DOB2","PAPPA") )


VlnPlot(AurN.combined_FBEC,pt.size = 0, features = c("ABCA9","ABCA6","ABCA10","TIMP1") )

VlnPlot(AurN.combined_FBEC, pt.size = 0, features = c("PDGF","PDRFRA","SPP1","IL1B","IL1R1","AREG"))
VlnPlot(AurN.combined_FBEC, pt.size = 0, features = c("EGFR","CCL2","CCR2","CSF1","CSF1R","IL34"))


# J Adv Res 2022 China AF snRNAseq
FB1 myofibroblast
VlnPlot(AurN.combined_FBEC, pt.size = 0, features = c("ACTG2","MYH11","ITGA8"))
FB2 
VlnPlot(AurN.combined_FBEC, pt.size = 0, features = c("POSTN","ELN"))
FB3 lipid trafficking
VlnPlot(AurN.combined_FBEC, pt.size = 0, features = c("NEGR1","ABCA8","ABCA6") )
FB4 cell adhesion and proliferaion
VlnPlot(AurN.combined_FBEC, pt.size = 0, features = c("KCND2","KCTD16","FGF7","FLRT2") )
FB5 macropharge
VlnPlot(AurN.combined_FBEC, pt.size = 0, features = c("CD163","MRC1") )
FB5 Tcell-like
VlnPlot(AurN.combined_FBEC, pt.size = 0, features = c("SKAPI","ITK") )


VlnPlot(AurN.combined_FBEC, pt.size = 0, features = c("FAP","THBS4","TSHZ2","FAM155A","COL1A2","MIR100HG","NOX4","COL22A1","COL1A1") )
VlnPlot(AurN.combined_FBEC, pt.size = 0, features = c("AEBP1","CHAD","COMP") )
VlnPlot(AurN.combined_FBEC, pt.size = 0.1, features = c("COL3A1","LUM","COMP") )
VlnPlot(AurN.combined_FBEC, pt.size = 0.1, features = c("ACTA2","SMA") )
VlnPlot(AurN.combined_FBEC, pt.size = 0, features = c("PDGF","PDRFRA","SPP1","IL1B","IL1R1","AREG"))
VlnPlot(AurN.combined_FBEC, pt.size = 0, features = c("EGFR","CCL2","CCR2","CSF1","CSF1R","IL34"))


VlnPlot(AurN.combined_FBEC,pt.size = 0,features = c("CCL2","CSF1","CSF2","CSF3","IFNG","IL4") )
VlnPlot(AurN.combined_FBEC,pt.size = 0,features = c("CD44","IL1B","AREG","EGFR","CSF1R","IL34") )
VlnPlot(AurN.combined_FBEC,pt.size = 0,features = c("CCR2","IL1B","AREG","EGFR","CSF1R","IL34") )
VlnPlot(AurN.combined_FBEC,pt.size = 0,features = c("PDGFA","PDGRRA","PDGFB","PDGFRB","PDGFC","PDGFRC") )
VlnPlot(AurN.combined_FBEC,pt.size = 0,features = c("PDGFRL","PDAP1","TGFB1","TGFB2","TGFBR1","TGFBR2") )
VlnPlot(AurN.combined_FBEC,pt.size = 0,features = c("TGFB3","TGFB3R","IL17A","IL17RA","IL17B","IL17RB") )
VlnPlot(AurN.combined_FBEC,pt.size = 0,features = c("IL17C","IL17RC","IL17D","IL17RD","IL17E","IL17RE") )
VlnPlot(AurN.combined_FBEC,pt.size = 0,features = c("IL17REL","IL17F","IL17RF","MMP12","MMP13","SLCA12") )
VlnPlot(AurN.combined_FBEC,pt.size = 0,features = c("HFE","CD63","ITGB1","CTRH1","IL6","IL11") )
VlnPlot(AurN.combined_FBEC,pt.size = 0,features = c("BRD4","CX3CR1","MEOX1","CXCL10","CCL19","IL11") )
VlnPlot(AurN.combined_FBEC,pt.size = 0,features = c("IL1R1","IL1R2","IL1RL1","IL1RAP","IL1RAP2","IL11") )
VlnPlot(AurN.combined_FBEC, pt.size = 0, features = c("TIMP1","TIMP2","TIMP3","MMP2","MMP9","IL10"))
VlnPlot(AurN.combined_FBEC2, pt.size = 0, features = c("C1QA","C1QB","C1QC","IL10","CD9","TREM2"))





AurN.combined_FBEC_trans <- SCTransform (AurN.combined_FBEC)


cd_genes <- c("DCN","COL1A1","PDGFRA","VWF","PECAM1","MYH11","ACTA2","ABCA6","ABCA8","CD9","FBLN1","FBLN2","TIMP1","MMP2","COL3A1","LUM","POSTN","CCL2","ELN","FBN1","PRG4","RGCC","CA4","NOTCH4","EGFL7","ITGA6","KDR","MKI67","NR2F2","ACKR1","VCAM1","CXCL2","TLR4","EGR1","SMOC1","NPR3","SEMA3G","EFNB2","MET","ADAMTS1","THBD")
DotPlot(AurN.combined_FBEC_trans,features = cd_genes)+RotatedAxis()+coord_flip()












#cardiomyocyteのみ抽出--------
AurN.combined_cardiomyocyte <- subset(AurN.combined_1, idents = c(1,2))
AurN.combined_cardiomyocyte <- ScaleData(AurN.combined_cardiomyocyte, verbose = FALSE)
AurN.combined_cardiomyocyte <- RunPCA(AurN.combined_cardiomyocyte, npcs = 30, verbose = FALSE)
AurN.combined_cardiomyocyte <- RunUMAP(AurN.combined_cardiomyocyte , reduction = "pca", dims = 1:30)
AurN.combined_cardiomyocyte <- FindNeighbors(AurN.combined_cardiomyocyte , reduction = "pca", dims = 1:30)
AurN.combined_cardiomyocyte <- FindClusters(AurN.combined_cardiomyocyte, resolution = 0.2)

DimPlot(AurN.combined_cardiomyocyte, reduction = "umap", label = FALSE, repel = TRUE)
DimPlot(AurN.combined_cardiomyocyte, reduction = "umap", split.by = "AFDuration")
DimPlot(AurN.combined_cardiomyocyte, reduction = "umap", split.by = "stim")

saveRDS(AurN.combined_cardiomyocyte, "~/Desktop/Human_LA_LAA/RDS/AurN_cluster_sn_LAA_cardiomyocyte.rds")
AurN.combined_cardiomyocyte<- readRDS ("~/Desktop/Human_LA_LAA/RDS/AurN_cluster_sn_LAA_cardiomyocyte.rds")

# ヒートマップの作成------------------------------------------------------------------------------

AurN.combined_cardiomyocyte.markers <- FindAllMarkers(AurN.combined_cardiomyocyte, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
AurN.combined_cardiomyocyte.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- AurN.combined_cardiomyocyte.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(AurN.combined_cardiomyocyte, features = top10$gene) 
write.table(top10, "table.txt", quote=F, col.names=F, append=T)



# Olink----
VlnPlot(AurN.combined_cardiomyocyte, pt.size = 0, features = c("IL17A","KCND2","RASGRF1","B4GALNT3","MTERF3", "FABP6") )
VlnPlot(AurN.combined_cardiomyocyte, pt.size = 0, features = c("RNF4","ARNT","BIRC3","B4GALNT3","MTERF3", "FABP6") )




# コンタミ探し
VlnPlot(AurN.combined_cardiomyocyte, pt.size = 0, features = c("NPPA","MYL7","MYL4","TTN","TNNT2","ACTC1"))
VlnPlot(AurN.combined_cardiomyocyte, pt.size = 0, features = c("RYR2","PLN","SLC8A1","DCN","GSN","PDGFRA"))
VlnPlot(AurN.combined_cardiomyocyte, pt.size = 0, features = c("DCN","GSN","PDGFRA","RGS5","ABCC9","KCNJ8"))
VlnPlot(AurN.combined_cardiomyocyte, pt.size = 0, features = c("MYH11","TAGLN","ACTA2","PLP1","NRXN1","NRXN3"))
VlnPlot(AurN.combined_cardiomyocyte, pt.size = 0, features = c("GPAM","FASN","LEP","MSLN","WT1","BNC1"))
VlnPlot(AurN.combined_cardiomyocyte, pt.size = 0, features = c("PTPRC","CD163","CD79B","CD14","IL7R","CD3E","CD4","CD8A","NKG7","GNLY","GZMB"))

VlnPlot(AurN.combined_cardiomyocyte, pt.size = 0, features = c("PTPRC","CD163","CD79B","CD14","CD68"))
VlnPlot(AurN.combined_cardiomyocyte, pt.size = 0, features = c("PTPRC","CD163","CD79B","CD14","CD68","VCAN","TIMP1"))

# Nature 2020 cells of the adult human heart
aCM1
VlnPlot(AurN.combined_cardiomyocyte, pt.size = 0, features = c("MYH6","NPPA","MYL4","ADGRL2","NFXL1","ROBO2"))

aCM2
VlnPlot(AurN.combined_cardiomyocyte, pt.size = 0, features = c("HAMP","SLIT3","ALDH1A2"))
aCM3
VlnPlot(AurN.combined_cardiomyocyte, pt.size = 0, features = c("CNN1","MYH9","DUSP27"))
aCM4
VlnPlot(AurN.combined_cardiomyocyte, pt.size = 0, features = c("CKM","COX41L","NDUFA4"))
aCM5
VlnPlot(AurN.combined_cardiomyocyte, pt.size = 0, features = c("DLC1","PLA2GS","MAML2"))
VlnPlot(AurN.combined_cardiomyocyte, pt.size = 0, features = c("SOX5","EBF1","KCNAB1"))

FeaturePlot (AurN.combined_cardiomyocyte, features = c("SORBS2"), min.cutoff=0.00005, max.cutoff='q90')
FeaturePlot(AurN.combined_cardiomyocyte, features=c('SORBS2'), min.cutoff=0.5, max.cutoff='q90', split.by = "AFDuration")
FeaturePlot(AurN.combined_cardiomyocyte, features=c('GJA1'), min.cutoff=0.5, max.cutoff='q90', split.by = "AFDuration")


VlnPlot(AurN.combined_cardiomyocyte, pt.size = 0, features = c("TNFAIP2"))
VlnPlot(AurN.combined_cardiomyocyte, pt.size = 0, features = c("ALDH1A2","ROR2","SYNPR","CASQ2"))
VlnPlot(AurN.combined_cardiomyocyte, pt.size = 0.1, features = c("SORBS2") , split.by = "AFDuration" , cols = c("EarlyAF" ="deepskyblue", "PermanentAF"="red"))


# J Adv Res 2022 China AF snRNAseq
CM1
VlnPlot(AurN.combined_cardiomyocyte, pt.size = 0, features = c("CACNA1D","HCN1"))
CM2
VlnPlot(AurN.combined_cardiomyocyte, pt.size = 0, features = c("COL1A2","FBN1","MECOM","SORBS2"))
CM3
VlnPlot(AurN.combined_cardiomyocyte, pt.size = 0, features = c("MYO18B","NPPA","NPPB","FHL2"))

VlnPlot(AurN.combined_cardiomyocyte, pt.size = 0, features = c("ACTA1","PDE3A","XIST","CEBPG","GLIS1"))
VlnPlot(AurN.combined_cardiomyocyte, pt.size = 0, features = c("MYH6","RYR2","TTN","MYH7"))
VlnPlot(AurN.combined_cardiomyocyte, pt.size = 0, features = c("GJA1"))


AurN.combined_cardiomyocyte_trans <- SCTransform (AurN.combined_cardiomyocyte)

cd_genes <- c("HIF3A","NPPA","MYH6","SORBS2","GJA1")
DotPlot(AurN.combined_cardiomyocyte_trans,features = cd_genes)+RotatedAxis()+coord_flip()

FeaturePlot(AurN.combined_cardiomyocyte, features=c('SORBS2'), min.cutoff=0.5, max.cutoff='q90')

FeaturePlot(AurN.combined_cardiomyocyte, features=c('SORBS2'), min.cutoff=1.5, max.cutoff='q90', split.by = "AFDuration")
VlnPlot(AurN.combined_cardiomyocyte, pt.size = 0, features = c("SORBS2"), split.by = "AFDuration")

#色変更
library(scales)
show_col(hue_pal()(3))

cols = c("0" = "#F8766D","1" = "#00BA38","2" = "#619CFF")

cols = c("0" = "#F8766D","1" = "#E58700","2" = "#C99800","3" = "#A3A500","4" = "#6BB100","5" = "#00BA38","6" = "#00BF7D","7" = "#00BF7D")




#IGF
VlnPlot(AurN.combined_cardiomyocyte,pt.size = 0,features = c("IGF1","IGF2","IGFL3","MYL7","MYL4","TTN") )
VlnPlot(AurN.combined_cardiomyocyte,pt.size = 0,features = c("IGF1R","IGF2R","ITGAV","ITGB3","ITGA6","ITGB4") )
VlnPlot(AurN.combined_cardiomyocyte,pt.size = 0,features = c("IGFLR1","IGF2R","ITGAV","ITGB3","ITGA6","ITGB4") )

VlnPlot(AurN.combined_cardiomyocyte,pt.size = 0,assay = "RNA",features = c("IGF1R","IGF2R","ITGAV","ITGB3","ITGA6","ITGB4") )

FeaturePlot(AurN.combined_cardiomyocyte, features=c("IGF1R","IGF2R"), min.cutoff=0.5, max.cutoff='q90')


AurN.combined_cardiomyocyte_prop. <- prop.table(table(Idents(AurN.combined_cardiomyocyte)))
write.table(AurN.combined_cardiomyocyte_prop., "AurN.combined_cardiomyocyte_prop.txt", quote=F, col.names=F, append=T)
write.csv(AurN.combined_cardiomyocyte_prop., "AurN.combined_cardiomyocyte_prop.csv")


AurN.combined_cardiomyocyte_prop. <- prop.table(table(Idents(AurN.combined_cardiomyocyte),AurN.combined_cardiomyocyte$AFDuration))
write.table(AurN.combined_cardiomyocyte_prop., "AurN.combined_cardiomyocyte_prop.AFDuration.txt", quote=F, col.names=F, append=T)
write.csv(AurN.combined_cardiomyocyte_prop., "AurN.combined_cardiomyocyte_prop.AFDuration.csv")



#CD45----
AurN.combined_CD45 <- subset(AurN.combined_1, idents = c(2,3,4))
AurN.combined_CD45 <- ScaleData(AurN.combined_CD45, verbose = FALSE)
AurN.combined_CD45<- RunPCA(AurN.combined_CD45, npcs = 30, verbose = FALSE)
AurN.combined_CD45 <- RunUMAP(AurN.combined_CD45 , reduction = "pca", dims = 1:30)
AurN.combined_CD45 <- FindNeighbors(AurN.combined_CD45 , reduction = "pca", dims = 1:30)
AurN.combined_CD45 <- FindClusters(AurN.combined_CD45, resolution = 0.5)



# Visualization 別々------
DimPlot(AurN.combined_CD45, reduction = "umap", label = TRUE, repel = FALSE)
DimPlot(AurN.combined_CD45, reduction = "umap", split.by = "stim")
DimPlot(AurN.combined_CD45, reduction = "umap", split.by = "AFDuration")

# RDS
saveRDS(AurN.combined_CD45, "~/Desktop/Human_LA_LAA/RDS/AurN_cluster_sn_LAA_CD45.rds")
AurN.combined_CD45 <- readRDS ("~/Desktop/Human_LA_LAA/RDS/AurN_cluster_sn_LAA_CD45.rds")

# コンタミ探し
VlnPlot(AurN.combined_CD45, pt.size = 0, features = c("NPPA","MYL7","MYL4","TTN","TNNT2","ACTC1"))
VlnPlot(AurN.combined_CD45, pt.size = 0, features = c("RYR2","PLN","SLC8A1","DCN","GSN","PDGFRA"))
VlnPlot(AurN.combined_CD45, pt.size = 0, features = c("DCN","GSN","PDGFRA","RGS5","ABCC9","KCNJ8"))
VlnPlot(AurN.combined_CD45, pt.size = 0, features = c("MYH11","TAGLN","ACTA2","PLP1","NRXN1","NRXN3"))
VlnPlot(AurN.combined_CD45, pt.size = 0, features = c("GPAM","FASN","LEP","MSLN","WT1","BNC1"))
VlnPlot(AurN.combined_CD45, pt.size = 0, features = c("PTPRC","CD163","CD79B","CD14","IL7R","CD3E","CD4","CD8A","NKG7","GNLY","GZMB"))

VlnPlot(AurN.combined_CD45, pt.size = 0, features = c("PTPRC","CD163","CD79B","CD14","CD68"))
VlnPlot(AurN.combined_CD45, pt.size = 0, features = c("PTPRC","CD163","CD79B","CD14","CD68","VCAN","TIMP1"))


DimPlot(AurN.combined_CD45, reduction = "umap", split.by = "stim")
DimPlot(AurN.combined_CD45, reduction = "umap", label = TRUE, split.by = "AFControl")
DimPlot(AurN.combined_CD45, reduction = "umap", label = TRUE,  repel = TRUE)


# Myeloid----
VlnPlot(AurN.combined_CD45, pt.size = 0, features = c("CD68","CD14","CSF1R","IL1B","S100A8", "S100A9") )
VlnPlot(AurN.combined_CD45, pt.size = 0, features = c("TREM2","CD9","APOE","CXCL3","C1QA","TNF","CXCL8") )
VlnPlot(AurN.combined_CD45, pt.size = 0, features = c("HLA-DPA1","HLA-DPB1","HLA-DQA1","CLEC10A","FCER1A","CD1C") )
VlnPlot(AurN.combined_CD45, pt.size = 0, features = c("PTPRC","HLA-DPB1","HLA-DQA1","CLEC10A","FCER1A","CD1C") )

# Bcells-----
VlnPlot(AurN.combined_CD45, pt.size = 0, features = c("CD79A","CD79B","FCER2","CD22","MYH11","CD34") )
VlnPlot(AurN.combined_CD45, pt.size = 0, features = c("LUM","TAGLN","CD34","CDH5","PECAM1","TYROBP") )

# Tcells-----
VlnPlot(AurN.combined_CD45, pt.size = 0, features = c("CD3E","CD4","CD8A","NKG7","GNLY","GZMB","IL7R") )
VlnPlot(AurN.combined_CD45, pt.size = 0, features = c("KLRB1","CCR6","CD1d1","MKI67","IL2RB","GZMB") )

VlnPlot(AurN.combined_CD45, pt.size = 0, features = c("TIMP1"), split.by = "stim" , cols = c("sa" ="deepskyblue", "Aur1Epi"="red"))

VlnPlot(AurN.combined_CD45, pt.size = 0, features = c("TIMP1") )
FeaturePlot (AurN.combined_CD45, features = c("TIMP1"), min.cutoff=0.00005, max.cutoff='q90')

AurN.combined_CD45_trans <- SCTransform (AurN.combined_CD45)

cd_genes <- c("CD3E","NKG7","CD8A","CD68","CD14","S100A8","HLA-DQA1","CD79A","CD79B","GNLY","TYROBP","KIT", "MYH11","PECAM1","CDH5","ACTC1","TNNC1","MKI67")
DotPlot(AurN.combined_CD45_trans,features = cd_genes)+RotatedAxis()+coord_flip()


AurN.combined_CD45.markers <- FindAllMarkers(AurN.combined_CD45, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
AurN.combined_CD45.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- AurN.combined_CD45.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(AurN.combined_CD45, features = top10$gene) 
write.table(top10, "table.txt", quote=F, col.names=F, append=T)




#CD45(myocardium除去)----
AurN.combined_CD45_1 <- subset(AurN.combined_CD45, idents = c(0,1,2,4))
AurN.combined_CD45_1 <- ScaleData(AurN.combined_CD45_1, verbose = FALSE)
AurN.combined_CD45_1<- RunPCA(AurN.combined_CD45_1, npcs = 30, verbose = FALSE)
AurN.combined_CD45_1 <- RunUMAP(AurN.combined_CD45_1 , reduction = "pca", dims = 1:30)
AurN.combined_CD45_1 <- FindNeighbors(AurN.combined_CD45_1 , reduction = "pca", dims = 1:30)
AurN.combined_CD45_1 <- FindClusters(AurN.combined_CD45_1, resolution = 0.5)



# Visualization 別々------
DimPlot(AurN.combined_CD45_1, reduction = "umap", label = FALSE, repel = FALSE)
DimPlot(AurN.combined_CD45_1, reduction = "umap", split.by = "stim")
DimPlot(AurN.combined_CD45_1, reduction = "umap", split.by = "AFDuration")

# RDS
saveRDS(AurN.combined_CD45_1, "~/Desktop/Human_LA_LAA/RDS/AurN_cluster_sn_LAA_CD45_1.rds")
AurN.combined_CD45_1 <- readRDS ("~/Desktop/Human_LA_LAA/RDS/AurN_cluster_sn_LAA_CD45_1.rds")

# コンタミ探し
VlnPlot(AurN.combined_CD45_1, pt.size = 0, features = c("NPPA","MYL7","MYL4","TTN","TNNT2","ACTC1"))
VlnPlot(AurN.combined_CD45_1, pt.size = 0, features = c("RYR2","PLN","SLC8A1","DCN","GSN","PDGFRA"))
VlnPlot(AurN.combined_CD45_1, pt.size = 0, features = c("DCN","GSN","PDGFRA","RGS5","ABCC9","KCNJ8"))
VlnPlot(AurN.combined_CD45_1, pt.size = 0, features = c("MYH11","TAGLN","ACTA2","PLP1","NRXN1","NRXN3"))
VlnPlot(AurN.combined_CD45_1, pt.size = 0, features = c("GPAM","FASN","LEP","MSLN","WT1","BNC1"))
VlnPlot(AurN.combined_CD45_1, pt.size = 0, features = c("PTPRC","CD163","CD79B","CD14","IL7R","CD3E","CD4","CD8A","NKG7","GNLY","GZMB"))

VlnPlot(AurN.combined_CD45_1, pt.size = 0, features = c("PTPRC","CD163","CD79B","CD14","CD68"))
VlnPlot(AurN.combined_CD45_1, pt.size = 0, features = c("PTPRC","CD163","CD79B","CD14","CD68","VCAN","TIMP1"))


DimPlot(AurN.combined_CD45_1, reduction = "umap", split.by = "stim")
DimPlot(AurN.combined_CD45_1, reduction = "umap", label = TRUE, split.by = "AFControl")
DimPlot(AurN.combined_CD45_1, reduction = "umap", label = TRUE,  repel = TRUE)


# Myeloid----
VlnPlot(AurN.combined_CD45_1, pt.size = 0, features = c("CD68","CD14","CSF1R","IL1B","S100A8", "S100A9") )
VlnPlot(AurN.combined_CD45_1, pt.size = 0, features = c("TREM2","CD9","APOE","CXCL3","C1QA","TNF","CXCL8") )
VlnPlot(AurN.combined_CD45_1, pt.size = 0, features = c("HLA-DPA1","HLA-DPB1","HLA-DQA1","CLEC10A","FCER1A","CD1C") )
VlnPlot(AurN.combined_CD45_1, pt.size = 0, features = c("PTPRC","HLA-DPB1","HLA-DQA1","CLEC10A","FCER1A","CD1C") )

# Bcells-----
VlnPlot(AurN.combined_CD45_1, pt.size = 0, features = c("CD79A","CD79B","FCER2","CD22","MYH11","CD34") )
VlnPlot(AurN.combined_CD45_1, pt.size = 0, features = c("LUM","TAGLN","CD34","CDH5","PECAM1","TYROBP") )

# Tcells-----
VlnPlot(AurN.combined_CD45_1, pt.size = 0, features = c("CD3E","CD4","CD8A","NKG7","GNLY","GZMB","IL7R") )
VlnPlot(AurN.combined_CD45_1, pt.size = 0, features = c("KLRB1","CCR6","CD1d1","MKI67","IL2RB","GZMB") )

VlnPlot(AurN.combined_CD45_1, pt.size = 0, features = c("TIMP1"), split.by = "stim" , cols = c("sa" ="deepskyblue", "Aur1Epi"="red"))

VlnPlot(AurN.combined_CD45_1, pt.size = 0, features = c("TIMP1") )
FeaturePlot (AurN.combined_CD45_1, features = c("TIMP1"), min.cutoff=0.00005, max.cutoff='q90')

AurN.combined_CD45_1_trans <- SCTransform (AurN.combined_CD45_1)

cd_genes <- c("CD3E","NKG7","CD8A","CD68","CD14","S100A8","HLA-DQA1","CD79A","CD79B","GNLY","TYROBP","KIT", "MYH11","PECAM1","CDH5","ACTC1","TNNC1","MKI67")
DotPlot(AurN.combined_CD45_1_trans,features = cd_genes)+RotatedAxis()+coord_flip()


AurN.combined_CD45_1.markers <- FindAllMarkers(AurN.combined_CD45_1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
AurN.combined_CD45_1.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- AurN.combined_CD45_1.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(AurN.combined_CD45_1, features = top10$gene) 
write.table(top10, "table.txt", quote=F, col.names=F, append=T)


#最終CD45(myocardium除去)----
AurN.combined_CD45_2 <- subset(AurN.combined_CD45_1, idents = c(0,1,2,3,4,5,6))
AurN.combined_CD45_2 <- ScaleData(AurN.combined_CD45_2, verbose = FALSE)
AurN.combined_CD45_2<- RunPCA(AurN.combined_CD45_2, npcs = 30, verbose = FALSE)
AurN.combined_CD45_2 <- RunUMAP(AurN.combined_CD45_2 , reduction = "pca", dims = 1:30)
AurN.combined_CD45_2 <- FindNeighbors(AurN.combined_CD45_2 , reduction = "pca", dims = 1:30)
AurN.combined_CD45_2 <- FindClusters(AurN.combined_CD45_2, resolution = 0.4)



# Visualization 別々------
DimPlot(AurN.combined_CD45_2, reduction = "umap", label = FALSE, repel = FALSE)
DimPlot(AurN.combined_CD45_2, reduction = "umap", split.by = "stim")
DimPlot(AurN.combined_CD45_2, reduction = "umap", split.by = "AFDuration")

# RDS
saveRDS(AurN.combined_CD45_2, "~/Desktop/Human_LA_LAA/RDS/AurN_cluster_sn_LAA_CD45_2.rds")
AurN.combined_CD45_2 <- readRDS ("~/Desktop/Human_LA_LAA/RDS/AurN_cluster_sn_LAA_CD45_2.rds")

# コンタミ探し
VlnPlot(AurN.combined_CD45_2, pt.size = 0, features = c("NPPA","MYL7","MYL4","TTN","TNNT2","ACTC1"))
VlnPlot(AurN.combined_CD45_2, pt.size = 0, features = c("RYR2","PLN","SLC8A1","DCN","GSN","PDGFRA"))
VlnPlot(AurN.combined_CD45_2, pt.size = 0, features = c("DCN","GSN","PDGFRA","RGS5","ABCC9","KCNJ8"))
VlnPlot(AurN.combined_CD45_2, pt.size = 0, features = c("MYH11","TAGLN","ACTA2","PLP1","NRXN1","NRXN3"))
VlnPlot(AurN.combined_CD45_2, pt.size = 0, features = c("GPAM","FASN","LEP","MSLN","WT1","BNC1"))
VlnPlot(AurN.combined_CD45_2, pt.size = 0, features = c("PTPRC","CD163","CD79B","CD14","IL7R","CD3E","CD4","CD8A","NKG7","GNLY","GZMB"))


# Myeloid----
VlnPlot(AurN.combined_CD45_2, pt.size = 0, features = c("CD68","CD14","CSF1R","IL1B","S100A8", "S100A9") )
VlnPlot(AurN.combined_CD45_2, pt.size = 0, features = c("TREM2","CD9","APOE","CXCL3","C1QA","TNF","CXCL8") )
VlnPlot(AurN.combined_CD45_2, pt.size = 0, features = c("HLA-DPA1","HLA-DPB1","HLA-DQA1","CLEC10A","FCER1A","CD1C") )
VlnPlot(AurN.combined_CD45_2, pt.size = 0, features = c("PTPRC","HLA-DPB1","HLA-DQA1","CLEC10A","FCER1A","CD1C") )

# Bcells-----
VlnPlot(AurN.combined_CD45_2, pt.size = 0, features = c("CD79A","CD79B","FCER2","CD22","MYH11","CD34") )
VlnPlot(AurN.combined_CD45_2, pt.size = 0, features = c("LUM","TAGLN","CD34","CDH5","PECAM1","TYROBP") )

# Tcells-----
VlnPlot(AurN.combined_CD45_2, pt.size = 0, features = c("CD3E","CD4","CD8A","NKG7","GNLY","GZMB","IL7R") )
VlnPlot(AurN.combined_CD45_2, pt.size = 0, features = c("KLRB1","CCR6","CD1d1","MKI67","IL2RB","GZMB") )

VlnPlot(AurN.combined_CD45_2, pt.size = 0, features = c("SPP1","TIMP1") )


AurN.combined_CD45_2.markers <- FindAllMarkers(AurN.combined_CD45_2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
AurN.combined_CD45_2.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top10 <- AurN.combined_CD45_2.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(AurN.combined_CD45_2, features = top10$gene) 
write.table(top10, "table.txt", quote=F, col.names=F, append=T)



