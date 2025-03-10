library(Seurat)
library(CellChat)
library(patchwork)
library(ggplot2)

sn_LAA_CM<- readRDS("~/Desktop/Human_LA_LAA/RDS/AurN_cluster_sn_LAA_cardiomyocyte.rds")
sc_LAA_imuune <- readRDS("~/Desktop/Human_LA_LAA/RDS/Aur_cluster_scMyeloid_final.rds")


DimPlot(sn_LAA_CM)
DimPlot(sc_LAA_imuune)


#sc_LAA_immune(myeloidのこと)から3Otherをぬく
sc_LAA_imuune <- subset(sc_LAA_imuune, subset = AFDuration1ex != c("3Other") )

#クラスタnameへんこう
new.cluster.ids <- c("My.10","My.11","My.12","My.13","My.14","My.15","My.16","My.17","My.18","My.19","My.20","My.21","My.22","My.23","My.24")
names(new.cluster.ids) <- levels(sc_LAA_imuune)
sc_LAA_imuune<- RenameIdents(sc_LAA_imuune, new.cluster.ids)

new.cluster.ids <- c("CM.0","CM.1","CM.2")
names(new.cluster.ids) <- levels(sn_LAA_CM)
sn_LAA_CM<- RenameIdents(sn_LAA_CM, new.cluster.ids)



#非免疫細胞全体+Myeloid

sc_sn_LAA <- merge(sn_LAA_CM, y = sc_LAA_imuune)



levels(sc_sn_LAA)

new.cluster.ids <- c("CM.0","CM.1","CM.2","My.0","My.1","My.2","My.3","My.4","My.5","My.6","My.7","My.8","My.9","My.10","My.11","My.12","My.13","My.14")
names(new.cluster.ids) <- levels(sc_sn_LAA)
sc_sn_LAA<- RenameIdents(sc_sn_LAA, new.cluster.ids)

DefaultAssay(sc_sn_LAA) <- "RNA"
cellchat <- createCellChat(object = sc_sn_LAA, group.by = "ident")



CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

CellChatDB[["interaction"]]$receptor[763] <- "IL1R1"
CellChatDB[["interaction"]]$receptor[763] 

CellChatDB[["interaction"]]$receptor[762] <- "IL1R1"
CellChatDB[["interaction"]]$receptor[762] 

cellchat@DB <- CellChatDB


#### 遺伝子発現データをcell-cell communication解析に使うための前処理 ####
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel

# 各クラスタで過剰発現したシグナル伝達遺伝子をobject@var.features に保存
cellchat <- identifyOverExpressedGenes(cellchat, thresh.p = 0.1)
# View(cellchat@var.features$features.info)
# 過剰発現しているリガンド-受容体相互作用を特定し、object@LRに保存
cellchat <- identifyOverExpressedInteractions(cellchat)
# View(cellchat@LR$LRsig)

#### コミュニケーションの確率の計算と、細胞間のコミュニケーションネットワークの推定####
# コミュニケーションの重要度を推定しobjrct@netに保存

cellchat <- computeCommunProb(cellchat, type="truncatedMean")

# 少数の細胞しか関与しないコミュニケーションはフィルターする
cellchat <- filterCommunication(cellchat, min.cells = 10)

#### シグナル伝達経路レベルでの細胞間情報伝達の推定 ####
# シグナル伝達経路の重要度を推定し object@netPに保存
cellchat <- computeCommunProbPathway(cellchat)
#write.csv(cellchat@netP$prob, "~/Desktop/JRAS/R-CellChat/output/netP.csv")
#### cell-cellコミュニケーションネットワークの集計を行う ####
# 集計された細胞間コミュニケーションの結果をobject@net$count, object@net$weight に保存
cellchat <- aggregateNet(cellchat)
# View(cellchat@net$count)
# View(cellchat@net$weight)

#sc_myeloid+sn_FB
saveRDS(cellchat, "~/Desktop/Human_LA_LAA/RDS/cellchat_result_Aur_scAll.Myeloid_sn.CM_merged_IL1new.rds")

cellchat <- readRDS("~/Desktop/Human_LA_LAA/RDS/cellchat_result_Aur_scAll.Myeloid_sn.CM_merged_IL1new.rds")



#---------------------------------ここから

groupSize <- as.numeric(table(cellchat@idents)) # クラスターごとの細胞数を保存
par(mfrow = c(1,2), xpd=TRUE) # 次のプロットを並べて表示するためのパラメータの設定
# 細胞間コミュニケーションの数とweight/strength をプロット
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# 各クラスタごとにバラバラにして細胞間コミュニケーションを可視化
mat <- cellchat@net$weight
#dev.off()
par(mfrow = c(4,6), xpd=TRUE) # 次のプロットを2行4列で表示するためのパラメータの設定
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i], arrow.size = 0.01)
}



#### 各シグナル伝達経路を可視化する（Hierarchy plot, Circle plot, Chord diagram）  ####
# 次のコマンドで表示されたpathwayから、可視化したいものを選ぶ（複数選択可）
cellchat@netP$pathways 
pathways.show <- c("EGF")  # 上記で選んだpathwayを記入

#print(cellchat@netP$pathways)
#print(cellchat@netP$prob[,,15])

# 細胞間コミュニケーションの可視化
# 各クラスタの並びを確認する
levels(cellchat@idents) 
# identsの中から左側に表示する細胞を番号でvertex.receiverに指定
vertex.receiver = c(4,5,6,7,8,9,10,11,12,13,14,15,16,17,18) # a numeric vector. myeloid

# # Hierarchy plot
# Circle plot
# クラスタを円上に並べてシグナル伝達を矢印で可視化
par(mfrow=c(1,1))
library(scales) # 色情報を取得するためのライブラリ
n_cluster <- 18 # cluster数を指定
palette <- hue_pal()(n_cluster)   # 色情報を取得
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle",color.use  = palette, show.legend = TRUE) # 色情報を指定


netVisual_aggregate(cellchat, signaling = pathways.show, layout="hierarchy", vertex.receiver = vertex.receiver)
# 左と右の部分はそれぞれ、選んだvertex.receiverへのオートクラインシグナルと、
# vertex.receiver以外へのパラクラインシグナルを強調している。
# 円と白抜き円はそれぞれソースとターゲットを表す。
# 円の大きさは各細胞群の細胞数に比例し、エッジの幅はコミュニケーション確率を表している。
# エッジの色はシグナル発信源と一致している。

# Circle plot
# クラスタを円上に並べてシグナル伝達を矢印で可視化
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Chord diagram
# クラスタを炎上に並べてChord diagramで可視化
par(mfrow=c(1,1))
pdf(file ="~/Desktop/JRAS/R-CellChat/cellchat.pdf", width = 20, height =16) # RStudioのキャンパスサイズだと、plotが小さくなるのでpdfに描画する
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
dev.off() #描画完了


#### ネットワーク中心性スコアの計算と可視化 ####
# 4つのネットワーク中心性指標に基づく各細胞グループの相対的重要性を示す
# https://github.com/sqjin/CellChat/issues/7 
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")  #スロット「netP」は、シグナル伝達経路の推定された細胞間通信ネットワーク
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

#### 複数の細胞タイプやシグナル伝達経路がどのように協調しているかを探るため、グローバルなコミュニケーションパターンを特定する ####
library(NMF)
library(ggalluvial)

# NMF Rパッケージに実装されているCopheneticとSilhouetteという2つのメトリクスに基づいて
# パターン数を推論する。
# どちらの指標も、コンセンサス行列の階層的クラスタリングに基づき、
# 特定のパターン数に対する安定性を測定する。
# パターン数の範囲では、CopheneticとSilhouetteの値が急に下がり始めるパターンが
# 適切なパターン数である。
selectK(cellchat, pattern = "outgoing") # 2つのプロットのが下がったNを調べる。

nPatterns = 6 # パターンの数を指定する
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)

netAnalysis_dot(cellchat, pattern = "outgoing")

# river plot
netAnalysis_river(cellchat, pattern = "outgoing")
#> Please make sure you have load `library(ggalluvial)` when running this function



#### 機能的な類似性に基づいてシグナル伝達グループを特定する ####
# ネットワークの類似性を計算
cellchat <- computeNetSimilarity(cellchat, type = "functional")
# シグナルネットワークに対してUMAPを作成
cellchat <- netEmbedding(cellchat, type = "functional", umap.method = "uwot")
# シグナルネットワークをUMAP上でクラスタリング
cellchat <- netClustering(cellchat, type = "functional")
# シグナル伝達ネットワークをUMAPで可視化
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)

View(CellChatDB[["interaction"]])

CellChatA = CellChatDB[["interaction"]]

View(CellChatA[CellChatA$pathway_name=="CD40",])

CellChatA[CellChatA$pathway_name=="CD40",]

write.csv(CellChatA[CellChatA$pathway_name=="CD40",],"~/Desktop/DCA_human/Results/230323_11_uapmi_CD40")

#View(CellChatDB["interaction$pathway_name" == "CD40"])




