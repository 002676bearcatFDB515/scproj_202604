library(Seurat)
library(scRepertoire)
library(devtools)
library(circlize)
library(scales)
#devtools::install_github("ncborcherding/scRepertoire@dev")
#manual: https://ncborcherding.github.io/vignettes/vignette.html
#not install from bioconductor,it is from github. 
#2.0 changes cloneType -> cloneSize compareClonaltypes -> clonalCompare highlightClonotypes -> highlightClones
# alluvialClonotypes -> alluvialClones
# split.by was combined into group.by
# expression2List was deprecated, but can be used with flag "force = TRUE"
# clonalCompare arg "numbers" is now "top.clones"
# StartracDiversity "sample" is not an arg any longer, replaced by "group.by"; "by" is deprecated but default is all 3 (expa, migr, tran)
# clonalNetwork "identity" is not an arg any longer, replaced by "group.by"

setwd("~/2026sc/scRepertoire")

G2_63 <- read.csv("~/2026sc/63/filtered_contig_annotations.csv")
G3_68 <- read.csv("~/2026sc/68/filtered_contig_annotations.csv")
G1_71 <- read.csv("~/2026sc/71/filtered_contig_annotations.csv")
G2_73 <- read.csv("~/2026sc/73/filtered_contig_annotations.csv")
G2_79 <- read.csv("~/2026sc/79/filtered_contig_annotations.csv")
G3_80 <- read.csv("~/2026sc/80/filtered_contig_annotations.csv")
G4_83 <- read.csv("~/2026sc/83/filtered_contig_annotations.csv")
G1_84 <- read.csv("~/2026sc/84/filtered_contig_annotations.csv")
G1_87 <- read.csv("~/2026sc/87/filtered_contig_annotations.csv")
G3_88 <- read.csv("~/2026sc/88/filtered_contig_annotations.csv")
G4_90 <- read.csv("~/2026sc/90/filtered_contig_annotations.csv")
G4_92 <- read.csv("~/2026sc/92/filtered_contig_annotations.csv")
contig_list <- list(G2_63, G3_68, G1_71, G2_73, G2_79, G3_80, G4_83, G1_84, G1_87, G3_88, G4_90, G4_92)
combined <- combineTCR(contig_list, 
                       samples = c("G2_63", "G3_68", "G1_71", "G2_73", "G2_79", "G3_80", "G4_83", "G1_84", "G1_87", "G3_88", "G4_90", "G4_92"), 
                                             )
# ID = c("G2", "G3", "G1", "G2", "G2", "G3", "G4", "G1", "G1", "G3", "G4", "G4")
clonalQuant(combined, cloneCall="strict", scale = T)

vizGenes(combined, 
         x.axis = "TRAV", 
         plot = "barplot", 
         order = "variance")
hms_cluster_id <- readRDS(file="~/2026sc/Only_t_cluster_id_test.rds")
seurat <- readRDS(file="~/2026sc/Only_t_cluster_id_test.rds")
DimPlot(seurat, label = T) + NoLegend()
table(Idents(seurat))
#add seurat
seurat <- combineExpression(combined, seurat, 
                            cloneCall= "gene", group.by = NULL, proportion = FALSE, 
                            cloneSize = c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

colorblind_vector <- colorRampPalette(rev(c("#0D0887FF", "#47039FFF", 
                                            "#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF",
                                            "#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))

names(seurat@meta.data)
head(seurat@meta.data)

DimPlot(seurat, group.by = "tech") +
  scale_color_manual(values=colorblind_vector(5), na.value="grey") + 
  theme(plot.title = element_blank())

#bottom references to cloneSize are invalid because we performed group.by = "tech" in the combineExpression step above
#cloneSize (or the v1.0 cloneType does not exist as a column for the data object)
slot(seurat, "meta.data")$cloneSize <- factor(slot(seurat, "meta.data")$cloneSize, 
                                              levels = c("Hyperexpanded (100 < X <= 500)", 
                                                         "Large (20 < X <= 100)", 
                                                         "Medium (5 < X <= 20)", 
                                                         "Small (1 < X <= 5)", 
                                                         "Single (0 < X <= 1)", NA))
DimPlot(seurat, group.by = "cloneSize") +
  scale_color_manual(values = colorblind_vector(5), na.value="grey") + 
  theme(plot.title = element_blank())

#added to help visualize proportions for top clone in each sample
clonalProportion(combined.TCR_filter, cloneCall = "aa", clonalSplit = c(1, 10, 25, 100, 500, 1000, 1e+05))
#只能做出blood和腫瘤旁邊相近的clonotype,不能做出全部
#add "exportTable = TRUE" to retrieve a table
clonalCompare(combined, top.clones = 3, samples = c("G2_63", "G3_68", "G1_71", "G2_73", "G2_79", "G3_80", "G4_83", "G1_84", "G1_87", "G3_88", "G4_90", "G4_92"),
                  cloneCall="aa", graph = "alluvial")
#top 20 across all samples: "CAVRMTTASLGKLQF_CASGGGPANSDYTF", "CAVSANSGTYQRF_CASSIKTGGYAEQFF", "CAVSNTGANTGKLTF_CASRRNERLFF", "CAVPTNAYKVIF_CASSQEPQGGINTGQLYF", "CALSRGSALGRLHF_CASSNTGGNERLFF", "CILRVPTASLGKLQF_CAWSPDTGGAETLYF", "CAAMATGGNNKLTF_CASSVHFNQAPLF", "CAASASNTGYQNFYF_CASGPRDRGRAEQFF", "CAVNYGGSGNKLIF_CASSQKLGGPTGQLYF", "CAAITGNYKYVF_CASRGDEQYF", "NA_CAWSPDTGGAETLYF", "CAVSLPGTGSNRLTF_CASSDRSTEVFF", "CAALLTGNTRKLIF;CAVRDLGGSNAKLTF_CGAREWDTSNERLFF", "CVLGITGNTRKLIF_CTCSADGQGNTGQLYF", "CAAIASSSFSKLVF_CASSDHRGANTEVFF", "CAVSANSGTYQRF;CVLGEPSGSWQLIF_CASSSRTGGYAEQFF", "CILRGGNEKITF_CASNTGAYEQYF", "CALGAPNAGAKLTF_CASSLDRVGEQYF", "CAALASSSFSKLVF_CGAREGGGGSAETLYF", "CAVSMPSGSWQLIF_CAWRGDNSPLYF"

library(ggraph)
#這裏的aa:CALSGVTSYDKVIF_CASSLLGGGNNEQFF是以前的結果
#saveRDS(seurat, file = "seurat.rds")
seurat <- highlightClones(seurat,  cloneCall= "aa", sequence = c("CAVRMTTASLGKLQF_CASGGGPANSDYTF", "CAVSANSGTYQRF_CASSIKTGGYAEQFF", 
                                                                 "CAVSNTGANTGKLTF_CASRRNERLFF", "CAVPTNAYKVIF_CASSQEPQGGINTGQLYF", 
                                                                 "CALSRGSALGRLHF_CASSNTGGNERLFF", "CILRVPTASLGKLQF_CAWSPDTGGAETLYF", 
                                                                 "CAAMATGGNNKLTF_CASSVHFNQAPLF", "CAASASNTGYQNFYF_CASGPRDRGRAEQFF", 
                                                                 "CAVNYGGSGNKLIF_CASSQKLGGPTGQLYF", "CAAITGNYKYVF_CASRGDEQYF"))

DimPlot(seurat, group.by = "highlight", pt.size = 0.5)
theme(plot.title = element_blank())
clonalOccupy(seurat, x.axis = "ident")

alluvialClones(seurat, cloneCall = "gene", 
                   y.axes = c("orig.ident", "ident", "tech"), 
                   color = c("TRAV3-3")) + 
  scale_fill_manual(values = c("grey", colorblind_vector(2)[2]))

alluvialClones(seurat, cloneCall = "gene", 
                   y.axes = c("orig.ident", "ident"), 
                   color = "ident") 
combined2 <- expression2List(seurat, group.by = "ident", force = TRUE)
length(combined2)
clonalDiversity(combined2, cloneCall = "nt")
clonalHomeostasis(combined2, cloneCall = "nt") + theme(axis.text.x = element_text(angle = 30, hjust=1))
clonalProportion(combined2, cloneCall = "nt") + theme(axis.text.x = element_text(angle = 30, hjust=1))
clonalOverlap(combined2, cloneCall="aa", method="overlap") + theme(axis.text.x = element_text(angle = 30, hjust=1))
clonalOccupy(seurat, x.axis = "ident") + theme(axis.text.x = element_text(angle = 30, hjust=1))

library(circlize)
library(scales)
circles <- getCirclize(seurat, 
                       group.by = "ident")
grid.cols <- scales::hue_pal()(length(unique(seurat@active.ident)))
names(grid.cols) <- levels(seurat@active.ident)
circlize::chordDiagram(circles,
                       self.link = 1, 
                       grid.col = grid.cols)

StartracDiversity(seurat, 
                  type = "tech", 
                  group.by = "celltype") + theme(axis.text.x = element_text(angle = 30, hjust=1))

library(ggraph)
clonalNetwork(seurat,       reduction = "umap", 
              group.by = "ident", filter.clones = NULL,
              filter.identity = NULL, cloneCall = "aa")

#only cd8+ T cells
Only_cd8<-subset(seurat, idents=c("proliferating cd8",
                                  "activated cd8",
                                  "effector cd8",
                                  "exhausted cd8"))
DimPlot(Only_cd8, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(Only_cd8, reduction = "umap", label = FALSE, pt.size = 0.5) 
saveRDS(Only_cd8, file = "Only_cd8_cluster_id_test.rds")
seurat <- readRDS(file="Only_cd8_cluster_id_test.rds")
DimPlot(seurat, label = T) + NoLegend()
DimPlot(seurat, label = F)
table(Idents(seurat))
seurat <- combineExpression(combined, seurat, 
                            cloneCall="gene", group.by = "sample", proportion = FALSE, 
                            cloneSize=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
circles <- getCirclize(seurat,  group.by = "ident")
grid.cols <- scales::hue_pal()(length(unique(seurat@active.ident)))
names(grid.cols) <- levels(seurat@active.ident)
circlize::chordDiagram(circles, self.link = 1, 
                       grid.col = grid.cols)
clonalNetwork(seurat, reduction = "umap", group.by = "ident",
              filter.clones = NULL,filter.identity = NULL,
              cloneCall = "aa")
StartracDiversity(seurat, type = "tech", 
                  group.by = "celltype")
alluvialClones(seurat, cloneCall = "gene", 
                   y.axes = c("celltype", "ident", "tech"),   color = "ident") 
combined2 <- expression2List(seurat, group.by = "ident")
length(combined2)
clonalDiversity(combined2, cloneCall = "nt")
clonalHomeostasis(combined2, cloneCall = "nt") + theme(axis.text.x = element_text(angle = 30, hjust=1))
clonalProportion(combined2, cloneCall = "nt") + theme(axis.text.x = element_text(angle = 30, hjust=1))
clonalOverlap(combined2, cloneCall="aa", method="overlap") + theme(axis.text.x = element_text(angle = 30, hjust=1))
clonalOccupy(seurat, x.axis = "ident")
clonalOccupy(seurat, label= FALSE,x.axis = "ident")

#only cd4+ T cells
hms_cluster_id <- readRDS(file="~/2026sc/Only_t_cluster_id_test.rds")
Only_cd4<-subset(hms_cluster_id, idents=c("naive cd4",
                                          "activated cd4",
                                          "cd4 Th2",
                                          "treg"))
DimPlot(Only_cd4, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(Only_cd4, reduction = "umap", label = FALSE, pt.size = 0.5) 
saveRDS(Only_cd4, file = "Only_cd4_cluster_id_test.rds")
seurat <- readRDS(file="Only_cd4_cluster_id_test.rds")
DimPlot(seurat, label = T) + NoLegend()
DimPlot(seurat, label = F) 
table(Idents(seurat))
seurat <- combineExpression(combined, seurat, 
                            cloneCall="gene", group.by = "sample", proportion = FALSE, 
                            cloneSize=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
circles <- getCirclize(seurat,  group.by = "ident")
grid.cols <- scales::hue_pal()(length(unique(seurat@active.ident)))
names(grid.cols) <- levels(seurat@active.ident)
circlize::chordDiagram(circles, self.link = 1, 
                       grid.col = grid.cols)
clonalNetwork(seurat, reduction = "umap", group.by = "ident",
              filter.clones = NULL,filter.identity = NULL,
              cloneCall = "aa")
StartracDiversity(seurat, type = "tech", 
                  group.by = "celltype")
alluvialClones(seurat, cloneCall = "gene", 
                   y.axes = c("celltype", "ident", "tech"),   color = "ident") 
combined2 <- expression2List(seurat, group.by = "ident", force = TRUE)
length(combined2)
clonalDiversity(combined2, cloneCall = "nt")
clonalHomeostasis(combined2, cloneCall = "nt")
clonalProportion(combined2, cloneCall = "nt")
clonalOverlap(combined2, cloneCall="aa", method="overlap")
clonalOccupy(seurat, x.axis = "ident")
clonalOccupy(seurat, label= FALSE,x.axis = "ident")
