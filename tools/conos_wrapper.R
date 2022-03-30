library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(conos)
library(pagoda2)
library(rhdf5)

options <- commandArgs(trailingOnly = TRUE)
exchange_dir <- options[1]
batch.key <- options[2]
space <- options[3]
n.hvg <- strtoi(options[4])
n.comps <- strtoi(options[5])

Convert(paste0(exchange_dir, "adata.h5ad"), dest = "h5seurat")
data <- LoadH5Seurat(paste0(exchange_dir, "adata.h5seurat"))
panel <- SplitObject(data, split.by = batch.key)

# Preprocessing
for (i in 1:length(panel)) {
    panel[[i]] = panel[[i]]$RNA@counts
}
panel <- lapply(
    panel,
    basicP2proc,
    n.cores = 30,
    min.cells.per.gene = 0,
    n.odgenes = n.hvg,
    get.largevis = FALSE,
    make.geneknn = FALSE,
    get.tsne = FALSE
)

# Conos
con <- Conos$new(panel)
con$buildGraph(
    k = 15,
    k.self = 5,
    space = space,
    ncomps = n.comps,
    n.odgenes = n.hvg,
    matching.method = "mNN",
    metric = "angular",
    verbose = TRUE
)
con$embedGraph(method = "UMAP")
saveRDS(con, file = paste0(exchange_dir, "conos.rds"))

# Saving for Scanpy
hdf5file <- "conos.h5"
saveConosForScanPy(
  con,
  output.path = exchange_dir,
  hdf5_filename = hdf5file,
  cm.norm = TRUE,
  pseudo.pca = TRUE,
  pca = TRUE,
  n.dims = 100,
  embedding = TRUE,
  alignment.graph = TRUE,
  verbose = TRUE
)

count.matrix <- h5read(paste0(exchange_dir, hdf5file), "count_matrix")
pseudo.pca <- h5read(paste0(exchange_dir, hdf5file), "pseudopca/pseudopca.df")
pca <- h5read(paste0(exchange_dir, hdf5file), "pca/pca.df")
embedding <- h5read(paste0(exchange_dir, hdf5file), "embedding/embedding.df")
graph.connectivities <- h5read(paste0(exchange_dir, hdf5file), "graph_connectivities")
graph.distances <- h5read(paste0(exchange_dir, hdf5file), "graph_distances")

write.table(embedding, paste0(exchange_dir, "embedding.csv"), sep = ",", quote = FALSE, row.names = FALSE)
write.table(pseudo.pca, paste0(exchange_dir, "pseudo_pca.csv"), sep = ",", quote = FALSE, row.names = FALSE)
write.table(pca, paste0(exchange_dir, "pca.csv"), sep = ",", quote = FALSE, row.names = FALSE)

count.matrix = sparseMatrix(
    x = as.numeric(count.matrix$data),  
    dims = as.numeric(c(count.matrix$shape[1], count.matrix$shape[2])), 
    p = as.numeric(count.matrix$indptr), 
    i = count.matrix$indices,
    index1 = FALSE
)
writeMM(count.matrix, paste0(exchange_dir, "count_matrix.mtx"))

graph.connectivities = sparseMatrix(
    x = as.numeric(graph.connectivities$data),  
    dims = as.numeric(c(graph.connectivities$shape[1], graph.connectivities$shape[2])), 
    p = as.numeric(graph.connectivities$indptr), 
    i = graph.connectivities$indices,
    index1 = FALSE
)
writeMM(graph.connectivities, paste0(exchange_dir, "graph_connectivities.mtx"))

graph.distances = sparseMatrix(
    x = as.numeric(graph.distances$data),  
    dims = as.numeric(c(graph.distances$shape[1], graph.distances$shape[2])), 
    p = as.numeric(graph.distances$indptr), 
    i = graph.distances$indices,
    index1 = FALSE
)
writeMM(graph.distances, paste0(exchange_dir, "graph_distances.mtx"))