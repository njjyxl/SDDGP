setwd("/home/test/Bcbio_project/Deconvolution")

options(warn = -1)
suppressMessages({
  library(Seurat); library(SeuratDisk); library(Matrix); library(CARD)
})

# ===== load data =====

# HD.Mouse.brain_cortex
scrna_path = 'data/HD.Mouse.Brain/allen_cortex_SC.RDS'
spatial_path = 'data/HD.Mouse.Brain/cortex_ST.RDS'
load("data/HD.Mouse.Brain/cortex_pos.Rdata")
sc_obj <- readRDS(scrna_path)
spatial_obj <- readRDS(spatial_path)
sc_count <- as.matrix(sc_obj@assays$RNA@counts)
spatial_count <- as.matrix(GetAssayData(spatial_obj, assay = "Spatial.016um", layer = "counts"))
cellType <- sc_obj@meta.data

gc(full = TRUE, reset = TRUE)

cat("\n========== benchmark START ==========\n")
t_start  <- Sys.time()
proc_start <- proc.time()

# CARD
CARD_obj = createCARDObject(
    sc_count = sc_count, sc_meta = cellType,
    spatial_count = spatial_count, spatial_location = pos,
    ct.varname = "celltype", ct.select = unique(cellType$celltype),
    sample.varname = "sample_id", minCountGene = 100,
    minCountSpot = 5) 
  
# Perform CARD deconvolution
CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
CARD_result <- CARD_obj@Proportion_CARD
save(CARD_result,file = "merfish_decon_result/CARD_result.Rdata")

proc_end <- proc.time()
t_end    <- Sys.time()
cat("========== benchmark END ==========\n\n")

# ===== Metrics recorded on the R side =====
wall_sec <- as.numeric(difftime(t_end, t_start, units = "secs"))
cpu_user <- (proc_end - proc_start)[["user.self"]]   + (proc_end - proc_start)[["user.child"]]
cpu_sys  <- (proc_end - proc_start)[["sys.self"]]    + (proc_end - proc_start)[["sys.child"]]

# R internal memory
gc_after <- gc(reset = FALSE)
r_mem_mb <- sum(gc_after[, "used"] * c(8, 8) / 1024) 
cat(sprintf("Wall-clock time    : %.2f s\n", wall_sec))
cat(sprintf("CPU time (user+sys): %.2f s (user=%.2f, sys=%.2f)\n",
            cpu_user + cpu_sys, cpu_user, cpu_sys))
cat(sprintf("R-side memory used : %.2f MB (lower bound; see /usr/bin/time -v for true Max RSS)\n",
            r_mem_mb))

