library(tools)
source(file_path_as_absolute("./.Rprofile"))
# options(renv.config.pak.enabled=TRUE)

bruceR::set.wd()
h5ad_dir <- file_path_as_absolute("./data/SCC/cached-results/H5ADs")
# h5ad_dir <- "./data/SCC/cached-results/H5ADs/patient_2_slice_0.h5ad"

library(Seurat)
# patient_2_slice_0 <- read.h5ad(paste0(h5ad_dir,"/aligned_patient_2_slice_0.h5ad")) |> 
#   create_SeuratObj_from_AnnData()
# patient_2_slice_1 <- read.h5ad(paste0(h5ad_dir,"/aligned_patient_2_slice_1.h5ad")) |>
#   create_SeuratObj_from_AnnData()
# patient_2_slice_2 <- read.h5ad(paste0(h5ad_dir,"/aligned_patient_2_slice_2.h5ad")) |>
#   create_SeuratObj_from_AnnData()

# saveRDS(patient_2_slice_0,file = ("./rds/patient_2_slice_0.Rds"))
# saveRDS(patient_2_slice_1,file = ("./rds/patient_2_slice_1.Rds"))
# saveRDS(patient_2_slice_2,file = ("./rds/patient_2_slice_2.Rds"))

# ##############################################################################
# read RDS files
# ##############################################################################
patient_2_slice_0 <- readRDS("./rds/patient_2_slice_0.Rds")
patient_2_slice_1 <- readRDS("./rds/patient_2_slice_1.Rds")
patient_2_slice_2 <- readRDS("./rds/patient_2_slice_2.Rds")

# ##############################################################################
# 3d plot
# ##############################################################################
# ##############################################################################
# TODO:
# 1. 绘制3d重建
# 2. 学习cellchat、spatialchat，跑细胞通讯的结果
# 3. 思考3d情况下怎么考虑细胞通讯
# ##############################################################################

