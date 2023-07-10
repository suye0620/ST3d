# set directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("./.Rprofile")
# options(renv.config.pak.enabled=TRUE)

h5ad_dir <- ("./data/SCC/cached-results/H5ADs/patient_2/")

library(Seurat)
# patient_2_slice_0 <- read.h5ad(paste0(h5ad_dir,"/aligned_patient_2_slice_0.h5ad")) |>
#   create_SeuratObj_from_AnnData()
# patient_2_slice_1 <- read.h5ad(paste0(h5ad_dir,"/aligned_patient_2_slice_1.h5ad")) |>
#   create_SeuratObj_from_AnnData()
# patient_2_slice_2 <- read.h5ad(paste0(h5ad_dir,"/aligned_patient_2_slice_2.h5ad")) |>
#   create_SeuratObj_from_AnnData()
# 
# saveRDS(patient_2_slice_0,file = ("./rds/patient_2_slice_0.Rds"))
# saveRDS(patient_2_slice_1,file = ("./rds/patient_2_slice_1.Rds"))
# saveRDS(patient_2_slice_2,file = ("./rds/patient_2_slice_2.Rds"))

# ##############################################################################
# read RDS files
# ##############################################################################
rds_dir <- "./rds/"
patient2_slices_dir <- paste0(rds_dir,"SCC-Patient2/slices/")
patient2_align_dir <- paste0(rds_dir,"SCC-Patient2/align/")

# load slices data
patient2_slices <- list()
for (rds_file in list.files(patient2_slices_dir)) {
  patient2_slices[[str_remove(rds_file,".Rds")]] <- readRDS(paste0(patient2_slices_dir,rds_file))
} 

# load mappings
patient2_mappings <- list()
for (rds_file in list.files(patient2_align_dir)) {
  patient2_mappings[[str_remove(rds_file,".Rds")]] <- readRDS(paste0(patient2_align_dir,rds_file))
} 
# patient_2_align01 <- readr::read_csv(paste0(h5ad_dir,"/pi01.csv"),col_names = F)
# patient_2_align12 <- readr::read_csv(paste0(h5ad_dir,"/pi12.csv"),col_names = F)
# saveRDS(patient_2_align01,file = ("./rds/patient_2_align01.Rds"))
# saveRDS(patient_2_align12,file = ("./rds/patient_2_align12.Rds"))

# ##############################################################################
# 3d plot
# TODO:
# 1. 绘制3d重建
# 2. 学习cellchat、spatialchat，跑细胞通讯的结果
# 3. 思考3d情况下怎么考虑细胞通讯
# ##############################################################################
library(plotly)
library(tidyverse)


get_nsamples_from_mat <- function(mat,n=200){
  selected_rows <- sample(nrow(mat), n)
  selected_matrix <- mat[selected_rows, ]
  return(selected_matrix)
}

Alignment_3Dplot <- function(list_slices,list_mappings,
                             z_scale=10){
  list_df <- vector("list",length = length(list_slices))
  forloop.counter <- 0
  
  for (slice in names(list_slices)) {
    forloop.counter <- forloop.counter+1
    df <- cbind(list_slices[[slice]]@images$image@coordinates,
          list_slices[[slice]]$original_clusters)
    colnames(df) <- c("x","y","clusters")
    df$layer <- slice
    df[[paste0("slice",(forloop.counter-1),"_index")]] <- 1:nrow(df)
    df$slice <- (forloop.counter-1)
    df$z <- z_scale*(forloop.counter-1)
    list_df[[forloop.counter]] <- df
  }
  df_plot <- dplyr::bind_rows(list_df)
  return(df_plot)
}

df_plot <- Alignment_3Dplot(list_slices = new_slices,z_scale = 20)

set.seed(666)
list_mappings <- list(
  "slice0_slice1"=get_positive_mapping(patient_2_align01) %>% get_nsamples_from_mat(),
  "slice1_slice2"=get_positive_mapping(patient_2_align12) %>% get_nsamples_from_mat()
)

p <- plot_ly(df_plot,type = 'scatter3d') %>% 
  add_markers(x = ~x, y = ~y, z = ~z,color=~clusters)

for (mapping in names(list_mappings)){
  slice_split <- str_split_1(mapping,"_")
  first_slice <- slice_split[1]
  second_slice <- slice_split[2]
  for (m in 1:nrow(list_mappings[[mapping]]) ) {
    p <- p %>% add_fun(function(plot){
      plot %>% filter(!!rlang::parse_expr(paste0(first_slice,"_index==",
                                                 list_mappings[[mapping]][m,1],
                                                 "|",
                                                 second_slice,"_index==",
                                                 list_mappings[[mapping]][m,2]))) %>%
        add_trace(x = ~x, y = ~y, z = ~z,
                  line = list(color = 'black', width = 3,dash='longdash'),
                  mode="lines",
                  name="",showlegend=F)
    })
  }
  
}

# ##############################################################################
# test code
# ##############################################################################
# 传入plot_ly是三个变量，所有的点被连在了一起
plot_ly(df_plot) %>% 
  add_markers(x = ~x, y = ~y, z = ~z,color=~clusters) %>% 
  add_fun(function(plot,index1=3,index2=5){
    plot %>% filter(!!rlang::parse_expr(paste0("slice0_index==",index1,
                                               "|slice1_index==",index2))) %>% 
    add_trace(x = ~x, y = ~y, z = ~z,
              line = list(color = 'black', width = 3,dash='longdash'),
              mode="lines",
              name="",showlegend=F)
  }) %>% 
  add_fun(function(plot,index1=6,index2=7){
    plot %>% filter(!!rlang::parse_expr(paste0("slice0_index==",index1,
                                               "|slice1_index==",index2))) %>% 
      add_trace(x = ~x, y = ~y, z = ~z,
                line = list(color = 'black', width = 3,dash='longdash'),
                mode="lines",
                name="",showlegend=F)
  })







