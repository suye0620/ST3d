# ##############################################################################
# set Python environment
# ##############################################################################
reticulate::use_condaenv("C:\\Users\\jinLab1\\.conda\\envs\\3dreconstruction\\")
library(reticulate)
options(bitmapType = "cairo")

# ##############################################################################
# convert Python object to R object
# ##############################################################################
.convert <- function(obj) {
  builtins <- reticulate::import_builtins()
  Mapping <- reticulate::import("typing")$Mapping
  DataFrame <- reticulate::import("pandas")$DataFrame
  issparse <- reticulate::import("scipy.sparse")$issparse
  isinstance <- builtins$isinstance
  
  if (issparse(obj)) {
    obj <- obj$tocoo()
    return(
      Matrix::sparseMatrix(
        i = as.vector(reticulate::py_to_r(obj$row)) + 1,
        j = as.vector(reticulate::py_to_r(obj$col)) + 1,
        x = as.vector(reticulate::py_to_r(obj$data)),
        dims = unlist(as.vector(reticulate::py_to_r(obj$shape)))
      )
    )
  }
  if (isinstance(obj, DataFrame)) {
    obj <- reticulate::py_to_r(obj)
    attr(obj, "pandas.index") <- NULL
    return(obj)
  }
  if (!isinstance(obj, Mapping)) {
    return(reticulate::py_to_r(obj))
  }
  ret <- list()
  for (key in builtins$list(obj$keys())) {
    ret[[key]] <- .convert(obj[[key]])
  }
  return(ret)
}


convert.h5ad <- function(adata) {
  X <- .convert(adata$X)
  layers <- .convert(adata$layers)
  obs <- .convert(adata$obs)
  var <- .convert(adata$var)
  obsm <- .convert(adata$obsm)
  varm <- .convert(adata$varm)
  obsp <- .convert(adata$obsp)
  varp <- .convert(adata$varp)
  uns <- .convert(adata$uns)
  rownames(X) <- rownames(obs)
  colnames(X) <- rownames(var)
  return(
    list(
      X = X,
      layers = layers,
      obs = obs,
      var = var,
      obsm = obsm,
      varm = varm,
      obsp = obsp,
      varp = varp,
      uns = uns
    )
  )
}


read.h5ad <- function(filename) {
  ad <- reticulate::import("anndata", convert = FALSE)
  adata <- ad$read_h5ad(filename)
  return(convert.h5ad(adata))
}


safe.sd <- function(x) {
  if (length(x) == 1)
    return(0.0)
  return(sd(x))
}

create_SeuratObj_from_AnnData <- function(list_anndata){
  # transpose
  mat <- Matrix::t(list_anndata$X)
  colnames(mat) <- rownames(list_anndata$obs)
  rownames(mat) <- rownames(list_anndata$var)
  
  new_SeuratObj <- Seurat::CreateSeuratObject(counts = mat, 
                                    meta.data = list_anndata$obs, 
                                    assay = "Spatial")
  coord.df <- data.frame(x = list_anndata$obsm$spatial[, 1], 
                           y = list_anndata$obsm$spatial[, 2], 
                           stringsAsFactors = FALSE)
  rownames(coord.df) <- Seurat::Cells(new_SeuratObj)
  
  new_SeuratObj@images$image <- new(
    Class = "SlideSeq", # can not change class name
    assay = "Spatial",
    key = "image_",
    coordinates = coord.df
  )
  
  return(new_SeuratObj)
}


ggplot.theme <- function(...) {
  return(
    ggplot2::theme(
      text = ggplot2::element_text(family = "ArialMT"),
      plot.background = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "#EEEEEE", linetype = "longdash"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "#FFFFFF"),
      legend.background = ggplot2::element_blank(),
      legend.box.background = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "#000000"),
      ...
    )
  )
}


ggplot.save <- function(filename, ...) {
  ggplot2::ggsave(filename, ..., dpi = 600, bg = "transparent")
}
