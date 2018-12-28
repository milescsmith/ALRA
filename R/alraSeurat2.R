#' @title alraSeurat2
#' @description Function to perform ALRA on a Seurat2 object
#'
#' @details This function performs ALRA on the "data" slot in a given Seurat object,
#' does the neccessary pre- and post-processing to fit with the object, places 
#' the imputed data in an assay slot named "ALRA", and returns that updated Seurat 
#' object.
#'
#' @param obj the Seurat object to run ALRA
#' @param assay.use The Seurat assay to impute. Default: "RNA"
#' @param slot.use The slot within the assay to impute. Default: "data"
#' @param ... parameters passing to alra()
#'
#' @importFrom glue glue
#' @importFrom Seurat GetAssayData SetAssayData
#'
#' @export
#' @return an updated Seurat object with the ALRA imputed data in the "data" 
#' slot of the ALRA assay slot
#'
alraSeurat2 <- function(obj, assay.use = "RNA", slot.use = "data", ...) {
  seurat_data <- GetAssayData(object = obj, 
                              assay = assay.use, 
                              slot = slot.use) %>% 
    as.matrix() %>%
    t()
  alra_data <- alra(seurat_data, 
                    ...) %>% 
    '[['(3) %>%
    t()
  colnames(alra_data) <- colnames(obj)
  rownames(alra_data) <- rownames(obj)
  alra_data <- Matrix(data_alra, sparse = T)
  obj[["alra"]] <- CreateAssayObject(data = alra_data)

  return(obj)
}
