options("experssion" = 500000)

calculate_limma <- function(data, class) {
  library(limma)
  data_mtx = data.matrix(data)
  class_mtx = data.matrix(class)
  design <- factor(class_mtx)
  design_matrix = model.matrix(~ design)
  fit <- lmFit(data_mtx, design_matrix)
  fit <- eBayes(fit)
  return(topTable(fit,sort="none",n=Inf))
}

