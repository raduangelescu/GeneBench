options("experssion" = 500000)

get_fake_data_mtx<- function() {
  set.seed(100)
  num_experiments = 6
  num_genes = 10
  #cols e experimente
  #rows e gene
  data<-matrix(rnorm(num_genes*num_experiments),ncol=num_experiments)
  
  dd<-sample(1:num_genes,size=num_genes)
  
  u<-matrix(2*rnorm(num_genes),ncol=num_experiments/2,nrow=num_genes)
  data[dd,(num_experiments/2+1):num_experiments]<-data[dd,(num_experiments/2+1):num_experiments]+u
  return(data)
}

calculate_sam <- function(data_mtx, class_mtx, gene_names) {
  library(DT)
  # source('http://bioconductor.org/biocLite.R') # Import biocLite() function into R environment biocLite('limma')
  library(samr)  # Load the library
  geneids <- as.character(1:nrow(data_mtx))
  
  data_param <- list(x = data_mtx, y = class_mtx, geneid=geneids,
               genenames=gene_names, logged2 = TRUE)
  samr.obj <- samr(data_param, resp.type = "Two class unpaired", assay.type=c("array"),testStatistic="wilcoxon",regression.method="ranks"
,nperms = 100)
  delta.table <- samr.compute.delta.table(samr.obj) 
  del <- 0
  siggenes.table<-samr.compute.siggenes.table(samr.obj, del, data_param, delta.table)
  data_frame_return = as.data.frame(do.call(rbind, siggenes.table))

  return(data_frame_return)
}

test_func <- function() {
  set.seed(100)
  num_experiments = 6
  num_genes = 10
  #cols e experimente
  #rows e gene
  data<-matrix(rnorm(num_genes*num_experiments),ncol=num_experiments)

  dd<-sample(1:num_genes,size=num_genes)
  
  u<-matrix(2*rnorm(num_genes),ncol=num_experiments/2,nrow=num_genes)
  data[dd,(num_experiments/2+1):num_experiments]<-data[dd,(num_experiments/2+1):num_experiments]+u
  
  y<-c(rep(1,num_experiments/2),rep(2,num_experiments/2))
  
  gene_names <- rbind(c('g1','g2','g3','g4','g5','g6','g7','g8','g9','g10'))
  result =calculate_sam(data, y, gene_names) 
  
  return(result)
}
#test_func()

