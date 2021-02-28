options("experssion" = 500000)

characteristic_direction<- function(data, class, genes) {
  library(GeoDE)
  class <- factor(data.matrix(class))
  uniV <- GeoDE::chdirAnalysis(data, class,gammas=1.0)
  return_res <- as.data.frame(uniV$results)

  return(return_res)
}
test<-function(){
  
  
  
  class <- rbind(c(1,1,1,2,2,2))
  genes <- rbind(c('g1','g2','g3','g4','g5','g6','g7','g8','g9','g10','g11','g12'))
  data <- cbind(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1),
                c(1,1,1,1,1,1,1,1,1,1,1,1,1,1),
                c(2,1,1,1,1,1,1,2,1,1,1,1,1,1),
                c(3,1,1,1,2,2,2,3,1,1,1,2,2,2),
                c(4,10,1,20,2,2,10,10,1,20,2,2,20),
                c(5,10,1,20,2,2,206,10,1,20,2,2,55),
                c(6,10,1,20,2,2,55,4,10,1,20,2,2,100)
                )
  
  do_data = as.data.frame(data)

  res <- chdirAnalysis(do_data,class
                ,CalculateSig=TRUE,nnull=10)
  res
}


