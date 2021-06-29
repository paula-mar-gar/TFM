plot_intervals <- function(n_total,n_interval,n_div,name_list){
  m <- 0
  list_tmp <- list()
  b <- n_interval-1
  for (a in 1:b){
    n <- m + 1
    m <- a * n_div
    list_tmp[[a]] <- c(n:m)
    cat(a,") ",n,":",m,"\n", sep = "")
  }
  n <- m + 1
  list_tmp[[n_interval]] <- c(n:n_total)
  cat(n_interval,") ",n,":",n_total,sep = "")
  assign(name_list,list_tmp,envir = parent.frame())
  
}
