foo <- function(data, just_msg = TRUE){
  
  out <- lapply(split(data, data$study.name), function(x){
    pos_constant <- length(unique(x$post)) == 1L
    if (!pos_constant) {
      lapply(split(x, x$outcome), function(x_sub){
        mps <- x_sub[x_sub$control==TRUE,"mpre"]
        sps <- x_sub[x_sub$control==TRUE,"sdpre"]
        mps_constant <- length(unique(mps)) %in% c(1L,0L)
        sps_constant <- length(unique(sps)) %in% c(1L,0L)
        r <- !mps_constant || !sps_constant
        if (r) {
          
          if(just_msg) {
            
            message(sprintf("Error: 'mpre' in '%s' for 'control' and 'outcome' rows are wrongly coded.", x[,'study.name'][1]))
            
          } else {
            
            stop(sprintf("'mpre' in '%s' for 'control' and 'outcome' rows are wrongly coded.", x[,'study.name'][1]), call. = FALSE)
            
          }
        } else if(just_msg) { 
          
          cat(sprintf("OK: No multi-outcome coding issues in '%s' detected.\n",x[,"study.name"][1])) 
          
        } 
        
        return(r)
      })
    } else if(just_msg) { 
      
      cat(sprintf("OK: No multi-outcome coding issues in '%s' detected.\n",x[,"study.name"][1])) 
      return(FALSE)
      
    }
  })

  if(just_msg) invisible(any(unlist(out)))
}

### EXAMPLE OF USE: ------------------------------------------------
dat <- read.csv("https://raw.githubusercontent.com/rnorouzian/m2/main/dat.csv")

foo(dat, just_msg = T)