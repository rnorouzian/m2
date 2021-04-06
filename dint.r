
Break = "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n"

notice = "   \"dint\", An R program for calculation of SMCD effect size \n     from complex longitudinal Second Language Research.
    Copyright (C) 2020-present Reza Norouzian, rnorouzian@gmail.com\n"

message(Break, notice, Break)


#================================================================================================================================


autoreg <- function(steps, r){
  
  steps <- if(length(steps) == 1) steps+1 else length(steps)+1
  x <- diag(steps)
  r <- data.frame(r^abs(row(x)-col(x)))
  rownames(r) <- colnames(r) <- c("pre", paste0("post", 1:(steps-1)))
  return(r)
} 


#===============================================================================================================================

t2d <- function(t, n1, n2 = NA){
  N <- ifelse(is.na(n2), n1, (n1 * n2)/(n1 + n2))
  t/sqrt(N)
}

#===============================================================================================================================

d2t <- function(d, n1, n2 = NA){
  N <- ifelse(is.na(n2), n1, (n1 * n2)/(n1 + n2))
  d*sqrt(N)
}

#===============================================================================================================================

trim <- function(X){
  X <- setNames(X, trimws(names(X)))
  y <- sapply(names(X), function(x) is.character(as.vector(X[[x]])))
  X[y] <- lapply(X[y], trimws)
  return(X)
}

#===============================================================================================================================


rm.allrowNA <- function(X) { 
  
  if(inherits(X, "list")){
    
    lapply(X, function(i) i[rowSums(is.na(i) | i == "") != ncol(i), , drop = FALSE])
    
  } else { X[rowSums(is.na(X) | X == "") != ncol(X), , drop = FALSE] }
}

#===============================================================================================================================

rm.allcolNA <- function(X) { 
  
  if(inherits(X, "list")){
    
    lapply(X, function(i) i[, colSums(is.na(i) | i == "") != nrow(i), drop = FALSE])
    
  } else { X[, colSums(is.na(X) | X == "") != nrow(X), drop = FALSE] }
}

#===============================================================================================================================

rm.colrowNA <- function(X){
  
  r <- rm.allrowNA(X)
  rm.allcolNA(r)  
  
}                                      

#================================================================================================================================

sdif <- function(n = NA, mpre = NA, mpos = NA, sdpre = NA, sdpos = NA, r = NA, t.pair = NA, df = NA,
                 sdp = NA){
  
  n <- ifelse(!is.na(n), n, ifelse(is.na(n) & !is.na(df), df + 1, NA))
  
  ifelse(!is.na(r) & !is.na(sdpre) & !is.na(sdpos), sqrt(sdpre^2+sdpos^2-2*r*sdpre*sdpos),
         ifelse(!is.na(n) & is.na(r) & !is.na(t.pair) & !is.na(mpre) & !is.na(mpos), sqrt((n*(mpos - mpre)^2)/t.pair^2), 
                ifelse(!is.na(r) & !is.na(sdp), sqrt(2*sdp^2*(1-r)), NA)))
}

#===============================================================================================================================

rdif <- function(n = NA, mpre = NA, mpos = NA, sdpre = NA, sdpos = NA, t.pair = NA, df = NA, sdif = NA, 
                 sdp = NA) {
  
  sdif <- ifelse(is.na(sdif), sdif(sdpre = sdpre, sdpos = sdpos, t.pair = t.pair, n = n, mpos = mpos, mpre = mpre, df = df, sdp = sdp), sdif)
  
  ifelse(!is.na(sdif) & is.na(sdp) & !is.na(sdpre) & !is.na(sdpos), (sdpre^2 + sdpos^2 - sdif^2)/(2*sdpre*sdpos), 
         ifelse(!is.na(sdp) & !is.na(sdif), 1 - (sdif^2/(2*sdp^2)), NA))
}

#===============================================================================================================================

cfactor <- function(df) exp(lgamma(df/2)-log(sqrt(df/2)) - lgamma((df-1)/2))

#===============================================================================================================================

reget <- function(List, what){
  
  s <- substitute(what)  
  
  if(!inherits(List, "list")) List <- list(List)  
  
  h <- lapply(List, function(x) do.call("subset", list(x, s)))
  
  res <- Filter(NROW, h)
  
  if(length(res) == 0) NULL else res
}

#===============================================================================================================================

dit <- Vectorize(function(dppc, dppt, nc, nt, rev.sign = FALSE){
  
  a <- dppc
  b <- dppt
  
  din <- b - a 
  
  test <- if(!rev.sign || rev.sign & b < 0 & a < 0 & abs(b) < abs(a)) FALSE else TRUE
  
  vc <- 1/nc + (1 - ( (nc-3) / ((nc-1)*cfactor(nc-1)^2) )) * dppc^2
  vt <- 1/nt + (1 - ( (nt-3) / ((nt-1)*cfactor(nt-1)^2) )) * dppt^2
  
  vi <- vt + vc
  
  din <- if(test) -din else din
  return(c(dint = din, vi = vi))
})  


#================================================================================================================================


handle_dint_errors <- function(res, just_msg = FALSE) {
  
  dat_obj <- data.frame(control = sapply(res[[1]][[1]], nrow), treatment = sapply(res[[2]][[1]], nrow))
  
  if(any(dat_obj == 0)){
    if(!any(rowSums(dat_obj == 0) == 1)){
      
      message("Note: ",dQuote(toString(names(res$tlist))),
                    " missing some posttests/outcomes/control groups.Check 'post','outcome','control' columns.")
      
      if(!just_msg){ 
        return(
          
          setNames(lapply(seq_along(res), function(i) Filter(NROW, res[[i]][[1]])), names(res))
        )}
      
    } else {  
      
      if(!just_msg){   
        
        stop(paste(dQuote(toString(names(res$tlist))),"has 'post','outcome','control' wrongly coded for its 'control' & 'treatment' rows."), call. = FALSE)
        
      } else {
        
        message("Error: ",dQuote(toString(names(res$tlist)))," has 'post','outcome','control' wrongly coded for its 'control' & 'treatment' rows.")  
      }
    }
  } else {
    
    v1 <- sapply(names(dat_obj), function(i) length(unique(dat_obj[[i]])))
    i1 <- names(which(v1 != 1))
    
    if(length(i1) == 1) {
      
      message("Note: ",dQuote(toString(names(res$tlist)))," missing some posttests/outcomes/controls perhaps in its ",i1, " row(s).")
      
      if(!just_msg){ 
        return(
          
          setNames(lapply(seq_along(res), function(i) Filter(NROW, res[[i]][[1]])), names(res))
        )}
      
    } else if(length(i1) > 1) {
      
      if(!just_msg){   
        
        stop(dQuote(toString(names(res$tlist)))," has posttests/outcomes/controls wrongly coded for its 'control' & 'treatment' rows.", call. = FALSE)
        
      } else {
        
        message("Error: ",dQuote(toString(names(res$tlist)))," has posttests/outcomes/controls wrongly coded for its 'control' & 'treatment' rows.")  
      }
      
    } else {
      
      if(!just_msg){ 
        return(
          
          setNames(lapply(seq_along(res), function(i) Filter(NROW, res[[i]][[1]])), names(res))
        ) 
      } else {
        
        cat(paste("OK: No dint coding issues in",dQuote(toString(names(res$tlist))),"detected.\n"))
      }
    }
  }
}


#=============================================================================================================================

ctlist_maker <- function(m, just_msg = FALSE){
  
  input <- setNames(lapply(m, function(i) 
    rev(expand.grid(outcome = seq_len(max(i$outcome, na.rm = TRUE)), 
                    post = seq_len(max(i$post, na.rm = TRUE))))), names(m))
  
  res <- setNames(lapply(1:0, function(i) lapply(input, function(inp) Map(function(p, o)
    do.call(rbind, lapply(m, function(x)
      x[x$control == i & x$post == p & x$outcome == o, , drop = FALSE])),
    inp$post, inp$outcome))), c("clist", "tlist"))
  
  handle_dint_errors(res, just_msg = just_msg)
}

#=============================================================================================================================


d.prepos <- function(d.pair = NA, study.name = NA, group.name = NA, n = NA, mdif.pair = NA, 
                     stder.pair = NA, mpre = NA, mpos = NA, sdpre = NA, sdpos = NA, r.prepos = NA, 
                     rev.sign = FALSE, rev.group = FALSE, autoreg = FALSE, t.pair = NA, 
                     df.pair = NA, sdif = NA, post = NA, control = NA, outcome = NA, time = NA, ...) 
{
  
  rev.sign <- ifelse(is.na(rev.sign), FALSE, rev.sign)
  rev.group <- ifelse(is.na(rev.group), FALSE, rev.group)
  autoreg <- ifelse(is.na(autoreg), FALSE, autoreg)
  control <- ifelse(is.na(control), FALSE, control)
  
  if(anyNA(post) || anyNA(outcome)) stop("'post' or 'outcome' missing in the EXCEL coding sheet.", call. = FALSE)
  
  d <- d.pair
  df <- df.pair
  mdif <- mdif.pair
  stder <- stder.pair
  r <- r.prepos
  
  r <- ifelse(autoreg == TRUE, autoreg(max(post, na.rm = TRUE), r)[,1][-1][post], r)
  
  n <- ifelse(!is.na(n), n, ifelse(is.na(n) & !is.na(df), df + 1, NA))
  mdif <- ifelse(!is.na(mdif), mdif, ifelse(!is.na(mpre) & !is.na(mpre) & is.na(mdif), mpos - mpre, NA))
  t.pair <- ifelse(!is.na(t.pair), t.pair, ifelse(is.na(t.pair) & !is.na(mdif) & !is.na(stder), mdif/stder, NA))
  d <- ifelse(!is.na(d), d, ifelse(!is.na(t.pair) & !is.na(n), t2d(t.pair, n), ifelse(!is.na(mdif) & !is.na(sdif), mdif/sdif, NA)))
  sdif <- ifelse(is.na(sdif), sdif(sdpre = sdpre, sdpos = sdpos, t.pair = t.pair, r = r, n = n, mpos = mpos, mpre = mpre), sdif)
  r <- ifelse(is.na(r), rdif(n = n, mpre = mpre, mpos = mpos, t.pair = t.pair, sdpre = sdpre, sdpos = sdpos, sdif = sdif), r)
  r <- ifelse(is.na(r) & is.na(sdif), .6, r)
  sdif <- sdif(sdpre = sdpre, sdpos = sdpos, t.pair = t.pair, r = r, n = n, mpos = mpos, mpre = mpre)
  d <- ifelse(!is.na(mdif) & is.na(d) & !is.na(sdif), mdif/sdif, d)*cfactor(n-1)
  d <- ifelse(rev.group, -d, d)
  
  out <- data.frame(d, group.name, n, sdif, r, rev.sign, post, control, outcome, time, ...)
  
  if(any(is.na(out$d))) stop("insufficient info. to calculate effect size(s).", call. = FALSE)
  
  return(out) 
}


#==========================================================================================================================================


pair2_ <- function(j, k){
  x1 <- expand.grid(d1 = j$d, d2 = k$d)
  split(as.matrix(x1), row(x1))
}

#==========================================================================================================================================


dinter2 <- function(clist, tlist, dot.names){
  
  nc <- clist$n 
  nt <- tlist$n
  dp_ <- pair2_(clist, tlist)
  dppc <- sapply(seq_along(dp_), function(i) dp_[[i]][1])
  dppt <- sapply(seq_along(dp_), function(i) dp_[[i]][2])
  rev.sign <- tlist$rev.sign
  mods <- tlist[dot.names]
  
  data.frame(t(dit(dppc = dppc, dppt = dppt, nc = nc, nt = nt, rev.sign = rev.sign)), mods)
}

#==========================================================================================================================================

get_dint <- function(m, dot.names){
  
  ctlists <- ctlist_maker(m)
  
  do.call(rbind, lapply(seq_len(unique(lengths(ctlists))), function(i) 
    dinter2(ctlists$clist[[i]], ctlists$tlist[[i]], dot.names)))
}

#==========================================================================================================================================

get_d_prepos <- function(data, m, ar, dot.names){
  
  args <- lapply(m, function(x) unclass(x[c(ar, dot.names)]))
  
  argsT <- setNames(lapply(names(args[[1]]), function(i) lapply(args, `[[`, i)), names(args[[1]]))
  
  do.call(Map, c(f = d.prepos, argsT))
}


#==========================================================================================================================================


d_prepo <- function(data = NULL, ar, dot.names)
{
  
  args <- unclass(data[c(ar, dot.names)])
  
  argsT <- setNames(args, names(args))
  
  do.call(d.prepos, argsT)
}


#==========================================================================================================================================


handle_prepos_errors <- function(data, ar, dot.names, check_sheet = TRUE){
  
  L <- split(data, data$study.name)  
  
  
  f1 <- function(number){
    
    ns <- names(L)[number]
    
    z <- try(d_prepo(L[[number]], ar, dot.names), silent = TRUE)
    
    
    if(inherits(z, "try-error")) { 
      
      Attrib <- TRUE
    message("Error: pre-post coding issues in ", toString(dQuote(ns)), " detected. Check descriptive columns ('n','mpre'...).")
    
      } else {
        
      Attrib <- FALSE
    cat(paste("Ok: No pre-post coding issues in", toString(dQuote(ns)), "detected.\n"))
      }
    
    y <- NA
    attr(y, "message") <- Attrib
    return(y)
  }
  
  
  f2 <- function(number){
    
    ns <- names(L)[number]
    
    z <- try(d_prepo(L[[number]], ar, dot.names), silent = TRUE)
    
    if(inherits(z, "try-error")) { 
      
      stop("pre-post coding issues in ", toString(dQuote(ns)), " detected. Check descriptive columns ('n','mpre'...).", call. = FALSE)
    }    
}  
  
  
  if(check_sheet){
  
  res <- invisible(lapply(seq_along(L), function(i) f1(i)))
  
  any(sapply(res, attr, "message"))
 }  
  else {
    
    invisible(lapply(seq_along(L), function(i) f2(i)))
  } 
  
}



#==========================================================================================================================================

check_sheet <- function(data, m, ar, dot.names){
  
  message(paste("\nError analysis of pre-post effects coding:\n"))
  
  bad_first_check <- handle_prepos_errors(data, ar, dot.names)
  
  if(!bad_first_check){
  
  L <- get_d_prepos(data, m, ar, dot.names)  
    
  message(paste("\nError analysis of dints effects coding:\n"))    
  
  invisible(lapply(names(L), function(i) ctlist_maker(L[i], just_msg = TRUE)))
  } else {
    
    message(paste("\nError analysis of dints effects coding stopped due to the 'Error' found above."))
  }
}


#==========================================================================================================================================

make_final_output <- function(L, dot.names){
  
  res <- setNames(lapply(names(L), function(i) get_dint(L[i], dot.names)), names(L))
  
  DF <- do.call(rbind, c(Map(cbind, studyID = names(res), res), make.row.names = FALSE))
  
  out <- cbind(esID = seq_len(nrow(DF)), DF)
  
  out[c("studyID","esID", setdiff(names(out), c("studyID","esID")))]
}


#==========================================================================================================================================


check_data_ <- function(data){
  
  data <- rm.allrowNA(trim(data))
  check <- "study.name" %in% names(data)
  if(!check) stop("Add a new column titled 'study.name'.", call. = FALSE)
  return(data)
}

#==========================================================================================================================================


dint <- function(data = NULL, check_sheet = FALSE){
  
  data <- check_data_(data)
  
  m <- split(data, data$study.name)
  
  if(is.null(reget(m, control))) stop("Required 'control/comparison' group not found.", call. = FALSE)
    
  ar <- formalArgs(d.prepos)[-c(3,21:22)]
  
  all_names_ <- names(data)
  
  dot.names <- all_names_[!all_names_ %in% ar]
  
  if(check_sheet){
    
   check_sheet(data, m, ar, dot.names)
  
  } else {
    
    handle_prepos_errors(data, ar, dot.names, check_sheet = FALSE)
    
    L <- get_d_prepos(data, m, ar, dot.names)
    
    make_final_output(L, dot.names)
  }
}
                         
