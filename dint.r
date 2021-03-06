
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

trim_ <- function(X){
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
  
  rm.allcolNA(rm.allrowNA(X))  

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

drop.col <- function(data, what){
  
    i1 <- !names(data) %in% what
    setNames(data[i1], names(data)[i1])
}           

#================================================================================================================================
           
           
handle_mpre_errors <- function(data, just_msg = TRUE){
  
  out <- lapply(split(data, data$study.name), function(x){
    pos_constant <- length(unique(x$post)) == 1L
    if (!pos_constant) {
      r_lst <- lapply(split(x, x$outcome), function(x_sub){
        mps <- x_sub[x_sub$control==TRUE,"mpre"]
        sps <- x_sub[x_sub$control==TRUE,"sdpre"]
        mps_constant <- length(unique(mps)) %in% c(1L,0L)
        sps_constant <- length(unique(sps)) %in% c(1L,0L)
        r <- !mps_constant || !sps_constant
        
        return(r)
      })
      r <- any(unlist(r_lst))
      if (r) {
        
        if(just_msg) {
          
          message(sprintf("Error: 'mpre' in '%s' for 'control' and 'outcome' rows are wrongly coded.", x[,'study.name'][1]))
          
        } else {
          
          stop(sprintf("'mpre' in '%s' for 'control' and 'outcome' rows are wrongly coded.", x[,'study.name'][1]), call. = FALSE)
          
        }
      } else if(just_msg) { 
        
        cat(sprintf("OK: No multi-outcome coding issues in '%s' detected.\n",x[,"study.name"][1])) 
        
      } 

    } else if(just_msg) { 
      
      cat(sprintf("OK: No multi-outcome coding issues in '%s' detected.\n",x[,"study.name"][1])) 
      return(FALSE)
      
    }
  })
  
  if(just_msg) invisible(any(unlist(out)))
}
           
           
#================================================================================================================================
           
handle_dint_errors <- function(res, just_msg = FALSE, smd = FALSE) {
  
  dat_obj <- data.frame(control = sapply(res[[1]][[1]], nrow), treatment = sapply(res[[2]][[1]], nrow))
  
  if(any(dat_obj == 0)){
    if(!any(rowSums(dat_obj == 0) == 1)){
      
      message("Note: ",dQuote(toString(names(res$tlist))),
              " missing some 'post','outcome' or 'control' row(s).")
      
      if(!just_msg){ 
        return(
          
          setNames(lapply(seq_along(res), function(i) Filter(NROW, res[[i]][[1]])), names(res))
        )}
      
    } else {  
      
      if(!just_msg){   
        
        stop(dQuote(toString(names(res$tlist)))," has 'post','outcome','control' wrongly coded for its 'control' & 'treatment' rows.", call. = FALSE)
        
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
        
        stop(dQuote(toString(names(res$tlist)))," has 'post','outcome','control' wrongly coded for its 'control' & 'treatment' rows.", call. = FALSE)
        
      } else {
        
        message("Error: ",dQuote(toString(names(res$tlist)))," has 'post','outcome','control' wrongly coded for its 'control' & 'treatment' rows.")  
      }
      
    } else {
      
      if(!just_msg){ 
        return(
          
          setNames(lapply(seq_along(res), function(i) Filter(NROW, res[[i]][[1]])), names(res))
        ) 
      } else {
        
        if(!smd){
        cat(paste("OK: No dint coding issues in",dQuote(toString(names(res$tlist))),"detected.\n"))
        } else {
          
        cat(paste("OK: No 'time', 'outcome' coding issues in",dQuote(toString(names(res$tlist))),"detected.\n"))
          
        }
      }
    }
  }
}


#=============================================================================================================================

ctlist_maker <- function(m, just_msg = FALSE, smd = FALSE){
  
  input <- setNames(lapply(m, function(i) 
    rev(expand.grid(outcome = seq_len(max(i$outcome, na.rm = TRUE)), 
                    post = seq_len(max(i$post, na.rm = TRUE))))), names(m))
  
  res <- setNames(lapply(1:0, function(i) lapply(input, function(inp) Map(function(p, o)
    do.call(rbind, lapply(m, function(x)
      x[x$control == i & x$post == p & x$outcome == o, , drop = FALSE])),
    inp$post, inp$outcome))), c("clist", "tlist"))
  
  handle_dint_errors(res, just_msg = just_msg, smd = smd)
}

#=============================================================================================================================


d.prepos <- function(d.pair = NA, study.name = NA, group.name = NA, n = NA, mdif.pair = NA, 
                     stder.pair = NA, mpre = NA, mpos = NA, sdpre = NA, sdpos = NA, r.prepos = NA, 
                     rev.sign = FALSE, rev.group = FALSE, autoreg = FALSE, t.pair = NA, 
                     df.pair = NA, sdif = NA, post = NA, control = NA, outcome = NA, ...) 
{
  
  rev.group <- ifelse(is.na(rev.group), FALSE, rev.group)
  autoreg <- ifelse(is.na(autoreg), FALSE, autoreg)
  rev.sign <- ifelse(is.na(rev.sign), FALSE, rev.sign)
    
  if(anyNA(post) || anyNA(outcome)) stop("'post' or 'outcome' missing in the EXCEL coding sheet.", call. = FALSE)
  
  d <- d.pair
  df <- df.pair
  mdif <- mdif.pair
  stder <- stder.pair
  r <- r.prepos
  
  r <- ifelse(autoreg == TRUE, autoreg(max(post, na.rm = TRUE), r)[,1][-1][post], r)
  
  n <- ifelse(!is.na(n), n, ifelse(is.na(n) & !is.na(df), df + 1, NA))
  mdif <- ifelse(!is.na(mdif), mdif, ifelse(!is.na(mpre) & is.na(mdif), mpos - mpre, NA))
  t.pair <- ifelse(!is.na(t.pair), t.pair, ifelse(is.na(t.pair) & !is.na(mdif) & !is.na(stder), mdif/stder, NA))
  d <- ifelse(!is.na(d), d, ifelse(!is.na(t.pair) & !is.na(n), t2d(t.pair, n), ifelse(!is.na(mdif) & !is.na(sdif), mdif/sdif, NA)))
  sdif <- ifelse(is.na(sdif), sdif(sdpre = sdpre, sdpos = sdpos, t.pair = t.pair, r = r, n = n, mpos = mpos, mpre = mpre), sdif)
  r <- ifelse(is.na(r), rdif(n = n, mpre = mpre, mpos = mpos, t.pair = t.pair, sdpre = sdpre, sdpos = sdpos, sdif = sdif), r)
  r <- ifelse(is.na(r) & is.na(sdif), .6, r)
  sdif <- sdif(sdpre = sdpre, sdpos = sdpos, t.pair = t.pair, r = r, n = n, mpos = mpos, mpre = mpre)
  d <- ifelse(!is.na(mdif) & is.na(d) & !is.na(sdif), mdif/sdif, d)*cfactor(n-1)
  d <- ifelse(rev.group, -d, d)
  
  out <- data.frame(d, group.name, n, sdif, r, rev.sign, post, control, outcome, ...)
  
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


handle_prepos_errors <- function(data, ar, dot.names, just_msg = TRUE, smd = FALSE){
  
  L <- split(data, data$study.name)  
  
  f1 <- function(number, smd = FALSE){
    
    ns <- names(L)[number]
    
    z <- try(d_prepo(L[[number]], ar, dot.names), silent = TRUE)
    
    
    if(inherits(z, "try-error")) { 
      
      Attrib <- TRUE
      
      if(!smd){
      message("Error: pre-post coding issues in ", toString(dQuote(ns)), " detected. Check all the columns ('n','mpre'...).")
      } else {
        
      message("Error: coding issues in ", toString(dQuote(ns)), " detected. Check all the columns ('n','mpre'...).")
      }
        
    } else {
      
      Attrib <- FALSE
      
      if(!smd){
      cat(paste("Ok: No pre-post coding issues in", toString(dQuote(ns)), "detected.\n"))
      } else { 
        
      cat(paste("Ok: No descriptive coding issues in", toString(dQuote(ns)), "detected.\n"))
        
        }
    }
    
    y <- NA
    attr(y, "message") <- Attrib
    return(y)
  }
  
  
  f2 <- function(number, smd = FALSE){
    
    ns <- names(L)[number]
    
    z <- try(d_prepo(L[[number]], ar, dot.names), silent = TRUE)
    
    if(inherits(z, "try-error")) { 
      
      if(!smd){
      stop("pre-post coding issues in ", toString(dQuote(ns)), " detected. Check all the columns ('n','mpre'...).", call. = FALSE)
      } else {
        
      stop("coding issues in ", toString(dQuote(ns)), " detected. Check all the columns ('n','mpre'...).", call. = FALSE)
      }
        }    
  }  
  
  if(just_msg){
    
    res <- invisible(lapply(seq_along(L), function(i) f1(i, smd = smd)))
    
    any(sapply(res, attr, "message"))
  }  
  else {
    
    invisible(lapply(seq_along(L), function(i) f2(i, smd = smd)))
  } 
  
}


#==========================================================================================================================================

                     
check_sheet <- function(data, m, ar, dot.names, smd = FALSE){
  
  message(paste("\nError analysis of multi-outcome study coding:\n"))
  
  bad_1st_check <- handle_mpre_errors(data)
  
  if(!bad_1st_check){
    
    message(paste("\nError analysis of pre-post effects coding:\n"))
    
    bad_2nd_check <- handle_prepos_errors(data, ar, dot.names, smd = smd)
    
  } else { bad_2nd_check <- TRUE
  
  message("\nError analysis of stage 2 coding stopped due to the 'Error' found above.")}
  
  
  if(!bad_2nd_check){
    
    L <- get_d_prepos(data, m, ar, dot.names)  
    
    if(!smd){
    message("\nError analysis of dints effects coding:\n")    
    } else { message("\nError analysis of 'outcome' and 'time' coding:\n")  }
    invisible(lapply(names(L), function(i) ctlist_maker(L[i], just_msg = TRUE, smd = smd)))
  } else {
    
    message("\nError analysis of stage 3 stopped due to the 'Error' found above.")
  }
}

#==========================================================================================================================================

make_final_output <- function(data, m, ar, dot.names, smd = FALSE, smd_raw = TRUE){
  
  handle_mpre_errors(data, just_msg = FALSE)
  
  handle_prepos_errors(data, ar, dot.names, just_msg = FALSE, smd = smd) 
  
if(!smd){ 
  
  L <- get_d_prepos(data, m, ar, dot.names)
  
  res <- setNames(lapply(names(L), function(i) get_dint(L[i], dot.names)), names(L))
  
  DF <- do.call(rbind, c(Map(cbind, studyID = names(res), res), make.row.names = FALSE))
  
  out <- cbind(esID = seq_len(nrow(DF)), DF)
  
  out2 <- out[c("studyID","esID", setdiff(names(out), c("studyID","esID")))]
  
  drop.col(out2, "sign_.")
  
} else {
  
  smd_maker(data, m, ar, dot.names, raw = smd_raw)
  
   }
}


#==========================================================================================================================================

check_data_ <- function(data, ar, smd = FALSE){
  
  data <- rm.allrowNA(trim_(data))
  idx <- ar %in% names(data)
  if(!all(idx)) stop("Column(s) ",toString(dQuote(ar[!idx])), " missing.", call. = FALSE)
  data <- transform(data, post_id = post, outcome_id = outcome, n_outcome = ave(outcome, study.name, FUN = max),
                    control = ifelse(is.na(control), FALSE, control),
                    sign_. = ifelse(is.na(rev.sign), FALSE, rev.sign))
  
  if(smd & anyNA(data$group.name)) stop("For SMD effect sizes, 'group.name' must be coded for.", call. = FALSE)
  return(data)
}


#==========================================================================================================================================                         
                         
Type_of_Input_ <- c(
  "numeric", "numeric/character",
  "numeric/character", "numeric",
  "numeric", "numeric","numeric",
  "numeric","numeric", "numeric",
  "numeric [0-1]", "logical","logical","logical",
  "numeric","numeric","numeric",
  "numeric","logical",
  "numeric", "numeric/character")        
                         
#============================================================================================================================================
                         
mask_ <- function(data, what, full = FALSE){
  
  data[] <- lapply(data, function(x) type.convert(as.character(x), as.is = TRUE))
  
  f1 <- function(x) as.numeric(factor(x, levels = unique(x)))
  f2 <- function(x) {
    temp <- substr(x, 1, 1)
    paste0(temp, ave(x, temp, FUN = function(y) match(y, unique(y))))
  }
  cols <- names(data)[sapply(data, is.numeric)]
  num.cols <- cols[cols %in% what]
  cols <- names(data)[sapply(data, is.character)]
  char.cols <- cols[cols %in% what]
  
  if(!full){
    
    if(length(num.cols))  data[num.cols] <- lapply(data[num.cols], f1)
    if(length(char.cols)) data[char.cols] <- lapply(data[char.cols], f2)
    
  }else{
    
    data[what] <- lapply(data[what], f1)
  }
  return(data)
}
                         
#================================================================================
                     
qcohen <- function(p, dbase = 0, n1, n2 = NA, lower.tail = TRUE, log.p = FALSE){
  
  q <- Vectorize(function(p, dbase, n1, n2, lower.tail, log.p){
    
    N <- ifelse(is.na(n2), n1, (n1 * n2)/(n1 + n2))
    df <- ifelse(is.na(n2), n1 - 1, (n1 + n2) - 2)
    ncp <- dbase*sqrt(N)
    
    qt(p, df, ncp, lower.tail = lower.tail, log.p = log.p)/sqrt(N)
  })
  q(p = p, dbase = dbase, n1 = n1, n2 = n2, lower.tail = lower.tail, log.p = log.p)
}

    
#=======================================================================================================================================


pcohen <- function(q, dbase = 0, n1, n2 = NA, lower.tail = TRUE, log.p = FALSE){
  
  p <- Vectorize(function(q, dbase, n1, n2, lower.tail, log.p){
  
  N <- ifelse(is.na(n2), n1, (n1 * n2)/(n1 + n2))
  df <- ifelse(is.na(n2), n1 - 1, (n1 + n2) - 2)
  ncp <- dbase*sqrt(N)
  
  pt(q*sqrt(N), df, ncp, lower.tail = lower.tail, log.p = log.p)
})
p(q = q, dbase = dbase, n1 = n1, n2 = n2, lower.tail = lower.tail, log.p = log.p)
}                     
                     
#=========================================================================================================================================
                     
smd_maker <- function(data, m, ar, dot.names, raw = FALSE){
  
vec.names <- c('mT','sdT','nT','mC','sdC','nC')

f <- function(m2, dot.names, vec.names){
  
  aa <- ctlist_maker(m2)
  
  b1 <- lapply(aa, function(x)  {
    y <- do.call(rbind,  x);
    y[order(y$group.name), c("mpre", "sdpre", "n", dot.names)] })
  
  cc1 <- do.call(cbind,rev(b1))
  cc1 <- cc1[grep(sprintf("clist.(%s)", paste0(dot.names, collapse="|")), names(cc1), invert = TRUE)]
  i1 <- grepl(sprintf("(%s)", paste0(dot.names, collapse="|")), names(cc1))
  cc1 <- cbind(cc1[!i1], cc1[i1])
  names(cc1)[1:6] <- vec.names
  cc_1 <- cc1[!duplicated(cc1[vec.names]),]
  cc_1 <- cbind(cc_1, time = 0)
  
  b2 <- lapply(aa, function(x)  {
    y <- do.call(rbind,  x);
    y[order(y$group.name), c("mpos", "sdpos", "n", dot.names)] })
  
  cc2 <- do.call(cbind,rev(b2))
  cc2 <- cc2[grep(sprintf("clist.(%s)", paste0(dot.names, collapse="|")), names(cc2), invert = TRUE)]
  i1 <- grepl(sprintf("(%s)", paste0(dot.names, collapse="|")), names(cc2))
  cc2 <- cbind(cc2[!i1], cc2[i1])
  names(cc2)[1:6] <- vec.names
  cc_2 <- cc2[!duplicated(cc2[vec.names]),]
  cc_2 <- transform(cc_2, time = tlist.post_id)
  
  out <- rbind(cc_1,cc_2, make.row.names = FALSE)
  
  out <- cbind(studyID = rep(names(m2), nrow(out)), out)
  
  fin <- setNames(out, sub("^tlist\\.", "", names(out)))
  
  fin[order(fin$time),]
}

LIST <- lapply(names(m), function(i) f(m2=m[i],dot.names=dot.names,vec.names=vec.names))

DF <- do.call(rbind, c(LIST, make.row.names = FALSE))

es_data <- cbind(esID = seq_len(nrow(DF)), DF)

es <- metafor::escalc("SMD", m1i = mT, m2i = mC, sd1i = sdT, sd2i = sdC, n1i = nT, n2i = nC, data = es_data)

output <- transform(es, yi = ifelse(sign_. == TRUE, -yi, yi)) 

output <- drop.col(output, if(!raw) c(vec.names,"sign_.","post_id") else "post_id")

output[c("studyID","esID","yi","vi", setdiff(names(output), c("studyID","esID","yi","vi")))]
}      
                   
#=========================================================================================================================================   
                   
make_final_output <- function(data, m, ar, dot.names, smd = FALSE, smd_raw = TRUE){
  
  handle_mpre_errors(data, just_msg = FALSE)
  
  handle_prepos_errors(data, ar, dot.names, just_msg = FALSE, smd = smd) 
  
if(!smd){ 
  
  L <- get_d_prepos(data, m, ar, dot.names)
  
  res <- setNames(lapply(names(L), function(i) get_dint(L[i], dot.names)), names(L))
  
  DF <- do.call(rbind, c(Map(cbind, studyID = names(res), res), make.row.names = FALSE))
  
  out <- cbind(esID = seq_len(nrow(DF)), DF)
  
  out2 <- out[c("studyID","esID", setdiff(names(out), c("studyID","esID")))]
  
  drop.col(out2, "sign_.")
  
} else {
  
  smd_maker(data, m, ar, dot.names, raw = smd_raw)
  
   }
}
                         
#=========================================================================================================================================                         
                         
dint <- function(data, check_sheet = FALSE, smd = FALSE, smd_raw = TRUE){
  
  ar <- formalArgs(d.prepos)[-c(3,21)] # [-c(3,21:22)] # 
  
  data <- check_data_(data, ar, smd = smd)
  
  m <- split(data, data$study.name)
  
  all_names_ <- names(data)
  
  dot.names <- all_names_[!all_names_ %in% ar]
  
  if(check_sheet){
    
    check_sheet(data, m, ar, dot.names, smd = smd)
    
  } else {
    
    make_final_output(data, m, ar, dot.names, smd = smd, smd_raw = smd_raw)
  }
}                       
