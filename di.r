Break = "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n"

notice = "   \"dint\", An R program for calculation of SMD effect size \n     from complex longitudinal Second Language Research.
    Copyright (C) 2020-present  Reza Norouzian, rnorouzian@gmail.com\n"

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

isFALSE <- function(x) identical(FALSE, x)              

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
           
rm.allrowNA2 <- function(X) { 
  
  if(inherits(X, "list")){
    
    lapply(X, function(i) if(NROW(i) != 0) Filter(NROW, i[rowSums(is.na(i) | i == "") != ncol(i), , drop = FALSE]) else Filter(NROW, i))
    
  } else { X[rowSums(is.na(X) | X == "") != ncol(X), , drop = FALSE] }
} 
           
#================================================================================================================================

roundi <- function(x, digits = 7){
  
  if(!inherits(x, "data.frame")) stop("Only used for a 'data.frame'.", call. = FALSE)
  
  num <- sapply(x, is.numeric)
  
  x[num] <- lapply(x[num], round, digits)
  
  return(x)
}                

#===============================================================================================================================

d_prepo <- function(data = NULL)
{
  
  if(is.null(reget(data, control))) stop("Required 'control' group not found.", call. = FALSE)
  
  ar <- formalArgs(d.prepos)[-c(21, 22)]
  
  dot.names <- names(data)[!names(data) %in% ar]
  
  args <- unclass(data[c(ar, dot.names)])
  
  argsT <- setNames(args, names(args))
  
  do.call(d.prepos, argsT)
}

#===============================================================================================================================

test.sheet <- function(data, metaset = FALSE){
  
  data <- rm.allrowNA(trim(data))
  
  check <- "study.name" %in% names(data)
  if(!check) stop("Add a new column named 'study.name'.", call. = FALSE)
  
  L <- split(data, data$study.name)         
  
  f <- function(number){
    
    ns <- names(L)[number]
    
    z <- try(d_prepo(L[[number]]), silent = TRUE)
    
    if(inherits(z, "try-error")) message("Error: pre-post coding problem in: *", toString(dQuote(ns)), "*") else message(if(metaset)"" else "Ok: ","No pre-post coding problem in ", toString(dQuote(ns)))
  }
  invisible(lapply(seq_along(L), function(i) f(i)))
}     

#===============================================================================================================================

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


get.uni <- function(data, what){
  
  data$study.name <- trimws(data$study.name)
  m <- split(data, data$study.name)
  m <- Filter(NROW, rm.allrowNA2(m)) 
  
  G <- substitute(what)
  E <- quote(x$x)
  E[[3]] <- G[[2]]
  G[[2]] <- E
  
  f <- sapply(m, function(x) sum(eval(G)) == nrow(x))
  
  h <- m[names(f)[f]]
  
  res <- Filter(NROW, h)
  
  if(length(res) == 0) NULL else res
}


#===============================================================================================================================


get.gen <- function(data, what){
  
  s <- substitute(what)  
  data$study.name <- trimws(data$study.name)
  m <- split(data, data$study.name)
  m <- Filter(NROW, rm.allrowNA2(m)) 
  
  h <- lapply(m, function(x) do.call("subset", list(x, s)))
  
  res <- Filter(NROW, h)
  
  if(length(res) == 0) NULL else res
}

#===============================================================================================================================

find.stud <- function(data, what, timevar = TRUE){
  
  s <- substitute(what)
  data$study.name <- trimws(data$study.name)
  if(!timevar) { unique(as.vector(subset(data, eval(s))$study.name))
    
  } else {
    chep <- sort(unique(na.omit(data$time)))
    G <- lapply(chep, function(x) bquote(.(s) & time == .(x)))
    setNames(lapply(seq_along(G), function(j) unique(as.vector(subset(data, eval(G[[j]]))$study.name))), as.character(chep))
  }
}       

#===============================================================================================================================

find.miss <- function(data, space = FALSE, all = FALSE){    
  
  res <- Filter(length, lapply(data, function(x) which(if(!space & !all) is.na(x) else if(space & !all) x == "" else is.na(x) | x == "")))
  if(length(res) == 0) NA else res
}                    

#===============================================================================================================================

mod.level <- function(data, what){
  
  sort(unique(na.omit(unlist(data[paste0(substitute(what))]))))
  
}

#===============================================================================================================================

mods.level <- function(data){
  
  f <- function(data, what) sort(unique(na.omit(unlist(data[what]))))
  
  ar <- c(formalArgs(d.prepos)[-c(21, 22)], c("dint", "SE", "id"))
  
  dot.names <- names(data)[!names(data) %in% ar]
  
  setNames(lapply(seq_along(dot.names), function(i) f(data = data, what = dot.names[i])), dot.names)
}                  


#===============================================================================================================================

pair_ <- function(j, k){
  lapply(seq_along(j), function(i) {
    x1 <- expand.grid(d1 = j[[i]]$d, d2 = k[[i]]$d);
    split(as.matrix(x1), row(x1))})
}                

#===============================================================================================================================


d.prepos <- function(d = NA, study.name = NA, group.name = NA, n = NA, mdif = NA, stder = NA, mpre = NA, mpos = NA, sdpre = NA, sdpos = NA, r = NA, rev.sign = FALSE, rev.group = FALSE, autoreg = FALSE, t.pair = NA, df = NA, sdif = NA, post = NA, control = NA, outcome = NA, time = NA, ...) 
{
  
  if(anyNA(control) || anyNA(post) || anyNA(outcome) || anyNA(time)) stop("'post', 'outcome', 'time', or 'control' missing in the EXCEL sheet.", call. = FALSE)
  
  rev.sign <- ifelse(is.na(rev.sign), FALSE, rev.sign)
  rev.group <- ifelse(is.na(rev.group), FALSE, rev.group)
  autoreg <- ifelse(is.na(autoreg), FALSE, autoreg)
  control <- ifelse(is.na(control), FALSE, control)
  
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
  
  out <- data.frame(d, n, sdif, r, rev.sign, post, control, outcome, time, ...)
  
  if(all(is.na(out$d))) stop("insufficient info. to calculate effect size(s).", call. = FALSE)
  
  return(out) 
}



#================================================================================================================================



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


dinter <- function(m, clist, tlist, pos, outcom){
    
if(!is.null(clist) & !is.null(tlist)){
  
nc <- m$n[m$control==TRUE & m$post == pos & m$outcome == outcom]
nt <- m$n[m$control==FALSE & m$post == pos & m$outcome == outcom]
dp_ <- pair_(clist, tlist)
dppc <- sapply(1:lengths(dp_), function(i) dp_[[1]][[i]][1])
dppt <- sapply(1:lengths(dp_), function(i) dp_[[1]][[i]][2])
rv <- m$rev.sign[m$control==FALSE & m$post == pos & m$outcom == outcom]

data.frame(t(dit(dppc = dppc, dppt = dppt, nc = nc, nt = nt, rev.sign = rv)))
 } else NULL
}
                  
                  
#================================================================================================================================                  
                  
                  
                  
dint <- function(data = NULL, impute_pre_post_r = FALSE)
{
  
  data <- rm.allrowNA(trim(data))
  check <- "study.name" %in% names(data)
  if(!check) stop("Add a new column named 'study.name'.", call. = FALSE) 
  
  
  m <- split(data, data$study.name)
  
  if(is.null(reget(m, control))) stop("Required 'control' group not found.", call. = FALSE)
  
  
  if(impute_pre_post_r) { 
    
    ar <- formalArgs(rdif)[c(-7, -9)]
    
    args <- lapply(m, function(x) unclass(x[ar]))
    
    argsT <- setNames(lapply(names(args[[1]]), function(i) lapply(args, `[[`, i)), names(args[[1]]))
    
    f <- do.call(Map, c(f = rdif, argsT))
    
    f <- lapply(f, zoo::na.locf0)
    
    m <- Map(function(x, y) transform(x, r = zoo::na.locf0(y, fromLast = TRUE)), m, f) 
  }     
  
  
  ar <- formalArgs(d.prepos)[-c(21, 22)]
  
  dot.names <- names(data)[!names(data) %in% ar]
  
  args <- lapply(m, function(x) unclass(x[c(ar, dot.names)]))
  
  argsT <- setNames(lapply(names(args[[1]]), function(i) lapply(args, `[[`, i)), names(args[[1]]))
  
  L <- do.call(Map, c(f = d.prepos, argsT))
  
  
  G <- function(m)
  {
    
    cdel3 <- reget(m, control & post == 4 & outcome == 1)
    cdel1 <- reget(m, control & post == 2 & outcome == 1)
    cdel2 <- reget(m, control & post == 3 & outcome == 1)
    cs <- reget(m, control & post == 1 & outcome == 1)
    
    tdel3 <- reget(m, !control & post == 4 & outcome == 1)
    tdel1 <- reget(m, !control & post == 2 & outcome == 1)
    tdel2 <- reget(m, !control & post == 3 & outcome == 1)
    ts <- reget(m, !control & post == 1 & outcome == 1) 
    
    if(all(sapply(list(cdel1, cdel2, cdel3, tdel1, tdel2, tdel3, ts, cs), is.null))) stop("Either 'control' or 'post' incorrectly coded.", call. = FALSE)
    
    short <- all(sapply(list(cs, ts), function(x) !is.null(x)))
    
    del1 <- all(sapply(list(cdel1, tdel1), function(x) !is.null(x)))
    
    del2 <- all(sapply(list(cdel2, tdel2), function(x) !is.null(x)))
    
    del3 <- all(sapply(list(cdel3, tdel3), function(x) !is.null(x)))
    
    cdel3..2 <- reget(m, control & post == 4 & outcome == 2)
    cdel1..2 <- reget(m, control & post == 2 & outcome == 2)
    cdel2..2 <- reget(m, control & post == 3 & outcome == 2)
    cs..2 <- reget(m, control & post == 1 & outcome == 2)
    
    tdel3..2 <- reget(m, !control & post == 4 & outcome == 2)
    tdel1..2 <- reget(m, !control & post == 2 & outcome == 2)
    tdel2..2 <- reget(m, !control & post == 3 & outcome == 2)
    ts..2 <- reget(m, !control & post == 1 & outcome == 2)
    
    
    short..2 <- all(sapply(list(cs..2, ts..2), function(x) !is.null(x)))
    
    del1..2 <- all(sapply(list(cdel1..2, tdel1..2), function(x) !is.null(x)))
    
    del2..2 <- all(sapply(list(cdel2..2, tdel2..2), function(x) !is.null(x)))
    
    del3..2 <- all(sapply(list(cdel3..2, tdel3..2), function(x) !is.null(x)))
    
    
    cdel3..3 <- reget(m, control & post == 4 & outcome == 3)
    cdel1..3 <- reget(m, control & post == 2 & outcome == 3)
    cdel2..3 <- reget(m, control & post == 3 & outcome == 3)
    cs..3 <- reget(m, control & post == 1 & outcome == 3)
    
    tdel3..3 <- reget(m, !control & post == 4 & outcome == 3)
    tdel1..3 <- reget(m, !control & post == 2 & outcome == 3)
    tdel2..3 <- reget(m, !control & post == 3 & outcome == 3)
    ts..3 <- reget(m, !control & post == 1 & outcome == 3)
    
    short..3 <- all(sapply(list(cs..3, ts..3), function(x) !is.null(x)))
    
    del1..3 <- all(sapply(list(cdel1..3, tdel1..3), function(x) !is.null(x)))
    
    del2..3 <- all(sapply(list(cdel2..3, tdel2..3), function(x) !is.null(x))) 
    
    del3..3 <- all(sapply(list(cdel3..3, tdel3..3), function(x) !is.null(x)))
    
    
    cdel3..4 <- reget(m, control & post == 4 & outcome == 4)
    cdel1..4 <- reget(m, control & post == 2 & outcome == 4)
    cdel2..4 <- reget(m, control & post == 3 & outcome == 4)
    cs..4 <- reget(m, control & post == 1 & outcome == 4)
    
    tdel3..4 <- reget(m, !control & post == 4 & outcome == 4)
    tdel1..4 <- reget(m, !control & post == 2 & outcome == 4)
    tdel2..4 <- reget(m, !control & post == 3 & outcome == 4)
    ts..4 <- reget(m, !control & post == 1 & outcome == 4)
    
    short..4 <- all(sapply(list(cs..4, ts..4), function(x) !is.null(x)))
    
    del1..4 <- all(sapply(list(cdel1..4, tdel1..4), function(x) !is.null(x)))
    
    del2..4 <- all(sapply(list(cdel2..4, tdel2..4), function(x) !is.null(x)))         
    
    del3..4 <- all(sapply(list(cdel3..4, tdel3..4), function(x) !is.null(x)))
    
    
   
    if(short){
      
    SHORT <- dinter(m, cs, ts, pos = 1, outcom = 1) 
    
    }
    
    
    if(short..2){
      
      SHORT..2 <- dinter(m, cs..2, ts..2, pos = 1, outcom = 2) 
    }
    
    
    if(short..3){
      
      SHORT..3 <- dinter(m, cs..3, ts..3, pos = 1, outcom = 3)
    }
    
    
    if(short..4){
      
      SHORT..4 <- dinter(m, cs..4, ts..4, pos = 1, outcom = 4)
    }
    
    
    if(del1){
      
      DEL1 <- dinter(m, cdel1, tdel1, pos = 2, outcom = 1)
    }
    
    
    if(del1..2){

      DEL1..2 <- dinter(m, cdel1..2, tdel1..2, pos = 2, outcom = 2)
    }
    
    
    if(del1..3){
      
      DEL1..3 <- dinter(m, cdel1..3, tdel1..3, pos = 2, outcom = 3)
    }
    
    
    if(del1..4){
      
      DEL1..4 <- dinter(m, cdel1..4, tdel1..4, pos = 2, outcom = 4)
    }
  
    
    if(del2){
      
      DEL2 <- dinter(m, cdel2, tdel2, pos = 3, outcom = 1)
    }
    
    
    if(del2..2){
      
      DEL2..2 <- dinter(m, cdel2..2, tdel2..2, pos = 3, outcom = 2)
    }
    
    
    if(del2..3){

      DEL2..3 <- dinter(m, cdel2..3, tdel2..3, pos = 3, outcom = 3)
    }
    
    
    if(del2..4){
      
      DEL2..4 <- dinter(m, cdel2..3, tdel2..3, pos = 3, outcom = 4)
    }
    
    
    if(del3){
      
      DEL3 <- dinter(m, cdel3, tdel3, pos = 4, outcom = 1)
    }
    
    
    if(del3..2){

      DEL3..2 <- dinter(m, cdel3..2, tdel3..2, pos = 4, outcom = 2)
    }
    
    
    if(del3..3){

      DEL3..3 <- dinter(m, cdel3..3, tdel3..3, pos = 4, outcom = 3)
    }
    
    
    if(del3..4){
      
      DEL3..4 <- dinter(m, cdel3..4, tdel3..4, pos = 4, outcom = 4)
    }
      
    
    list(SHORT = if(short) SHORT else NULL, SHORT..2 = if(short..2) SHORT..2 else NULL, SHORT..3 = if(short..3) SHORT..3 else NULL, SHORT..4 = if(short..4) SHORT..4 else NULL,
         DEL1 = if(del1) DEL1 else NULL, DEL1..2 = if(del1..2) DEL1..2 else NULL, DEL1..3 = if(del1..3) DEL1..3 else NULL, DEL1..4 = if(del1..4) DEL1..4 else NULL, 
         DEL2 = if(del2) DEL2 else NULL, DEL2..2 = if(del2..2) DEL2..2 else NULL, DEL2..3 = if(del2..3) DEL2..3 else NULL, DEL2..4 = if(del2..4) DEL2..4 else NULL,
         DEL3 = if(del3) DEL3 else NULL, DEL3..2 = if(del3..2) DEL3..2 else NULL, DEL3..3 = if(del3..3) DEL3..3 else NULL, DEL3..4 = if(del3..4) DEL3..4 else NULL) 
  }
  
  setNames(lapply(L, G), names(L))
}


#================================================================================================================================================================



meta.set <- function(data, file.name = NULL, na = "", impute_pre_post_r = FALSE){
  
  data <- rm.allrowNA(trim(data))
  
  L <- dint(data, impute_pre_post_r = impute_pre_post_r)
  
  d <- do.call(rbind, 
               Map(cbind, pp <- Filter(Negate(is.null), lapply(L, function(x) 
                 do.call(rbind, x))), 
                 studyID = seq_along(pp)))
  
  h <- cbind(study.name = sub("(.*)\\.(SHORT|DEL(1|2|3))(\\.+\\d.*)?", "\\1", rownames(d)), 
             studyID = d$studyID, dintID = seq_len(nrow(d)), d[-3])
  
  rownames(h) <- NULL
  
  n <- split(h, h$study.name)  
  hh <- setNames(lapply(n, NROW), names(n))
  
  D <- data
  
  m <- split(D, D$study.name)
  
  ff <- setNames(lapply(seq_along(m), function(i) NROW(subset(m[[i]], !control))), names(m))
  
  eq <- identical(ff, hh)
  
  if(eq){
    
    message(paste0("'dints' successfully computed and reshaped into 'long form'."))
    
    ar <- formalArgs(d.prepos)[-(20:22)]
    
    mod.names <- names(D)[!names(D) %in% ar]
    
    mods <- subset(D[order(D$study.name), ], !control, select = mod.names)
    
    H <- roundi(cbind(h, mods))
    
    if(!is.null(file.name)) { 
      file.name <- trimws(file.name)
      nm <- paste0(file.name, ".csv")
      ur <- try(write.csv(H, nm, row.names = FALSE, na = na), silent = TRUE)
      if(inherits(ur, "try-error")) stop(paste0("\nClose the Excel file '", nm, "' and try again OR pick another file name."), call. = FALSE)
      message(paste0("\nNote: Check folder '", basename(getwd()),"' for the Excel file '", nm, "'.\n"))
    }
    
    return(H)
    
  } else {
    
    message(paste0("\nProblem in coding sheet detected. See error analysis below:\n"))
    ur <- try(test.sheet(D, metaset = TRUE), silent = TRUE)
    
    if(inherits(ur, "try-error")) { stop("\nIncomplete data: check descrptive columns ('n', 'mpre' etc.)", call. = FALSE) 
      
    } else {
      
      message("\nSo, pre-post effect sizes can be calculated but not 'dints'... why?")    
      
      ur <- try(dint(D, impute_pre_post_r = impute_pre_post_r), silent = TRUE)
      
      if(inherits(ur, "try-error")) { stop("\nDue to incorrect/incomplete data: Check columns 'post', 'control', & 'outcome'.", call. = FALSE) 
        
      } else if(is.null(unlist(ur))) {
        
        message("Due probably to column 'post' having inconsistent values.")
        
      } else { message("\nDue to coding inconsistencies/errors (ex. part of a study coded or entered into this program but not all of it).")}
    }
  }
}     


#================================================================================================================================================================


need54_ <- c("zoo") #, "bayesmeta", "distr", "robumeta", "ellipse", "lavaan", "semPlot", "tidyverse", "weightr", "nlme", "lme4","effects"
not.have54_ <- need54_[!(need54_ %in% installed.packages()[,"Package"])]
if(length(not.have54_)) install.packages(not.have54_)


suppressWarnings(                                         
  suppressMessages({ 
    
    for(i in need54_){
      library(i, character.only = TRUE)
    }
  }))      


# Example of use:
# w7 <- read.csv('https://raw.githubusercontent.com/rnorouzian/m/master/w7.csv')

# meta.set(w7, impute_pre_post_r = FALSE)

 
