
####################################################################################
####			               interpeval                                       			####
####################################################################################

#' Function to evaluate the accuracy of the interpolation.
#' 
#' 
#' @param fit output of saemix.
#' @param rgrid individual grids and predictions.
#' @param simulations a list of simulations.
#' @param output a list with fields
#' \itemize{
#'   \item \code{name} path of the output
#'   \item \code{scaleFactor} ratio between the unit of data and the unit of model output (default=1)
#'   }
#' @param id a vector of selected id's (default: all the id's are used)
#' @param index a vector of selected indexes (default: all the indexes are used).
#' @param extreme a boolean, compare the predictions obtained with extreme parameter values (default = TRUE)
#' @return A list.
#' @author Marc Lavielle.
#' 
#' @examples
#' 
#' r <- interpeval(fit.saemix, rgrid, simulations)
#' 
#' plot(r$plot[[1]])  
#' plot(r$plot[[2]])
#' plot(r$plot[[3]])
#' head(r$pred.extreme)
#' 
#' r <- interpeval(fit.saemix, rgrid, simulations, index=c(3,6,9))
#' r <- interpeval(fit.saemix, rgrid, simulations, id=c("id-03","id-06","id-09"))
#'
#' @export interpeval

interpeval <- function(fit, rgrid, simulations, output=NULL, index=NULL, id=NULL, extreme=T) {
  
  model <- fit@model@modelPKSim
  data <- fit@data@data %>% rename(id=fit@data@name.group, time=fit@data@name.X, y=fit@data@name.response) %>%
    mutate(f.approx=fit@results@predictions$ipred)
  
  list.index <- seq_len(fit@data@N)
  list.id <- levels(data$id)
  if (!is.null(id)) {
    list.index <- sort(list.index[match(id, list.id)])
  } else if (!is.null(index)) {
    list.index <- index
  }
  list.id <- list.id[list.index]
  simulations <- simulations[list.id]
  rgrid[1:3] <- lapply(rgrid[1:3], function(x) x[list.index])
  data <- data %>% dplyr::filter(index %in% list.index) %>% droplevels()
  
  psif <- psi(fit)
  d <- ncol(psif)
  psi.ebe <- as.matrix(psif[list.index,], ncol=d)
  if (d==1) {
    psi.ebe <- as.matrix(psi.ebe, ncol=1)
    colnames(psi.ebe) <- colnames(psif) 
  }
  N <- length(simulations)
  pred <- NULL
  cat(" // Predictions using EBEs //\n")
  pb <- progress::progress_bar$new(format = "   [:bar] :percent // id :what  ", total = N, width = 60)
  pb$tick(0)
  Sys.sleep(0.5)
  for (i in seq_len(N)) {
    pb$tick(tokens = list(what = i))
    pred <- c(pred, model(psi.ebe[i,], simulations[[i]], output))
  }
  
  Qebe <- data %>% mutate(f.exact=pred) %>% select(id, time, y, f.exact, f.approx)
  Ef <- Qebe %>% reframe(rmse=sqrt(mean(((f.approx-f.exact)/f.exact)^2)), maxe=max(abs((f.approx-f.exact)/f.exact)), .by=id)
  
  txt <- NULL
  for (j in seq_len(ncol(psi.ebe)))
    txt <- paste0(txt, paste0(colnames(psi.ebe)[j], " = ", signif(psi.ebe[,j],3),"\n"))
  Ttxt <- data.frame(id=factor(levels(data$id), levels=levels(data$id)), text=txt)
  pl1 <- ggplot(Qebe)  +
    geom_line(aes(time, f.approx, color="f.approx"), linewidth=0.75)  + geom_line(aes(time, f.exact, color="f.exact")) + 
    geom_point(aes(time, f.approx, color="f.approx")) + geom_point(aes(time, f.exact, color="f.exact"), size=1) +  
    facet_wrap(~id, scales="free") + ylab("Predictions using EBEs") + scale_color_manual(values = c("f.approx"="red", "f.exact"="blue")) + theme(legend.title = element_blank())
  pl1 <- pl1 + geom_text(data=Ttxt,  mapping=aes(x=Inf, y=Inf, label=text), hjust = 1.1, vjust = 1.1) 
  
  if (!extreme) {
    pl <- pl1
    Q <- wm <- NULL
  } else {
    W <- rgrid[[1]]
    L <- rgrid[[2]]
    f <- rgrid[[3]]
    j1 <- which(unlist(lapply(W[[1]], function(x) {length(x)>1})))
    d <- length(j1)
    mi <- 2*d
    wpm <- list()
    for (i in 1:N) {
      wpm[[i]] <- matrix(unlist(lapply(W[[i]], function(x) {mean(x)})), nrow=1)
      W[[i]] <- W[[i]][j1]
      L[[i]] <- array(drop((L[[i]])), dim=dim(L[[i]])[j1])
    }
    r.ext <- extreme.psi(L)
    pname <- rgrid$parameter
    Q <- wm <- NULL
    cat(" // Predictions using extreme parameter values //\n")
    pb <- progress::progress_bar$new(format = "   [:bar] :percent // id :what  ", total = N, width = 60)
    pb$tick(0)
    Sys.sleep(0.5)
    for (i in seq_len(N)) {
      pb$tick(tokens = list(what = i))
      idi <- levels(data$id)[[i]]
      ti <- (data %>% dplyr::filter(id==idi))[["time"]]
      nti <- length(ti)
      
      ri=r.ext[[i]]
      Wi <- W[[i]]
      Li <- L[[i]]
      fi <- f[[i]]
      predi0 <- predim <- wim <- NULL
      wmi <- wpm[[i]]
      for (j in seq_len(mi)) {
        if (d==1)
          xij <- matrix(ri[((j-1)*2^d+1):(j*2^d),], ncol=ncol(ri))
        else
          xij <- ri[((j-1)*2^d+1):(j*2^d),]
        wijm <- matrix(1, ncol=d)
        for (k in seq_len(ncol(xij)))
          wijm[k] <- mean(Wi[[k]][xij[,k]])
        uij <- NULL
        for (k in seq_len(2^d))
          uij <- c(uij, do.call(`[`, c(list(Li), xij[k,])))
        wmi[j1] <- wijm
        predi0 <- c(predi0, model(wmi, simulations[[i]], output))
        predim <- c(predim, rowMeans(matrix(fi[, uij], nrow=nrow(fi))))
        wim <- rbind(wim, data.frame(wmi))
      }
      names(wim) <- pname
      wm <- rbind(wm, wim %>% mutate(id=idi, q=1:mi))
      Qi <- data.frame(id=idi, q=rep(1:mi, each=nti), time=rep(ti,mi), f.exact=predi0, f.approx=predim)
      Q <- rbind(Q,  Qi)
    }
    Q <- Q %>% mutate(id=factor(id, levels=levels(data$id)), q=factor(q)) # %>% relocate(id, .before=1) %>% relocate(q, .after=1)
    wm <- wm %>% mutate(id=factor(id, levels=levels(data$id)), q=factor(q)) %>% relocate(id, .before=1) %>% relocate(q, .after=1)
    pl2 <- Q %>% select(time, f.exact, f.approx, id, q) %>% pivot_longer(-c(id, q, time)) %>% 
      mutate(name=factor(name, levels=c("f.exact", "f.approx"))) %>%
      ggplot() + geom_line(aes(time,value,color=q, linetype=name), linewidth=0.75) + 
      facet_wrap(~id, scales="free") + ylab("Predictions using exteme parameter values") + labs(linetype = NULL) 
    
    psi.ebe <- as.data.frame(psi.ebe) %>% select(pname[j1])
    plg <- gridplot(psi.ebe, list(W,L,parameter=pname[j1]), list.id)
    pl <- append(list(pl1, pl2), plg)
  }
  return(list(plot=pl, pred.EBE=Qebe, errEBE=Ef, pred.extreme=Q, psi.extreme=wm))
}

# -----------------------------------------------------------------
# -----------------------------------------------------------------
# -----------------------------------------------------------------

gridplot <- function(psi.ebe, rgrid, ids) {
  W <- rgrid[[1]]
  L <- rgrid[[2]]
  p <- rgrid$parameter
  N <- length(L)
  d <- length(W[[1]])
  
  D <- M <- NULL
  for (i in seq_len(N)) {
    Wi <- W[[i]]
    Li <- array(drop((L[[i]])), dim=dim(L[[i]]))
    u <- extreme.psi(Li)[[1]]
    n <- max(Li)
    if (d==1)
      s <- matrix(sapply(1:n, function(x) which(c(Li)==x, arr.ind = TRUE)), ncol=1)
    else
      s <- t(sapply(1:n, function(x) which(Li==x, arr.ind = TRUE)))
    Di = data.frame(id=ids[i], k=1:n)
    Mi = data.frame(id=ids[i], q=factor(rep(seq_len(2*d), each=2^d)))
    for (j in seq_len(length(Wi))) {
      Di[p[j]] <- Wi[[j]][s[,j]]
      Mi[p[j]] <- Wi[[j]][u[,j]]
    }
    D <- rbind(D, Di)
    M <- rbind(M, Mi)
  }
  D <- D %>% mutate(id=factor(id, levels=ids))
  M <- M %>% mutate(id=factor(id, levels=ids))
  psi.ebe <- psi.ebe %>% mutate(id=factor(ids, levels=ids))
  if (d == 1) {
    pl <- list(ggplot(D) + geom_point(aes(.data[[p]], 0), color="grey70") + 
                 geom_point(data=M, aes(.data[[p]], 0, color=q), size=2.5) + geom_point(data=psi.ebe, aes(.data[[p]], 0), shape=4, size=1.5, stroke=2) +
                 facet_wrap(~id, scales="free") + theme(axis.ticks.y=element_blank(),axis.text.y=element_blank(), legend.position="none") + ylab(NULL)) 
  } else {
    pl <- list()
    k <- 0
    for (i1 in 1:(d-1))
      for (i2 in (i1+1):d) {
        k <- k+1
        pl[[k]] <- ggplot(D) + geom_point(aes(.data[[p[i1]]], .data[[p[i2]]]), color="grey70") + 
          geom_point(data=M, aes(.data[[p[i1]]], .data[[p[i2]]], color=q), size=2.5) + 
          geom_point(data=psi.ebe, aes(.data[[p[i1]]], .data[[p[i2]]]), shape=4, size=2, stroke=2) + facet_wrap(~id, scales="free")
      }
  }
  return(pl)
}

# -----------------------------------------------------------------
# -----------------------------------------------------------------
# -----------------------------------------------------------------

extreme.psi <- function(L) {
  if (!is.list(L))  L <- list(L)
  N <- length(L)
  d <- sum(dim(L[[1]])>1)
  
  Rf <- list()
  for (i in seq_len(N)) {
    Li <- L[[i]]
    X <- which(Li>0, arr.ind = TRUE)
    Res <- NULL
    if (d==1) {
      Res <- matrix(c(min(X), min(X)+1, max(X)-1, max(X)), ncol=1)
    } else {
      for (j in seq_len(d)) {
        imin <- which(X[,j]==min(X[,j]))
        v <- X[imin[which.min(rowMeans(X[imin,]))],]
        names(v) <- NULL
        w <- list()
        for (m in seq_len(d))
          w <- append(w, list(c(v[m], v[m]+1)))
        Res <- rbind(Res, expand.grid(w))
        imax <- which(X[,j]==max(X[,j]))
        v <- X[imax[which.max(rowMeans(X[imax,]))],]
        names(v) <- NULL
        w <- list()
        for (m in seq_len(d))
          w <- append(w, list(c(v[m], v[m]-1)))
        Res <- rbind(Res, expand.grid(w))
      }
    }
    Rf[[i]] <- Res
  }
  return(Rf)
}

# -----------------------------------------------------------------
# -----------------------------------------------------------------
# -----------------------------------------------------------------

extreme.bak <- function(L) {
  if (!is.list(L))  L <- list(L)
  N <- length(L)
  L1 <- L[[1]]
  dL <- dim(L1)
  d <- length(dL)
  perm <- permutations(d)
  
  M <- matrix(1:d, nrow=1)
  for (k in d:1) {
    A <- M
    A[,k:d] <- A[,k:d] + max(M) + 1 - A[1,k]
    M <- rbind(M, A)
  }
  Dc <- dim(L1)/2
  
  Rf <- list()
  for (i in seq_len(N)) {
    Li <- L[[i]]
    R.IM <- array(dim=c(2^d, 2^d, factorial(d)))
    R.DM <- array(dim=c(2^d, factorial(d)))
    for (jp in seq_len(nrow(perm))) {
      Lp <- aperm(Li, perm[jp,])
      X <- which(Lp>0, arr.ind=TRUE)
      operm <- order(perm[jp,])
      
      jj <<- NULL
      mm <<- NULL
      r <- test(X, d)
      S <- sign(matrix(r$jj[M], ncol=d))
      So <- S[, operm]
      IM <- matrix(nrow=2^d, ncol=2^d)
      Djp <- matrix(r$mm[M], ncol=d)
      DM <- vector(length=2^d)
      for (j in 1 : nrow(M)) {
        Q <- matrix(r$mm[M[j,]], nrow=1)
        for (k in d:1) {
          A <- Q
          A[,k] <- A[,k]  - S[j,k] 
          Q <- rbind(Q, A)
        }
        jm <- NULL
        for (k in 1:nrow(Q))
          jm <- c(jm, do.call(`[`, c(list(Lp), Q[k,])))
        cS <- which(apply(S, 1, identical, So[j,]))
        IM[, cS] <- jm
        DM[cS] <- sum((Djp[j,] - Dc[perm[jp,]])^2)
      }
      colnames(IM) <- NULL
      R.IM[,,jp] <- IM
      R.DM[,jp] <- DM
    }
    rf <- NULL
    for (j in seq_len(nrow(R.DM))) {
      jM <- which.max(R.DM[j,])
      rf <- cbind(rf , R.IM[,j,jM])
    }
    Rf[[i]] <- t(unique(t(apply(rf, 2, sort))))
  }
  return(Rf)
}

permutations <- function(n){
  if(n==1){
    return(matrix(1))
  } else {
    sp <- permutations(n-1)
    p <- nrow(sp)
    A <- matrix(nrow=n*p,ncol=n)
    for(i in 1:n)
      A[(i-1)*p+1:p,] <- cbind(i,sp+(sp>=i))
    return(A)
  }
}

test <- function(X, d, op=NULL) {
  if (is.null(op)) {
    jj <<- c(jj, -(d+1 - ncol(X)))
    test(X, d, 1)
    jj <<- c(jj, (d+1 - ncol(X)))
    test(X, d, 2)
  } else {
    if (op==1) {
      i1 <- which(X[,1]==min(X[,1]))
      mm <<- c(mm, min(X[,1]))
    } else {
      i1 <- which(X[,1]==max(X[,1]))
      mm <<- c(mm, max(X[,1]))
    }
    if (ncol(X) > 1) {
      jj <<- c(jj, -(d+2 - ncol(X)))
      test(matrix(X[i1, -1], nrow=length(i1)), d, 1)
      jj <<- c(jj, (d+2 - ncol(X)))
      test(matrix(X[i1, -1], nrow=length(i1)), d, 2)
    }
    else
      return(list(mm=mm, jj=jj))
  }
}

