
####################################################################################
####			               gridfuncd                                        			####
####################################################################################

#' Function to create the grid used by saemixPBPK.
#' 
#' 
#' @param simulations a list of simulations.
#' @param data a data frame with the data.
#' @param output a list with fields
#' \itemize{
#'   \item \code{name} path of the output
#'   \item \code{scaleFactor} ratio between the unit of data and the unit of model output (default=1)
#'   }
#' @param psi0 a vector of initial values.
#' @param psi.min a vector of minimum parameter values (default = 0.05*psi0).
#' @param psi.max a vector of maximum parameter values (default = 10*psi0).
#' @param e.obj the maximum relative error between interpolations and model predictions (default=0.05).
#' @param h.obj the minimum RMSE between model predictions and observations (default=0.25). Used to eliminate unrealistic values in the grid.
#' @param rh a coefficient between 0 and 1 (default = 0.9).
#' @param nw0 number of points in the initial grid (default = 11).
#' @param alpha power used to define the initial grid (default = 0.5).
#' @param N0 number of id's used to define the initial grid (default = 6).
#' @param id a vector of selected id's (default = NULL)
#' @param index a vector of selected indexes (default = NULL).
#' @param model model function (default = modelBatch).
#' @param seed seed used to randomly select N0 id's for the initial grid (default = 12345).
#' @return A list of d grids.
#' @author Marc Lavielle.
#' 
#' @examples
#' 
#' x.grid05 <- gridfuncd(simulations=theod1.simulations, data=theo.data, psi0=c(Cl=0.005))
#' x.grid02 <- gridfuncd(simulations=theod2.simulations, data=theo.data, psi0=c(Cl=0.005, lipo=1), nw0=7, N.select=5, e.obj=0.02)
#
#'
#' @export gridfuncd

gridfuncd <- function(simulations=NULL, data=NULL, output=NULL, psi0=NULL, psi.min=NULL, psi.max=NULL, e.obj=0.05, 
                      h.obj=0.25, rh=0.9, nw0=11, alpha=0.5, N0=6, model=modelBatch, seed=12345,
                      index=NULL, id=NULL) {
  
  model <<- model
  pred.mse <- function(y, list.sim, p, output) {
    pred <- NULL
    for (sim in list.sim) 
      pred <-  c(pred, model(p, sim, output))
    h <- 1 - mean((y-pred)^2)/mean(y^2)
    return(h)
  }
  
  # browser()
  # 
  # if (identical(model,modelBatch)) {
  #   
  # }
  
  names(data) <- tolower(names(data))
  if (!("id" %in% names(data)) )
    stop("'id' column (identifiers) is missing in the data file", call. = FALSE)
  if (!("y" %in% names(data)) )
    stop("'y' column (observations) is missing in the data file", call. = FALSE)
  if (!is.null(psi.min)) {
    if (length(psi.min) != length(psi0))
      stop("psi.min and psi0 should have the same length", call. = FALSE)
    if (min(psi0-psi.min) <= 0)
      stop("psi0 should be greater than psi.min", call. = FALSE)
  }
  if (!is.null(psi.max)) {
    if (length(psi.max) != length(psi0))
      stop("psi.max and psi0 should have the same length", call. = FALSE)
    if (min(psi.max-psi0) <= 0)
      stop("psi.max should be greater than psi0", call. = FALSE)
  }
  if (e.obj<=0 | e.obj>=1) 
    stop("e.obj should be strictly between 0 and 1", call. = FALSE)
  if (h.obj<=0 | h.obj>=1) 
    stop("h.obj should be strictly between 0 and 1", call. = FALSE)
  if (rh<=0 | rh>=1) 
    stop("rh should be strictly between 0 and 1", call. = FALSE)
  if (alpha<=0) 
    stop("alpha should be strictly positive", call. = FALSE)
  
  cat(" // Compute grid //\n")
  set.seed(seed)
  N <- length(simulations)
  if (N0>N) 
    stop(paste0("N0 should be less than or equal to ", N), call. = FALSE)
  
  data <- data %>% select(id, y) %>% mutate(id=factor(id))
  list.id <- levels(data$id)
  if (!is.null(id)) {
    index <- match(id, list.id)
  if (anyNA(index))
    stop("input list of id's does not match with the id's of the data file", call. = FALSE)
    index <- sort(index)
  }
  if (!is.null(index)) {
    if (min(index)<1 | max(index)>N)
      stop(paste0("input list of indexes should be between 0 and ", N), call. = FALSE)
    select.index <- index
  } else
    select.index <- sort(sample(N, N0))
  select.sim <- simulations[select.index]
  
  sim1 <- select.sim[[1]]$simBatch
  if (!"SimulationBatch" %in% class(sim1))
    stop("simulation must be a OSP SimulationBatch", call. = FALSE)
  # if (length(sim1$getVariableParameters())!=length(psi0))
  #   stop("length of psi0 does not match with the number of parameters required for the simulation ", call. = FALSE)
  
  select.data <- data %>% dplyr::filter(id %in% list.id[select.index]) %>% droplevels()
  y <- select.data$y
  
  d <- length(psi0)
  psi0.names <- names(psi0)
  names(psi0) <- NULL
  
  if (is.null(psi.min))
    psi.min <- 0.05*psi0 
  if (is.null(psi.max))
    psi.max <- 10*psi0 
  
  for (jp in seq_len(d)) {
    w.min <- w.max <- psi0
    
    w.min[jp] <- psi.min[jp]
    h.min <- try(pred.mse(y, select.sim, w.min, output), silent = T)
    if (inherits(h.min, "try-error"))  h.min <- 0
    k.min <- 0
    while (h.min < h.obj) {
      k.min <- k.min + 1
      if (k.min == 50)
        stop("Unable to find a reasonable fit... try another parameter value", call. = FALSE)
      w.min[jp] <- (1-rh)*psi0[jp] + rh*w.min[jp]
      h.min <- try(pred.mse(y, select.sim, w.min, output), silent = T)
      if (inherits(h.min, "try-error"))  h.min <- 0
    }
    psi.min[jp] <- w.min[jp]
    
    w.max[jp] <- psi.max[jp]
    h.max <- try(pred.mse(y, select.sim, w.max, output), silent = T)
    if (inherits(h.max, "try-error"))  h.max <- 0
    k.max <- 0
    while (h.max < h.obj) {
      k.max <- k.max + 1
      if (k.max == 50)
        stop("Unable to find a reasonable fit... try another parameter value", call. = FALSE)
      w.max[jp] <- (1-rh)*psi0[jp] + rh*w.max[jp]
      h.max <- try(pred.mse(y, select.sim, w.max, output), silent = T)
      if (inherits(h.max, "try-error"))  h.max <- 0
    }
    psi.max[jp] <- w.max[jp]
  }
  
  xg <- list()
  #  pb <- progress_bar$new(format = paste0("   [:bar] :percent // ", psi0.names[jp], " : w0 =  :what  "), total = length(lw0), width = 60)
  pb <- progress::progress_bar$new(format = paste0("   [:bar] :percent //  :param =  :what  "), total = d*nw0, width = 60)
  pb$tick(0)
  Sys.sleep(0.5)
  # pb <- progress::progress_bar$new(format = paste0("   [:bar] :percent //  :param =  :what  "), total = NA, width = 60)
  for (jp in seq_len(d)) {
    wd <- wm <- psi0
    lwa <- seq(psi.min[jp]^alpha, psi.max[jp]^alpha, length=nw0)^(1/alpha)
    qwa <- rep(TRUE, nw0)
    lw1 <- lw0 <- NULL
    test.w <- F
    # w0 <- lw0 <- lwa[1]
    jw <- 0
    # while (!test.w) {
    #  for (jw in seq_len(nw0)) {
    for (w0 in lwa) {
      pb$tick(tokens = list(param=psi0.names[jp], what = signif(w0, 4)))
      jw <- jw+1
      if (qwa[jw]) {
        wd[jp] <- w0
        wb <- wd
        pred0 <- NULL
        for (sim in select.sim) 
          pred0 <-  c(pred0, model(wd, sim, output))
        w.max <- 2*w0
        wb[jp] <- w.max
        wm[jp] <- (w0 + w.max)/2
        e <- rmse(wb, select.sim, wm, pred0, output)
        while (e<e.obj) {
          w.max <- 2*w.max
          wb[jp] <- w.max
          wm[jp] <- (w0 + w.max)/2
          e <- rmse(wb, select.sim, wm, pred0, output)
        }
        r <- my.optimize(select.sim, jp, wd, w.max, pred0, e.obj, output)
        lw1 <- c(lw1, r)
        #  browser()
        # if (w0==max(lwa))
        #   test.w <- T
        # else {
        ja <- which(lwa < r)
        qwa[1:(length(ja)-1)] <- FALSE
        # lwa <- lwa[-(1:(length(ja)-1))]
        # w0 <- lwa[1]
        lw0 <- c(lw0, w0)
        # }
      } else {
        Sys.sleep(10 / 100)
      }
    }
    
    w <- lw0[1]
    x.grid <- NULL
    while (!is.na(w)) {
      x.grid <- c(x.grid, w)
      w <- approx(lw0, lw1, w)$y
    }
    
    xg[[jp]] <- x.grid
  }
  return(xg)
}

test.grid <- function(w0, w1, select.sim, output) {
  K <- 11
  wx <- seq(w0, w1, length=11)
  pred <- NULL
  for (w in wx) {
    predw <- NULL
    for (sim in select.sim) 
      predw <- c(predw, model(w, sim, output))
    pred <- cbind(pred, predw)
  }
  app <- NULL
  for (k in 1:K) {
    app <- cbind(app, ((K-k)*pred[,1]+(k-1)*pred[,K])/(K-1) )
  }
  e <- apply((abs((pred-app)/app)), MARGIN = 2, max)
  return(data.frame(x=wx, e=e))
}


# -------------------------------------------------------------------------------

rmse <- function(w1, select.sim, wm, pred0, output) { 
  pred1 <- predm <- NULL
  for (sim in select.sim) {
    pred1 <- c(pred1, model(w1, sim, output))
    predm <- c(predm, model(wm, sim, output))
  }
  predm[predm==0] <- 0.000001
  preda <- (pred0 + pred1)/2
  e <- max(abs((predm - preda)/predm))
  return(e)
}

# -------------------------------------------------------------------------------

my.optimize <- function(select.sim, j, wd, w.max, pred0, e.obj, output) {
  e.tol <- e.obj/10
  wa <- w0 <- wd[j]
  wb <- w.max
  wdk <- wdm <- wd
  wdk[j] <- wb
  wdm[j] <- (w0 + wb)/2
  ea <- 0
  eb <- rmse(wdk, select.sim, wdm, pred0, output)
  
  de <- 0.3
  pa <- (eb - e.obj)^de/((e.obj - ea)^de + (eb - e.obj)^de)
  wk <- pa * wa + (1 - pa) * wb
  wdk[j] <- wk
  wdm[j] <- (w0 + wk)/2
  ek <- rmse(wdk, select.sim, wdm, pred0, output)
  n_iter <-1 
  while (abs(ek - e.obj) > e.tol) {
    if (ek > e.obj) {
      wb <- wk
      eb <- ek
    } else {
      wa <- wk
      ea <- ek
    }
    pa <- (eb - e.obj)^de/((e.obj - ea)^de + (eb - e.obj)^de)
    wk <- pa * wa + (1 - pa) * wb
    wdk[j] <- wk
    wdm[j] <- (w0 + wk)/2
    ek <- rmse(wdk, select.sim, wdm, pred0, output)
    n_iter <- n_iter+1
    if(n_iter==100) 
      stop("Unable to find a reasonable fit... try to decrease e.obj value", call. = FALSE)
    
  }
  return(wk)
}


