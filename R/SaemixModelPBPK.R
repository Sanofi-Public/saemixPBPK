####################################################################################
####			SaemixModelPBPK class - User-level function			                      ####
####################################################################################

#' Function to create a SaemixModelPBPK object
#' 
#' 
#' @param model name of the function used to compute the structural model (default=modelBatch).
#' @param psi0 a matrix with a number of columns equal to the number d of
#' parameters in the model, and one (when no covariates are available) or two
#' (when covariates enter the model) giving the initial estimates for the fixed
#' effects. The column names of the matrix should be the names of the
#' parameters in the model, and will be used in the plots and the summaries.
#' When only the estimates of the mean parameters are given, psi0 may be a
#' named vector.
#' @param x.grid a list of d grids.
#' @param output a list with fields
#' \itemize{
#'   \item \code{name} path of the output
#'   \item \code{scaleFactor} ratio between the unit of data and the unit of model output (default=1)
#'   }
#' @param file file name with a list rgrid  that contains individual predictions previously computed).
#' @param simulations a list of simulations.
#' @param model.interpolate a boolean, whether interpolation should be used or not (default=T).
#' @param error.model type of residual error model. Valid types are constant,
#' proportional, combined and exponential (default=combined).
#' @param transform.par the distribution for each parameter (0=normal,
#' 1=log-normal, 2=probit, 3=logit). Defaults to a vector of 1s (all parameters
#' have a log-normal distribution)
#' @param fixed.estim whether parameters should be estimated (1) or fixed to
#' their initial estimate (0). Defaults to a vector of 1s
#' @param covariate.model a matrix giving the covariate model. Defaults to no
#' covariate in the model
#' @param covariance.model a square matrix of size equal to the number of
#' parameters in the model, giving the variance-covariance matrix of the model:
#' 1s correspond to estimated variances (in the diagonal) or covariances
#' (off-diagonal elements). Defaults to the identity matrix
#' @param omega.init a square matrix of size equal to the number of parameters
#' in the model, giving the initial estimate for the variance-covariance matrix
#' of the model.
#' @param error.init a vector of size 2 giving the initial value of a and b in
#' the error model. Defaults to 1 for each estimated parameter in the error
#' model
#' @return A SaemixModelPBPK object.
#' @author Marc Lavielle' 
#' @examples
#' 
#' output <- list(name = "Organism|PeripheralVenousBlood|Theophylline|Plasma (Peripheral Venous Blood)",
#'                scaleFactor = 180.17/1000)
#'
#' 
#' saemix.modeld <- saemixModelPBPKd(psi0=p0, 
#'                                   x.grid=x.grid05, 
#'                                   simulations=theophylline.simulations, 
#'                                   output=output)

#' 
#' @export saemixModelPBPK

saemixModelPBPK <-function(psi0=NULL, transform.par=NULL, x.grid=NULL, model=modelBatch,
                           simulations=NULL, output=NULL,
                           fixed.estim=NULL, omega.init=NULL, covariance.model=NULL, covariate.model=NULL,
                           error.model="combined", error.init=NULL, file=NULL, model.interpolate=TRUE) {
  
  if (!is.null(output) && is.null(output[['scaleFactor']]))
    output$scaleFactor <- 1
  
  sim1 <- simulations[[1]]$simBatch
  if (!"SimulationBatch" %in% class(sim1))
    stop("simulations[[i]][['simBatch']] must be a OSP SimulationBatch", call. = FALSE)
  if (is.null(names(psi0)))
    stop("names of psi0 are missing", call. = FALSE)
  if (is.matrix(psi0)) {
    d <- ncol(psi0)
  } else {
    d <- length(psi0)
    psi0 <- matrix(psi0, ncol=d, dimnames=list(NULL, names(psi0)))
  }
  # if (length(sim1$getVariableParameters())!=d)
  #   stop("length of psi0 does not match with the number of parameters required for the simulation ", call. = FALSE)
  if (!is.null(x.grid) && length(x.grid)!=d)
    stop("length of psi0 does not match with the dimension of the initial grid", call. = FALSE)
  if (!is.logical(model.interpolate))
    stop("model.interpolate should be logical", call. = FALSE)
  
  if (is.null(transform.par))
    transform.par <- rep(1, d)
  trans.parj <- transform.par
  if (is.null(fixed.estim))
    fixed.estim <- rep(1, d)
  if (is.null(error.init)) {
    if (grepl("combined", error.model))
      error.init <- c(0.1, 0.2)
    else
      error.init <- 0.3
  }
  if (is.null(omega.init))
    omega.init <- diag(rep(1, d))
  if (is.null(covariance.model))
    covariance.model <- diag(rep(1, d))
  
  N <- length(simulations)
  if (is.null(file) ) {
    if (model.interpolate) {
      xdim <- unlist(lapply(x.grid, length))
      if (d==1 & !is.list(x.grid))  x.grid <- list(x.grid)
      rgrid <- list(x=replicate(N, x.grid, FALSE), b=replicate(N, array(0, xdim), FALSE), f=replicate(N, NULL, FALSE), parameter=colnames(psi0))
    } else {
      rgrid <- FALSE
    }
  } else {
    if (!file.exists(file))
      stop(paste0("file ",file," does not exists"), call. = FALSE)
    load(file)
    if (!exists("rgrid") || length(rgrid) != 4)
      stop(paste0("file ",file," should contain a valid list rgrid obtained from a previous run"), call. = FALSE)
  }
  assign("rgrid",rgrid, .GlobalEnv)
  assign("model.interpolate",model.interpolate, .GlobalEnv)
  assign("output",output, .GlobalEnv)
  
  
  model.approx <-function(psi, id, xidep) {
    
    rgrid <- get('rgrid',.GlobalEnv)
    model.interpolate <- get('model.interpolate',.GlobalEnv)
    output <- get('output',.GlobalEnv)
    
    alg <- 0.5
    if (is.list(psi) & !is.data.frame(psi))   psi <- psi[[1]]
    if (!is.matrix(psi))   psi=matrix(psi, nrow=1)
    ypred <- NULL
    N <- nrow(psi)
    d <- ncol(psi)
    st <- 0
    for (k in seq_len(N)) {
      j <- id[k]
      psii <- psi[k,]
      if (model.interpolate) {
        Wi <- rgrid$x[[j]]
        lWi <- rgrid$b[[j]]
        fi <- rgrid$f[[j]]
        ii <- pii <- vector(mode="list", length=d)
        for (l in seq_len(d)) {
          Wil <- Wi[[l]]
          if (length(Wil)>1) {
            i2l <- detect_index(Wil, function(x) x>psii[l])
            if (i2l==0) {
              ni <- length(Wil)
              m <- 0
              while (psii[l] > Wil[ni]) {
                if (trans.parj[l]==0) 
                  Wil <- c( Wil, Wil[ni]+(Wil[ni] - Wil[ni-1]))
                else
                  Wil <- c(Wil, ((2*Wil[ni]^alg) - Wil[ni-1]^alg)^(1/alg))
                ni <- length(Wil)
                m <- m+1
              }
              di <- dim(lWi)
              di[l] <- m
              lWi <- abind::abind(lWi, array(0,dim=di), along=l)
              i2l <- length(Wil)
              Wi[[l]] <- Wil
            } else if (i2l==1) {
              m <- 0
              while (psii[l] < Wil[1]) {
                if (trans.parj[l]==0) 
                  Wil <- c(Wil[1]-(Wil[2] - Wil[1]), Wil)
                else
                  Wil <- c((Wil[1]^2)/Wil[2], Wil)
                m <- m+1
              }
              if (trans.parj[l]==1 && Wil[1]<=0)  
                Wil[1] <- psii[l]/5
              di <- dim(lWi)
              di[l] <- m
              lWi <- abind::abind(array(0,dim=di), lWi, along=l)
              i2l <- 2
              Wi[[l]] <- Wil
            }
            ii[[l]] <- c(i2l-1, i2l)
            pil <- (Wil[i2l] - psii[l])/(Wil[i2l] - Wil[i2l-1])
            pii[[l]] <- c(pil, 1-pil)
          } else {
            ii[[l]] <- pii[[l]] <- 1
          }
        }
        ig <- as.matrix(expand.grid(ii))
        colnames(ig) <- NULL
        whi <- apply(as.matrix(expand.grid(pii)), 1, prod)
        iw <- NULL
        sk <- 0
        for (l in 1:nrow(ig)) {
          igl <- ig[l,]
          iwl <- do.call(`[`, c(list(lWi), igl)) # = lWi[igl[1], igl[2], ...]
          if (iwl == 0) {
            sk <- sk+1
            m.psi <- Wi[[1]][igl[1]]
            if (d >= 2)
              for (m in 2:d)
                m.psi <- c(m.psi, Wi[[m]][igl[m]])
            f.m <- model(m.psi, simulations[[j]], output)
            fi <- cbind(fi, f.m)
            iwl <- ncol(fi)
            lWi <- do.call(`[<-`, c(list(lWi), igl, iwl))  # <=>  lWi[igl[1], igl[2], ...] <- iwl
          }
          iw <- c(iw, iwl)
        }
        st <- st+sk
        if (sk>0) {
          rgrid$f[[j]] <- fi
          rgrid$b[[j]] <- lWi
          rgrid$x[[j]] <- Wi
        }
        ypredi <- fi[,iw]%*%whi
      } else {
        ypredi <- model(psii, simulations[[j]], output)
      }
      ypred <- c(ypred, ypredi)
    }
    assign("rgrid",rgrid, .GlobalEnv)
    return(ypred)
  }
  saemix.model <-saemixModel(model=model.approx,
                             modelPKSim=model,
                             psi0=psi0,
                             fixed.estim=fixed.estim,
                             transform.par=transform.par,
                             omega.init=omega.init,
                             covariance.model=covariance.model,
                             covariate.model=covariate.model,
                             error.model=error.model,  
                             error.init=error.init)
  
  return(saemix.model)
}




#' Function to run simulation batches with given parameter values
#' 
#' 
#' @param psi a vector or a matrix with a number of elements, resp. with the number of columns, equal to the number d of
#' parameters in the model
#' @param simulation a list with fields
#' \itemize{
#'   \item \code{simBatch} a SimulationBatch
#'   \item \code{indexTimeRemove} indexes of times to remove (optional)
#'   }
#' @param output a list with fields (optional)
#' \itemize{
#'   \item \code{name} path of the output
#'   \item \code{scaleFactor} ratio between the unit of data and the unit of model output (default=1)
#'   }
#' @return A vector of output values.
#' @author Marc Lavielle 
#' @export modelBatch
#' 
modelBatch <- function(psi, simulation, output=NULL) {
  simBatch <- simulation$simBatch
  indexTimeRemove <- simulation$indexTimeRemove
  if (!is.matrix(psi)) psi <- as.matrix(matrix(psi, nrow=1))
  ypred <- NULL
  for (i in seq_len(nrow(psi))) {
    ids <- simBatch$addRunValues(parameterValues = psi[i,])
    simResults <- runSimulationBatches(simulationBatches = simBatch)
    simulatedData <- getOutputValues(simulationResults = simResults[[1]][[1]])$data
    if (is.null(output)) simulatedValues <- simulatedData[,3]
    else simulatedValues <- simulatedData[[output$name]]*output$scaleFactor
    if (length(indexTimeRemove)>0) simulatedValues <- simulatedValues[-indexTimeRemove]
    ypred <- c(ypred, simulatedValues)
  }
  return(ypred)
}





