library(av)
library(tuneR)
library(seewave)
library(soundgen)
library(fda)
library(plotly)
library(tidyverse)
library(quadprog)

library(Rcpp)
library(RcppParallel)
library(RcppArmadillo)

rm(list = ls())

# sourceCpp("src/kron_helper.cpp")
sourceCpp("src/bvar_kron.cpp")

# A <- matrix(rnorm(4*6), 4, 6)
# B <- matrix(rnorm(2*3), 2, 3)
# C <- matrix(rnorm(6*2), 6, 2)
# res <- KroneckerProdByBlocks(A, B, C)
# 
# 
# library(microbenchmark)
# microbenchmark(
#   serial = KroneckerProdByBlocks(A, B, C),
#   parallel = KroneckerProdByBlocksParallel(A, B, C),
#   times = 10
# )



# >> Data Reading ---------------------------------
# WARNING: maybe change the input name
load("data/data_1.RData")

# subsample gufi due to computational issues
gufi_ids_exclude = read.table("data/gufi_ids_exclude.txt",
                              header = T)

gufi_meanspec_amps = gufi_meanspec_amps[,-gufi_ids_exclude$x]
gufi = gufi[-gufi_ids_exclude$x,]

# remove the zero frequency peak since it doesn't help 
# with the inference
# the peaks are at the first position only

to_remove_indexes = c(1,2,3)

gufi_meanspec_amps = gufi_meanspec_amps[-to_remove_indexes,]
gabbiani_meanspec_amps = gabbiani_meanspec_amps[-to_remove_indexes,]

gufi_meanspec_freqs = gufi_meanspec_freqs[-to_remove_indexes]
gabbiani_meanspec_freqs = gabbiani_meanspec_freqs[-to_remove_indexes]

load("results/prima_parte/outputs/basis_selection_work_space.RData")


# save(
#   falchi,
#   falchi_meanspec_amps,
#   falchi_meanspec_freqs,
#   falchi_fdanova_result,
#   falchi_fda_anova,
#   gufi,
#   gufi_meanspec_amps,
#   gufi_meanspec_freqs,
#   gufi_fdanova_result,
#   gufi_fda_anova,
#   gabbiani,
#   gabbiani_meanspec_amps,
#   gabbiani_meanspec_freqs,
#   gabbiani_fdanova_result,
#   gabbiani_fda_anova,
#   file = "data_1_bis.RData"
# )


# >> Constraint basis fit ------------------------------------


# fit a single observation
# return coeffiecients
# constraint_grid_npoints ideally has to be high
# to guarantee the constraint on a dense points grid
ConstraintSplinesFit = function(x_grid,
                                y_vals, # vector
                                basis, # fda object
                                box_constraints = c(0,1),
                                constraint_grid_npoints = 1000){ # to optimize eval the basis once
  
  # this for loss
  B <- eval.basis(x_grid, basis)
  
  # this for constraint
  x_range = range(x_grid)
  constranint_x_grid = seq(x_range[1], x_range[2],
                           length = constraint_grid_npoints)
  
  B_const <- eval.basis(constranint_x_grid, basis)
  
  # QP components
  Dmat <- t(B) %*% B
  dvec <- t(B) %*% y_vals
  
  # Constraints: B_check %*% coef <= 1 and >= 0
  Amat <- t(rbind(-B_const, diag(1, ncol(B_const))))
  
  bvec <- c(- rep(box_constraints[2], nrow(B_const)),
            rep(box_constraints[1], ncol(B_const)))
  
  # Solve QP
  fit <- solve.QP(Dmat, dvec, Amat, bvec, meq = 0)
  
  
  return(list("basis" = basis,
              "coef" = fit$solution))
  
}


# add linear differential operator penalty
# second derivative
# this accept just a basis expansion already evaluated on a grid

ConstraintSplinesFitGeneral = function(y_vals, # vector
                                B_observed_y, # basis evaluated on observed y values
                                B_constraint, # basis evaluated on dense grid: used in constraint
                                Dmat, # quadratic penalty component
                                Amat, # constraint matrix,
                                bvec){ # constraint vector
  
  dvec <- t(B_observed_y) %*% y_vals
  
  # Solve QP
  fit <- solve.QP(Dmat, dvec, Amat, bvec, meq = 0)
  
  # return coef
  return(fit$solution)
  
}


# Leave One Out Cross Validation for constraint Splines
# varying basis number
LOOCVConstraintSplinesInt = function(x_grid,
                                    y_matrix, # each column is an observation
                                    basis_num_seq, # sequence of integers
                                    box_constraints = c(0,1),
                                    constraint_grid_npoints = 1000){
  
  error_grid = rep(NA, length(basis_num_seq))
  
  
  for (b in 1:length(basis_num_seq)){
    
    
    cat("nbasis: ")
    cat(basis_num_seq[b])
    cat("\n")
    
    # for each observation
    # compute once quantities
    
    # this for loss
    basis <- create.bspline.basis(range(x_grid),basis_num_seq[b])
    B_temp_obs <- eval.basis(x_grid, basis)
    
    # this for constraint
    x_range = range(x_grid)
    constranint_x_grid = seq(x_range[1], x_range[2],
                             length = constraint_grid_npoints)
    
    B_const <- eval.basis(constranint_x_grid, basis)
    
    
    # Constraints: B_check %*% coef <= 1 and >= 0
    Amat <- t(rbind(-B_const, diag(1, ncol(B_const))))
    
    bvec <- c(- rep(box_constraints[2], nrow(B_const)),
              rep(box_constraints[1], ncol(B_const)))
    
    
    temp_cols_error = rep(NA, NCOL(y_matrix))
    
    # for each column
    for(colj in 1:NCOL(y_matrix)){
      
      temp_err = rep(NA, length(colj))
      
      # for each point
      for(point_index in 1:length(colj)){
        # QP components
        Dmat <- t(B_temp_obs[-point_index,]) %*% B_temp_obs[-point_index,]
        dvec <- t(B_temp_obs[-point_index,]) %*% y_matrix[-point_index,colj]
        
        temp_coef <- solve.QP(Dmat, dvec, Amat, bvec, meq = 0)$solution
        
        # prediction error
        temp_err[point_index] = (y_matrix[point_index,colj] - 
                                   B_const[point_index,] %*% temp_coef)^2
        
      }
      
      temp_cols_error[colj] = mean(temp_err)
      
    }
    
    error_grid[b] = mean(temp_cols_error)
  }
  
  return(list("basis_num_seq" = basis_num_seq,
              "loocv_err" = error_grid,
              "basis_min" = basis_num_seq[which.min(error_grid)]))
  
}


# Leave One Out Cross Validation for constraint Splines
# fixed basis number
# using second derivative linear differential operator integral
# penalty, on a lambda grid
LOOCVConstraintSplinesDiff = function(x_grid,
                                     y_matrix, # each column is an observation
                                     basis_num, # integer
                                     lambda_grid,
                                     box_constraints = c(0,1),
                                     constraint_grid_npoints = 1000){
  
  error_grid = rep(NA, length(lambda_grid))
  
  # once for all
  
  x_range = range(x_grid)
  
  basis <- create.bspline.basis(x_range, basis_num)
  B_temp_obs <- eval.basis(x_grid, basis)
  
  constranint_x_grid = seq(x_range[1], x_range[2],
                           length = constraint_grid_npoints)
  
  B_const <- eval.basis(constranint_x_grid, basis)
  
  # get Penalty matrix
  # same for all folds
  
  fdPar_obj = fdPar(
    fdobj = basis,
    Lfdobj = int2Lfd(2),
    lambda = 1 # actual penalty is considered below
  )
  
  penmat <-smooth.basis(
    x_grid,
    y_matrix,
    fdPar_obj
  )$penmat
  
  
  
  for (i in 1:length(lambda_grid)){
    
    
    cat("iter lambda: ")
    cat(i)
    cat("\n")
    
    # for each observation
    # compute once quantities
    
    # this for loss
    
    
    # this for constraint
    
    # Constraints: B_check %*% coef <= 1 and >= 0
    Amat <- t(rbind(-B_const, diag(1, ncol(B_const))))
    
    bvec <- c(- rep(box_constraints[2], nrow(B_const)),
              rep(box_constraints[1], ncol(B_const)))
    
    
    temp_cols_error = rep(NA, NCOL(y_matrix))
    
    # for each column
    for(colj in 1:NCOL(y_matrix)){
      
      temp_err = rep(NA, length(colj))
      
      # for each point
      for(point_index in 1:length(colj)){
        # QP components
        Dmat <- t(B_temp_obs[-point_index,]) %*% B_temp_obs[-point_index,] + lambda_grid[i] * penmat
        dvec <- t(B_temp_obs[-point_index,]) %*% y_matrix[-point_index,colj]
        
        temp_coef <- solve.QP(Dmat, dvec, Amat, bvec, meq = 0)$solution
        
        # prediction error
        temp_err[point_index] = (y_matrix[point_index,colj] - 
                                   B_const[point_index,] %*% temp_coef)^2
        
      }
      
      temp_cols_error[colj] = mean(temp_err)
      
    }
    
    error_grid[i] = mean(temp_cols_error)
  }
  
  return(list("lambda_grid" = lambda_grid,
              "loocv_err" = error_grid,
              "lambda_min" = lambda_grid[which.min(error_grid)],
              "basis_num" = basis_num))
  
}

# return an fd object with constraint splines
# (unpenalized)
ToFdConstraintSplinesInt = function(x_grid,
                                    y_matrix, # each column is an observation
                                    basis_num, # integer
                                    box_constraints = c(0,1),
                                    constraint_grid_npoints = 1000){
  
  
  basis <- create.bspline.basis(range(x_grid),basis_num)
  B = eval.basis(x_grid, basis)
  
  # allocate, to fill
  BasisCoefMatr = matrix(NA,
                         nrow = basis_num,
                         ncol = NCOL(y_matrix))
  
  
  for(i in 1:NCOL(y_matrix)){
    temp_coef = ConstraintSplinesFit(x_grid = x_grid,
                                     y_vals = y_matrix[,i],
                                     basis = basis,
                                     box_constraints = c(0,1),
                                     constraint_grid_npoints = constraint_grid_npoints)$coef
    BasisCoefMatr[,i] = temp_coef
    
  }
  
  
  return(fd(coef = BasisCoefMatr,
            basisobj = basis,
            fdnames = list("frequency" = x_grid,
                           "reps" = rep("y", NCOL(y_matrix)),
                           "values" = "Amplitude")))
  
  
}

# return an fd object with constraint splines
# (integral squared second derivated penalized)
ToFdConstraintSplinesDiff = function(x_grid,
                                    y_matrix, # each column is an observation
                                    basis_num, # integer
                                    my.lambda,
                                    box_constraints = c(0,1),
                                    constraint_grid_npoints = 1000){
  
  
  basis <- create.bspline.basis(range(x_grid),basis_num)
  B = eval.basis(x_grid, basis)
  
  x_range = range(x_grid)
  constranint_x_grid = seq(x_range[1], x_range[2],
                           length = constraint_grid_npoints)
  
  B_const <- eval.basis(constranint_x_grid, basis)
  
  fdPar_obj = fdPar(
    fdobj = basis,
    Lfdobj = int2Lfd(2),
    lambda = 1 # actual penalty is considered below
  )
  
  penmat <-smooth.basis(
    x_grid,
    y_matrix,
    fdPar_obj
  )$penmat
  
  
  Amat <- t(rbind(-B_const, diag(1, ncol(B_const))))
  
  bvec <- c(- rep(box_constraints[2], nrow(B_const)),
            rep(box_constraints[1], ncol(B_const)))
  
  Dmat <- t(B) %*% B + my.lambda * penmat
  
  
  # allocate, to fill
  BasisCoefMatr = matrix(NA,
                         nrow = basis_num,
                         ncol = NCOL(y_matrix))
  
  
  for(i in 1:NCOL(y_matrix)){
  
    dvec <- t(B) %*% y_matrix[,i]
    
    BasisCoefMatr[,i] <- solve.QP(Dmat, dvec, Amat, bvec, meq = 0)$solution
    
  }
  
  
  return(fd(coef = BasisCoefMatr,
            basisobj = basis,
            fdnames = list("frequency" = x_grid,
                           "reps" = rep("y", NCOL(y_matrix)),
                           "values" = "Amplitude")))
  
  
}





# >> GCV Smoothing functions --------------------------------

# find best number of non penalized basis according to gcv
GCV_IntegerSmoothBasis <- function(
  n_basis_seq = seq(20, 120), # basis number sequence
  basis_order = 4, # polynomial order
  basis_rangeval, # evaluation range
  x_grid, # domain values (vector)
  response_matrix
) {
  # each column is an observation

  basis_gcv_vals = rep(NA, length(n_basis_seq))

  for (i in 1:length(n_basis_seq)) {
    bspl_base = create.bspline.basis(
      rangeval = basis_rangeval,
      norder = basis_order,
      nbasis = n_basis_seq[i]
    )

    basis_gcv_vals[i] <- mean(smooth.basis(
      x_grid,
      response_matrix,
      bspl_base,
    )$gcv)
  }

  return(list(
    "n_basis_seq" = n_basis_seq,
    "gcv_vals" = basis_gcv_vals,
    "best_n_basis" = n_basis_seq[which.min(basis_gcv_vals)]
  ))
}


# find best lambda of non penalized basis according to gcv
# using a
GCV_DifferentialSmoothBasis <- function(
  n_basis = 120, # basis number
  lambda_grid = 10^seq(-4, 4, length = 100),
  basis_order = 4, # polynomial order
  basis_rangeval, # evaluation range
  x_grid, # domain values (vector)
  response_matrix
) {
  # each column is an observation

  basis_gcv_vals = rep(NA, length(lambda_grid))

  bspl_base = create.bspline.basis(
    rangeval = basis_rangeval,
    norder = basis_order,
    nbasis = n_basis
  )

  for (i in 1:length(lambda_grid)) {
    fdpar_obj <- fdPar(
      fdobj = bspl_base,
      Lfdobj = int2Lfd(2),
      lambda = lambda_grid[i]
    )

    basis_gcv_vals[i] <- mean(smooth.basis(
      argvals = x_grid,
      y = response_matrix,
      fdParobj = fdpar_obj
    )$gcv)
  }

  return(list(
    "lambda_grid" = lambda_grid,
    "gcv_vals" = basis_gcv_vals,
    "best_lambda" = lambda_grid[which.min(basis_gcv_vals)]
  ))
}

# >> ANOVA functions -------------------------------------


# fit functional ANOVA
# with group indexes j = 1,,,J
# with constraint on zero sum for beta_j(t) for each t
# i.e beta_1(t) + beta_2(t) + ... + beta_j(t) + .. + beta_J(t) = 0
# for each t

FitANOVAConstrast = function(factor,
                             X,
                             dom,
                             basis_y,
                             coef_y,
                             y_names,
                             basis_beta,
                             lambda){
  factor_matrix = cbind(1, model.matrix(~ -1 + factor))
  yt = cbind(X, 0)
  Xt = rbind(factor_matrix, c(0, rep(1, length(levels(factor)))))
  
  
  B = cbind(coef_y, 0) # constraint
  yfd = fd(B, basis_y, y_names)
  
  xlist = lapply(
    1:NCOL(factor_matrix),
    function(x) c(factor_matrix[, x], 1)
  )
  xlist[[1]][length(xlist[[1]])] = 0 # constraint
  
  betalist = lapply(
    1:(length(levels(factor)) + 1),
    function(x) fdPar(basis_beta, Lfdobj = int2Lfd(2), lambda = lambda)
  )
  # functional intercept
  betalist[[1]] = fdPar(basis_beta, lambda = 0)
  
  return(fRegress(yfd, xlist, betalist))
}


# helper: use more time but less resources
# avoiding saving explicitly the kronecker(B, C)
# compute A %*% (kronecker(B, C))


# Funzione per scegliere lambda ottimo per basi di beta

# copied from fda library
# this does Leave One Out CV for each observation index in CVobs
fRegress.CV2 = function(
  y,
  xfdlist,
  betalist,
  wt = NULL,
  CVobs = 1:N, # observation on which to compute the error
  returnMatrix = FALSE,
  ...
) {
  argList <- fda:::fRegressArgCheck(y, xfdlist, betalist, wt)
  yfdobj <- argList$yfd
  xfdlist <- argList$xfdlist
  betalist <- argList$betalist
  wt <- argList$wt
  p <- length(xfdlist)
  N <- dim(xfdlist[[1]]$coef)[2]
  M <- length(CVobs)
  yfd <- yfdobj
  SSE.CV <- 0
  errcoefs <- c()
  for (m in 1:length(CVobs)) {
    out_ids <- CVobs[[m]]
    wti <- wt[-out_ids]
    txfdlist <- xfdlist
    for (k in 1:p) {
      txfdlist[[k]] <- xfdlist[[k]][-out_ids]
    }
    yfdi <- yfd[-out_ids]
    tres <- fRegress(yfdi, txfdlist, betalist, wti)
    betaestlisti <- tres$betaestlist
    for (i in out_ids) {
      yhatfdi <- 0
      for (k in 1:p) {
        betafdPark = betaestlisti[[k]]
        betafdk = betafdPark$fd
        xfdk = xfdlist[[k]]
        xfdik = xfdk[i]
        tempfd = xfdik * betafdk
        yhatfdi <- yhatfdi + tempfd
      }
      errfdi <- yfd[i] - yhatfdi
      SSE.CV <- SSE.CV + inprod(errfdi, errfdi) # integration here
      errcoefs <- cbind(errcoefs, errfdi$coefs)
    }
  }
  errfd <- fd(errcoefs, errfdi$basis)
  names(errfd$fdnames)[[3]] <- "Xval Errors"
  return(list(SSE.CV = SSE.CV, errfd.cv = errfd))
}

#' @description
#' Make a vector of sets balanced for a factor variable
#' 
#' @param factor (vector of factor indicator): length equal to observation number
#' @return vector with fold indicators
MakeBalancedFolds = function(factor,
                             nfold = 5){
  folds = c()
  for (f in levels(factor)) {
    idx = cut(sample(1:length(which(factor == f))), breaks = nfold)
    levels(idx) = 1:nfold
    folds = c(folds, idx)
  }
  
  return(folds)
}

#' @description
#' Compute Cross validation integrated error
#' for functional ANOVA
#' @param factor (vector of factor indicator): length equal to observation number
#' @param X (matrix): function observation matrix, each column is an observation
#' each row correspond to a value on the domain grid
#' @param dom (vector of num): grid of domain points
#' where the function is evaluated (for example times)
CvFunctionalANOVA = function(
  factor,
  X,
  dom,
  basis_y,
  coef_y,
  y_names,
  basis_beta,
  lambda_grid,
  nfold = 5
) {
  cv_err = rep(NA, length(lambda_grid))
  
  
  folds = MakeBalancedFolds(factor = factor,
                            nfold = nfold)
  
  factor_matrix = cbind(1, model.matrix(~ -1 + factor))
  fdobj = smooth.basis(dom, X, basis_y)$fd
  
  
  B = cbind(coef_y, 0) # constraint
  yfd = fd(B, basis_y, y_names)
  
  xlist = lapply(
    1:NCOL(factor_matrix),
    function(x) c(factor_matrix[, x], 1)
  )
  xlist[[1]][length(xlist[[1]])] = 0 # constraint
  # 0 1 1 1 | 0
  
  for (i in 1:length(lambda_grid)) {
    cat(i, "")
    

    betalist = lapply(
      1:(length(levels(factor)) + 1),
      function(x)
        fdPar(basis_beta, Lfdobj = int2Lfd(2), lambda = lambda_grid[i])
    )

    # functional intercept
    # no differential penalty
    betalist[[1]] = fdPar(basis_beta, lambda = 0)

    folds_list = lapply(1:nfold, function(f) which(folds == f))
    
    
    # NOTE: since the constraint observation is NOT
    # in the folds -> it is in each fitting fold
    # because CVobs are the test observations
    # so the contraint is always respected
    m = fRegress.CV2(yfd, xlist, betalist, CVobs = folds_list)
    cv_err[i] = m$SSE.CV
  }
  
  
  return(list("lambda_grid" = lambda_grid,
              "cv_error" = cv_err,
              "lambda_min" = lambda_grid[which.min(cv_err)]))
}


# taken from fda 
My.fRegress.stderr = function (y, y2cMap, SigmaE, returnMatrix = FALSE, ...) 
{
  yfdobj <- y$yfdobj
  xfdlist <- y$xfdlist
  betalist <- y$betalist
  betaestlist <- y$betaestlist
  yhatfdobj <- y$yhatfdobj
  Cmat <- y$Cmat
  Dmat <- y$Dmat
  Cmatinv <- y$Cmatinv
  wt <- y$wt
  df <- y$df
  betastderrlist <- y$betastderrlist
  YhatStderr <- y$YhatStderr
  Bvar <- y$Bvar
  c2bMap <- y$c2bMap
  p <- length(xfdlist)
  ncoef <- 0
  for (j in 1:p) {
    betaParfdj <- betalist[[j]]
    ncoefj <- betaParfdj$fd$basis$nbasis
    ncoef <- ncoef + ncoefj
  }
  if (inherits(yfdobj, "fdPar") || inherits(yfdobj, "fd")) {
    if (inherits(yfdobj, "fdPar")) 
      yfdobj <- yfdobj$fd
    N <- dim(yfdobj$coefs)[2]
    ybasisobj <- yfdobj$basis
    rangeval <- ybasisobj$rangeval
    ynbasis <- ybasisobj$nbasis
    ninteg <- max(501, 10 * ynbasis + 1)
    tinteg <- seq(rangeval[1], rangeval[2], len = ninteg)
    deltat <- tinteg[2] - tinteg[1]
    ybasismat <- eval.basis(tinteg, ybasisobj, 0, returnMatrix)
    basisprodmat <- matrix(0, ncoef, ynbasis * N)
    mj2 <- 0
    for (j in 1:p) {
      betafdParj <- betalist[[j]]
      betabasisj <- betafdParj$fd$basis
      ncoefj <- betabasisj$nbasis
      bbasismatj <- eval.basis(tinteg, betabasisj, 0, 
                               returnMatrix)
      xfdj <- xfdlist[[j]]
      tempj <- eval.fd(tinteg, xfdj, 0, returnMatrix)
      mj1 <- mj2 + 1
      mj2 <- mj2 + ncoefj
      indexj <- mj1:mj2
      mk2 <- 0
      for (k in 1:ynbasis) {
        mk1 <- mk2 + 1
        mk2 <- mk2 + N
        indexk <- mk1:mk2
        tempk <- bbasismatj * ybasismat[, k]
        basisprodmat[indexj, indexk] <- deltat * crossprod(tempk, 
                                                           tempj)
      }
    }
    y2cdim <- dim(y2cMap)
    if (y2cdim[1] != ynbasis || y2cdim[2] != dim(SigmaE)[1]) 
      stop("Dimensions of Y2CMAP not correct.")
    Varc <- y2cMap %*% SigmaE %*% t(y2cMap)
    print("[DEBUG]: dim(Varc) = ")
    print(dim(Varc))
    
    C2BMap <- Cmatinv %*% basisprodmat
    
    print("[DEBUG]: dim(C2BMap) = ")
    print(dim(C2BMap))
    
    # this causes memory problems, it's to big
    # CVar <- kronecker(Varc, diag(rep(1, N)))
    # 

    # Bvar <- C2BMap %*% CVar %*% t(C2BMap)
    
    # Bvar <- KroneckerProdByBlocksParallel(C2BMap, Varc, diag(rep(1, N))) %*% t(C2BMap)
    
    # Bvar <- KroneckerProdByBlocks(C2BMap, Varc, diag(rep(1, N))) %*% t(C2BMap)
    
    Bvar <- ComputeBvar(C2BMap, Varc, diag(rep(1, N)))
    
    print("[DEBUG]: Bvar with Kronecker computed!")
    
    nplot <- max(51, 10 * ynbasis + 1)
    tplot <- seq(rangeval[1], rangeval[2], len = nplot)
    betastderrlist <- vector("list", p)
    PsiMatlist <- vector("list", p)
    mj2 <- 0
    for (j in 1:p) {
      betafdParj <- betalist[[j]]
      betabasisj <- betafdParj$fd$basis
      ncoefj <- betabasisj$nbasis
      mj1 <- mj2 + 1
      mj2 <- mj2 + ncoefj
      indexj <- mj1:mj2
      bbasismat <- eval.basis(tplot, betabasisj, 0, returnMatrix)
      PsiMatlist <- bbasismat
      bvarj <- Bvar[indexj, indexj]
      bstderrj <- sqrt(diag(bbasismat %*% bvarj %*% t(bbasismat)))
      bstderrfdj <- smooth.basis(tplot, bstderrj, betabasisj)$fd
      betastderrlist[[j]] <- bstderrfdj
    }
    YhatStderr <- matrix(0, nplot, N)
    B2YhatList <- vector("list", p)
    for (iplot in 1:nplot) {
      YhatVari <- matrix(0, N, N)
      tval <- tplot[iplot]
      for (j in 1:p) {
        Zmat <- eval.fd(tval, xfdlist[[j]])
        betabasisj <- betalist[[j]]$fd$basis
        PsiMatj <- eval.basis(tval, betabasisj)
        B2YhatMapij <- t(Zmat) %*% PsiMatj
        B2YhatList[[j]] <- B2YhatMapij
      }
      m2j <- 0
      for (j in 1:p) {
        m1j <- m2j + 1
        m2j <- m2j + betalist[[j]]$fd$basis$nbasis
        B2YhatMapij <- B2YhatList[[j]]
        m2k <- 0
        for (k in 1:p) {
          m1k <- m2k + 1
          m2k <- m2k + betalist[[k]]$fd$basis$nbasis
          B2YhatMapik <- B2YhatList[[k]]
          YhatVari <- YhatVari + B2YhatMapij %*% Bvar[m1j:m2j, 
                                                      m1k:m2k] %*% t(B2YhatMapik)
        }
      }
      YhatStderr[iplot, ] <- matrix(sqrt(diag(YhatVari)), 
                                    1, N)
    }
  }
  else {
    ymat <- as.matrix(yfdobj)
    N <- dim(ymat)[1]
    Zmat <- NULL
    for (j in 1:p) {
      xfdj <- xfdlist[[j]]
      if (inherits(xfdj, "fd")) {
        xcoef <- xfdj$coefs
        xbasis <- xfdj$basis
        betafdParj <- betalist[[j]]
        bbasis <- betafdParj$fd$basis
        Jpsithetaj <- inprod(xbasis, bbasis)
        Zmat <- cbind(Zmat, t(xcoef) %*% Jpsithetaj)
      }
      else if (inherits(xfdj, "numeric")) {
        Zmatj <- xfdj
        Zmat <- cbind(Zmat, Zmatj)
      }
    }
    c2bMap <- Cmatinv %*% t(Zmat)
    y2bmap <- c2bMap
    Bvar <- y2bmap %*% as.matrix(SigmaE) %*% t(y2bmap)
    betastderrlist <- vector("list", p)
    mj2 <- 0
    for (j in 1:p) {
      betafdParj <- betalist[[j]]
      betabasisj <- betafdParj$fd$basis
      ncoefj <- betabasisj$nbasis
      mj1 <- mj2 + 1
      mj2 <- mj2 + ncoefj
      indexj <- mj1:mj2
      bvarj <- Bvar[indexj, indexj]
      xfdj <- xfdlist[[j]]
      if (inherits(xfdj, "fd")) {
        betarng <- betabasisj$rangeval
        ninteg <- max(c(501, 10 * ncoefj + 1))
        tinteg <- seq(betarng[1], betarng[2], len = ninteg)
        bbasismat <- eval.basis(tinteg, betabasisj, 
                                0, returnMatrix)
        bstderrj <- sqrt(diag(bbasismat %*% bvarj %*% 
                                t(bbasismat)))
        bstderrfdj <- smooth.basis(tinteg, bstderrj, 
                                   betabasisj)$fd
      }
      else {
        bsterrj <- sqrt(diag(bvarj))
        onebasis <- create.constant.basis(betabasisj$rangeval)
        bstderrfdj <- fd(t(bstderrj), onebasis)
      }
      betastderrlist[[j]] <- bstderrfdj
    }
    B2YhatList <- vector("list", p)
    YhatVari <- matrix(0, N, N)
    for (j in 1:p) {
      betabasisj <- betalist[[j]]$fd$basis
      Xfdj <- xfdlist[[j]]
      B2YhatMapij <- inprod(Xfdj, betabasisj)
      B2YhatList[[j]] <- B2YhatMapij
    }
    m2j <- 0
    for (j in 1:p) {
      m1j <- m2j + 1
      m2j <- m2j + betalist[[j]]$fd$basis$nbasis
      B2YhatMapij <- B2YhatList[[j]]
      m2k <- 0
      for (k in 1:p) {
        m1k <- m2k + 1
        m2k <- m2k + betalist[[k]]$fd$basis$nbasis
        B2YhatMapik <- B2YhatList[[k]]
        YhatVari <- YhatVari + B2YhatMapij %*% Bvar[m1j:m2j, 
                                                    m1k:m2k] %*% t(B2YhatMapik)
      }
    }
    YhatStderr <- matrix(sqrt(diag(YhatVari)), N, 1)
  }
  fRegressList <- list(yfdobj = y$yfdobj, xfdlist = y$xfdlist, 
                       betalist = y$betalist, betaestlist = y$betaestlist, 
                       yhatfdobj = y$yhatfdobj, Cmat = y$Cmat, Dmat = y$Dmat, 
                       Cmatinv = y$Cmatinv, wt = y$wt, df = y$df, y2cMap = y2cMap, 
                       SigmaE = SigmaE, betastderrlist = betastderrlist, YhatStderr = YhatStderr, 
                       Bvar = Bvar, c2bMap = c2bMap)
  class(fRegressList) = "fRegress"
  return(fRegressList)
}

ComputeBetaSd = function(factor,
                         X,
                         dom,
                         basis_y,
                         coef_y,
                         y_names,
                         basis_beta,
                         lambda){
  
  factor_matrix = cbind(1, model.matrix(~ -1 + factor))
  yt = cbind(X, 0)
  Xt = rbind(factor_matrix, c(0, rep(1, length(levels(factor)))))
  
  
  B = cbind(coef_y, 0) # constraint
  yfd = fd(B, basis_y, y_names)
  
  xlist = lapply(
    1:NCOL(factor_matrix),
    function(x) c(factor_matrix[, x], 1)
  )
  xlist[[1]][length(xlist[[1]])] = 0 # constraint
  
  betalist = lapply(
    1:(length(levels(factor)) + 1),
    function(x) fdPar(basis_beta, Lfdobj = int2Lfd(2), lambda = lambda)
  )
  # functional intercept
  betalist[[1]] = fdPar(basis_beta, lambda = 0)
  
  m1 = fRegress(yfd, xlist, betalist)
  
  
  # compute functional standard errors
  # NOTE: change yfd with another basis system for errors 
  SigmaE = diag(apply(eval.fd(dom, yfd), 1, var))
  
  
  return(list("model" = m1,
              "beta_se" = My.fRegress.stderr(
                m1,
                smooth.basis(dom, X, basis_y)$y2cMap,
                SigmaE)
              )
    )
  
}

# non parametric bootstrap beta IC
# for each (stratified by group) bootstrap sample
# estimate beta_j(f) for each group j,
# discretize its values on a grid of {f_g} g = 1,,,,G
# when all bootstrap iterations are finished
# for each point on the grid take the wanted quantiles (ex. 0.025, 0.975)
# and smooth them to get the functional beta 
BootBetaIC = function(B, # bootstrap replicates
                      factor,
                      X,
                      dom,
                      basis_y,
                      coef_y,
                      y_names,
                      basis_beta,
                      lambda,
                      n_discrete = 500, # number of discrete points, where to evaluate each beta
                      my.quantiles = c(0.025, 0.975)){
  
  # domain points grid
  # not used
  eval_points = seq(min(dom), max(dom), length = n_discrete)
  
  # store the discretize points
  beta_boot_array = array(NA,
                          dim = c(basis_beta[["nbasis"]], # basis
                                  B,
                                  length(unique(factor)) + 1)) # intercept
  
  print(dim(beta_boot_array))
  
  for(b in 1:B){
    
    if((b %% 100) == 0){
      cat(b)
      cat("\n")
    }
    
    
    # computational inefficient, but not much problematic
    b_sample_indexes = c()
    
    # get stratified bootstrap sample
    for(fact in unique(factor)){
      
      factor_candidates_indexes = which(factor == fact)
      
      b_sample_indexes = c(b_sample_indexes,
                           sample(factor_candidates_indexes,
                                  size = length(factor_candidates_indexes),
                                  replace = T))
    }
    
    
    temp_betaestlist = FitANOVAConstrast(factor = factor[b_sample_indexes],
                      X = X[,b_sample_indexes],
                      dom = dom,
                      basis_y = basis_y,
                      coef_y = coef_y[, b_sample_indexes],
                      y_names = y_names,
                      basis_beta = basis_beta,
                      lambda = lambda)$betaestlist
    
    
    # save beta values on the grid
    # first intercept, then betas
    for (i in 1:dim(beta_boot_array)[3]){
      
      beta_boot_array[,b,i] = temp_betaestlist[[i]]$fd$coefs
      
      
    }
    
  }
  
  
  # (non-bootstrapped) fit
  original_fit_beta = FitANOVAConstrast(factor = factor,
                                        X = X,
                                        dom = dom,
                                        basis_y = basis_y,
                                        coef_y = coef_y,
                                        y_names = y_names,
                                        basis_beta = basis_beta,
                                        lambda = lambda)$betaestlist
  
  
  # 2. save quantiles as functional object to plotting
  quantile_betas = list()
  
  for(i in 1:dim(beta_boot_array)[3]){  # over intercept and factors
    quantile_betas[[i]] = list()
    for(q in my.quantiles){
      
      # splines coefs quantile for each basis function
      quant_coefs = apply(beta_boot_array[,,i], 1, quantile, probs = q)
      quant_coefs = matrix(quant_coefs, nrow = basis_beta$nbasis, ncol = 1)
      
      quantile_betas[[i]][[as.character(q)]] = fd(coef = quant_coefs,
                                                  basisobj = basis_beta)
    }
  }
  
  names(quantile_betas) = c("intercept", levels(factor))
  
  return(list(
    "original_fit_beta" = original_fit_beta,
    "quantile_betas" = quantile_betas
  ))
  
}


PermutFANOVA = function(factor,
                          X,
                          dom,
                          basis_y,
                          coef_y,
                          y_names,
                          basis_beta,
                          lambda,
                          my.quantile = 0.95,
                          n_perm = 100,
                          seed = 123){
  
  set.seed(seed)
  
  factor_matrix = cbind(1, model.matrix(~ -1 + factor))
  yt = cbind(X, 0)
  Xt = rbind(factor_matrix, c(0, rep(1, length(levels(factor)))))
  
  B = cbind(coef_y, 0) # vincolo
  yfd = fd(B, basis_y, y_names)
  
  xlist = lapply(
    1:NCOL(factor_matrix),
    function(x) c(factor_matrix[, x], 1)
  )
  xlist[[1]][length(xlist[[1]])] = 0 # vincolo
  
  betalist = lapply(
    1:(length(levels(factor)) + 1),
    function(x) fdPar(basis_beta, Lfdobj = int2Lfd(2), lambda = lambda)
  )
  # functional intercept
  betalist[[1]] = fdPar(basis_beta, lambda = 0)
  
  m1 = fRegress(yfd, xlist, betalist)
  
  
  # Permutations 
  
  # including intercept
  FtNew = function(yh, my.mean, Y){
    # tmean = t(my.mean)
    num = apply(yh, 1, function(row) (row - my.mean)^2)
    den = Y - yh
    return(apply(yh, 2, mean) / apply(den^2, 2, mean))
  }
  
  
  ytmp = yfd
  Fmatr = matrix(NA, n_perm, length(dom))
  
  
  for (b in 1:n_perm) {
    id = sample(1:ncol(X))
    ytmp$coefs = ytmp$coefs[, c(id, ncol(X) + 1)]
    mtmp = fRegress(ytmp, xlist, betalist)
    
    Y = t(eval.fd(dom, ytmp)[, 1:ncol(X)])
    intercept = t(eval.fd(dom, mtmp$betaestlist[[1]]$fd))
    Yh = t(eval.fd(dom, mtmp$yhatfdobj)[, 1:ncol(X)])
    
    # Fmatr[b, ] = Ft(Yh, Y)
    
    Fmatr[b, ] = FtNew(yh = Yh, Y = Y, my.mean = intercept)
    
    
    cat(b, "\n")
  }
  
  Y = t(eval.fd(dom, yfd)[, 1:ncol(X)])
  Yh = t(eval.fd(dom, m1$yhatfdobj)[, 1:ncol(X)])
  intercept = t(eval.fd(dom, mtmp$betaestlist[[1]]$fd))
  
  fobs = FtNew(yh = Yh, Y = Y, my.mean = intercept)
  
  return(list(
    "fobs" = fobs,
    "qF" = apply(Fmatr, 2, quantile, my.quantile),
    "qFmax" = quantile(apply(Fmatr, 1, max), ,my.quantile)
  ))
  
}


# >> Plotting Functions ---------------------------------
FunctionalMeanBandPlot = function(fd_means,
                                  fd_sds,
                                  my.main,
                                  my.ylim,
                                  my.lwd = 2,
                                  my.xlab = "Frequenza",
                                  my.ylab = "Ampiezza",
                                  my.levels.variable,
                                  save_path = NULL,
                                  my.width = 600,
                                  my.height = 350){
  
  if(!is.null(save_path)){
    png(save_path,
         width = my.width,
         height = my.height)
  }
  
  plot(
    fd_means[[1]],
    main = my.main,
    ylim = my.ylim,
    lwd = my.lwd,
    xlab = my.xlab,
    ylab = my.ylab
  )
  for (i in 2:length(fd_means))
    lines(fd_means[[i]], col = i, lwd = my.lwd)
  
  for (i in 1:length(fd_means)) {
    lines(
      fd_means[[i]] + fd_sds[[i]],
      col = i,
      lwd = my.lwd,
      lty = 2
    )
    lines(
      fd_means[[i]] - fd_sds[[i]],
      col = i,
      lwd = my.lwd,
      lty = 2
    )
  }
  
  abline(h = 0, lty = 3)
  
  
  legend("topright",
         legend = levels(my.levels.variable),
         col = 1:length(unique(my.levels.variable)),
         lwd = my.lwd,
         bty = "n")
  
  if(!is.null(save_path)){
    dev.off()
  }
}


PlotBetaWithQuantiles <- function(original_fit,
                                  quantile_betas,
                                  my.name,
                                  n_points = 500,
                                  col_main = "black",
                                  col_quant = "black",
                                  lwd_main = 2,
                                  lwd_quant = 2,
                                  lty_quant = 2,
                                  save_path = NULL,
                                  my.width = 600,
                                  my.height = 350,
                                  my.layout.matr = cbind(matrix(1, 2, 2),
                                                         matrix(2:5, 2, 2))) {
  
  
  if(!is.null(save_path)){
    png(save_path,
        width = my.width,
        height = my.height)
  }
  
  graphics::layout(mat = my.layout.matr)
  
  # common basis
  basis_beta = quantile_betas[[1]][[1]][[2]]
  
  eval_points <- seq(basis_beta$rangeval[1],
                     basis_beta$rangeval[2],
                     length.out = n_points)
  
  # Loop over all beta functions (intercept and factors)
  for (i in seq_along(original_fit)) {
    
    
    # Evaluate original beta
    main_vals <- eval.fd(eval_points, original_fit[[i]]$fd)
    
    # Evaluate quantiles
    quant_list <- quantile_betas[[i]]
    lower_vals <- eval.fd(eval_points, quant_list[[1]])  # e.g. 0.025
    upper_vals <- eval.fd(eval_points, quant_list[[2]])  # e.g. 0.975
    
    # Y-axis range
    y_min <- min(lower_vals)
    y_max <- max(upper_vals)
    
    # Set plot title
    beta_name <- names(quantile_betas)[i]
    
    # Plot
    plot(eval_points, main_vals, type = "l",
         ylim = c(y_min, y_max),
         lwd = lwd_main,
         col = col_main,
         main = paste(my.name, " : ", beta_name),
         xlab = "Frequency", ylab = "Beta(f)")
    
    # Add quantile bounds
    lines(eval_points, lower_vals, col = col_quant, lty = lty_quant,
          lwd = lwd_quant)
    lines(eval_points, upper_vals, col = col_quant, lty = lty_quant,
          lwd = lwd_quant)
    
    abline(h = 0, lty = 2)
    
    # intercept
    if(i == 1){
      legend("topright", legend = c("Mean", "Quantiles"),
             col = c(col_main, col_quant),
             lwd = c(lwd_main, 1),
             lty = c(1, lty_quant),
             bty = "n")
    }
    }
    

  
  if(!is.null(save_path)){
    dev.off()
  }
  
  par(mfrow = c(1, 1))
}



fANOVABetaSdPlot = function(my.betaestlist,
                            my.betastderrlist,
                            my.factor,
                            my.name,
                            save_path = NULL,
                            my.width = 600,
                            my.height = 350,
                            my.layout.matr = cbind(matrix(1, 2, 2),
                                                   matrix(2:5, 2, 2))){
  
  
  if(!is.null(save_path)){
    png(save_path,
         width = my.width,
         height = my.height)
  }

  graphics::layout(mat = my.layout.matr)
  
  plot(
    my.betaestlist[[1]]$fd,
    main = paste0(c(my.name, "Intercetta funzionale"), collapse = " - "),
    ylim = c(0, 1)
  ) # mu(t)
  lines(
    my.betaestlist[[1]]$fd +
      2 * my.betastderrlist[[1]],
    lty = 3
  )
  lines(
    my.betaestlist[[1]]$fd -
      2 * my.betastderrlist[[1]],
    lty = 3
  )
  sapply(
    2:(length(levels(my.factor)) + 1),
    function(x) {
      plot(
        my.betaestlist[[x]],
        main = levels(my.factor)[x - 1],
        ylim = c(-0.1, 0.1)
      )
      abline(h = 0, lty = 2)
      lines(
        my.betaestlist[[x]]$fd +
          2 * my.betastderrlist[[x]],
        lty = 3
      )
      lines(
        my.betaestlist[[x]]$fd -
          2 * my.betastderrlist[[x]],
        lty = 3
      )
    } # alpha_j(t)
  )
  
  if(!is.null(save_path)){
    dev.off()
  }
  
  par(mfrow = c(1, 1))
}

# >> Constants ------------------------------------

# if TRUE save RData, at the risk of override existing ones 
DO_SAVE_RDATA = TRUE



# save Plotting
MY.WIDTH = 2000
MY.HEIGHT = 2000




# >>Regularize fit via GCV -------------------


# ╭──────╮
# │Falchi│------------------------------------
# ╰──────╯


# Changing basis number
n_basis_seq = seq(20, nrow(falchi_meanspec_amps) - 2)
NORDER = 4
rangeval_falchi <- range(falchi_meanspec_freqs)


falchi_nbasis_gcv <- GCV_IntegerSmoothBasis(
  n_basis_seq = 20:(nrow(falchi_meanspec_amps) - 2),
  basis_order = NORDER,
  basis_rangeval = range(falchi_meanspec_freqs),
  x_grid = falchi_meanspec_freqs,
  response_matrix = falchi_meanspec_amps
)

plot(falchi_nbasis_gcv$n_basis_seq, falchi_nbasis_gcv$gcv_vals, type = "l")


# by eye 64 basis seem enough
falchi_basis_int = create.bspline.basis(
  rangeval = range(falchi_meanspec_freqs),
  norder = NORDER,
  nbasis = falchi_nbasis_gcv$best_n_basis
)

falchi_meanspec_fd_int = smooth.basis(
  falchi_meanspec_freqs,
  falchi_meanspec_amps,
  falchi_basis_int
)$fd

falchi_meanspec_fd_int$fdnames = list("frequency" = falchi_meanspec_freqs,
                                      "reps" = rep("y", length(falchi_meanspec_fd_int$fdnames$reps)),
                                      "values" = "Amplitude")

# Changing lambda penalty

falchi_LD_gcv <- GCV_DifferentialSmoothBasis(
  n_basis = nrow(falchi_meanspec_amps),
  lambda_grid = 10^seq(-10, -5, length = 100),
  basis_order = NORDER,
  basis_rangeval = range(falchi_meanspec_freqs),
  x_grid = falchi_meanspec_freqs,
  response_matrix = falchi_meanspec_amps
)

plot(
  log(falchi_LD_gcv$lambda_grid, base = 10),
  falchi_LD_gcv$gcv_vals,
  type = "l"
)

# compare gcv values
min(falchi_nbasis_gcv$gcv_vals)
min(falchi_LD_gcv$gcv_vals)

# compare fits

falchi_basis_max = create.bspline.basis(
  rangeval = range(falchi_meanspec_freqs),
  norder = NORDER,
  nbasis = nrow(falchi_meanspec_amps)
)


falchi_fdPar = fdPar(
  fdobj = falchi_basis_max,
  Lfdobj = int2Lfd(2),
  lambda = falchi_LD_gcv$best_lambda
)

falchi_basis_diff <- smooth.basis(
  falchi_meanspec_freqs,
  falchi_meanspec_amps,
  falchi_fdPar
)

falchi_basis_diff$fd$fdnames = list("frequency" = falchi_meanspec_freqs,
                                      "reps" = rep("y", length(falchi_basis_diff$fd$fdnames$reps)),
                                      "values" = "Amplitude")


falchi_meanspec_fd_diff = falchi_basis_diff$fd



# ╭────╮
# │Gufi│ ----------------------------------------------------------
# ╰────╯

dim(gufi_meanspec_amps)

# Changing basis number
n_basis_seq = seq(20, nrow(gufi_meanspec_amps) - 2)
NORDER = 4
rangeval_gufi <- range(gufi_meanspec_freqs)


gufi_nbasis_gcv <- GCV_IntegerSmoothBasis(
  n_basis_seq = 20:(nrow(gufi_meanspec_amps) - 2),
  basis_order = NORDER,
  basis_rangeval = range(gufi_meanspec_freqs),
  x_grid = gufi_meanspec_freqs,
  response_matrix = gufi_meanspec_amps
)

plot(gufi_nbasis_gcv$n_basis_seq, gufi_nbasis_gcv$gcv_vals, type = "l")
abline(v = gufi_nbasis_gcv$best_n_basis, col = "red")


# by eye 90 basis seem enough
gufi_basis_int = create.bspline.basis(
  rangeval = range(gufi_meanspec_freqs),
  norder = NORDER,
  nbasis = gufi_nbasis_gcv$best_n_basis
) # too much variability

gufi_meanspec_fd_int = smooth.basis(
  gufi_meanspec_freqs,
  gufi_meanspec_amps,
  gufi_basis_int
)$fd

gufi_meanspec_fd_int$fdnames = list("frequency" = gufi_meanspec_freqs,
                                      "reps" = rep("y", length(gufi_meanspec_fd_int$fdnames$reps)),
                                      "values" = "Amplitude")


# Changing lambda penalty

gufi_LD_gcv <- GCV_DifferentialSmoothBasis(
  n_basis = nrow(gufi_meanspec_amps),
  lambda_grid = 10^seq(-10, -6, length = 100),
  basis_order = NORDER,
  basis_rangeval = range(gufi_meanspec_freqs),
  x_grid = gufi_meanspec_freqs,
  response_matrix = gufi_meanspec_amps
)

plot(log(gufi_LD_gcv$lambda_grid, base = 10), gufi_LD_gcv$gcv_vals, type = "l")

# compare gcv values
min(gufi_nbasis_gcv$gcv_vals) # much smaller
min(gufi_LD_gcv$gcv_vals)

# compare fits

gufi_basis_max = create.bspline.basis(
  rangeval = range(gufi_meanspec_freqs),
  norder = NORDER,
  nbasis = nrow(gufi_meanspec_amps)
)
gufi_fdPar = fdPar(
  fdobj = gufi_basis_max,
  Lfdobj = int2Lfd(2),
  lambda = gufi_LD_gcv$best_lambda
)

gufi_meanspec_fd_diff <- smooth.basis(
  gufi_meanspec_freqs,
  gufi_meanspec_amps,
  gufi_fdPar
)$fd

gufi_meanspec_fd_diff$fdnames = list("frequency" = gufi_meanspec_freqs,
                                    "reps" = rep("y", length(gufi_meanspec_fd_diff$fdnames$reps)),
                                    "values" = "Amplitude")

# compare fitting
gufi_meanspec_fd = gufi_meanspec_fd_diff

# ╭────────╮
# │Gabbiani│ ------------------------------------------------------------------
# ╰────────╯

dim(gabbiani_meanspec_amps)

n_basis_seq = seq(20, nrow(gabbiani_meanspec_amps) - 2)
NORDER = 4
rangeval_gabbiani <- range(gabbiani_meanspec_freqs)


gabbiani_nbasis_gcv <- GCV_IntegerSmoothBasis(
  n_basis_seq = 20:(nrow(gabbiani_meanspec_amps) - 2),
  basis_order = NORDER,
  basis_rangeval = range(gabbiani_meanspec_freqs),
  x_grid = gabbiani_meanspec_freqs,
  response_matrix = gabbiani_meanspec_amps
)

plot(gabbiani_nbasis_gcv$n_basis_seq, gabbiani_nbasis_gcv$gcv_vals, type = "l")


# by eye 100 basis seem enough
gabbiani_basis_int = create.bspline.basis(
  rangeval = range(gabbiani_meanspec_freqs),
  norder = NORDER,
  nbasis = gabbiani_nbasis_gcv$best_n_basis
)

gabbiani_meanspec_fd_int = smooth.basis(
  gabbiani_meanspec_freqs,
  gabbiani_meanspec_amps,
  gabbiani_basis_int
)$fd


gabbiani_meanspec_fd_int$fdnames = list("frequency" = gabbiani_meanspec_freqs,
                                         "reps" = rep("y", length(gabbiani_meanspec_fd_int$fdnames$reps)),
                                         "values" = "Amplitude")

# Changing lambda penalty

gabbiani_LD_gcv <- GCV_DifferentialSmoothBasis(
  n_basis = (nrow(gabbiani_meanspec_amps)),
  lambda_grid = 10^seq(-8, -6, length = 100),
  basis_order = NORDER,
  basis_rangeval = range(gabbiani_meanspec_freqs),
  x_grid = gabbiani_meanspec_freqs,
  response_matrix = gabbiani_meanspec_amps
)

plot(
  log(gabbiani_LD_gcv$lambda_grid, base = 10),
  gabbiani_LD_gcv$gcv_vals,
  type = "l"
)

# compare gcv values
min(gabbiani_nbasis_gcv$gcv_vals)
min(gabbiani_LD_gcv$gcv_vals)

# compare fits

gabbiani_basis_max = create.bspline.basis(
  rangeval = range(gabbiani_meanspec_freqs),
  norder = NORDER,
  nbasis = (nrow(gabbiani_meanspec_amps))
)
gabbiani_fdPar = fdPar(
  fdobj = gabbiani_basis_max,
  Lfdobj = int2Lfd(2),
  lambda = gabbiani_LD_gcv$best_lambda
)

gabbiani_meanspec_fd_diff <- smooth.basis(
  gabbiani_meanspec_freqs,
  gabbiani_meanspec_amps,
  gabbiani_fdPar
)$fd

gabbiani_meanspec_fd_diff$fdnames = list("frequency" = gabbiani_meanspec_freqs,
                                     "reps" = rep("y", length(gabbiani_meanspec_fd_diff$fdnames$reps)),
                                     "values" = "Amplitude")


# >> Regularize fit via Constraint Splines -----------------------

# ╭──────╮
# │Falchi│------------------------------------
# ╰──────╯

falchi_loocv_pen_int = LOOCVConstraintSplinesInt(x_grid = falchi_meanspec_freqs,
                                                 y_matrix = falchi_meanspec_amps,
                                                 basis_num_seq = 20:50,
                                                 box_constraints = c(0,1))


falchi_loocv_pen_diff = LOOCVConstraintSplinesDiff(x_grid = falchi_meanspec_freqs,
                                                   y_matrix = falchi_meanspec_amps,
                                                   basis_num = 110,
                                                   lambda_grid = 10^(seq(-6, -3, length = 20)),
                                                   box_constraints = c(0,1))

falchi_meanspec_fd_con_int = ToFdConstraintSplinesInt(x_grid = falchi_meanspec_freqs,
                                                      y_matrix = falchi_meanspec_amps,
                                                      basis_num = falchi_loocv_pen_int$basis_min)

falchi_meanspec_fd_con_diff = ToFdConstraintSplinesDiff(x_grid = falchi_meanspec_freqs,
                                                        y_matrix = falchi_meanspec_amps,
                                                        basis_num = falchi_loocv_pen_diff$basis_num,
                                                        my.lambda = falchi_loocv_pen_diff$lambda_min)

falchi_meanspec_fd = falchi_meanspec_fd_con_diff

par(mfrow = c(1,2))
plot(falchi_loocv_pen_int$basis_num_seq,
     falchi_loocv_pen_int$loocv_err, type = "b", pch = 16,
     xlab = "n basis",
     ylab = "LOOCV error")

plot(log(falchi_loocv_pen_diff$lambda_grid, base = 10),
     falchi_loocv_pen_diff$loocv_err, type = "b", pch = 16,
     xlab = "log(lambda, base = 10)",
     ylab = "LOOCV error")
par(mfrow = c(1,1))





# ╭────╮
# │Gufi│ ----------------------------------------------------------
# ╰────╯

gufi_loocv_pen_int = LOOCVConstraintSplinesInt(x_grid = gufi_meanspec_freqs,
                                                 y_matrix = gufi_meanspec_amps,
                                                 basis_num_seq = 20:50,
                                                 box_constraints = c(0,1))


gufi_loocv_pen_diff = LOOCVConstraintSplinesDiff(x_grid = gufi_meanspec_freqs,
                                                   y_matrix = gufi_meanspec_amps,
                                                   basis_num = 110,
                                                   lambda_grid = 10^(seq(-10, -3, length = 20)),
                                                   box_constraints = c(0,1))

gufi_meanspec_fd_con_int = ToFdConstraintSplinesInt(x_grid = gufi_meanspec_freqs,
                                                      y_matrix = gufi_meanspec_amps,
                                                      basis_num = gufi_loocv_pen_int$basis_min)

gufi_meanspec_fd_con_diff = ToFdConstraintSplinesDiff(x_grid = gufi_meanspec_freqs,
                                                        y_matrix = gufi_meanspec_amps,
                                                        basis_num = gufi_loocv_pen_diff$basis_num,
                                                        my.lambda = gufi_loocv_pen_diff$lambda_min)
gufi_meanspec_fd = gufi_meanspec_fd_con_diff

par(mfrow = c(1,2))
plot(gufi_loocv_pen_int$basis_num_seq,
     gufi_loocv_pen_int$loocv_err, type = "b", pch = 16,
     xlab = "n basis",
     ylab = "LOOCV error")

plot(log(gufi_loocv_pen_diff$lambda_grid, base = 10),
     gufi_loocv_pen_diff$loocv_err, type = "b", pch = 16,
     xlab = "log(lambda, base = 10)",
     ylab = "LOOCV error")
par(mfrow = c(1,1))





# ╭────────╮
# │Gabbiani│ ------------------------------------------------------------------
# ╰────────╯

gabbiani_loocv_pen_int = LOOCVConstraintSplinesInt(x_grid = gabbiani_meanspec_freqs,
                                               y_matrix = gabbiani_meanspec_amps,
                                               basis_num_seq = 20:50,
                                               box_constraints = c(0,1))


gabbiani_loocv_pen_diff = LOOCVConstraintSplinesDiff(x_grid = gabbiani_meanspec_freqs,
                                                 y_matrix = gabbiani_meanspec_amps,
                                                 basis_num = 110,
                                                 lambda_grid = 10^(seq(-10, -3, length = 20)),
                                                 box_constraints = c(0,1))

gabbiani_meanspec_fd_con_int = ToFdConstraintSplinesInt(x_grid = gabbiani_meanspec_freqs,
                                                    y_matrix = gabbiani_meanspec_amps,
                                                    basis_num = gabbiani_loocv_pen_int$basis_min)

gabbiani_meanspec_fd_con_diff = ToFdConstraintSplinesDiff(x_grid = gabbiani_meanspec_freqs,
                                                      y_matrix = gabbiani_meanspec_amps,
                                                      basis_num = gabbiani_loocv_pen_diff$basis_num,
                                                      my.lambda = gabbiani_loocv_pen_diff$lambda_min)
gabbiani_meanspec_fd = gabbiani_meanspec_fd_con_diff

par(mfrow = c(1,2))
plot(gabbiani_loocv_pen_int$basis_num_seq,
     gabbiani_loocv_pen_int$loocv_err, type = "b", pch = 16,
     xlab = "n basis",
     ylab = "LOOCV error")

plot(log(gabbiani_loocv_pen_diff$lambda_grid, base = 10),
     gabbiani_loocv_pen_diff$loocv_err, type = "b", pch = 16,
     xlab = "log(lambda, base = 10)",
     ylab = "LOOCV error")
par(mfrow = c(1,1))

# .. parameters table ------------------------------------

representation_selection_df = data.frame("animal" = c(rep("falchi", 4), rep("gufi", 4), rep("gabbiani", 4)),
                                         "constraint" = rep(c(FALSE, FALSE, TRUE, TRUE), 3) ,
                                         "penalty type" = rep(c("INT", "DIFF", "INT", "DIFF"), 3) ,
                                         "min error parameter" = c(falchi_nbasis_gcv$best_n_basis,
                                                                   falchi_LD_gcv$best_lambda,
                                                                   falchi_loocv_pen_int$basis_min,
                                                                   falchi_loocv_pen_diff$lambda_min,
                                                                   gufi_nbasis_gcv$best_n_basis,
                                                                   gufi_LD_gcv$best_lambda,
                                                                   gufi_loocv_pen_int$basis_min,
                                                                   gufi_loocv_pen_diff$lambda_min,
                                                                   gabbiani_nbasis_gcv$best_n_basis,
                                                                   gabbiani_LD_gcv$best_lambda,
                                                                   gabbiani_loocv_pen_int$basis_min,
                                                                   gabbiani_loocv_pen_diff$lambda_min),
                                         "domain unique points" = c(rep(nrow(falchi_meanspec_amps), 4),
                                                    rep(nrow(gufi_meanspec_amps), 4),
                                                    rep(nrow(gabbiani_meanspec_amps), 4)))

representation_selection_df

save(representation_selection_df,
     file = "results/prima_parte/outputs/representation_selection_df.RData")

# .. joint plot -------------------------------------------
# here join is meant same species but different criterions

# ╭──────╮
# │Falchi│------------------------------------
# ╰──────╯

png("results/prima_parte/images/falchi_fits_crit.png",
    width = MY.WIDTH, height = MY.HEIGHT)

# compare fitting
par(mfrow = c(2, 2))

plot(falchi_meanspec_fd_int, main = paste0("Falchi GCV Int - nbasis: ",
                                           falchi_meanspec_fd_int$basis$nbasis, collapse = "")) # border problem
plot(falchi_meanspec_fd_diff, main = paste0("Falchi GCV Diff - log(lambda): ",
                                            round(log(falchi_LD_gcv$best_lambda, base = 10),2), collapse = "")) # choose this one


plot(falchi_meanspec_fd_con_int,
     main = paste0("Falchi Constraint LOOCV Int - nbasis: ",
                   falchi_loocv_pen_int$basis_min, collapse = ""))
plot(falchi_meanspec_fd_con_diff,
     main = paste0("Falchi Constraint LOOCV Diff - log(lambda): ",
                   round(log(falchi_loocv_pen_diff$lambda_min, base = 10),2), collapse = ""))
par(mfrow = c(1,1))

dev.off()

# ╭────╮
# │Gufi│ ----------------------------------------------------------
# ╰────╯

png("results/prima_parte/images/gufi_fits_crit.png",
    width = MY.WIDTH, height = MY.HEIGHT)

par(mfrow = c(2, 2))

plot(gufi_meanspec_fd_int, main = paste0("Gufi GCV Int - nbasis: ",
                                           gufi_meanspec_fd_int$basis$nbasis, collapse = "")) # border problem
plot(gufi_meanspec_fd_diff, main = paste0("Gufi GCV Diff - log(lambda): ",
                                            round(log(gufi_LD_gcv$best_lambda, base = 10),2), collapse = "")) # choose this one


plot(gufi_meanspec_fd_con_int,
     main = paste0("Gufi Constraint LOOCV Int - nbasis: ",
                   gufi_loocv_pen_int$basis_min, collapse = ""))
plot(gufi_meanspec_fd_con_diff,
     main = paste0("Gufi Constraint LOOCV Diff - log(lambda): ",
                   round(log(gufi_loocv_pen_diff$lambda_min, base = 10),2), collapse = ""))
par(mfrow = c(1,1))

dev.off()

# ╭────────╮
# │Gabbiani│ ------------------------------------------------------------------
# ╰────────╯

png("results/prima_parte/images/gabbiani_fits_crit.png",
    width = MY.WIDTH, height = MY.HEIGHT)

par(mfrow = c(2, 2))

plot(gabbiani_meanspec_fd_int, main = paste0("Gabbiani GCV Int - nbasis: ",
                                           gabbiani_meanspec_fd_int$basis$nbasis, collapse = "")) # border problem
plot(gabbiani_meanspec_fd_diff, main = paste0("Gabbiani GCV Diff - log(lambda): ",
                                            round(log(gabbiani_LD_gcv$best_lambda, base = 10),2), collapse = "")) # choose this one


plot(gabbiani_meanspec_fd_con_int,
     main = paste0("Gabbiani Constraint LOOCV Int - nbasis: ",
                   gabbiani_loocv_pen_int$basis_min, collapse = ""))
plot(gabbiani_meanspec_fd_con_diff,
     main = paste0("Gabbiani Constraint LOOCV Diff - log(lambda): ",
                   round(log(gabbiani_loocv_pen_diff$lambda_min, base = 10),2), collapse = ""))
par(mfrow = c(1,1))

dev.off()

# .. Save image ----------------------------
gc()
save.image(file = "results/prima_parte/outputs/basis_selection_work_space.RData")

# >>f Means ---------------------------------

# ╭──────╮
# │Falchi│------------------------------------
# ╰──────╯


falchi_meanspec_fd_mean = lapply(
  levels(falchi$Climate_zone),
  function(i) mean.fd(falchi_meanspec_fd[falchi$Climate_zone == i])
)

falchi_meanspec_fd_sd = lapply(
  levels(falchi$Climate_zone),
  function(i) sd.fd(falchi_meanspec_fd[falchi$Climate_zone == i])
)


# ╭────╮
# │Gufi│ ----------------------------------------------------------
# ╰────╯

gufi_meanspec_fd_mean = lapply(
  levels(gufi$Climate_zone),
  function(i) mean.fd(gufi_meanspec_fd[gufi$Climate_zone == i])
)

gufi_meanspec_fd_sd = lapply(
  levels(gufi$Climate_zone),
  function(i) sd.fd(gufi_meanspec_fd[gufi$Climate_zone == i])
)



# ╭────────╮
# │Gabbiani│ ------------------------------------------------------------------
# ╰────────╯


gabbiani_meanspec_fd_mean = lapply(
  levels(gabbiani$Cluster),
  function(i) mean.fd(gabbiani_meanspec_fd[gabbiani$Cluster == i])
)

gabbiani_meanspec_fd_sd = lapply(
  levels(gabbiani$Cluster),
  function(i) sd.fd(gabbiani_meanspec_fd[gabbiani$Cluster == i])
)

# ..joint plot ---------------


png("results/prima_parte/images/f_mean_sd.png",
     width = MY.WIDTH, height = MY.HEIGHT)

par(mfrow = c(1,3))

FunctionalMeanBandPlot(fd_means = falchi_meanspec_fd_mean,
                       fd_sds = falchi_meanspec_fd_sd,
                       my.main = "Medie Funzionali ± Sd Falchi",
                       my.ylim = c(0,0.8),
                       my.lwd = 2,
                       my.xlab = "Frequenza",
                       my.ylab = "Ampiezza",
                       my.levels.variable = falchi$Climate_zone,
                       save_path = NULL,
                       my.width = 600,
                       my.height = 350)

FunctionalMeanBandPlot(fd_means = gufi_meanspec_fd_mean,
                       fd_sds = gufi_meanspec_fd_sd,
                       my.main = "Medie Funzionali ± Sd Gufi",
                       my.ylim = c(0,1),
                       my.lwd = 2,
                       my.xlab = "Frequenza",
                       my.ylab = "Ampiezza",
                       my.levels.variable = gufi$Climate_zone,
                       save_path = NULL,
                       my.width = 600,
                       my.height = 350)

FunctionalMeanBandPlot(fd_means = gabbiani_meanspec_fd_mean,
                       fd_sds = gabbiani_meanspec_fd_sd,
                       my.main = "Medie Funzionali ± Sd Gabbiani",
                       my.ylim = c(0,1),
                       my.lwd = 2,
                       my.xlab = "Frequenza",
                       my.ylab = "Ampiezza",
                       my.levels.variable = gabbiani$Cluster,
                       save_path = NULL,
                       my.width = 600,
                       my.height = 350)

par(mfrow = c(1,1))

dev.off()


# >>fPCA ---------------------------------------------------------

# 10 components
falchi_pcf = pca.fd(falchi_meanspec_fd, nharm = 10)
gufi_pcf = pca.fd(gufi_meanspec_fd, nharm = 10)
gabbiani_pcf = pca.fd(gabbiani_meanspec_fd, nharm = 10)


png("results/prima_parte/images/f_pca_explained_var.png",
     width = MY.WIDTH, height = MY.HEIGHT)

par(mfrow = c(3, 1))
plot(cumsum(falchi_pcf$varprop), type = "b",
     pch = 16, main = "fPCA varianza spiegata cumulata Falchi")
plot(cumsum(gufi_pcf$varprop), type = "b",
     pch = 16, main = "fPCA varianza spiegata cumulata Gufi")
plot(cumsum(gabbiani_pcf$varprop), type = "b",
     pch = 16, main = "fPCA varianza spiegata cumulata Gabbiani")

dev.off()

par(mfrow = c(1, 1))


# 3 components -> see mean +- eps*harmonic
falchi_pcf = pca.fd(falchi_meanspec_fd, nharm = 3)
gufi_pcf = pca.fd(gufi_meanspec_fd, nharm = 3)
gabbiani_pcf = pca.fd(gabbiani_meanspec_fd, nharm = 3)


png("results/prima_parte/images/f_pca_harmonics.png",
     width = MY.WIDTH, height = MY.HEIGHT)

par(mfrow = c(3, 3), mar = rep(0, 4) + 1)
plot(falchi_pcf)
plot(gufi_pcf)
plot(gabbiani_pcf)
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))

dev.off()

# 2 components -> see scores
falchi_pcf = pca.fd(falchi_meanspec_fd, nharm = 2)
gufi_pcf = pca.fd(gufi_meanspec_fd, nharm = 2)
gabbiani_pcf = pca.fd(gabbiani_meanspec_fd, nharm = 2)


png("results/prima_parte/images/f_pca_scores.png",
     width = MY.WIDTH, height = MY.HEIGHT)

par(mfrow = c(3, 1))
plot(falchi_pcf$scores, type = "p", col = falchi$Climate_zone, pch = 16,
     main = "fPCA punteggi Falchi")
legend("topright", legend = levels(falchi$Climate_zone), col = 1:3, pch = 16,
       bty = "n")
plot(gufi_pcf$scores, type = "p", col = gufi$Climate_zone, pch = 16,
     main = "fPCA punteggi Gufi")
legend("topright", legend = levels(gufi$Climate_zone), col = 1:3, pch = 16,
       bty = "n")
plot(gabbiani_pcf$scores, type = "p", col = gabbiani$Cluster, pch = 16,
     main = "fPCA punteggi Gabbiani")
legend("topright", legend = levels(gabbiani$Cluster), col = 1:4, pch = 16,
       bty = "n")

par(mfrow = c(1, 1))

dev.off()


# >>f Clustering --------------------------------------

# >>f ANOVA ---------------------------------

# ╭──────╮
# │Falchi│------------------------------------
# ╰──────╯

load("results/prima_parte/outputs/cv_fanova_res_falchi.RData")

set.seed(123)

# select lambda 
cv_fanova_res_falchi = CvFunctionalANOVA(factor = falchi$Climate_zone,
                  X = falchi_meanspec_amps,
                  dom = falchi_meanspec_freqs,
                  basis_y = falchi_meanspec_fd$basis,
                  coef_y = falchi_meanspec_fd$coefs,
                  y_names = falchi_meanspec_fd_diff$fdnames,
                  basis_beta = falchi_meanspec_fd_diff$basis,
                  lambda_grid = 10^seq(-10, 2, by = 0.5),
                  nfold = 5)

save(cv_fanova_res_falchi,
     file = "results/prima_parte/outputs/cv_fanova_res_falchi.RData")


falchi_bs_beta = BootBetaIC(B = 1000,
                            factor = falchi$Climate_zone,
                            X = falchi_meanspec_amps,
                            dom = falchi_meanspec_freqs,
                            basis_y = falchi_meanspec_fd$basis,
                            coef_y = falchi_meanspec_fd$coefs,
                            y_names = falchi_meanspec_fd_diff$fdnames,
                            basis_beta = falchi_meanspec_fd_diff$basis,
                            lambda = 0.01,
                            n_discrete = 500)

save(boot_fanova_beta_falchi,
     file = "results/prima_parte/outputs/boot_fanova_beta_falchi.RData")


# fit model + beta se
# falchi_anova_model = ComputeBetaSd(factor = falchi$Climate_zone,
#                                X = falchi_meanspec_amps,
#                                dom = falchi_meanspec_freqs,
#                                basis_y = falchi_meanspec_fd$basis,
#                                coef_y = falchi_meanspec_fd$coefs,
#                                y_names = falchi_meanspec_fd_diff$fdnames,
#                                basis_beta = falchi_meanspec_fd_diff$basis,
#                                lambda = cv_fanova_res_falchi$lambda_min)

gc()

# ftest
perm_fanova_res_falchi = PermutFANOVA(factor = falchi$Climate_zone,
                                      X = falchi_meanspec_amps,
                                      dom = falchi_meanspec_freqs,
                                      basis_y = falchi_meanspec_fd$basis,
                                      coef_y = falchi_meanspec_fd$coefs,
                                      y_names = falchi_meanspec_fd_diff$fdnames,
                                      basis_beta = falchi_meanspec_fd_diff$basis,
                                      lambda = 0.01,
                                      n_perm = 1000,
                                      seed = 123)

gc()


# fANOVABetaSdPlot(my.betaestlist = falchi_anova_model$model$betaestlist,
#                  my.betastderrlist = falchi_anova_model$beta_se$betastderrlist,
#                  my.factor = falchi$Climate_zone,
#                  my.name = "Falchi",
#                  save_path = "results/prima_parte/images/f_beta_falchi.png",
#                  my.width = MY.WIDTH,
#                  my.height = MY.HEIGHT,
#                  my.layout.matr = cbind(matrix(1, 2, 2),
#                                         matrix(2:5, 2, 2)))



PlotBetaWithQuantiles(original_fit = falchi_bs_beta$original_fit_beta,
                      quantile_betas = falchi_bs_beta$quantile_betas,
                      my.name = "Falchi",
                      save_path = "results/prima_parte/images/f_beta_quant_falchi.png",
                      my.width = MY.WIDTH,
                      my.height = MY.HEIGHT,
                      my.layout.matr = cbind(matrix(1, 2, 2),
                                             matrix(2:5, 2, 2)))

# ╭────╮
# │Gufi│ ----------------------------------------------------------
# ╰────╯

load("results/prima_parte/outputs/cv_fanova_res_gufi.RData")

cv_fanova_res_gufi = CvFunctionalANOVA(factor = gufi$Climate_zone,
                                         X = gufi_meanspec_amps,
                                         dom = gufi_meanspec_freqs,
                                         basis_y = gufi_meanspec_fd_diff$basis,
                                         coef_y = gufi_meanspec_fd_diff$coefs,
                                         y_names = gufi_meanspec_fd_diff$fdnames,
                                         basis_beta = gufi_meanspec_fd_diff$basis,
                                         lambda_grid = 10^seq(-10, 2, by = 0.5),
                                         nfold = 5)

save(cv_fanova_res_gufi,
     file = "results/prima_parte/outputs/cv_fanova_res_gufi.RData")


# fit model + beta se
# gufi_anova_model = ComputeBetaSd(factor = gufi$Climate_zone,
#                                    X = gufi_meanspec_amps,
#                                    dom = gufi_meanspec_freqs,
#                                    basis_y = gufi_meanspec_fd_diff$basis,
#                                    coef_y = gufi_meanspec_fd_diff$coefs,
#                                    y_names = gufi_meanspec_fd_diff$fdnames,
#                                    basis_beta = gufi_meanspec_fd_diff$basis,
#                                    lambda = cv_fanova_res_gufi$lambda_min)

gc()

boot_fanova_beta_gufi = BootBetaIC(B = 1000,
                            factor = gufi$Climate_zone,
                            X = gufi_meanspec_amps,
                            dom = gufi_meanspec_freqs,
                            basis_y = gufi_meanspec_fd$basis,
                            coef_y = gufi_meanspec_fd$coefs,
                            y_names = gufi_meanspec_fd_diff$fdnames,
                            basis_beta = gufi_meanspec_fd_diff$basis,
                            lambda = cv_fanova_res_gufi$lambda_min,
                            n_discrete = 500)

save(boot_fanova_beta_gufi,
     file = "results/prima_parte/outputs/boot_fanova_beta_gufi.RData")

# ftest
perm_fanova_res_gufi = PermutFANOVA(factor = gufi$Climate_zone,
                                      X = gufi_meanspec_amps,
                                      dom = gufi_meanspec_freqs,
                                      basis_y = gufi_meanspec_fd_diff$basis,
                                      coef_y = gufi_meanspec_fd_diff$coefs,
                                      y_names = gufi_meanspec_fd_diff$fdnames,
                                      basis_beta = gufi_meanspec_fd_diff$basis,
                                      lambda = cv_fanova_res_gufi$lambda_min,
                                      n_perm = 1000,
                                      seed = 123)

gc()


# fANOVABetaSdPlot(my.betaestlist = gufi_anova_model$model$betaestlist,
#                  my.betastderrlist = gufi_anova_model$beta_se$betastderrlist,
#                  my.factor = gufi$Climate_zone,
#                  my.name = "Gufi",
#                  save_path = "results/prima_parte/images/f_beta_gufi.png",
#                  my.width = MY.WIDTH,
#                  my.height = MY.HEIGHT,
#                  my.layout.matr = cbind(matrix(1, 2, 2),
#                                         matrix(2:5, 2, 2)))


PlotBetaWithQuantiles(original_fit = gufi_bs_beta$original_fit_beta,
                      quantile_betas = gufi_bs_beta$quantile_betas,
                      my.name = "Gufi",
                      save_path = "results/prima_parte/images/f_beta_quant_gufi.png",
                      my.width = MY.WIDTH,
                      my.height = MY.HEIGHT,
                      my.layout.matr = cbind(matrix(1, 2, 2),
                                             matrix(2:5, 2, 2)))



# ╭────────╮
# │Gabbiani│ ------------------------------------------------------------------
# ╰────────╯

load("results/prima_parte/outputs/cv_fanova_res_gabbiani.RData")

cv_fanova_res_gabbiani = CvFunctionalANOVA(factor = gabbiani$Cluster,
                                       X = gabbiani_meanspec_amps,
                                       dom = gabbiani_meanspec_freqs,
                                       basis_y = gabbiani_meanspec_fd_diff$basis,
                                       coef_y = gabbiani_meanspec_fd_diff$coefs,
                                       y_names = gabbiani_meanspec_fd_diff$fdnames,
                                       basis_beta = gabbiani_meanspec_fd_diff$basis,
                                       lambda_grid = 10^seq(-10, 2, by = 0.5),
                                       nfold = 5)

save(cv_fanova_res_gabbiani,
     file = "results/prima_parte/outputs/cv_fanova_res_gabbiani.RData")


gc()


boot_fanova_beta_gabbiani = BootBetaIC(B = 1000,
                          factor = gabbiani$Cluster,
                          X = gabbiani_meanspec_amps,
                          dom = gabbiani_meanspec_freqs,
                          basis_y = gabbiani_meanspec_fd$basis,
                          coef_y = gabbiani_meanspec_fd$coefs,
                          y_names = gabbiani_meanspec_fd_diff$fdnames,
                          basis_beta = gabbiani_meanspec_fd_diff$basis,
                          lambda = cv_fanova_res_gabbiani$lambda_min,
                          n_discrete = 500)

save(boot_fanova_beta_gabbiani,
     file = "results/prima_parte/outputs/boot_fanova_beta_gabbiani.RData")

# fit model + beta se
# gabbiani_anova_model = ComputeBetaSd(factor = gabbiani$Cluster,
#                                  X = gabbiani_meanspec_amps,
#                                  dom = gabbiani_meanspec_freqs,
#                                  basis_y = gabbiani_meanspec_fd_diff$basis,
#                                  coef_y = gabbiani_meanspec_fd_diff$coefs,
#                                  y_names = gabbiani_meanspec_fd_diff$fdnames,
#                                  basis_beta = gabbiani_meanspec_fd_diff$basis,
#                                  lambda = cv_fanova_res_gabbiani$lambda_min)

# ftest
perm_fanova_res_gabbiani = PermutFANOVA(factor = gabbiani$Cluster,
                                    X = gabbiani_meanspec_amps,
                                    dom = gabbiani_meanspec_freqs,
                                    basis_y = gabbiani_meanspec_fd_diff$basis,
                                    coef_y = gabbiani_meanspec_fd_diff$coefs,
                                    y_names = gabbiani_meanspec_fd_diff$fdnames,
                                    basis_beta = gabbiani_meanspec_fd_diff$basis,
                                    lambda = cv_fanova_res_gabbiani$lambda_min,
                                    n_perm = 1000,
                                    seed = 123)

gc()


# fANOVABetaSdPlot(my.betaestlist = gabbiani_anova_model$model$betaestlist,
#                  my.betastderrlist = gabbiani_anova_model$beta_se$betastderrlist,
#                  my.factor = gabbiani$Cluster,
#                  my.name = "gabbiani",
#                  save_path = "results/prima_parte/images/f_beta_gabbiani.png",
#                  my.width = MY.WIDTH,
#                  my.height = MY.HEIGHT,
#                  my.layout.matr = cbind(matrix(1, 2, 2),
#                                         matrix(2:5, 2, 2)))


PlotBetaWithQuantiles(original_fit = boot_fanova_beta_gabbiani$original_fit_beta,
                      quantile_betas = boot_fanova_beta_gabbiani$quantile_betas,
                      my.name = "Gabbiani",
                      save_path = "results/prima_parte/images/f_beta_quant_gabbiani.png",
                      my.width = MY.WIDTH,
                      my.height = MY.HEIGHT,
                      my.layout.matr = cbind(matrix(1, 2, 2),
                                             matrix(2:5, 2, 2)))

# .. joint plot -----------------------------

#.. cv err ------


png("results/prima_parte/images/f_anova_cv_err.png",
     width = MY.WIDTH, height = MY.HEIGHT)

par(mfrow = c(3, 1))

plot(cv_fanova_res_falchi$lambda_grid,
     cv_fanova_res_falchi$cv_error,
     pch = 16, type = "b",
     xlab = "lambda",
     ylab = "err",
     main = "Errore integrato di Cv - Falchi")

plot(cv_fanova_res_gufi$lambda_grid,
     cv_fanova_res_gufi$cv_error,
     pch = 16, type = "b",
     xlab = "lambda",
     ylab = "err",
     main = "Errore integrato di Cv - Gufi")


plot(cv_fanova_res_gabbiani$lambda_grid,
     cv_fanova_res_gabbiani$cv_error,
     pch = 16, type = "b",
     xlab = "lambda",
     ylab = "err",
     main = "Errore integrato di Cv - Gabbiani")

dev.off()

# .. f test --------
par(mfrow = c(3, 1))

png("results/prima_parte/images/f_anova_f_test.png",
     width = MY.WIDTH, height = MY.HEIGHT)

# falchi
plot(
  falchi_meanspec_freqs,
  perm_fanova_res_falchi$fobs,
  type = 'l',
  lwd = 2,
  main = "Ftest Falchi",
  xlab = "Frequenza",
  ylab = "F"
)
lines(falchi_meanspec_freqs, perm_fanova_res_falchi$qF, lty = 3)
abline(h = perm_fanova_res_falchi$qFmax, lwd = 2, col = 2)

# gufi
plot(
  gufi_meanspec_freqs,
  perm_fanova_res_gufi$fobs,
  type = 'l',
  lwd = 2,
  ylim = c(0, 1.1),
  main = "Ftest Gufi",
  xlab = "Frequenza",
  ylab = "F"
)
lines(gufi_meanspec_freqs, perm_fanova_res_gufi$qF, lty = 3)
abline(h = perm_fanova_res_gufi$qFmax, lwd = 2, col = 2)


# gabbiani
plot(
  gabbiani_meanspec_freqs,
  perm_fanova_res_gabbiani$fobs,
  type = 'l',
  lwd = 2,
  ylim = c(0, 1.1),
  main = "Ftest gabbiani",
  xlab = "Frequenza",
  ylab = "F"
)
lines(gabbiani_meanspec_freqs, perm_fanova_res_gabbiani$qF, lty = 3)
abline(h = perm_fanova_res_gabbiani$qFmax, lwd = 2, col = 2)

dev.off()


par(mfrow = c(1, 1))


