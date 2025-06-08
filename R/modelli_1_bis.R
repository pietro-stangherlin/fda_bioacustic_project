library(av)
library(tuneR)
library(seewave)
library(soundgen)
library(fda)
library(plotly)
library(tidyverse)
library(quadprog)

rm(list = ls())


# >> Data Reading ---------------------------------
# WARNING: maybe change the input name
load("data/data_1.RData")

# subsample gufi due to computational issues
gufi_ids_exclude = read.table("data/gufi_ids_exclude.txt",
                              header = T)

gufi_meanspec_amps = gufi_meanspec_amps[,-gufi_ids_exclude$x]
gufi = gufi[-gufi_ids_exclude$x,]

# sagufi_meanspec_amps# save(
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
# on a single curve
# (leave one point at a time)
LOOCVConstraintSPlinesSingleCurve = function(y_vals, # vector
                                            B_observed_y, # basis evaluated on observed y values
                                            B_constraint, # basis evaluated on dense grid: used in constraint
                                            Dmat, # quadratic penalty component
                                            Amat, # constraint matrix,
                                            bvec, # constraint vector
                                            box_constraints = c(0,1)){
  
  dvec_all <- t(B_observed_y) %*% y_vals
  
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
  
  return(error_grid)
  
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
    lambda = 1
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
  
  return(error_grid)
  
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
              "beta_se" = fRegress.stderr(
                m1,
                smooth.basis(dom, X, basis_y)$y2cMap,
                SigmaE)
              )
    )
  
}

PermutFANOVA = function(factor,
                          X,
                          dom,
                          basis_y,
                          coef_y,
                          y_names,
                          basis_beta,
                          lambda,
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
    "qF" = apply(Fmatr, 2, quantile, .95),
    "qFmax" = quantile(apply(Fmatr, 1, max), .95)
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
MY.WIDTH = 1500
MY.HEIGHT = 1000



# >> Regularize fit via Constraint Splines -----------------------

# ╭──────╮
# │Falchi│------------------------------------
# ╰──────╯


falchi_basis_seq = 20:50

falchi_loocv_pen_int = LOOCVConstraintSplinesInt(x_grid = falchi_meanspec_freqs,
                          y_matrix = falchi_meanspec_amps,
                          basis_num_seq = falchi_basis_seq,
                          box_constraints = c(0,1))

plot(falchi_basis_seq, falchi_loocv_pen_int, type = "l")

falchi_lambda_grid = 10^(seq(-6, -3, length = 20))

falchi_loocv_pen_diff = LOOCVConstraintSplinesDiff(x_grid = falchi_meanspec_freqs,
                                                 y_matrix = falchi_meanspec_amps,
                                                 basis_num = 110, # increment
                                                 lambda_grid = falchi_lambda_grid,
                                                 box_constraints = c(0,1))


plot(log(falchi_lambda_grid, base = 10), falchi_loocv_pen_diff,
     type = "b", pch = 16)

nbasis = 80

basis <- create.bspline.basis(range(falchi_meanspec_freqs),nbasis)
B = eval.basis(falchi_meanspec_freqs, basis)

BasisCoefMatr = matrix(NA,
                       nrow = nbasis,
                       ncol = NCOL(falchi_meanspec_amps))


temp_coef = ConstraintSplinesFit(x_grid = falchi_meanspec_freqs,
                                 y_vals = falchi_meanspec_amps[,1],
                                 basis = basis,
                                 box_constraints = c(0,1))$coef

BasisCoefMatr[,1] = temp_coef



x_fine <- seq(min(falchi_meanspec_freqs),
              max(falchi_meanspec_freqs), length.out = 500)
B_fine <- eval.basis(x_fine, basis)
y_fine <- B_fine %*% temp_coef

plot(x_fine, y_fine, col = 1,
     type = "l",
     ylim = c(0,1.1))

max_y = max(y_fine)

for(i in 2:NCOL(falchi_meanspec_amps)){
  temp_coef = ConstraintSplinesFit(x_grid = falchi_meanspec_freqs,
                                   y_vals = falchi_meanspec_amps[,i],
                                   basis = basis,
                                   box_constraints = c(0,1))$coef
  
  
  BasisCoefMatr[,i] = temp_coef
  
  y_fine <- B_fine %*% as.matrix(temp_coef)
  max_y = max(c(max_y, y_fine))
  
  lines(x_fine, y_fine, col = i)
}

max_y
abline(h = 1, lwd = 1)

falchi_meanspec_fd = fd(coef = BasisCoefMatr,
                        basisobj = basis)

# ╭────╮
# │Gufi│ ----------------------------------------------------------
# ╰────╯

basis <- create.bspline.basis(range(gufi_meanspec_freqs),110)
B = eval.basis(gufi_meanspec_freqs, basis)

temp_coef = ConstraintSplinesFit(x_grid = gufi_meanspec_freqs,
                                 y_vals = gufi_meanspec_amps[,1],
                                 basis = basis,
                                 box_constraints = c(0,1))$coef


x_fine <- seq(min(gufi_meanspec_freqs),
              max(gufi_meanspec_freqs), length.out = 1000)
B_fine <- eval.basis(x_fine, basis)


y_fine <- B_fine %*% as.matrix(temp_coef)

plot(x_fine, y_fine, col = 1,
     type = "l",
     ylim = c(0,1.1))

max_y = max(y_fine)

for(i in 2:NCOL(gufi_meanspec_amps)){
  temp_coef = ConstraintSplinesFit(x_grid = gufi_meanspec_freqs,
                                   y_vals = gufi_meanspec_amps[,i],
                                   basis = basis,
                                   box_constraints = c(0,1))$coef
  
  y_fine <- B_fine %*% temp_coef
  max_y = max(c(max_y, y_fine))
  
  lines(x_fine, y_fine, col = i)
}

# abline(h = 1, lwd = 1)

max_y



# ╭────────╮
# │Gabbiani│ ------------------------------------------------------------------
# ╰────────╯

basis <- create.bspline.basis(range(gabbiani_meanspec_freqs),110)
B = eval.basis(gabbiani_meanspec_freqs, basis)

temp_coef = ConstraintSplinesFit(x_grid = gabbiani_meanspec_freqs,
                                 y_vals = gabbiani_meanspec_amps[,1],
                                 basis = basis,
                                 box_constraints = c(0,1))$coef


x_fine <- seq(min(gabbiani_meanspec_freqs),
              max(gabbiani_meanspec_freqs), length.out = 1000)
B_fine <- eval.basis(x_fine, basis)


y_fine <- B_fine %*% as.matrix(temp_coef)

plot(x_fine, y_fine, col = 1,
     type = "l",
     ylim = c(0,1.1))

max_y = max(y_fine)

for(i in 2:NCOL(gabbiani_meanspec_amps)){
  temp_coef = ConstraintSplinesFit(x_grid = gabbiani_meanspec_freqs,
                                   y_vals = gabbiani_meanspec_amps[,i],
                                   basis = basis,
                                   box_constraints = c(0,1))$coef
  
  y_fine <- B_fine %*% temp_coef
  max_y = max(c(max_y, y_fine))
  
  lines(x_fine, y_fine, col = i)
}

abline(h = 1, lwd = 2)

max_y


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
  nbasis = 64
)

falchi_meanspec_fd_int = smooth.basis(
  falchi_meanspec_freqs,
  falchi_meanspec_amps,
  falchi_basis_int
)$fd


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

falchi_meanspec_fd_diff <- smooth.basis(
  falchi_meanspec_freqs,
  falchi_meanspec_amps,
  falchi_fdPar
)$fd

falchi_basis_diff <- smooth.basis(
  falchi_meanspec_freqs,
  falchi_meanspec_amps,
  falchi_fdPar
)

diag(falchi_basis_diff$penmat)


# compare fitting
par(mfrow = c(1, 2))

plot(falchi_meanspec_fd_int, main = "Int penalty") # border problem
plot(falchi_meanspec_fd_diff, main = "Diff penalty") # choose this one

par(mfrow = c(1, 1))

falchi_meanspec_fd = falchi_meanspec_fd_diff


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
  nbasis = 90
) # too much variability

gufi_meanspec_fd_int = smooth.basis(
  gufi_meanspec_freqs,
  gufi_meanspec_amps,
  gufi_basis_int
)$fd


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
  nbasis = row(gufi_meanspec_amps)
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

# compare fitting
par(mfrow = c(1, 2))

plot(gufi_meanspec_fd_int, main = "Int penalty") # border problem
plot(gufi_meanspec_fd_diff, main = "Diff penalty") # choose this one

par(mfrow = c(1, 1))

gufi_meanspec_fd = gufi_meanspec_fd_diff

# ╭────────╮
# │Gabbiani│ ------------------------------------------------------------------
# ╰────────╯

dim(gabbiani_meanspec_amps)

# test
gabbiani_meanspec_amps = gabbiani_meanspec_amps[,which(gabbiani_meanspec_amps[1,] > 0.1)]

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
  nbasis = 100
)

gabbiani_meanspec_fd_int = smooth.basis(
  gabbiani_meanspec_freqs,
  gabbiani_meanspec_amps,
  gabbiani_basis_int
)$fd


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

# compare fitting
par(mfrow = c(1, 2))
# here the problematic is inverted
plot(gabbiani_meanspec_fd_int, main = "Int penalty") # choose this one
plot(gabbiani_meanspec_fd_diff, main = "Diff penalty") # border problem

par(mfrow = c(1, 1))

gabbiani_meanspec_fd = gabbiani_meanspec_fd_int


# >>f Means ---------------------------------

# ╭──────╮
# │Falchi│------------------------------------
# ╰──────╯


falchi_meanspec_fd_mean = lapply(
  levels(falchi$Climate_zone),
  function(i) mean.fd(falchi_meanspec_fd[falchi$Climate_zone == i])
)

# Deviazioni standard funzionali
falchi_freqs_grid = seq(
  rangeval_falchi[1],
  rangeval_falchi[2],
  length.out = 201
)
falchi_sd_fd_matrix = eval.fd(falchi_freqs_grid, falchi_meanspec_fd)


falchi_meanspec_fd_sd = lapply(levels(falchi$Climate_zone), function(i) {
  Data2fd(
    argvals = falchi_freqs_grid,
    y = apply(falchi_sd_fd_matrix[, falchi$Climate_zone == i], 1, sd),
    basisobj = falchi_fdPar$fd$basis,
    nderiv = falchi_fdPar$Lfd$nderiv,
    lambda = falchi_fdPar$lambda
  )
})


# ╭────╮
# │Gufi│ ----------------------------------------------------------
# ╰────╯

gufi_meanspec_fd_mean = lapply(
  levels(gufi$Climate_zone),
  function(i) mean.fd(gufi_meanspec_fd[gufi$Climate_zone == i])
)

# Deviazioni standard funzionali
gufi_freqs_grid = seq(
  rangeval_gufi[1],
  rangeval_gufi[2],
  length.out = 201
)

gufi_sd_fd_matrix = eval.fd(gufi_freqs_grid, gufi_meanspec_fd)

gufi_meanspec_fd_sd = lapply(levels(gufi$Climate_zone), function(i) {
  Data2fd(
    argvals = gufi_freqs_grid,
    y = apply(gufi_sd_fd_matrix[, gufi$Climate_zone == i], 1, sd),
    basisobj = gufi_fdPar$fd$basis,
    nderiv = gufi_fdPar$Lfd$nderiv,
    lambda = gufi_fdPar$lambda
  )
})



# ╭────────╮
# │Gabbiani│ ------------------------------------------------------------------
# ╰────────╯


gabbiani_meanspec_fd_mean = lapply(
  levels(gabbiani$Cluster),
  function(i) mean.fd(gabbiani_meanspec_fd[gabbiani$Cluster == i])
)

# Deviazioni standard funzionali
gabbiani_freqs_grid = seq(
  rangeval_gabbiani[1],
  rangeval_gabbiani[2],
  length.out = 201
)
gabbiani_sd_fd_matrix = eval.fd(gabbiani_freqs_grid, gabbiani_meanspec_fd)


gabbiani_meanspec_fd_sd = lapply(levels(gabbiani$Cluster), function(i) {
  Data2fd(
    argvals = gabbiani_freqs_grid,
    y = apply(gabbiani_sd_fd_matrix[, gabbiani$Cluster == i], 1, sd),
    basisobj = gabbiani_basis_int,
  )
})


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
                       my.ylim = c(0,0.8),
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
                       my.ylim = c(0,0.8),
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

set.seed(123)

# select lambda 
cv_fanova_res_falchi = CvFunctionalANOVA(factor = falchi$Climate_zone,
                  X = falchi_meanspec_amps,
                  dom = falchi_meanspec_freqs,
                  basis_y = falchi_meanspec_fd$basis,
                  coef_y = falchi_meanspec_fd$coefs,
                  y_names = falchi_meanspec_fd_diff$fdnames,
                  basis_beta = falchi_meanspec_fd_diff$basis,
                  lambda_grid = 1:5,
                  nfold = 5)

# fit model + beta se
falchi_anova_model = ComputeBetaSd(factor = falchi$Climate_zone,
                               X = falchi_meanspec_amps,
                               dom = falchi_meanspec_freqs,
                               basis_y = falchi_meanspec_fd$basis,
                               coef_y = falchi_meanspec_fd$coefs,
                               y_names = falchi_meanspec_fd_diff$fdnames,
                               basis_beta = falchi_meanspec_fd_diff$basis,
                               lambda = cv_fanova_res_falchi$lambda_min)

gc()

# ftest
perm_fanova_res_falchi = PermutFANOVA(factor = falchi$Climate_zone,
                                      X = falchi_meanspec_amps,
                                      dom = falchi_meanspec_freqs,
                                      basis_y = falchi_meanspec_fd$basis,
                                      coef_y = falchi_meanspec_fd$coefs,
                                      y_names = falchi_meanspec_fd_diff$fdnames,
                                      basis_beta = falchi_meanspec_fd_diff$basis,
                                      lambda = cv_fanova_res_falchi$lambda_min,
                                      n_perm = 100,
                                      seed = 123)

gc()


fANOVABetaSdPlot(my.betaestlist = falchi_anova_model$model$betaestlist,
                 my.betastderrlist = falchi_anova_model$beta_se$betastderrlist,
                 my.factor = falchi$Climate_zone,
                 my.name = "Falchi",
                 save_path = "results/prima_parte/images/f_beta_falchi.png",
                 my.width = MY.WIDTH,
                 my.height = MY.HEIGHT,
                 my.layout.matr = cbind(matrix(1, 2, 2),
                                        matrix(2:5, 2, 2)))



# ╭────╮
# │Gufi│ ----------------------------------------------------------
# ╰────╯


cv_fanova_res_gufi = CvFunctionalANOVA(factor = gufi$Climate_zone,
                                         X = gufi_meanspec_amps,
                                         dom = gufi_meanspec_freqs,
                                         basis_y = gufi_meanspec_fd_diff$basis,
                                         coef_y = gufi_meanspec_fd_diff$coefs,
                                         y_names = gufi_meanspec_fd_diff$fdnames,
                                         basis_beta = gufi_meanspec_fd_diff$basis,
                                         lambda_grid = 10^seq(-4, 0, by = 0.5),
                                         nfold = 5)

# fit model + beta se
gufi_anova_model = ComputeBetaSd(factor = gufi$Climate_zone,
                                   X = gufi_meanspec_amps,
                                   dom = gufi_meanspec_freqs,
                                   basis_y = gufi_meanspec_fd_diff$basis,
                                   coef_y = gufi_meanspec_fd_diff$coefs,
                                   y_names = gufi_meanspec_fd_diff$fdnames,
                                   basis_beta = gufi_meanspec_fd_diff$basis,
                                   lambda = cv_fanova_res_gufi$lambda_min)

gc()

# ftest
perm_fanova_res_gufi = PermutFANOVA(factor = gufi$Climate_zone,
                                      X = gufi_meanspec_amps,
                                      dom = gufi_meanspec_freqs,
                                      basis_y = gufi_meanspec_fd_diff$basis,
                                      coef_y = gufi_meanspec_fd_diff$coefs,
                                      y_names = gufi_meanspec_fd_diff$fdnames,
                                      basis_beta = gufi_meanspec_fd_diff$basis,
                                      lambda = cv_fanova_res_gufi$lambda_min,
                                      n_perm = 100,
                                      seed = 123)

gc()


fANOVABetaSdPlot(my.betaestlist = gufi_anova_model$model$betaestlist,
                 my.betastderrlist = gufi_anova_model$beta_se$betastderrlist,
                 my.factor = gufi$Climate_zone,
                 my.name = "Gufi",
                 save_path = "results/prima_parte/images/f_beta_gufi.png",
                 my.width = MY.WIDTH,
                 my.height = MY.HEIGHT,
                 my.layout.matr = cbind(matrix(1, 2, 2),
                                        matrix(2:5, 2, 2)))




# ╭────────╮
# │Gabbiani│ ------------------------------------------------------------------
# ╰────────╯

cv_fanova_res_gabbiani = CvFunctionalANOVA(factor = gabbiani$Cluster,
                                       X = gabbiani_meanspec_amps,
                                       dom = gabbiani_meanspec_freqs,
                                       basis_y = gabbiani_meanspec_fd_diff$basis,
                                       coef_y = gabbiani_meanspec_fd_diff$coefs,
                                       y_names = gabbiani_meanspec_fd_diff$fdnames,
                                       basis_beta = gabbiani_meanspec_fd_diff$basis,
                                       lambda_grid = 10^seq(-4, 0, by = 0.5),
                                       nfold = 5)

gc()

# fit model + beta se
gabbiani_anova_model = ComputeBetaSd(factor = gabbiani$Cluster,
                                 X = gabbiani_meanspec_amps,
                                 dom = gabbiani_meanspec_freqs,
                                 basis_y = gabbiani_meanspec_fd_diff$basis,
                                 coef_y = gabbiani_meanspec_fd_diff$coefs,
                                 y_names = gabbiani_meanspec_fd_diff$fdnames,
                                 basis_beta = gabbiani_meanspec_fd_diff$basis,
                                 lambda = cv_fanova_res_gabbiani$lambda_min)

# ftest
perm_fanova_res_gabbiani = PermutFANOVA(factor = gabbiani$Cluster,
                                    X = gabbiani_meanspec_amps,
                                    dom = gabbiani_meanspec_freqs,
                                    basis_y = gabbiani_meanspec_fd_diff$basis,
                                    coef_y = gabbiani_meanspec_fd_diff$coefs,
                                    y_names = gabbiani_meanspec_fd_diff$fdnames,
                                    basis_beta = gabbiani_meanspec_fd_diff$basis,
                                    lambda = cv_fanova_res_gabbiani$lambda_min,
                                    n_perm = 100,
                                    seed = 123)

gc()


fANOVABetaSdPlot(my.betaestlist = gabbiani_anova_model$model$betaestlist,
                 my.betastderrlist = gabbiani_anova_model$beta_se$betastderrlist,
                 my.factor = gabbiani$Cluster,
                 my.name = "gabbiani",
                 save_path = "results/prima_parte/images/f_beta_gabbiani.png",
                 my.width = MY.WIDTH,
                 my.height = MY.HEIGHT,
                 my.layout.matr = cbind(matrix(1, 2, 2),
                                        matrix(2:5, 2, 2)))

# .. joint plot -----------------------------

#.. cv err ------
par(mfrow = c(3, 1))


png("results/prima_parte/images/f_anova_cv_err.png",
     width = MY.WIDTH, height = MY.HEIGHT)

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


