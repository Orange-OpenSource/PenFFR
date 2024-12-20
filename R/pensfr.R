# Software Name : PenFFR
# SPDX-FileCopyrightText: Copyright (c) Orange SA
# SPDX-License-Identifier:  GPL-2.0-or-later

# This software is distributed under the GNU General Public License v2.0 or later
# see the "LICENSE.md" file for more details or https://spdx.org/licenses/GPL-2.0-or-later.html

# Authors: see CONTRIBUTORS.md
# Software description: Linear regression for fully functional data.
my_funmat4 <- function(X_fd.list, basis.beta, t.mat, Y.mat){

  source("R/func-to-mat1.R")

  d <- length(X_fd.list)
  n <- ncol(X_fd.list[[1]][["coefs"]])
  L.beta <- length(basis.beta[["names"]])

  ## Compute the basis function integral
  tmp.b <- lapply(1:d, function(l){
    basis.integral(t0 = min(t.mat),
                   tf = max(t.mat),
                   basis.beta = basis.beta,
                   basis.X = X_fd.list[[l]][["basis"]])
  })

  ## R.mat at a point t for a covariate l and subject i ----
  R.mat <- function(t, l, i){

    coefs.i <- as.vector(X_fd.list[[l]][["coefs"]][,i])

    b <- eval.basis(evalarg = t, basisobj = basis.beta)
    B1 <- t(diag(rep(b, length(b))))

    B <- eval.basis(evalarg = t, basisobj = X_fd.list[[l]][["basis"]])

    obs.grid <- as.vector(na.omit(as.vector(unlist(t.mat[i,]))))

    j <- which(obs.grid == t)
    B2 <- matrix(0, nrow = L.beta**2, ncol = length(coefs.i))
    for (k in 1:j) {
      tmp.B2 <- readRDS(file = paste(paste("BasisIntegral", L.beta, sep = ""),
                                     "/mat",l, k,".rds", sep = ""))
      B2 <- B2 + as.matrix(tmp.B2)
    }

    ### Vecteur R ----
    return(coefs.i %*% t(B2) %*% B1)
  }

  ## Build the design matrix ----
  R <- lapply(1:n, function(i){
    tmp.time <- data$time[data$id == i]

    t(sapply(tmp.time, function(tj) {
      B0 <- eval.basis(evalarg = tj, basisobj = basis.beta)
      c(B0, as.vector(sapply(1:d, function(l){
        R.mat(tj, l, i)
      })))
    }))

  })

  new.R <- c()
  for (i in 1:n) {
    new.R <- rbind(new.R, R[[i]])
  }
  rm(R)

  data <- cbind.data.frame(data, data.frame(new.R))
  colnames(data) <- c("output","id", "time",
                      paste("X", 1:(ncol(data)-3), sep = "."))
  data$id <- as.factor(data$id)

  return(data)
}

get.data3 <- function(X_fd.list, Y.mat, t.mat, nbasis, n.order = 4){
  source("R/func-to-mat1.R")
  basis.beta <- create.bspline.basis(rangeval = range(t.mat, na.rm = T),
                                     nbasis = nbasis, norder = n.order)
  data <- my_funmat4(X_fd.list = X_fd.list,
                     basis.beta = basis.beta,
                     Y.mat = Y.mat,
                     t.mat = t.mat)
  return(data)
}

pensfr <- function(Y.mat, t.mat = "all.ok", X_fd.list, X.scal = data.frame(), nbasis, n.order=4, pen = T){

  '
  Y.mat is matrix response with of n rows and m = max m_i columns with na values after the end of the observation for each sample
  t.mat is a matrix of time with the corresponded values of Y.mat
  X_fd.list is the list of functional predictors
  X.scal is the data frame of the non functional predictors
  nbasis is an integer for the number of basis for functional parameter
  Pen is the boolean value for the penalization
  '

  # Check if the t.mat is not provided ----
  n <- nrow(Y.mat)
  m <- ncol(Y.mat)
  d <- length(X_fd.list)
  if (all(t.mat == "all.ok")) {
    obs.grid <- ((1:m)-1)/(m-1)
    t.mat <- matrix(rep(obs.grid, n), nrow = n, byrow = T)
  }

  # Check if nbasis is provided ----
  if (all(nbasis == "all.ok")){
    nbasis = min(20, floor(((n*m)/(d+1))**0.5))
  }

  # Build the dataframe and penalty matrix for the model ----
  source("R/func-to-mat2.R")
  data <- get.data2(X_fd.list = X_fd.list,
                    Y.mat = Y.mat, nbasis = nbasis,
                    t.mat = t.mat, n.order = n.order)

  if (ncol(X.scal) > 0){
    # Formatting the non-functional data ----
    X.scal.new <- data.frame()
    for (i in 1:n) {
      ni <- length(data$time[data$id == i])
      tmp <- matrix(rep(as.vector(unlist(X.scal[i,])), each = ni),
                    nrow = ni, byrow = F)
      X.scal.new <- rbind(X.scal.new, tmp)
    }
    X.scal.new <- data.frame(X.scal.new)
    colnames(X.scal.new) <- paste("Z", 1:ncol(X.scal.new), sep = ".")
    data <- cbind.data.frame(data, X.scal.new)
  }

  if (pen == F){

    ## Run the model ----
    my_formula <- as.formula(paste("output ~ 0", paste(colnames(data)[-(1:3)],
                                                       collapse = " + "), sep = " + "))
    model <- lm(my_formula, data = cbind.data.frame(data))

    ## get the functional parameters ----
    params <- as.vector(model$coefficients)
    params[is.na(params)] <- 0
    beta.fd <- lapply(1:(d+1), function(l){
      if (l == 1){
        fd(coef = params[1:nbasis],
           basisobj = create.bspline.basis(rangeval = range(obs.grid),
                                           nbasis = nbasis, norder = n.order),
           fdnames = list(main = TeX("$\\beta_0(t)"),
                          xlab = "Time",
                          ylab = TeX("$\\beta_0(t)")))
      } else{
        bifd(coef = matrix(params[nbasis+(((nbasis**2)*(l-2)+1):((nbasis**2)*(l-1)))],
                           nbasis, nbasis, byrow = F),
             sbasisobj = create.bspline.basis(rangeval = range(obs.grid),
                                              nbasis = nbasis, norder = n.order),
             tbasisobj = create.bspline.basis(rangeval = range(obs.grid),
                                              nbasis = nbasis, norder = n.order))
      }
    })

    if (ncol(X.scal) > 0) {
      ## get the non functional parameters ----
      beta.scal <- params[-c(1:(nbasis + (nbasis**2)*d))]
      ## results ----
      return(list(model = model, beta.fd = beta.fd, beta.scal = beta.scal, n.grid = m))
    } else {
      ## results ----
      return(list(model = model, beta.fd = beta.fd, n.grid = m))
    }

  } else {

    source("R/func-to-mat2.R")
    ## Build the penalty matrix ----
    vals <- my_penmat2(obs.grid = t.mat, d = length(X_fd.list)+1,
                       L.beta = nbasis)

    # Run the functional model ----
    source("R/Pensim2.R")
    d <- length(X_fd.list)
    n.lam <- min(10, floor(exp(log(100)/(d+1))))
    fofreg <- Pensim2(data = data, lams = seq(0, 5, length.out = n.lam),
                      vals = vals, d = d, d.z = ncol(X.scal),
                      t.mat = t.mat)

    model <- fofreg$model

    ## get the functional parameters ----
    params <- model$coefficients
    params[is.na(params)] <- 0
    beta.fd <- lapply(1:(d+1), function(l){
      if (l == 1){
        fd(coef = params[1:nbasis],
           basisobj = create.bspline.basis(rangeval = range(obs.grid),
                                           nbasis = nbasis, norder = 4),
           fdnames = list(main = TeX("$\\beta_0(t)"),
                          xlab = "Time",
                          ylab = TeX("$\\beta_0(t)")))
      } else{
        bifd(coef = matrix(params[nbasis+(((nbasis**2)*(l-2)+1):((nbasis**2)*(l-1)))],
                           nbasis, nbasis, byrow = F),
             sbasisobj = create.bspline.basis(rangeval = range(obs.grid),
                                              nbasis = nbasis, norder = 4),
             tbasisobj = create.bspline.basis(rangeval = range(obs.grid),
                                              nbasis = nbasis, norder = 4))
      }
    })

    if (ncol(X.scal) > 0) {
      ## get the non functional parameters ----
      beta.scal <- params[-c(1:(nbasis + (nbasis**2)*d))]
      ## results ----
      return(list(model = model, beta.fd = beta.fd, n.grid = m,
                  beta.scal = beta.scal, lambda = fofreg$lambda))
    } else {
      ## results ----
      return(list(model = model, beta.fd = beta.fd,
                  lambda = fofreg$lambda, n.grid = m))
    }
  }
}


pred.penffr2 <- function(model, newX_fd.list = list(), newX.scal = data.frame(), t.mat = "all.ok"){

  '
  t.mat must be the time matrix of predictions with of n rows and m columns
  newX_fd.list is the list of functional predictors
  newX.scal is the data frame of the non functional predictors
  model is the penffr model
  '

  # Check the providing of new functional predictors ----
  if (length(newX_fd.list) == 0) {
    if (nrow(newX.scal) == 0) {
      pred <- predict(model$model)
    } else {
      return(NA)
    }

  } else {

    # observation grid of predictions ----
    n <- ncol(newX_fd.list[[1]][["coefs"]])
    obs.grid <- seq(newX_fd.list[[1]][["basis"]][["rangeval"]][1],
                    newX_fd.list[[1]][["basis"]][["rangeval"]][2],
                    length.out = model$n.grid)
    t.mat <- matrix(rep(obs.grid, n), nrow = n, byrow = T)

    # Build the dataframe and penalty matrix for the model ----
    source("R/func-to-mat1.R")
    nbasis <- length(as.vector(model$beta.fd[[1]][["coefs"]]))
    data <- get.data2(X_fd.list = newX_fd.list,
                      Y.mat = t.mat, nbasis = nbasis,
                      t.mat = t.mat)

    if (ncol(newX.scal) > 0){
      # Formatting the non-functional data ----
      X.scal.new <- data.frame()
      for (i in 1:n) {
        ni <- length(data$time[data$id == i])
        tmp <- matrix(rep(as.vector(unlist(newX.scal[i,])), each = ni),
                      nrow = ni, byrow = F)
        X.scal.new <- rbind(X.scal.new, tmp)
      }
      X.scal.new <- data.frame(X.scal.new)
      colnames(X.scal.new) <- paste("Z", 1:ncol(X.scal.new), sep = ".")
      data <- cbind.data.frame(data, X.scal.new)
    }

    # prediction ----
    pred <- predict(model$model, newdata = data, type = "response")
  }

  # result ----
  return(pred)

}
