# Software Name : PenFFR
# SPDX-FileCopyrightText: Copyright (c) Orange SA
# SPDX-License-Identifier:  GPL-2.0-or-later

# This software is distributed under the GNU General Public License v2.0 or later
# see the "LICENSE.md" file for more details or https://spdx.org/licenses/GPL-2.0-or-later.html

# Authors: see CONTRIBUTORS.md
# Software description: Linear regression for fully functional data.
penffr.pert <- function(Y.mat, t.mat = "all.ok", X_fd.list, X.scal = data.frame(), N = 100,
                          n.order = 4, nbasis = "all.ok", pen = F, mod.type = "concurrent"){

  '
  Y.mat       : matrix response with of n rows and m = max m_i columns with na values after the end of the observation for each sample
  t.mat       : a matrix of time with the corresponded values of Y.mat
  X_fd.list   : list of functional predictors
  X.scal      : data frame of the non functional predictors
  nbasis      : integer for the number of basis for functional parameter
  Pen         : boolean value for the penalization
  mod.type    : Type of model "concurrent" or "integral"
  tau         : level of desired quantile
  '

  # Test of the type of model
  if (mod.type == "concurrent") {

    ## compute the design matrix for the model
    n <- nrow(Y.mat)
    m <- ncol(Y.mat)
    d <- length(X_fd.list)
    if (all(t.mat == "all.ok")) {
      obs.grid <- ((1:m)-1)/(m-1)
      t.mat <- matrix(rep(obs.grid, n), nrow = n, byrow = T)
    }
    if (all(nbasis == "all.ok")){
      nbasis = min(20, floor((n*m)/(d+1)))
    }
    source("R/func-to-mat1.R")
    data <- get.data1(X_fd.list = X_fd.list,
                      Y.mat = Y.mat, nbasis = nbasis,
                      t.mat = t.mat, n.order = n.order)

    if (ncol(X.scal) > 0){
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

      ## Run the weighted models
      mods <- lapply(1:N, function(w){

        ### generate the weights
        set.seed(w)
        weights <- rexp(n = length(na.omit(as.vector(Y.mat))), rate = 0.05)

        ### models
        my_formula <- as.formula(paste("output ~ 0", paste(colnames(data)[-(1:3)],
                                                           collapse = " + "), sep = " + "))
        model <- lm(my_formula, data = cbind.data.frame(data), weights = weights)
        params <- model$coefficients
        params[is.na(params)] <- 0
        beta.fd <- lapply(1:(d+1), function(l){
          if (l != 1){
            fd(coef = params[(nbasis*(l-1)+1):(nbasis*l)],
               basisobj = create.bspline.basis(rangeval = range(t.mat, na.rm = T),
                                               nbasis = nbasis, norder = 4),
               fdnames = list(main = paste("beta_",(l-1), sep = ""),
                              xlab = "Time (millisec.)",
                              ylab = paste("beta_",(l-1),sep = "")))
          } else{
            fd(coef = params[(nbasis*(l-1)+1):(nbasis*l)],
               basisobj = create.bspline.basis(rangeval = range(t.mat, na.rm = T),
                                               nbasis = nbasis, norder = 4),
               fdnames = list(main = paste("beta_",(l-1),sep = ""),
                              xlab = "Time (millisec.)",
                              ylab = paste("beta_",(l-1),sep = "")))
          }
        })

        if (ncol(X.scal) > 0) {
          beta.scal <- model$coefficients[-c(1:(nbasis*(d+1)))]
          tmp.mod <- list(model = model, beta.fd = beta.fd, beta.scal = beta.scal,
                       n.grid = m)
        } else {
          tmp.mod <- list(model = model, beta.fd = beta.fd, n.grid = m)
        }
        tmp.mod
      })

    } else {

      ## Run the weighted models
      vals <- my_penmat1(obs.grid = t.mat, d = length(X_fd.list)+1,
                         nbasis = nbasis, n.order = n.order, deg = 2)
      n.lam <- min(10, floor(exp(log(100)/(d+1))))

      mods <- lapply(1:N, function(w){

        ### generate the weights
        set.seed(w)
        weights <- rexp(n = length(na.omit(as.vector(Y.mat))) + nbasis*(d+1),
                        rate = 0.05)

        ### models
        source("R/Pensim1.R")
        fofreg <- Pensim1(data = data, lams = seq(0.01, 5.0,length.out = n.lam),
                          vals = vals, d = d, d.z = ncol(X.scal),
                          t.mat = t.mat, weights = weights)
        model <- fofreg$model
        lams <- fofreg$lambda
        params <- model$coefficients
        params[is.na(params)] <- 0
        beta.fd <- lapply(1:(d+1), function(l){
          if (l != 1){
            fd(coef = params[(nbasis*(l-1)+1):(nbasis*l)],
               basisobj = create.bspline.basis(rangeval = range(t.mat, na.rm = T),
                                               nbasis = nbasis, norder = n.order),
               fdnames = list(main = paste("beta_",(l-1), sep = ""),
                              xlab = "Time (millisec.)",
                              ylab = paste("beta_",(l-1),sep = "")))
          } else{
            fd(coef = params[(nbasis*(l-1)+1):(nbasis*l)],
               basisobj = create.bspline.basis(rangeval = range(t.mat, na.rm = T),
                                               nbasis = nbasis, norder = n.order),
               fdnames = list(main = paste("beta_",(l-1),sep = ""),
                              xlab = "Time (millisec.)",
                              ylab = paste("beta_",(l-1),sep = "")))
          }
        })
        if (ncol(X.scal) > 0) {
          beta.scal <- params[-c(1:(nbasis*(d+1)))]
          tmp.mod <- list(model = model, beta.fd = beta.fd,
                       beta.scal = beta.scal, lambda = lams, n.grid = m)
        } else {
          tmp.mod <- list(model = model, beta.fd = beta.fd, lambda = lams,
                       n.grid = m)
        }
        tmp.mod
      })
    }

    return(mods)

  } else if (mod.type == "integral"){

    ## compute the design matrix for the model
    n <- nrow(Y.mat)
    m <- ncol(Y.mat)
    d <- length(X_fd.list)
    if (all(t.mat == "all.ok")) {
      obs.grid <- ((1:m)-1)/(m-1)
      t.mat <- matrix(rep(obs.grid, n), nrow = n, byrow = T)
    }
    if (all(nbasis == "all.ok")){
      nbasis = min(20, floor(((n*m)/(d+1))**0.5))
    }
    source("R/func-to-mat2.R")
    data <- get.data2(X_fd.list = X_fd.list,
                      Y.mat = Y.mat, nbasis = nbasis,
                      t.mat = t.mat, n.order = n.order)
    if (ncol(X.scal) > 0){
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

      mods <- lapply(1:N, function(w){

        ### generate the weights
        set.seed(w)
        weights <- rexp(n = length(na.omit(as.vector(Y.mat))), rate = 0.05)

        ### model
        my_formula <- as.formula(paste("output ~ 0", paste(colnames(data)[-(1:3)],
                                                           collapse = " + "), sep = " + "))
        model <- lm(my_formula, data = cbind.data.frame(data), weights = weights)

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
          tmp.mod <- list(model = model, beta.fd = beta.fd,
                          beta.scal = beta.scal, n.grid = m)
        } else {
          ## results ----
          tmp.mod <- list(model = model, beta.fd = beta.fd, n.grid = m)
        }
        tmp.mod
      })

    } else {

      source("R/func-to-mat2.R")
      vals <- my_penmat2(obs.grid = t.mat, d = length(X_fd.list)+1,
                         L.beta = nbasis)
      source("R/Pensim2.R")
      d <- length(X_fd.list)
      n.lam <- min(10, floor(exp(log(100)/(d+1))))

      mods <- lapply(1:N, function(w){

        ### generate the weights
        set.seed(w)
        weights <- rexp(n = length(na.omit(as.vector(Y.mat))) + nbasis*(1+d*nbasis),
                        rate = 0.05)

        ### model
        fofreg <- Pensim2(data = data, lams = seq(0, 5, length.out = n.lam),
                          vals = vals, d = d, d.z = ncol(X.scal),
                          t.mat = t.mat, weights = weights)
        model <- fofreg$model
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
          tmp.mod <- list(model = model, beta.fd = beta.fd, n.grid = m,
                          beta.scal = beta.scal, lambda = fofreg$lambda)
        } else {
          ## results ----
          tmp.mod <- list(model = model, beta.fd = beta.fd,
                          lambda = fofreg$lambda, n.grid = m)
        }
        tmp.mod
      })
    }

    return(mods)

  } else {
    msg = "No good type of model. It should be 'concurrent' or 'integral'"
    return(msg)
  }
}

pred.penffr.pert <- function(models, mod.type = "concurrent", newX_fd.list = list(),
                             newX.scal = data.frame(), t.mat = "all.ok", tau = 0.5){
  N <- length(models)
  library(jubilee)

  if (mod.type == "concurrent") {

    if (length(newX_fd.list) == 0) {
      if (nrow(newX.scal) == 0) {
        preds <- jubilee.mcsapply(1:N, function(w){
          predict(models[[w]]$model)
        })
      } else {
        return(NA)
      }

    } else {

      # Compute the design matrix
      n <- ncol(newX_fd.list[[1]][["coefs"]])
      if (all(t.mat == "all.ok")) {
        obs.grid <- seq(newX_fd.list[[1]][["basis"]][["rangeval"]][1],
                        newX_fd.list[[1]][["basis"]][["rangeval"]][2],
                        length.out = models[[1]]$n.grid)
        t.mat <- matrix(rep(obs.grid, n), nrow = n, byrow = T)
      }

      d <- length(newX_fd.list)

      # Build the dataframe and penalty matrix for the model ----
      source("R/func-to-mat1.R")
      nbasis <- length(as.vector(models[[1]]$beta.fd[[1]][["coefs"]]))
      n.order <- nbasis - length(models[[1]]$beta.fd[[1]][["basis"]][["params"]])
      data <- get.data1(X_fd.list = newX_fd.list,
                        Y.mat = t.mat, nbasis = nbasis,
                        t.mat = t.mat, n.order = n.order)

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
      preds <- jubilee.mcsapply(1:N, function(w){
        predict(models[[w]]$model, newdata = data, type = "response")
      })
    }
    return(preds)

  } else if (mod.type == "integral"){

    if (length(newX_fd.list) == 0) {
      if (nrow(newX.scal) == 0) {
        preds <- jubilee.mcsapply(1:N, function(w){
          predict(models[[w]]$model)
        })
      } else {
        return(NA)
      }

    } else {
      # Compute the design matrix
      n <- ncol(newX_fd.list[[1]][["coefs"]])
      obs.grid <- seq(newX_fd.list[[1]][["basis"]][["rangeval"]][1],
                      newX_fd.list[[1]][["basis"]][["rangeval"]][2],
                      length.out = models[[1]]$n.grid)
      t.mat <- matrix(rep(obs.grid, n), nrow = n, byrow = T)

      # Build the dataframe and penalty matrix for the model ----
      source("R/func-to-mat1.R")
      nbasis <- length(as.vector(models[[1]]$beta.fd[[1]][["coefs"]]))
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

      # prediction
      preds <- jubilee.mcsapply(1:N, function(w){
        predict(models[[w]]$model, newdata = data, type = "response")
      })
    }
    return(preds)

  } else {
    msg = "No good type of model. It should be -concurrent- or -integral-"
    return(msg)

  }
}



