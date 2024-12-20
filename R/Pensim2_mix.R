# Software Name : PenFFR
# SPDX-FileCopyrightText: Copyright (c) Orange SA
# SPDX-License-Identifier:  GPL-2.0-or-later

# This software is distributed under the GNU General Public License v2.0 or later
# see the "LICENSE.md" file for more details or https://spdx.org/licenses/GPL-2.0-or-later.html

# Authors: see CONTRIBUTORS.md
# Software description: Linear regression for fully functional data.
Pensim2_mix <- function(data, data.gated, ind, lams, vals, t.mat, d, d.z, K.comp = "BIC"){

  '
  data is the data frame of functional data expand as matrix
  data.gated is the data frame for the gated features
  ind is the index of non pertinent features
  lams is a vector of possible value of penalty parameters
  vals is the list of penalty matrix for every parameters
  t.mat is a n*m matrix of observation
  d is the number of functional predictors
  d.z is the number of scalar predictors
  '

  ## grid search of the optimum lambdas ----
  grid_search <- as.matrix(expand.grid(rep(list(lams), d+1)))
  cvals <- lapply(1:length(vals), function(l){
    tryCatch({
      chol(vals[[l]])
    },
    error = function(e){
      vals[[l]]
      return(vals[[l]])
    })
  })
  nbasis <- nrow(vals[[1]])

  source("R/utils.R")
  mspe.vect <- sapply(1:nrow(grid_search), function(v){

    ### Build the covariates and response for penalization ----
    datapen <- grid_search[v,1] * cvals[[1]]
    for (l in 2:(d+1)) {
      tmp <- grid_search[v,l] * cvals[[l]]
      datapen <- rbind(cbind(datapen, matrix(0, nrow = nrow(datapen), ncol = ncol(tmp))),
                       cbind(matrix(0, nrow = nrow(tmp), ncol = ncol(datapen)), tmp))
    }
    datapen.gated <- datapen
    datapen.gated <- datapen.gated[,c(1:nbasis)]
    if (d>1) {
      for (l in 2:d) {
        datapen.gated <- rbind(datapen.gated, datapen[,c(1:nbasis)])
      }
    }

    for (i in 1:2) {
      datapen <- cbind.data.frame(0, datapen)
    }
    if (d.z > 0) {
      for (i in 1:d.z) {
        datapen <- cbind.data.frame(datapen, 0)
        datapen.gated <- cbind.data.frame(datapen.gated, 0)
      }
    }
    tmp <- (n+1):(n+nrow(datapen))
    datapen <- cbind.data.frame(id = tmp, datapen)
    datapen$id <- as.factor(datapen$id)
    colnames(datapen) <- colnames(data)
    colnames(datapen.gated) <- colnames(data.gated)
    newdat <- rbind.data.frame(data, datapen)
    newdat.gated <- rbind.data.frame(data.gated, datapen.gated)

    ### model estimation ----
    K.max = 5
    n.repeat = 20
    my_formula <- as.formula(paste(" ~ 0", paste(colnames(newdat)[-c(1:3, ind)],
                                                 collapse = " + "),
                                   sep = " + "))
    tmp.model <- FLXMRglmfix(my_formula)
    my_formula2 <- as.formula(paste(" ~ ", paste(colnames(newdat.gated),
                                                 collapse = " + "),
                                    sep = ""))
    model2 <- stepFlexmix(output ~ 0 | id,
                          model = tmp.model,
                          concomitant = FLXPmultinom(my_formula2),
                          data = cbind(newdat, newdat.gated),
                          nrep = n.repeat, k = 1:K.max)
    if (K.comp == "BIC") {
      m.BIC <- getModel(model2, "BIC")
    } else {
      m.BIC <- model2@models[[K.comp]]
    }

    ### Compute the BIC
    BIC(m.BIC)

  })

  ind <- which.min(mspe.vect)
  lams <- as.vector(grid_search[ind,])

  for (l in 2:(d+1)) {
    tmp <- grid_search[v,l] * cvals[[l]]
    datapen <- rbind(cbind(datapen, matrix(0, nrow = nrow(datapen), ncol = ncol(tmp))),
                     cbind(matrix(0, nrow = nrow(tmp), ncol = ncol(datapen)), tmp))
  }
  datapen.gated <- datapen
  datapen.gated <- datapen.gated[,-c(1:nbasis)]
  for (i in 1:2) {
    datapen <- cbind.data.frame(0, datapen)
  }
  if (d.z > 0) {
    for (i in 1:d.z) {
      datapen <- cbind.data.frame(datapen, 0)
      datapen.gated <- cbind.data.frame(datapen.gated, 0)
    }
  }
  tmp <- (n+1):(n+nrow(datapen))
  datapen <- cbind.data.frame(id = tmp, datapen)
  datapen$id <- as.factor(datapen$id)
  colnames(datapen) <- colnames(data)
  colnames(datapen.gated) <- colnames(data.gated)
  newdat <- rbind.data.frame(data, datapen)
  newdat.gated <- rbind.data.frame(data.gated, datapen.gated)

  ## Model estimation ----
  K.max = 5
  n.repeat = 20
  my_formula <- as.formula(paste(" ~ 0", paste(colnames(newdat)[-(1:3)],
                                               collapse = " + "),
                                 sep = " + "))
  tmp.model <- FLXMRglmfix(my_formula)
  my_formula2 <- as.formula(paste(" ~ ", paste(colnames(newdat.gated),
                                               collapse = " + "),
                                  sep = ""))
  if (K.comp == "BIC") {
    model2 <- stepFlexmix(output ~ 0 | id,
                          model = tmp.model,
                          concomitant = FLXPmultinom(my_formula2),
                          data = cbind(data, data.gated),
                          nrep = n.repeat, k = 1:K.max)
    m.BIC <- getModel(model2, "BIC")
  } else {
    m.BIC <- flexmix(output ~ 0 | id,
                     model = tmp.model,
                     concomitant = FLXPmultinom(my_formula2),
                     data = cbind(data, data.gated),
                     k = K.comp)
  }

  return(list(lambda = lams, model = m.BIC))
}
