# Software Name : PenFFR
# SPDX-FileCopyrightText: Copyright (c) Orange SA
# SPDX-License-Identifier:  GPL-2.0-or-later

# This software is distributed under the GNU General Public License v2.0 or later
# see the "LICENSE.md" file for more details or https://spdx.org/licenses/GPL-2.0-or-later.html

# Authors: see CONTRIBUTORS.md
# Software description: Linear regression for fully functional data.

Pensim1_mix <- function(data, data.gated, X.scal, lams, vals, t.mat, K.min, K.max, n, m, d, d.z, K.comp, seed){

  '
  data is the data frame of functional data expand as matrix
  data.gated is the data frame for the gated features
  lams is a vector of possible value of penalty parameters
  vals is the list of penalty matrix for every parameters
  t.mat is a n*m matrix of observation
  d is the number of functional predictors
  d.z is the number of scalar predictors
  '
  if (length(lams) > 5) {
    set.seed(seed)
    ind <- sample(2:length(lams), 5)
    lams <- c(lams[1], lams[ind])
  }

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
  vals2 <- cvals[[1]]
  for (l in 1:d) {
    vals2 <- rbind(cbind(vals2, matrix(0, nrow = nrow(vals2), ncol = ncol(cvals[[l]]))),
                   cbind(matrix(0, nrow = nrow(vals[[l]]), ncol = ncol(vals2)), cvals[[l]]))
  }

  source("R/utils.R")
  cv.bic <- lapply(1:length(lams), function(v){

    ### v=1 case
    if (v == 1) {
      ## Model estimation ----
      # tmp.formula <- paste(colnames(data)[-(1:3)], collapse = " + ")
      # my_formula <- as.formula(paste(paste("output ~ 0 ", tmp.formula, sep = " + "),
      #                                "id", sep = "|"))
      # my_formula2 <- as.formula(paste("~", paste(colnames(data.gated),
      #                                            collapse = " + "),
      #                                 sep = ""))
      # tmp.model <- FLXMRlmm(random = ~ 1, lm.fit = "lm.wfit")
      # control <- list(verbose = -0, iter.max = 500, minprior = 0.001)

      tmp.formula <- paste(colnames(data)[-(1:3)], collapse = " + ")
      my_formula <- as.formula(paste(paste("output ~ 0 ", tmp.formula, sep = " + "),
                                     "id", sep = "|"))
      my_formula <- as.formula(paste("output ~ 0 ", tmp.formula, sep = " + "))
      my_formula2 <- as.formula(paste("~", paste(colnames(data.gated),
                                                 collapse = " + "),
                                      sep = ""))
      #tmp.model <- FLXMRglm(my_formula)
      tmp.model <- FLXMRlmm(random = ~ 1, lm.fit = "lm.wfit")
      control <- list(verbose = -0, iter.max = 2000, minprior = 0.1)

      if (K.comp == "BIC") {
        set.seed(seed)
        model <- initFlexmix(my_formula,
                             model = tmp.model, #FLXMRlmm(random = ~ 1),
                             concomitant = FLXPmultinom(my_formula2),
                             data = cbind(data, data.gated),
                             nrep = 5, k = K.min:K.max,
                             control = control, verbose = F)
        #print(model)
        m.BIC <- getModel(model, "BIC")

      } else {
        set.seed(seed)
        m.BIC <- flexmix(my_formula,
                         model = tmp.model,
                         k = K.comp,
                         concomitant = FLXPmultinom(my_formula2),
                         data = cbind(data, data.gated),
                         control = control)
      }
    } else {
      ### Build the covariates and response for penalization ----
      datapen <- (lams[v]**2)*vals2
      #datapen <- matrix(rep(as.vector(t(datapen)), each = m), ncol = ncol(datapen), byrow = F)
      datapen.gated <- datapen[,-c(1:nbasis)]
      datapen <- cbind.data.frame(((1:nrow(datapen))-1)/(nrow(datapen)-1), datapen)
      if (d.z > 0) {
        for (i in 1:d.z) {
          datapen <- cbind.data.frame(datapen, mean(X.scal[,i]))
          datapen.gated <- cbind.data.frame(datapen.gated, 0)
        }
      }
      tmp <- (n+1):(n+nrow(datapen))
      datapen <- cbind.data.frame(id = tmp, datapen)
      datapen <- cbind.data.frame(0, datapen)
      colnames(datapen) <- colnames(data)
      colnames(datapen.gated) <- colnames(data.gated)
      datapen$id <- as.factor(datapen$id)
      newdat <- rbind.data.frame(data, datapen)
      newdat.gated <- rbind.data.frame(data.gated, datapen.gated)

      ## Model estimation ----
      # tmp.formula <- paste(colnames(newdat)[-(1:3)], collapse = " + ")
      # my_formula <- as.formula(paste(paste("output ~ 0", tmp.formula, sep = " + "),
      #                                "id", sep = "|"))
      # my_formula2 <- as.formula(paste("~", paste(colnames(newdat.gated),
      #                                            collapse = " + "),
      #                                 sep = ""))
      # tmp.model <- FLXMRlmm(random = ~ 1)
      # control <- list(verbose = -0, iter.max = 500, minprior = 0.001)

      tmp.formula <- paste(colnames(newdat)[-(1:3)], collapse = " + ")
      my_formula <- as.formula(paste(paste("output ~ 0 ", tmp.formula, sep = " + "),
                                     "id", sep = "|"))
      my_formula <- as.formula(paste("output ~ 0 ", tmp.formula, sep = " + "))
      my_formula2 <- as.formula(paste("~", paste(colnames(newdat.gated),
                                                 collapse = " + "),
                                      sep = ""))
      #tmp.model <- FLXMRglm(my_formula)
      tmp.model <- FLXMRlmm(random = ~ 1, lm.fit = "lm.wfit")
      control <- list(verbose = -0, iter.max = 2000, minprior = 0.1)
      if (K.comp == "BIC") {
        set.seed(seed)
        model2 <- initFlexmix(my_formula,
                              model = tmp.model,
                              concomitant = FLXPmultinom(my_formula2),
                              data = cbind(newdat, newdat.gated),
                              nrep = 5, k = K.min:K.max,
                              control = control, verbose = F,
                              init = list(name = "tol.em"))
        m.BIC <- getModel(model2, "BIC")

      } else {
        set.seed(seed)
        m.BIC <- flexmix(my_formula,
                         model = tmp.model,
                         k = K.comp,
                         concomitant = FLXPmultinom(my_formula2),
                         data = cbind(newdat, newdat.gated),
                         control = control)
      }
    }


    ### Compute the BIC
    bic <- BIC(m.BIC)

    list(BIC = bic, model = m.BIC)

  })

  mspe.vect <- sapply(1:length(cv.bic), function(v){
    cv.bic[[v]]$BIC
  })

  v <- which.min(mspe.vect)


  return(list(lambda = lams[v], model = cv.bic[[v]]$model))
}
