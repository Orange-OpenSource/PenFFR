# Software Name : PenFFR
# SPDX-FileCopyrightText: Copyright (c) Orange SA
# SPDX-License-Identifier:  GPL-2.0-or-later

# This software is distributed under the GNU General Public License v2.0 or later
# see the "LICENSE.md" file for more details or https://spdx.org/licenses/GPL-2.0-or-later.html

# Authors: see CONTRIBUTORS.md
# Software description: Linear regression for fully functional data.

#' Function that build the design matrix of the gated network model
#' @param X_fd.list the list of functional predictors
#' @param t.mat the matrix of time with the corresponded values of Y.mat
#' @param nbasis the number of basis for functional parameters
#' @param n.order the order of the splines basis
#' @import fda
#' @importFrom stats integrate
#' @return Dataframe of the design matrix of the gated network model
#' @export
gated.features <- function(t.mat, nbasis, n.order, X_fd.list){

  ## Parameters ----
  d <- length(X_fd.list)
  basis.beta <- create.bspline.basis(rangeval = range(t.mat, na.rm = T),
                                     nbasis = nbasis, norder = n.order)

  L.X <- sum(sapply(1:d, function(l){
    length(X_fd.list[[l]][["basis"]][["names"]])
  }))
  L.beta <- length(basis.beta[["names"]])


  ## Compute the integral ----
  vals <- sapply(1:L.X, function(j){
    t(sapply(1:(d*L.beta), function(i){

      f <- function(x){
        tmp.mat1 <- lapply(1:d, function(l){
          eval.basis(evalarg = x,
                     basisobj = X_fd.list[[l]][["basis"]])
        })
        b1 <- matrix(tmp.mat1[[1]], nrow = 1)
        if (d>1) {
          for (l in 2:d) {
            tmp <- matrix(tmp.mat1[[l]], nrow = 1)
            #cpt <- length(tmp)
            b1 <- c(b1, tmp)
            # rbind(cbind(b1, matrix(0, nrow = nrow(b1), ncol = cpt)),
            #             c(rep(0,ncol(b1)), tmp))
          }
        }

        tmp.mat2 <- eval.basis(evalarg = x,
                               basisobj = basis.beta)
        b2 <- matrix(tmp.mat2, nrow = 1)
        if (d>1) {
          for (l in 2:d) {
            tmp <- matrix(tmp.mat2, nrow = 1)
            #cpt <- length(tmp)
            b2 <- c(b2, tmp)
            # rbind(cbind(b2, matrix(0, nrow = nrow(b2), ncol = cpt)),
            #             c(rep(0,ncol(b2)), tmp))
          }
        }
        tmp.mat <- b1[j] * b2[i]

        return(tmp.mat)
      }

      stats::integrate(Vectorize(f), lower = range(t.mat, na.rm = T)[1],
                       upper = range(t.mat, na.rm = T)[2],
                       rel.tol = 1e-15)[["value"]]

    }))
  })


  ## Functional coefficients of covariates ----
  coefs <- c()
  for (l in 1:d) {
    coefs <- rbind(coefs, X_fd.list[[l]][["coefs"]])
  }


  ## Gated features matrix ----
  r <- t(vals %*% coefs)

  # Build the gated data ----
  data.gated <- data.frame()
  for (i in 1:nrow(t.mat)) {
    ki <- length(na.omit(as.vector(unlist(t.mat[i,]))))
    tmp.B <- matrix(rep(r[i,],ki), nrow = ki, byrow = T)
    data.gated <- rbind(data.gated, data.frame(tmp.B))
  }

  return(data.gated)

}
