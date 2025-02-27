
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PenFFR

<!-- badges: start -->

The goal of PenFFR is to design a penalized funtion-on-funtion linear
regression with multiples functional and scalar covariates. This package
allows you to build two types of functional linear models: - The
concurrent linear model whose equation is given by :
$$\mathrm{y}_i(t) = \beta_0(t) + \sum_{\ell=1}^p\beta_\ell(t)\mathrm{x}_i^\ell(t) + \varepsilon_i(t)$$ -
The integral linear model given by :
$$\mathrm{y}_i(t) = \beta_0(t) + \sum_{\ell=1}^p\int_0^t\beta_\ell(s,t)\mathrm{x}_i^\ell(s)\,ds + \varepsilon_i(t)$$

Parameters are estimated using a spline-based expansion with a number of
basis functions that you can set.

To handle heterogeneous data, we also provide a Mixture-of-Experts (MoE)
version of the linear concurrent model. The conditional density of
$\mathrm{Y}(t)$ according to the function-on-function MoE (FFMoE) model
is $$\begin{eqnarray}
    f(\mathrm{Y}(t)|\mathrm{X}(t),\Psi(t)) &=& \sum_{k=1}^\mathrm{K}\pi_k(\mathrm{X}(t),\alpha_k(t)) \Phi(\mathrm{Y}(t);\mathrm{X}(t)\beta_k(t),\sigma^2_k),\label{ME2}
\end{eqnarray}$$ with $$\begin{itemize}
    \item $\pi_k(\mathrm{X}(t),\ \alpha_k(t))$ the mixture proportion of group $k$, also called the $k^{\text{th}}$ gated network function, depending on the covariate $\mathrm{X}(t)$ through group specific functional parameter $\alpha_k(t)$. More details will be provided in the next section;
    \item $\Psi_k(t) = (\beta_k(t),\alpha_k(t))$ are the functional parameters;
    \item $\Phi(\mathrm{Y}(t); \mathrm{X}_i(t)\beta_k(t), \sigma^2_k)$ is the Gaussian density probability function of mean $\mathrm{X}(t)\beta_k(t)$ and variance $\sigma^2_k$.
\end{itemize}$$

## Download

You can download the package

- [Version
  0.0.1](https://gitlab.tech.orange/rlsoftwarenet/penffr/-/archive/main/penffr-main.zip)

Or clone from Gitlab :

    $ git clone https://gitlab.tech.orange/rlsoftwarenet/penffr.git

## Cite

If you use this software, please cite the following work :

- Jean Steve Tamo Tchomgui, Julien Jacques, Guillaume Fraysse, Vincent
  Barriac, Stéphane Chrétien. A mixture of experts regression model for
  functional response with functional covariates. Statistics and
  Computing, 2024, 34 (154),
  [s11222-024-10455-z](https://link.springer.com/article/10.1007/s11222-024-10455-z).

- Jean Steve Tamo Tchomgui, Julien Jacques, Vincent Barriac, Guillaume
  Fraysse, Stéphane Chrétien. A Penalized Spline Estimator for
  Functional Linear Regression with Functional Response. 2024.
  [hal-04120709v3](https://hal.science/hal-04120709v3).

For any questions, Please contact Jean Steve Tamo Tchomgui at
<jeanstevetamo at yahoo dot fr> or Guillaume Fraysse at
<guillaume dot fraysse at orange dot com>.

## License

Copyright (c) 2021 — 2025 Orange

This code is released under the GPL2 license. See the `LICENSE.md` file for
more information.

## Contact

- Homepage: [opensource.orange.com](http://opensource.orange.com/)
- e-mail: <guillaume.fraysse@orange.com>

## Link documents

- [Contributors list](CONTRIBUTORS.md)

<!-- badges: end -->
