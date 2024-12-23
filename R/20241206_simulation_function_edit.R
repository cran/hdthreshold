#' @import stats
NULL

#Consider adding importFrom("stats", "approx", "fft", "median", "pnorm", "qnorm", "rexp", "rnorm", "runif", "setNames", "var") to your NAMESPACE file.

#' Uniform kernel function
#'
#' @param x a vector
#'
#' @return a vector of values
#' @export
#'
#' @examples
#' K(0)
K = function(x){ifelse((abs(x)) < 1, 1 / 2, 0)}

#' Simulate an MA infinity process with algrebraic decay
#'
#' @param N sample size
#' @param beta algebraic decay parameter
#'
#' @return simulated MA infinity process
#' @export
#'
#' @examples
#' x = MAinf_normal(100, 1.5)
MAinf_normal = function(N, beta) {
  NN = 2 * N
  a = (1:NN) ^ (-beta)
  u = rnorm(NN)
  X = Re(fft(fft(a) * fft(u), inverse = TRUE)) / NN
  return(X[1:N])
}

#' Simulate an example data frame
#'
#' @param N cross-sectional dimension
#' @param TL time series length
#' @param p fraction of non-zero coefficients
#' @param gamma value of non-zero coefficients
#'
#' @return simulated data frame
#' @export
#'
#' @examples
#'d = threshold.example(10, 200, 0.1, 2)
threshold.example = function(N, TL, p, gamma) {
  x = t(replicate(N, MAinf_normal(TL, 1.5)))
  e = t(replicate(N, MAinf_normal(TL, 1.5)))
  tau = rep(0, N)
  tau[0:floor(N * p)] = gamma
  y = x + ifelse(x > 0, tau, 0) + e
  id = t(replicate(TL, 1:N))
  d = data.frame(as.vector(t(y)), as.vector(t(x)), as.vector(id))
  colnames(d) = c("y", "x", "id")
  return(d)
}


#' Monte Carlo simulation for existence of threshold effects under known threshold location
#'
#' Monte Carlo simulation to study the size and power properties of the uniform test
#' for existence of threshold effects under known threshold locations. Provides the Monte Carlo
#' distribution of the test statistic and empirical rejection probabilities at 10%, 5% and 1% level.
#'
#' @param N cross-sectional dimension
#' @param TL time series length
#' @param p fraction of non-zero coefficients
#' @param M number of Monte Carlo runs
#' @param epsilon specification of error term. If \code{"iid"} is selected the error term is iid standard normal.
#' If \code{"factor"} is selected, the error term follows a factor model with strong cross-sectional and weak
#' temporal dependence.
#' @param running specification of running variable. If \code{"iid"} is selected the running variable is iid uniformly
#' distributed. If \code{"factor"} is selected, the running variable follows a factor model with strong cross-sectional and weak
#' temporal dependence.
#' @param hetero if \code{hetero=1} the error term is heteroskedastic, if \code{hetero=0} the error term is homoskedastic.
#'
#' @return A list containing the value of the test statistic for each Monte Carlo run and the
#' empirical rejection rate for a 10%, 5% and 1% confidence level.
#' @export
#'
#' @examples
#' result_threshold = simulation.threshold(10, 200, 0, 10, epsilon = "iid",
#'                    running = "iid", hetero = 0)
simulation.threshold = function(N,
                                TL,
                                p,
                                M,
                                epsilon = c("iid", "factor"),
                                running = c("iid", "factor"),
                                hetero = c(0, 1)) {
  test90 = test95 = test99 = TAU = vector(length = M)
  for (m in 1:M) {
    message(paste("Iteration ", m))
    if (epsilon == "iid") {
      eps = matrix(rnorm(N * TL, 0, 1), nrow = N)
    }
    if (epsilon == "factor") {
      beta = 1.5
      eps = t(replicate(N, MAinf_normal(TL, beta)))
      f = MAinf_normal(TL, beta)
      gamma = rnorm(N)
      nu = 0
      eps = gamma %*% t(f) + nu + eps
    }
    if (running == "iid") {
      X = matrix(runif(N * TL, min = -1, max = 1), nrow = N)
    }
    if (running == "factor") {
      beta = 1.5
      X = t(replicate(N, MAinf_normal(TL, beta)))
      f_X = MAinf_normal(TL, beta)
      gamma_X = rnorm(N)
      X = gamma_X %*% t(f_X) + X
      X = X / 4
    }
    tau = rep(0, N)
    tau[0:floor(N * p)] = TL ^ (-2 / 5) * sqrt(log(N)) * runif(floor(N *
                                                                      p), 2, 10)
    u = matrix(runif(N * TL, min = -1, max = 1), nrow = N)
    u = X * u
    if (hetero == 1) {
      sigma = abs(3 / 2 - abs(X)) * (3 / 2) ^ (2 * u)
      sigma = matrix(1, N, TL) + 0.25 * sigma
    }
    if (hetero == 0) {
      sigma = 1
    }
    Y = cos(X) + sin(u) + ifelse(X > 0, 1, 0) * tau + sigma * eps
    tauhat = vartau = vector(length = N)
    for (i in 1:N) {
      h = rdrobust::rdbwselect(Y[i, ], X[i, ], kernel = "uniform")$bws[1]
      S0p = sum(ifelse(X[i, ] > 0, K(X[i, ] / h), 0))
      S1p = sum(ifelse(X[i, ] > 0, X[i, ] * K(X[i, ] / h), 0))
      S2p = sum(ifelse(X[i, ] > 0, X[i, ] ^ 2 * K(X[i, ] / h), 0))
      S0n = sum(ifelse(X[i, ] <= 0, K(X[i, ] / h), 0))
      S1n = sum(ifelse(X[i, ] <= 0, X[i, ] * K(X[i, ] / h), 0))
      S2n = sum(ifelse(X[i, ] <= 0, X[i, ] ^ 2 * K(X[i, ] / h), 0))
      wp = ifelse(X[i, ] > 0, K(X[i, ] / h), 0) * (S2p - X[i, ] * S1p) / (S2p *
                                                                            S0p - S1p ^ 2)
      wn = ifelse(X[i, ] <= 0, K(X[i, ] / h), 0) * (S2n - X[i, ] * S1n) / (S2n *
                                                                             S0n - S1n ^ 2)
      tauhat[i] = sum((wp - wn) * Y[i, ])
      fit = KernSmooth::locpoly(X[i, ],
                    Y[i, ] - ifelse(X[i, ] > 0, 1, 0) * tauhat[i],
                    bandwidth = h,
                    truncate = FALSE)
      e = Y[i, ] - ifelse(X[i, ] > 0, 1, 0) * tauhat[i] - approx(fit$x, fit$y, xout =
                                                                   X[i, ])$y
      wpos = ifelse(X[i, ] > 0, ifelse(X[i, ] < h, 1, 0), 0)
      wneg = ifelse(X[i, ] < 0, ifelse(X[i, ] > -h, 1, 0), 0)
      vartau[i] = var(e[wpos + wneg != 0]) * sum(wp ^ 2 + wn ^ 2)
    }
    TAU[m] = max(abs(tauhat) / sqrt(vartau))
    test90[m] = as.numeric(ifelse(max(abs(tauhat) / sqrt(vartau)) > fdrtool::qhalfnorm((0.9) ^
                                                                                (1 / N)), 1, 0))
    test95[m] = as.numeric(ifelse(max(abs(tauhat) / sqrt(vartau)) > fdrtool::qhalfnorm((0.95) ^
                                                                                (1 / N)), 1, 0))
    test99[m] = as.numeric(ifelse(max(abs(tauhat) / sqrt(vartau)) > fdrtool::qhalfnorm((0.99) ^
                                                                                (1 / N)), 1, 0))
  }
  list = list(TAU, c(mean(test90), mean(test95), mean(test99)))
}

#' Monte Carlo simulation for heterogeneity of threshold effects under known threshold location
#'
#' Monte Carlo simulation to study the size and power properties of the uniform test
#' for heterogeneity of threshold effects under known threshold locations. Provides the Monte Carlo
#' distribution of the test statistic and empirical rejection probabilities at 10%, 5% and 1% level.
#'
#' @param N cross-sectional dimension
#' @param TL time series length
#' @param p fraction of non-zero coefficients
#' @param M number of Monte Carlo runs
#' @param epsilon specification of error term. If \code{"iid"} is selected the error term is iid standard normal.
#' If \code{"factor"} is selected, the error term follows a factor model with strong cross-sectional and weak
#' temporal dependence.
#' @param running specification of running variable. If \code{"iid"} is selected the running variable is iid uniformly
#' distributed. If \code{"factor"} is selected, the running variable follows a factor model with strong cross-sectional and weak
#' temporal dependence.
#' @param hetero if \code{hetero=1} the error term is heteroskedastic, if \code{hetero=0} the error term is homoskedastic.
#'
#' @return A list containing the value of the test statistic for each Monte Carlo run and the
#' empirical rejection rate for a 10%, 5% and 1% confidence level.
#' @export
#'
#' @examples
#' result_hetero = simulation.hetero(10, 200, 0, 10, epsilon = "iid",
#'                 running = "iid", hetero = 0)
simulation.hetero = function(N,
                             TL,
                             p,
                             M,
                             epsilon = c("iid", "factor"),
                             running = c("iid", "factor"),
                             hetero = c(0, 1)) {
  test90 = test95 = test99 = TAU = vector(length = M)
  for (m in 1:M) {
    message(paste("Iteration ", m))
    if (epsilon == "iid") {
      eps = matrix(rnorm(N * TL, 0, 1), nrow = N)
    }
    if (epsilon == "factor") {
      beta = 1.5
      eps = t(replicate(N, MAinf_normal(TL, beta)))
      f = MAinf_normal(TL, beta)
      gamma = rnorm(N)
      nu = 0
      eps = gamma %*% t(f) + nu + eps
    }
    if (running == "iid") {
      X = matrix(runif(N * TL, min = -1, max = 1), nrow = N)
    }
    if (running == "factor") {
      beta = 1.5
      X = t(replicate(N, MAinf_normal(TL, beta)))
      f_X = MAinf_normal(TL, beta)
      gamma_X = rnorm(N)
      X = gamma_X %*% t(f_X) + X
      X = X / 4
    }
    tau = rep(0, N)
    tau[0:floor(N * p)] = TL ^ (-2 / 5) * sqrt(log(N)) * runif(floor(N *
                                                                      p), 2, 10)
    u = matrix(runif(N * TL, min = -1, max = 1), nrow = N)
    u = X * u
    if (hetero == 1) {
      sigma = abs(3 / 2 - abs(X)) * (3 / 2) ^ (2 * u)
      sigma = matrix(1, N, TL) + 0.25 * sigma
    }
    if (hetero == 0) {
      sigma = 1
    }
    Y = cos(X) + sin(u) + ifelse(X > 0, 1, 0) * tau + sigma * eps
    tauhat = vartau = vector(length = N)
    for (i in 1:N) {
      h = rdrobust::rdbwselect(Y[i, ], X[i, ], kernel = "uniform")$bws[1]
      S0p = sum(ifelse(X[i, ] > 0, K(X[i, ] / h), 0))
      S1p = sum(ifelse(X[i, ] > 0, X[i, ] * K(X[i, ] / h), 0))
      S2p = sum(ifelse(X[i, ] > 0, X[i, ] ^ 2 * K(X[i, ] / h), 0))
      S0n = sum(ifelse(X[i, ] <= 0, K(X[i, ] / h), 0))
      S1n = sum(ifelse(X[i, ] <= 0, X[i, ] * K(X[i, ] / h), 0))
      S2n = sum(ifelse(X[i, ] <= 0, X[i, ] ^ 2 * K(X[i, ] / h), 0))
      wp = ifelse(X[i, ] > 0, K(X[i, ] / h), 0) * (S2p - X[i, ] * S1p) / (S2p *
                                                                            S0p - S1p ^ 2)
      wn = ifelse(X[i, ] <= 0, K(X[i, ] / h), 0) * (S2n - X[i, ] * S1n) / (S2n *
                                                                             S0n - S1n ^ 2)
      tauhat[i] = sum((wp - wn) * Y[i, ])
      fit = KernSmooth::locpoly(X[i, ],
                    Y[i, ] - ifelse(X[i, ] > 0, 1, 0) * tauhat[i],
                    bandwidth = h,
                    truncate = FALSE)
      e = Y[i, ] - ifelse(X[i, ] > 0, 1, 0) * tauhat[i] - approx(fit$x, fit$y, xout =
                                                                   X[i, ])$y
      wpos = ifelse(X[i, ] > 0, ifelse(X[i, ] < h, 1, 0), 0)
      wneg = ifelse(X[i, ] < 0, ifelse(X[i, ] > -h, 1, 0), 0)
      vartau[i] = var(e[wpos + wneg != 0]) * sum(wp ^ 2 + wn ^ 2)
    }
    vartautilde = vector(length = N)
    for (i in 1:N) {
      vartautilde[i] = vartau[i] * (1 - 1 / N) ^ 2 + sum(vartau[-i]) / N ^ 2
    }
    taubar = mean(tauhat)
    test90[m] = as.numeric(ifelse(max(
      abs(tauhat - taubar) / sqrt(vartautilde)
    ) > fdrtool::qhalfnorm((0.9) ^ (1 / N)), 1, 0))
    test95[m] = as.numeric(ifelse(max(
      abs(tauhat - taubar) / sqrt(vartautilde)
    ) > fdrtool::qhalfnorm((0.95) ^ (1 / N)), 1, 0))
    test99[m] = as.numeric(ifelse(max(
      abs(tauhat - taubar) / sqrt(vartautilde)
    ) > fdrtool::qhalfnorm((0.99) ^ (1 / N)), 1, 0))
    TAU[m] = max(abs(tauhat - taubar) / sqrt(vartautilde))
  }
  list = list(TAU, c(mean(test90), mean(test95), mean(test99)))
}

#' Monte Carlo simulation for existence of derivative threshold effects under known threshold location
#'
#' Monte Carlo simulation to study the size and power properties of the uniform test
#' for existence of threshold effects in the first derivative under known threshold locations. Provides the Monte Carlo
#' distribution of the test statistic and empirical rejection probabilities at 10%, 5% and 1% level.
#'
#' @param N cross-sectional dimension
#' @param TL time series length
#' @param p fraction of non-zero coefficients
#' @param M number of Monte Carlo runs
#' @param epsilon specification of error term. If \code{"iid"} is selected the error term is iid standard normal.
#' If \code{"factor"} is selected, the error term follows a factor model with strong cross-sectional and weak
#' temporal dependence.
#' @param running specification of running variable. If \code{"iid"} is selected the running variable is iid uniformly
#' distributed. If \code{"factor"} is selected, the running variable follows a factor model with strong cross-sectional and weak
#' temporal dependence.
#' @param hetero if \code{hetero=1} the error term is heteroskedastic, if \code{hetero=0} the error term is homoskedastic.
#'
#' @return A list containing the value of the test statistic for each Monte Carlo run and the
#' empirical rejection rate for a 10%, 5% and 1% confidence level.
#' @export
#'
#' @examples
#' result_derivative = simulation.derivative(10, 200, 0, 10, epsilon = "iid",
#'                     running = "iid", hetero = 0)
simulation.derivative = function(N,
                                 TL,
                                 p,
                                 M,
                                 epsilon = c("iid", "factor"),
                                 running = c("iid", "factor"),
                                 hetero = c(0, 1)) {
  test90 = test95 = test99 = TAU = vector(length = M)
  for (m in 1:M) {
    message(paste("Iteration ", m))
    if (epsilon == "iid") {
      eps = matrix(rnorm(N * TL, 0, 1), nrow = N)
    }
    if (epsilon == "factor") {
      beta = 1.5
      eps = t(replicate(N, MAinf_normal(TL, beta)))
      f = MAinf_normal(TL, beta)
      gamma = rnorm(N)
      nu = 0
      eps = gamma %*% t(f) + nu + eps
    }
    if (running == "iid") {
      X = matrix(runif(N * TL, min = -1, max = 1), nrow = N)
    }
    if (running == "factor") {
      beta = 1.5
      X = t(replicate(N, MAinf_normal(TL, beta)))
      f_X = MAinf_normal(TL, beta)
      gamma_X = rnorm(N)
      X = gamma_X %*% t(f_X) + X
      X = X / 4
    }
    tau = rep(0, N)
    tau[0:floor(N * p)] = 1
    u = matrix(runif(N * TL, min = -1, max = 1), nrow = N)
    u = X * u
    if (hetero == 1) {
      sigma = abs(3 / 2 - abs(X)) * (3 / 2) ^ (2 * u)
      sigma = matrix(1, N, TL) + 0.25 * sigma
    }
    if (hetero == 0) {
      sigma = 1
    }
    Y = p * abs(X) * 5 + (1 - p) * X * 5 + sigma * eps
    tauhat = tauhat1 = vartau1 = vector(length = N)
    for (i in 1:N) {
      h = rdrobust::rdbwselect(Y[i, ], X[i, ], kernel = "uniform")$bws[1] * sqrt(3)
      S0p = sum(ifelse(X[i, ] > 0, K(X[i, ] / h), 0))
      S1p = sum(ifelse(X[i, ] > 0, X[i, ] * K(X[i, ] / h), 0))
      S2p = sum(ifelse(X[i, ] > 0, X[i, ] ^ 2 * K(X[i, ] / h), 0))
      S0n = sum(ifelse(X[i, ] <= 0, K(X[i, ] / h), 0))
      S1n = sum(ifelse(X[i, ] <= 0, X[i, ] * K(X[i, ] / h), 0))
      S2n = sum(ifelse(X[i, ] <= 0, X[i, ] ^ 2 * K(X[i, ] / h), 0))
      wp = ifelse(X[i, ] > 0, K(X[i, ] / h), 0) * (S2p - X[i, ] * S1p) / (S2p *
                                                                            S0p - S1p ^ 2)
      wn = ifelse(X[i, ] <= 0, K(X[i, ] / h), 0) * (S2n - X[i, ] * S1n) / (S2n *
                                                                             S0n - S1n ^ 2)
      wp1 = ifelse(X[i, ] > 0, K((X[i, ] - 0) / h), 0) * (S1p - (X[i, ] -
                                                                   0) * S0p) / (-S2p * S0p + S1p ^ 2)
      wn1 = ifelse(X[i, ] <= 0, K((X[i, ] - 0) / h), 0) * (S1n - (X[i, ] -
                                                                    0) * S0n) / (-S2n * S0n + S1n ^ 2)
      tauhat[i] = sum((wp - wn) * Y[i, ])
      fit = KernSmooth::locpoly(X[i, ],
                    Y[i, ] - ifelse(X[i, ] > 0, 1, 0) * tauhat[i],
                    bandwidth = h,
                    truncate = FALSE)
      e = Y[i, ] - ifelse(X[i, ] > 0, 1, 0) * tauhat[i] - approx(fit$x, fit$y, xout =
                                                                   X[i, ])$y
      wpos = ifelse(X[i, ] > 0, ifelse(X[i, ] < h, 1, 0), 0)
      wneg = ifelse(X[i, ] < 0, ifelse(X[i, ] > -h, 1, 0), 0)
      tauhat1[i] = sum((wp1 - wn1) * Y[i, ])
      vartau1[i] = var(e[wpos + wneg != 0]) * sum(wp1 ^ 2 + wn1 ^ 2)
    }
    TAU[m] = max(abs(tauhat) / sqrt(vartau1))
    test90[m] = as.numeric(ifelse(max(abs(tauhat1) / sqrt(vartau1)) > fdrtool::qhalfnorm((0.9) ^
                                                                                  (1 / N)), 1, 0))
    test95[m] = as.numeric(ifelse(max(abs(tauhat1) / sqrt(vartau1)) > fdrtool::qhalfnorm((0.95) ^
                                                                                  (1 / N)), 1, 0))
    test99[m] = as.numeric(ifelse(max(abs(tauhat1) / sqrt(vartau1)) > fdrtool::qhalfnorm((0.99) ^
                                                                                  (1 / N)), 1, 0))
  }
  list = list(TAU, c(mean(test90), mean(test95), mean(test99)))
}

#' Monte Carlo simulation for uniform test of existence of threshold effects under unknown threshold location
#'
#' Monte Carlo simulation to study the size and power properties of the uniform test
#' for existence of threshold effects under unknown threshold locations. Provides the Monte Carlo
#' distribution of the test statistic and empirical rejection probabilities at 10%, 5% and 1% level.
#'
#' @param N cross-sectional dimension
#' @param TL time series length
#' @param p fraction of non-zero coefficients
#' @param M number of Monte Carlo runs
#' @param epsilon specification of error term. If \code{"iid"} is selected the error term is iid standard normal.
#' If \code{"factor"} is selected, the error term follows a factor model with strong cross-sectional and weak
#' temporal dependence.
#' @param running specification of running variable. If \code{"iid"} is selected the running variable is iid uniformly
#' distributed. If \code{"factor"} is selected, the running variable follows a factor model with strong cross-sectional and weak
#' temporal dependence.
#' @param hetero if \code{hetero=1} the error term is heteroskedastic, if \code{hetero=0} the error term is homoskedastic.
#' @param threshold specifies the distribution for the non-zero threshold coefficients, possible values are \code{"normal"} for the
#' standard normel, \code{"exponential"} for an exponential distribution with parameter 1, or \code{"uniform"} for a uniform distribution.
#'
#' @return A list containing the value of the test statistic for each Monte Carlo run and the
#' empirical rejection rate for a 10%, 5% and 1% confidence level.
#' @export
#'
#' @examples
#' result_unknown = simulation.unknown(2, 800, 0, 10, epsilon = "iid", running = "iid",
#'                  hetero = 0, threshold = "gaussian")
simulation.unknown = function(N,
                              TL,
                              p,
                              M,
                              epsilon = c("iid", "factor"),
                              running = c("iid", "factor"),
                              hetero = c(0, 1),
                              threshold = c("uniform", "exponential", "gaussian")) {
  test90 = test95 = test99 = TAU = vector(length = M)
  for (m in 1:M) {
    message(paste("Iteration ", m))
    if (epsilon == "iid") {
      eps = matrix(rnorm(N * TL, 0, 1), nrow = N)
    }
    if (epsilon == "factor") {
      beta = 1.5
      eps = t(replicate(N, MAinf_normal(TL, beta)))
      f = MAinf_normal(TL, beta)
      gamma = rnorm(N)
      nu = 0
      eps = gamma %*% t(f) + nu + eps
    }
    if (running == "iid") {
      X = matrix(runif(N * TL, min = -1, max = 1), nrow = N)
    }
    if (running == "factor") {
      beta = 1.5
      X = t(replicate(N, MAinf_normal(TL, beta)))
      f_X = MAinf_normal(TL, beta)
      gamma_X = rnorm(N)
      X = gamma_X %*% t(f_X) + X
      X = X / 4
    }
    tau = rep(0, N)
    if (threshold == "uniform") {
      tau[0:floor(N * p)] = TL ^ (-2 / 5) * sqrt(log(N)) * runif(floor(N * p), 2, 10)
    }
    if (threshold == "exponential") {
      tau[0:floor(N * p)] = rexp(floor(N * p))
    }
    if (threshold == "gaussian") {
      tau[0:floor(N * p)] = rnorm(floor(N * p))
    }
    u = matrix(runif(N * TL, min = -1, max = 1), nrow = N)
    u = X * u
    if (hetero == 1) {
      sigma = abs(3 / 2 - abs(X)) * (3 / 2) ^ (2 * u)
      sigma = matrix(1, N, TL) + 0.25 * sigma
    }
    if (hetero == 0) {
      sigma = 1
    }
    Y = cos(X) + sin(u) + ifelse(X > 0, 1, 0) * tau + sigma * eps
    if (running == "factor") {
      C = c(-0.3, -0.15, 0, 0.15, 0.3)
    }
    if (running == "iid") {
      C = c(-0.6, -0.3, 0, 0.3, 0.6)
    }
    tauhat = vartau = matrix(nrow = N, ncol = length(C))
    for (j in 1:length(C)) {
      c = C[j]
      for (i in 1:N) {
        h = rdrobust::rdbwselect(Y[i, ], X[i, ], kernel = "uniform", c = c)$bws[1]
        S0p = sum(ifelse(X[i, ] - c > 0, K((X[i, ] - c) / h), 0))
        S1p = sum(ifelse(X[i, ] - c > 0, (X[i, ] - c) * K((X[i, ] - c) /
                                                            h), 0))
        S2p = sum(ifelse(X[i, ] - c > 0, (X[i, ] - c) ^ 2 * K((X[i, ] -
                                                                 c) / h), 0))
        S0n = sum(ifelse(X[i, ] - c <= 0, K((X[i, ] - c) / h), 0))
        S1n = sum(ifelse(X[i, ] - c <= 0, (X[i, ] - c) * K((X[i, ] - c) /
                                                             h), 0))
        S2n = sum(ifelse(X[i, ] - c <= 0, (X[i, ] - c) ^ 2 * K((X[i, ] -
                                                                  c) / h), 0))
        wp = ifelse(X[i, ] > c, K((X[i, ] - c) / h), 0) * (S2p - (X[i, ] -
                                                                    c) * S1p) / (S2p * S0p - S1p ^ 2)
        wn = ifelse(X[i, ] <= c, K((X[i, ] - c) / h), 0) * (S2n - (X[i, ] -
                                                                     c) * S1n) / (S2n * S0n - S1n ^ 2)
        tauhat[i, j] = sum((wp - wn) * Y[i, ])
        fit = KernSmooth::locpoly(
          X[i, ],
          Y[i, ] - ifelse(X[i, ] > c, 1, 0) * tauhat[i, j],
          bandwidth = h,
          truncate = TRUE
        )
        e = Y[i, ] - ifelse(X[i, ] > c, 1, 0) * tauhat[i, j] - approx(fit$x, fit$y, xout =
                                                                        X[i, ])$y
        wpos = ifelse(X[i, ] > c, ifelse((X[i, ] - c) < h, 1, 0), 0)
        wneg = ifelse(X[i, ] <= c, ifelse((X[i, ] - c) > -h, 1, 0), 0)
        vartau[i, j] = var(e[wpos + wneg != 0]) * sum(wp ^ 2 + wn ^ 2)
      }
    }
    TAU[m] = max(abs(tauhat) / sqrt(vartau))
    test90[m] = as.numeric(ifelse(max(abs(tauhat) / sqrt(vartau)) > fdrtool::qhalfnorm((0.9) ^
                                                                                (1 / (N * length(C)))), 1, 0))
    test95[m] = as.numeric(ifelse(max(abs(tauhat) / sqrt(vartau)) > fdrtool::qhalfnorm((0.95) ^
                                                                                (1 / (N * length(C)))), 1, 0))
    test99[m] = as.numeric(ifelse(max(abs(tauhat) / sqrt(vartau)) > fdrtool::qhalfnorm((0.99) ^
                                                                                (1 / (N * length(C)))), 1, 0))
  }
  list = list(TAU, c(mean(test90), mean(test95), mean(test99)))
}

#' Monte Carlo simulation for pooled test of existence of threshold effects
#'
#' Monte Carlo simulation to study the size and power properties of the pooled test
#' for existence of threshold effects under unknown threshold locations. The pooled test can
#' be based on a nonparametric regression model or a linear panel threshold regression model.
#' Provides the Monte Carlo distribution of the test statistic and empirical rejection probabilities
#' at 10%, 5% and 1% level.
#'
#' @param N cross-sectional dimension
#' @param TL time series length
#' @param p fraction of non-zero coefficients
#' @param M number of Monte Carlo runs
#' @param epsilon specification of error term. If \code{"iid"} is selected the error term is iid standard normal.
#' If \code{"factor"} is selected, the error term follows a factor model with strong cross-sectional and weak
#' temporal dependence.
#' @param running specification of running variable. If \code{"iid"} is selected the running variable is iid uniformly
#' distributed. If \code{"factor"} is selected, the running variable follows a factor model with strong cross-sectional and weak
#' temporal dependence.
#' @param hetero if \code{hetero=1} the error term is heteroskedastic, if \code{hetero=0} the error term is homoskedastic.
#' @param threshold specifies the distribution for the non-zero threshold coefficients, possible values are \code{"normal"} for the
#' standard normel, \code{"exponential"} for an exponential distribution with parameter 1, or \code{"uniform"} for a uniform distribution.
#' @param method method of estimation (\code{"nonparametric"} vs. \code{"parametric"})
#'
#' @return A list containing the value of the test statistic for each Monte Carlo run and the
#' empirical rejection rate for a 10%, 5% and 1% confidence level.
#' @export
#'
#' @examples
#' result_pooled = simulation.pooled(5, 400, 0, 10, epsilon = "iid", running = "iid",
#'                 hetero = 0, threshold = "gaussian", method = "nonparametric")
simulation.pooled = function(N,
                             TL,
                             p,
                             M,
                             epsilon = c("iid", "factor"),
                             running = c("iid", "factor"),
                             hetero = c(0, 1),
                             threshold = c("uniform", "exponential", "gaussian"),
                             method = c("parametric", "nonparametric")) {
  test90 = test95 = test99 = TAU = vector(length = M)
  for (m in 1:M) {
    message(paste("Iteration ", m))
    if (epsilon == "iid") {
      eps = matrix(rnorm(N * TL, 0, 1), nrow = N)
    }
    if (epsilon == "factor") {
      beta = 1.5
      eps = t(replicate(N, MAinf_normal(TL, beta)))
      f = MAinf_normal(TL, beta)
      gamma = rnorm(N)
      nu = 0
      eps = gamma %*% t(f) + nu + eps
    }
    if (running == "iid") {
      X = matrix(runif(N * TL, min = -1, max = 1), nrow = N)
    }
    if (running == "factor") {
      beta = 1.5
      X = t(replicate(N, MAinf_normal(TL, beta)))
      f_X = MAinf_normal(TL, beta)
      gamma_X = rnorm(N)
      X = gamma_X %*% t(f_X) + X
      X = X / 4
    }
    tau = rep(0, N)
    if (threshold == "uniform") {
      tau[0:floor(N * p)] = TL ^ (-2 / 5) * sqrt(log(N)) * runif(floor(N * p), 2, 10)
    }
    if (threshold == "exponential") {
      tau[0:floor(N * p)] = rexp(floor(N * p))
    }
    if (threshold == "gaussian") {
      tau[0:floor(N * p)] = rnorm(floor(N * p))
    }
    u = matrix(runif(N * TL, min = -1, max = 1), nrow = N)
    u = X * u
    if (hetero == 1) {
      sigma = abs(3 / 2 - abs(X)) * (3 / 2) ^ (2 * u)
      sigma = matrix(1, N, TL) + 0.25 * sigma
    }
    if (hetero == 0) {
      sigma = 1
    }
    Y = cos(X) + sin(u) + ifelse(X > 0, 1, 0) * tau + sigma * eps
    tauhat = vartau = vector(length = length(C))
    X = as.vector(X)
    Y = as.vector(Y)
    if (running == "factor") {
      C = c(-0.3, -0.15, 0, 0.15, 0.3)
    }
    if (running == "iid") {
      C = c(-0.6, -0.3, 0, 0.3, 0.6)
    }
    for (j in 1:length(C)) {
      c = C[j]
      if (method == "parametric") {
        h = 20
      }
      if (method == "nonparametric") {
        h = rdrobust::rdbwselect(Y, X, kernel = "uniform", c = c)$bws[1]
      }
      S0p = sum(ifelse(X - c > 0, K((X - c) / h), 0))
      S1p = sum(ifelse(X - c > 0, (X - c) * K((X - c) / h), 0))
      S2p = sum(ifelse(X - c > 0, (X - c) ^ 2 * K((X - c) / h), 0))
      S0n = sum(ifelse(X - c <= 0, K((X - c) / h), 0))
      S1n = sum(ifelse(X - c <= 0, (X - c) * K((X - c) / h), 0))
      S2n = sum(ifelse(X - c <= 0, (X - c) ^ 2 * K((X - c) / h), 0))
      wp = ifelse(X > c, K((X - c) / h), 0) * (S2p - (X - c) * S1p) / (S2p *
                                                                         S0p - S1p ^ 2)
      wn = ifelse(X <= c, K((X - c) / h), 0) * (S2n - (X - c) * S1n) / (S2n *
                                                                          S0n - S1n ^ 2)
      tauhat[j] = sum((wp - wn) * Y)
      fit = KernSmooth::locpoly(X,
                    Y - ifelse(X > c, 1, 0) * tauhat[j],
                    bandwidth = h,
                    truncate = TRUE)
      e = Y - ifelse(X > c, 1, 0) * tauhat[j] - approx(fit$x, fit$y, xout =
                                                         X)$y
      wpos = ifelse(X > c, ifelse((X - c) < h, 1, 0), 0)
      wneg = ifelse(X <= c, ifelse((X - c) > -h, 1, 0), 0)
      vartau[j] = var(e[wpos + wneg != 0]) * sum(wp ^ 2 + wn ^ 2)
    }
    TAU[m] = max(abs(tauhat) / sqrt(vartau))
    test90[m] = as.numeric(ifelse(max(abs(tauhat) / sqrt(vartau)) > fdrtool::qhalfnorm((0.9) ^
                                                                                (1 / (length(C)))), 1, 0))
    test95[m] = as.numeric(ifelse(max(abs(tauhat) / sqrt(vartau)) > fdrtool::qhalfnorm((0.95) ^
                                                                                (1 / (length(C)))), 1, 0))
    test99[m] = as.numeric(ifelse(max(abs(tauhat) / sqrt(vartau)) > fdrtool::qhalfnorm((0.99) ^
                                                                                (1 / (length(C)))), 1, 0))
  }
  list = list(TAU, c(mean(test90), mean(test95), mean(test99)))
}



