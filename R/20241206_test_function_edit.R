
#' Uniform test for existence of threshold effects
#'
#' Uniform test for existence of threshold effects in a nonparametric panel regression. Both the
#' known and unknown threshold location case are covered. Apart from the uniform test statistic and
#' the corresponding p-value, a table for the results of the individual threshold estimates and test statistics
#' is provided.
#'
#' @param data a data frame containing the response, running and id variables
#' @param response name of the dependent variable (aka response variable)
#' @param running name of the running variable (aka forcing variable)
#' @param id name of the id variable
#' @param bw an optional scalar bandwidth parameter for the local linear estimation. If not specified, the bandwidth
#' is selected by the command [rdrobust::rdbwselect()].
#' @param C a scalar value for the true threshold location (for the known case) or a grid of candidate threshold locations
#' (for the unknown case)
#' @param alpha specifies a threshold to determine which and how many individual-specific threshold effects and test
#' statistics are displayed in the output table. Only individuals which are significant at the alpha confidence level are selected.
#' @param alternative specifies whether we consider a two-sided alternative (default) or left-/right-sided alternative.
#'
#' @seealso [threshold.example()], [rdrobust::rdbwselect()].
#'
#' @return A list containing:\tabular{ll}{
#'    \code{I_hat} \tab the value of the uniform test statistic \cr
#'    \tab \cr
#'    \code{p_value} \tab the corresponding p-value \cr
#'    \tab \cr
#'    \code{N} \tab the cross-sectional dimension \cr
#'    \tab \cr
#'    \code{Critical_values} \tab critical values at 10%, 5%, 1%, and 0.1% confidence level \cr
#'    \tab \cr
#'    \code{Table} \tab a table displaying the estimation result for a selection of individuals, including the id variable, the threshold location,
#'    the estimated coefficient, the estimated standard error, and the individual test statistic. \cr
#' }
#' @export
#'
#' @examples
#' d = threshold.example(10, 200, 0.1, 2)
#' threshold.test(data = d, response = "y", running = "x", id = "id", C = 0)
threshold.test = function(data,
                          response,
                          running,
                          id,
                          bw = NULL,
                          C = 0,
                          alpha = NULL,
                          alternative = "two") {
  names = unique(data[, id])
  N = length(names)
  result = data.frame(matrix(ncol = 5, nrow = 0))
  colnames(result) = c("ID", "c_0j", "gamma_hat_j", "v_hat_j", "I_hat_j")
  result_j = data.frame(matrix(nrow = length(C), ncol = 5))
  colnames(result_j) = c("ID", "c_0j", "gamma_hat_j", "v_hat_j", "I_hat_j")
  for (j in 1:N) {
    message(paste("Individual ", j))
    y = data[data[, id] == names[j], response]
    x = data[data[, id] == names[j], running]
    result_j[, 1] = names[j]
    result_j[, 2] = C
    for (i in 1:length(C)) {
      c = C[i]
      h = ifelse(is.null(bw),
                 rdrobust::rdbwselect(y, x, kernel = "uniform", c = c)$bws[1],
                 bw)
      K = function(x) {
        ifelse((abs(x)) <= 1, 1 / 2, 0)
      }
      S0p = sum(ifelse(x - c > 0, K((x - c) / h), 0))
      S1p = sum(ifelse(x - c > 0, (x - c) * K((x - c) / h), 0))
      S2p = sum(ifelse(x - c > 0, (x - c) ^ 2 * K((x - c) / h), 0))
      S0n = sum(ifelse(x - c <= 0, K((x - c) / h), 0))
      S1n = sum(ifelse(x - c <= 0, (x - c) * K((x - c) / h), 0))
      S2n = sum(ifelse(x - c <= 0, (x - c) ^ 2 * K((x - c) / h), 0))
      wp = ifelse(x - c > 0, K((x - c) / h), 0) * (S2p - (x - c) * S1p) / (S2p *
                                                                             S0p - S1p ^ 2)
      wn = ifelse(x - c <= 0, K((x - c) / h), 0) * (S2n - (x - c) * S1n) / (S2n *
                                                                              S0n - S1n ^ 2)
      tauhat = sum((wp - wn) * y)
      result_j[i, 3] = tauhat
      fit1 = KernSmooth::locpoly(x,
                     y - ifelse(x > c, 1, 0) * tauhat,
                     bandwidth = h,
                     truncate = FALSE)
      e = y - ifelse((x - c) > 0, 1, 0) * tauhat - approx(fit1$x, fit1$y, xout =
                                                            x)$y
      wpos = ifelse(x - c > 0, ifelse((x - c) < h, 1, 0), 0)
      wneg = ifelse(x - c < 0, ifelse((x - c) > -h, 1, 0), 0)
      result_j[i, 4] = sqrt(var(e[wpos + wneg != 0]) * sum(wp ^ 2 + wn ^
                                                             2))
    }
    result_j[, 5] = as.numeric(result_j[, 3]) / as.numeric(result_j[, 4])
    if (alternative == "two") {
      result_j_max = result_j[abs(result_j[, 5]) == max(abs(result_j[, 5])), ]
    }
    if (alternative == "right") {
      result_j_max = result_j[result_j[, 5] == max(result_j[, 5]), ]
    }
    if (alternative == "left") {
      result_j_max = result_j[result_j[, 5] == min(result_j[, 5]), ]
    }
    result = rbind(result, result_j_max)
  }
  I_hat = ifelse(alternative == "left", min(result[, 5]), max(abs(result[, 5])))
  if (alternative == "two") {
    p_value = 1 - fdrtool::phalfnorm(I_hat) ^ (N * length(C))
  }
  if (alternative == "right") {
    p_value = 1 - pnorm(I_hat) ^ (N * length(C))
  }
  if (alternative == "left") {
    p_value = 1 - pnorm(-I_hat) ^ (N * length(C))
  }
  alpha = ifelse(is.null(alpha), 1, alpha)
  if (alternative == "two") {
    result_1 = subset(result, abs(result$I_hat_j) > fdrtool::qhalfnorm((1 - alpha) ^ (1 / (N *
                                                                             length(
                                                                               C
                                                                             )))))
  }
  else {
    result_1 = subset(result, abs(result$I_hat_j) > qnorm((1 - alpha) ^ (1 / (N * length(
      C
    )))))
  }
  result_1 = result_1[order(result_1$c_0j), ]
  result_1[, 2:5] = round(result_1[, 2:5], digits = 3)
  rownames(result_1) <- NULL
  if (alternative == "two") {
    critical = c(fdrtool::qhalfnorm((0.90) ^ (1 / (N * length(
      C
    )))),
    fdrtool::qhalfnorm((0.95) ^ (1 / (N * length(
      C
    )))),
    ###5
    fdrtool::qhalfnorm((0.99) ^ (1 / (N * length(
      C
    )))),
    fdrtool::qhalfnorm((0.999) ^ (1 / (N * length(
      C
    )))))
  }
  if (alternative == "right") {
    critical = c(qnorm((0.90) ^ (1 / (N * length(
      C
    )))), qnorm((0.95) ^ (1 / (N * length(
      C
    )))),
    qnorm((0.99) ^ (1 / (N * length(
      C
    )))), qnorm((0.999) ^ (1 / (N * length(
      C
    )))))
  }
  if (alternative == "left") {
    critical = -c(qnorm((0.90) ^ (1 / (N * length(
      C
    )))), qnorm((0.95) ^ (1 / (N * length(
      C
    )))),
    qnorm((0.99) ^ (1 / (N * length(
      C
    )))), qnorm((0.999) ^ (1 / (N * length(
      C
    )))))
  }
  critical = setNames(critical, c("90%", "95%", "99%", "99.9%"))

  list = list(I_hat, p_value, N, critical, result_1)
  list = setNames(list, c("I_hat", "p_value", "N", "Critical_values", "Table"))
  return(list)
}


#' Uniform test for heterogeneity of threshold effects
#'
#' Uniform test for heterogeneity of threshold effects in a nonparametric panel regression under known threshold locations.
#' Apart from the uniform test statistic and the corresponding p-value, a table for the results of the individual threshold
#' estimates and test statistics is provided.
#'
#' @param data a data frame containing the response, running and id variables
#' @param response name of the dependent variable (aka response variable)
#' @param running name of the running variable (aka forcing variable)
#' @param id name of the id variable
#' @param bw an optional scalar bandwidth parameter for the local linear estimation. If not specified, the bandwidth
#' is selected by the command [rdrobust::rdbwselect()].
#' @param c a scalar value for the true threshold location
#' @param alpha specifies a threshold to determine which and how many individual-specific threshold effects and test
#' statistics are displayed in the output table. Only individuals which are significant at the alpha confidence level are selected.
#' @param alternative specifies whether we consider a two-sided alternative (default) or left-/right-sided alternative.
#' @param use.median if \code{TRUE}, the median replaces the mean as a robust alternative in the test for heterogeneity
#'
#' @seealso [threshold.example()], [rdrobust::rdbwselect()].
#'
#' @return A list containing:\tabular{ll}{
#'    \code{Q_hat} \tab the value of the uniform test statistic \cr
#'    \tab \cr
#'    \code{p_value} \tab the corresponding p-value \cr
#'    \tab \cr
#'    \code{N} \tab the cross-sectional dimension \cr
#'    \tab \cr
#'    \code{Critical_values} \tab critical values at 10%, 5%, 1%, and 0.1% confidence level \cr
#'    \tab \cr
#'    \code{Table} \tab a table displaying the estimation result for a selection of individuals, including the id variable, the threshold location,
#'    the estimated coefficient, the estimated standard error, and the individual test statistic. \cr
#' }
#' @export
#'
#' @examples
#' d = threshold.example(10, 200, 0.1, 2)
#' threshold.heterogeneity.test(data = d, response = "y", running = "x", id = "id", c = 0)
threshold.heterogeneity.test = function(data,
                                        response,
                                        running,
                                        id,
                                        bw = NULL,
                                        c = 0,
                                        alpha = NULL,
                                        alternative = "two",
                                        use.median = FALSE) {
  names = unique(data[, id])
  N = length(names)
  result = data.frame(matrix(ncol = 5, nrow = 0))
  colnames(result) = c("ID",
                       "c_0j",
                       "gamma_hat_j-gamma_bar",
                       "v_tilde_j",
                       "Q_hat_j")
  result_j = data.frame(matrix(nrow = 1, ncol = 5))
  colnames(result_j) = c("ID",
                         "c_0j",
                         "gamma_hat_j-gamma_bar",
                         "v_tilde_j",
                         "Q_hat_j")
  for (j in 1:N) {
    message(paste("Individual ", j))
    y = data[data[, id] == names[j], response]
    x = data[data[, id] == names[j], running]
    result_j[, 1] = names[j]
    result_j[, 2] = c
    h = ifelse(is.null(bw),
               rdrobust::rdbwselect(y, x, kernel = "uniform", c = c)$bws[1],
               bw)
    K = function(x) {
      ifelse((abs(x)) <= 1, 1 / 2, 0)
    }
    S0p = sum(ifelse(x - c > 0, K((x - c) / h), 0))
    S1p = sum(ifelse(x - c > 0, (x - c) * K((x - c) / h), 0))
    S2p = sum(ifelse(x - c > 0, (x - c) ^ 2 * K((x - c) / h), 0))
    S0n = sum(ifelse(x - c <= 0, K((x - c) / h), 0))
    S1n = sum(ifelse(x - c <= 0, (x - c) * K((x - c) / h), 0))
    S2n = sum(ifelse(x - c <= 0, (x - c) ^ 2 * K((x - c) / h), 0))
    wp = ifelse(x - c > 0, K((x - c) / h), 0) * (S2p - (x - c) * S1p) / (S2p *
                                                                           S0p - S1p ^ 2)
    wn = ifelse(x - c <= 0, K((x - c) / h), 0) * (S2n - (x - c) * S1n) / (S2n *
                                                                            S0n - S1n ^ 2)
    tauhat = sum((wp - wn) * y)
    result_j[, 3] = tauhat
    fit1 = KernSmooth::locpoly(x,
                   y - ifelse(x > c, 1, 0) * tauhat,
                   bandwidth = h,
                   truncate = FALSE)
    e = y - ifelse((x - c) > 0, 1, 0) * tauhat - approx(fit1$x, fit1$y, xout =
                                                          x)$y
    wpos = ifelse(x - c > 0, ifelse((x - c) < h, 1, 0), 0)
    wneg = ifelse(x - c < 0, ifelse((x - c) > -h, 1, 0), 0)
    result_j[, 4] = sqrt(var(e[wpos + wneg != 0]) * sum(wp ^ 2 + wn ^ 2))
    result_j[, 5] = as.numeric(result_j[, 3]) / as.numeric(result_j[, 4])
    result = rbind(result, result_j)
  }
  if (use.median == TRUE) {
    result[, 3] = result[, 3] - median(result[, 3])
  }
  else{
    result[, 3] = result[, 3] - mean(result[, 3])
  }
  for (j in 1:N) {
    result[j, 4] = result[j, 4] * (1 - 1 / N) ^ 2 + sum(result[-j, 4]) / N ^
      2
  }
  result[, 5] = as.numeric(result[, 3]) / as.numeric(result[, 4])
  Q_hat = ifelse(alternative == "left", min(result$Q_hat_j), max(abs(result$Q_hat_j)))
  if (alternative == "two") {
    p_value = 1 - fdrtool::phalfnorm(Q_hat) ^ N
  }
  if (alternative == "right") {
    p_value = 1 - pnorm(Q_hat) ^ N
  }
  if (alternative == "left") {
    p_value = 1 - pnorm(-Q_hat) ^ N
  }
  alpha = ifelse(is.null(alpha), 1, alpha)
  if (alternative == "two") {
    result_1 = subset(result, abs(result$Q_hat_j) > fdrtool::qhalfnorm((1 - alpha) ^ (1 / N)))
  }
  else {
    result_1 = subset(result, abs(result$Q_hat_j) > qnorm((1 - alpha) ^ (1 / N)))
  }
  result_1 = result_1[order(result_1$c_0j), ]
  result_1[, 2:5] = round(result_1[, 2:5], digits = 3)
  rownames(result_1) <- NULL

  if (alternative == "two") {
    critical = c(fdrtool::qhalfnorm((0.90) ^ (1 / N)),
                 fdrtool::qhalfnorm((0.95) ^ (1 / N)),
                 fdrtool::qhalfnorm((0.99) ^ (1 / N)),
                 fdrtool::qhalfnorm((0.999) ^ (1 / N)))
  }
  if (alternative == "right") {
    critical = c(qnorm((0.90) ^ (1 / N)), qnorm((0.95) ^ (1 / N)),
                 qnorm((0.99) ^ (1 / N)), qnorm((0.999) ^ (1 / N)))
  }
  if (alternative == "left") {
    critical = -c(qnorm((0.90) ^ (1 / N)), qnorm((0.95) ^ (1 / N)),
                  qnorm((0.99) ^ (1 / N)), qnorm((0.999) ^ (1 / N)))
  }
  critical = setNames(critical, c("90%", "95%", "99%", "99.9%"))
  list = list(Q_hat, p_value, N, critical, result_1)
  list = setNames(list, c("Q_hat", "p_value", "N", "Critical_values", "Table"))
  return(list)
}


#' Uniform test for existence of derivative threshold effects
#'
#' Uniform test for existence of threshold effects in the first derivative for nonparametric panel regressions. Both the
#' known and unknown threshold location case are covered. Apart from the uniform test statistic and
#' the corresponding p-value, a table for the results of the individual threshold estimates and test statistics
#' is provided.
#'
#' @param data a data frame containing the response, running and id variables
#' @param response name of the dependent variable (aka response variable)
#' @param running name of the running variable (aka forcing variable)
#' @param id name of the id variable
#' @param bw an optional scalar bandwidth parameter for the local linear estimation. If not specified, the bandwidth
#' is selected by the command [rdrobust::rdbwselect()].
#' @param C a scalar value for the true threshold location (for the known case) or a grid of candidate threshold locations
#' (for the unknown case)
#' @param alpha specifies a threshold to determine which and how many individual-specific threshold effects and test
#' statistics are displayed in the output table. Only individuals which are significant at the alpha confidence level are selected.
#' @param alternative specifies whether we consider a two-sided alternative (default) or left-/right-sided alternative.
#'
#' @seealso [threshold.example()], [rdrobust::rdbwselect()].
#'
#' @return A list containing:\tabular{ll}{
#'    \code{I_hat} \tab the value of the uniform test statistic \cr
#'    \tab \cr
#'    \code{p_value} \tab the corresponding p-value \cr
#'    \tab \cr
#'    \code{N} \tab the cross-sectional dimension \cr
#'    \tab \cr
#'    \code{Critical_values} \tab critical values at 10%, 5%, 1%, and 0.1% confidence level \cr
#'    \tab \cr
#'    \code{Table} \tab a table displaying the estimation result for a selection of individuals, including the id variable, the threshold location,
#'    the estimated coefficient, the estimated standard error, and the individual test statistic. \cr
#' }
#' @export
#'
#' @examples
#' d = threshold.example(10, 200, 0.1, 2)
#' threshold.derivative.test(data = d, response = "y", running = "x", id = "id", C = 0)
threshold.derivative.test = function(data,
                                     response,
                                     running,
                                     id,
                                     bw = NULL,
                                     C = 0,
                                     alpha = NULL,
                                     alternative = "two") {
  names = unique(data[, id])
  N = length(names)
  result = data.frame(matrix(ncol = 5, nrow = 0))
  colnames(result) = c("ID", "c_0j", "gamma_hat_j", "v_hat_j", "I_hat_j")
  result_j = data.frame(matrix(nrow = length(C), ncol = 5))
  colnames(result_j) = c("ID", "c_0j", "gamma_hat_j", "v_hat_j", "I_hat_j")
  for (j in 1:N) {
    message(paste("Individual ", j))
    y = data[data[, id] == names[j], response]
    x = data[data[, id] == names[j], running]
    result_j[, 1] = names[j]
    result_j[, 2] = C
    for (i in 1:length(C)) {
      c = C[i]
      h = ifelse(is.null(bw),
                 rdrobust::rdbwselect(y, x, kernel = "uniform", c = c)$bws[1],
                 bw)
      K = function(x) {
        ifelse((abs(x)) <= 1, 1 / 2, 0)
      }
      S0p = sum(ifelse(x - c > 0, K((x - c) / h), 0))
      S1p = sum(ifelse(x - c > 0, (x - c) * K((x - c) / h), 0))
      S2p = sum(ifelse(x - c > 0, (x - c) ^ 2 * K((x - c) / h), 0))
      S0n = sum(ifelse(x - c <= 0, K((x - c) / h), 0))
      S1n = sum(ifelse(x - c <= 0, (x - c) * K((x - c) / h), 0))
      S2n = sum(ifelse(x - c <= 0, (x - c) ^ 2 * K((x - c) / h), 0))
      wp = ifelse(x - c > 0, K((x - c) / h), 0) * (S2p - (x - c) * S1p) / (S2p *
                                                                             S0p - S1p ^ 2)
      wn = ifelse(x - c <= 0, K((x - c) / h), 0) * (S2n - (x - c) * S1n) / (S2n *
                                                                              S0n - S1n ^ 2)
      wp1 = ifelse(x > c, K((x - c) / h), 0) * (S1p - (x - c) * S0p) / (-S2p *
                                                                          S0p + S1p ^ 2)
      wn1 = ifelse(x <= c, K((x - c) / h), 0) * (S1n - (x - c) * S0n) / (-S2n *
                                                                           S0n + S1n ^ 2)
      tauhat = sum((wp - wn) * y)
      result_j[i, 3] = sum((wp1 - wn1) * y)
      fit1 = KernSmooth::locpoly(x,
                     y - ifelse(x > c, 1, 0) * tauhat,
                     bandwidth = h,
                     truncate = FALSE)
      e = y - ifelse((x - c) > 0, 1, 0) * tauhat - approx(fit1$x, fit1$y, xout =
                                                            x)$y
      wpos = ifelse(x - c > 0, ifelse((x - c) < h, 1, 0), 0)
      wneg = ifelse(x - c < 0, ifelse((x - c) > -h, 1, 0), 0)
      result_j[i, 4] = sqrt(var(e[wpos + wneg != 0]) * sum(wp1 ^ 2 + wn1 ^
                                                             2))
    }
    result_j[, 5] = as.numeric(result_j[, 3]) / as.numeric(result_j[, 4])
    if (alternative == "two") {
      result_j_max = result_j[abs(result_j$I_hat_j) == max(abs(result_j$I_hat_j)), ]
    }
    if (alternative == "right") {
      result_j_max = result_j[result_j$I_hat_j == max(result_j$I_hat_j), ]
    }
    if (alternative == "left") {
      result_j_max = result_j[result_j$I_hat_j == min(result_j$I_hat_j), ]
    }
    result = rbind(result, result_j_max)
  }
  I_hat = ifelse(alternative == "left", min(result$I_hat_j), max(abs(result$I_hat_j)))
  if (alternative == "two") {
    p_value = 1 - fdrtool::phalfnorm(I_hat) ^ (N * length(C))
  }
  if (alternative == "right") {
    p_value = 1 - pnorm(I_hat) ^ (N * length(C))
  }
  if (alternative == "left") {
    p_value = 1 - pnorm(-I_hat) ^ (N * length(C))
  }
  alpha = ifelse(is.null(alpha), 1, alpha)
  if (alternative == "two") {
    result_1 = subset(result, abs(result$I_hat_j) > fdrtool::qhalfnorm((1 - alpha) ^ (1 / (N *
                                                        length(C)))))
  }
  else {
    result_1 = subset(result, abs(result$I_hat_j) > qnorm((1 - alpha) ^ (1 / (N *
                                                        length(C)))))
  }
  result_1 = result_1[order(result_1$c_0j), ]
  result_1[, 2:5] = round(result_1[, 2:5], digits = 3)
  rownames(result_1) <- NULL
  if (alternative == "two") {
    critical = c(fdrtool::qhalfnorm((0.90) ^ (1 / (N * length(C)))),
                 fdrtool::qhalfnorm((0.95) ^ (1 / (N * length(C)))),
                 fdrtool::qhalfnorm((0.99) ^ (1 / (N * length(C)))),
                 fdrtool::qhalfnorm((0.999) ^ (1 / (N * length(C)))))
  }
  if (alternative == "right") {
    critical = c(qnorm((0.90) ^ (1 / (N * length(C)))),
                 qnorm((0.95) ^ (1 / (N * length(C)))),
                 qnorm((0.99) ^ (1 / (N * length(C)))),
                 qnorm((0.999) ^ (1 / (N * length(C)))))
  }
  if (alternative == "left") {
    critical = -c(qnorm((0.90) ^ (1 / (N * length(C)))),
                 qnorm((0.95) ^ (1 / (N * length(C)))),
                 qnorm((0.99) ^ (1 / (N * length(C)))),
                 qnorm((0.999) ^ (1 / (N * length(C)))))
  }
  critical = setNames(critical, c("90%", "95%", "99%", "99.9%"))
  list = list(I_hat, p_value, N, critical, result_1)
  list = setNames(list, c("I_hat", "p_value", "N", "Critical_values", "Table"))
  return(list)
}
