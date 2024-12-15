# 定义函数
#' @importFrom stats model.matrix
#' @title Weighted GLM estimation for poisson distribution
#' @description
#' A short description...
#'
#' @param x A numeric vector of independent variable (predictor).
#' @param y A numeric vector of dependent variable (response).
#' @param w A numeric vector of weights for each observation.
#'
#'
#' @return A list with two components:
#' \describe{
#'   \item{m}{The number of iterations taken to converge.}
#'   \item{b}{The estimated coefficients of the model (intercept and slope).}
#' }
#' @export
#'
#'
weighted_glm <- function(x, y, w) {
  # 创建数据框
  dt <- data.frame(y = y, x = x)

  # 初始化参数
  X <- model.matrix(~ x, data = dt) # 设计矩阵
  b0 <- rep(0, ncol(X))
  m <- 0

  repeat {
    m <- m + 1
    b <- b0

    # 计算线性预测值 eta 8*1
    eta <- X %*% b

    # 计算均值估计 mu 8*1
    mu <- w * exp(eta)

    # 计算工作变量 z 8*1
    z <- eta + (dt$y - mu) / mu

    # 计算权重矩阵 W
    W <- diag(as.vector(mu) * w)

    # 使用加权最小二乘法更新参数
    b0 <- solve(t(X) %*% W %*% X) %*% (t(X) %*% W %*% z)

    # 检查收敛条件
    if (max(abs(b - b0)) < 1e-8) break
  }

  return(list(m = m, b = b0))
}
