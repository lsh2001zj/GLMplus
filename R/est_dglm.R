#' @title an example for dglm
#'
#' @param gamma gamma
#' @param beta beta
#' @param n n
#'
#' @return two estimations
#' @export
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom statmod rinvgauss
est_dglm <- function(gamma, beta, n) {
  # 随机生成自变量 x1 和 x2
  x1 <- rnorm(n, mean = 0, sd = 1)
  x2 <- runif(n, min = 0, max = 1)
  # 计算均值参数 mu 和离散参数 phi
  mu <- exp(beta[1] + beta[2]*x1 + beta[3]*x2 + beta[4]*x1*x2)  # GLM1 模型中的期望
  phi <- exp(gamma[1] + gamma[2]*x1 + gamma[3]*x2)  # GLM2 模型中的方差
  # 生成响应变量 Y
  Y <- mapply(function(mu_val, phi_val) {
    rinvgauss(1, mean = mu_val,dispersion = phi_val)
  }, mu, phi)
  X_mu <- cbind(1, x1, x2, x1*x2)  # 设计矩阵
  X_phi <- cbind(1, x1, x2)
  m <- 0                 # 迭代计数
  j <- 0                 # 迭代计数
  # 初始值
  mu0 <- mean(Y)
  z0 <-  log(mu0) + (Y - mu0) / mu0
  W0 <- diag(c(rep(1/mu0,n)))
  b0 <- solve(t(X_mu) %*% W0 %*% X_mu) %*% (t(X_mu) %*% W0 %*% z0)
  repeat {
    m <- m + 1
    b <- b0
    eta <- X_mu %*% b       # 计算线性预测值
    mu <- exp(eta)       # 计算均值估计
    z <- eta + (Y - mu) / mu  # 计算工作变量
    W <- diag(c(1/mu))
    b0 <- solve(t(X_mu) %*% W %*% X_mu) %*% (t(X_mu) %*% W %*% z)  # 更新参数
    if (max(abs(b - b0)) < 1e-8) break  # 检查收敛条件
  }
  d0 <- (Y - mu)^2 / (mu^2 * Y)
  phi0 <- mean(d0)
  z0 <-  log(phi0) + (d0 - phi0) / phi0
  W0 <- diag(c(rep(1/phi0,n)))
  g0 <- solve(t(X_phi) %*% W0 %*% X_phi) %*% (t(X_phi) %*% W0 %*% z0)

  repeat {
    j <- j + 1
    b <- b0
    # 更新 GLM2 的参数 γ
    d <- (Y - mu)^2 / (mu^2 * Y)  # 残差平方
    repeat {
      g <- g0
      eta2 <- X_phi %*% g
      phi <- exp(eta2)  # 计算方差估计
      z2 <- eta2 + (d - phi) / phi  # 伽马回归的工作变量
      W2 <- diag(c(rep(1 / 2, n)))  # 权重矩阵
      g0 <- solve(t(X_phi) %*% W2 %*% X_phi) %*% (t(X_phi) %*% W2 %*% z2)  # 更新γ
      if (max(abs(g - g0)) < 1e-10) break
    }
    eta <- X_mu %*% b
    mu <- exp(eta)
    z <- eta + (Y - mu) / mu
    W <- diag(1 / c(exp(X_phi %*% g0) * mu))
    b0 <- solve(t(X_mu) %*% W %*% X_mu) %*% (t(X_mu) %*% W %*% z)
    if (max(abs(b - b0)) < 1e-10) break
  }
  return(list(params1 = b0, params2 = g0))
}
