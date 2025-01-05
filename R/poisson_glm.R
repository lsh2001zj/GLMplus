# 定义函数
#' @importFrom stats glm
#' @importFrom stats rpois
#' @importFrom stats predict
#' @importFrom stats as.formula
#' @importFrom stats poisson
#' @title Weighted GLM estimation for experiment
#' @description
#' A short description...
#'
#' @param data A dataframe
#' @param dep_var_name response variable
#' @param explanatory_vars explanatory variables
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
poisson_glm <- function(data, dep_var_name, explanatory_vars) {
  dep_var <- data[[dep_var_name]]
  indep_vars <- data[, explanatory_vars, drop = FALSE]
  for (col in names(indep_vars)) {
    if (is.factor(indep_vars[[col]])) {
      indep_vars <- cbind(indep_vars, model.matrix(~ indep_vars[[col]] - 1))
      indep_vars[[col]] <- NULL  # 删除原始的分类变量
    }
  }
  X <- cbind(1, indep_vars)
  formula <- as.formula(paste(dep_var_name, "~ ."))
  offset_var <- "Exposure"
  if (!is.null(offset_var)) {
    offset_term <- log(data[[offset_var]])
  } else {
    offset_term <- NULL
  }
  model <- glm(formula, data = data, family = poisson(), offset = offset_term)
  result <- summary(model)
  coefficients <- result$coefficients
  output <- data.frame(Variable = rownames(coefficients),
                       Estimate = coefficients[, 1],
                       StdError = coefficients[, 2],
                       zValue = coefficients[, 3],
                       PrValue = coefficients[, 4])
  fitted_values <- predict(model, newdata = data, type = "response")
  data$fitGLM1 <- fitted_values
  output_test <- data.frame(fitGLM1 = data$fitGLM1,
                            ClaimNb = data[[dep_var_name]],
                            Exposure = data$Exposure)
  return(list(coefficients = output, test_data = output_test))
}
