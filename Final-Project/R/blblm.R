#' @import purrr
#' @import future
#' @import furrr
#' @import stats
#' @importFrom magrittr %>%
#' @importFrom utils capture.output
#' @aliases NULL
#' @details
#' Linear Regression with Little Bag of Bootstraps
"_PACKAGE"

## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("."))

#' main function to run linear regression with mini bag bootstrap
#'
#' @param formula an object with class "formula": a symbolic description of the model to be fitted.
#' @param data a data frame or a vector of filenames
#' @param m the dataset is split into m files.
#' @param B bootstrap iteration times
#' @param cl number of cores we use in parallel computing
#' @param type type of the model, either lm for linear regression or glm for logisic regression
#'
#' @return blblm object
#'
#' @export
#'
blblm <- function(formula, data, m = 10, B = 5000, cl = 1, type = "lm") {
  if (is.data.frame(data)) {
      data_list <- split_data(data, m)
      if (cl == 1) {
        if (type == "lm"){
        estimates <- map(data_list,
                         ~ lm_each_subsample(
                           formula = formula,
                           data = .,
                           n = nrow(data),
                           B = B
                         ))
        } else {
          estimates <- map(data_list,
                           ~ glm_each_subsample(
                             formula = formula,
                             data = .,
                             n = nrow(data),
                             B = B
                           ))
        }
      } else{
      suppressWarnings(plan(multiprocess, workers = cl))
      estimates <- future_map(data_list,
                              ~ lm_each_subsample(
                                formula = formula,
                                data = .,
                                n = nrow(data),
                                B = B
                              ))
    }
  } else{
    if (cl > 1) {
      suppressWarnings(plan(multiprocess, workers = cl))
      N = data %>%
        map(~ {
          df <- read.csv(.,)
          nrow(df)
        }) %>% reduce(`+`)
      estimates = data %>%
        future_map(~ {
          df <- read.csv(.,)
          lm_each_subsample(
            formula = formula,
            data = df,
            n = N,
            B = B
          )
        })
      } else {
        N = data %>%
          map( ~ {
            df <- read.csv(., )
            nrow(df)
          }) %>% reduce(`+`)
        estimates = data %>%
          map( ~ {
            df <- read.csv(., )
            lm_each_subsample(
              formula = formula,
              data = df,
              n = N,
              B = B
            )
          })
      }
  }
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}

#' This function take a data frame and output a list of data frame.
#' The list has length m. Each element in the list has the same dimension
#' as the original data but they could be repeated because they are
#' formed through bootstrap resampling.
#'
#' @param data a data frame
#' @param m the dataset is split into m files.
#'
#' @export
split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}

#' replicate linear regression fitting
#'
#' @param formula a formula class object
#' @param data data frame for analysis
#' @param n sample size
#' @param B bootstrap iteration times
#' @export
lm_each_subsample <- function(formula, data, n, B) {
  replicate(B, lm_each_boot(formula, data, n), simplify = FALSE)
}

#' replicate glm fitting
#'
#' @param formula a formula class object
#' @param data data frame for analysis
#' @param n sample size
#' @param B bootstrap iteration times
#' @export
glm_each_subsample <- function(formula, data, n, B) {
  replicate(B, glm_each_boot(formula, data, n), simplify = FALSE)
}

#' linear regression on one bootstrap
#'
#' @param formula an object of class "formula"
#' @param data data frame
#' @param n total sample size
#'
#' @export
lm_each_boot <- function(formula, data, n) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  lm1(formula, data, freqs)
}

#' logistic regression on one bootstrap
#'
#' @param formula an object of class "formula"
#' @param data data frame
#' @param n total sample size
#'
#' @export
glm_each_boot <- function(formula, data, n) {
  environment(formula) <- environment()
  indx <- rep(1:nrow(data), size = n, replace = TRUE)
  newdata <- data[indx,]
  fit <- suppressWarnings(glm(formula, newdata, family = "binomial"))
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}

#' fit regression model
#'
#' @param formula an object of class "formula"
#' @param data an data frame
#' @param freqs freqency of each observations
#'
#' @export
lm1 <- function(formula, data, freqs) {
  environment(formula) <- environment()
  fit <- lm(formula, data, weights = freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}

#' obtain coefficients
#'
#' @param fit fitted regression
#'
#' @export
blbcoef <- function(fit) {
  coef(fit)
}

#' obtain blblm sigma
#'
#' @param fit fitted regression
#'
#' @export
blbsigma <- function(fit) {
  p <- fit$rank
  y <- model.extract(fit$model, "response")
  e <- fitted(fit) - y
  w <- fit$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}

#' print for blblm class
#'
#' @param x blblm object
#' @param ... dditional arguments to be passed to the low level regression fitting functions
#'
#' @export
#' @method print blblm
print.blblm <- function(x, ...) {
  cat("blblm model:", capture.output(x$formula))
  cat("\n")
}

##################

#' obtain sigma
#'
#' @param object an object of class "blblm"
#' @param confidence need confidence interval or not
#' @param level the confidence level required.
#' @param ... dditional arguments to be passed to the low level regression fitting functions
#'
#' @export
#' @method sigma blblm
sigma.blblm <- function(object, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha <- 1 - 0.95
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
}

#' obtain coefficients
#'
#' @param object blblm object
#' @param ... dditional arguments to be passed to the low level regression fitting functions
#'
#' @export
#' @method coef blblm
coef.blblm <- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}

#' obtain confidence interval
#'
#' @param object an object with class "blblm
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level user specified confidence level, with default value 95%
#' @param ... dditional arguments to be passed to the low level regression fitting functions
#'
#' @export
#' @method confint blblm
confint.blblm <- function(object, parm = NULL, level = 0.95, ...) {
  if (is.null(parm)) {
    parm <- attr(terms(object$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates
  out <- map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
  })
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}

#' prediciton for blbm class
#'
#' @param object an object of class "blblm"
#' @param new_data new data frame for prejection
#' @param confidence draw CI or not
#' @param level user specified confidence level, with default value 95%
#' @param ... additional arguments to be passed to the low level regression fitting functions
#'
#' @export
#' @method predict blblm
predict.blblm <- function(object, new_data, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms(object$formula), "term.labels")), new_data)
  if (confidence) {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>%
      apply(1, mean_lwr_upr, level = level) %>%
      t())
  } else {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans())
  }
}

mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}

map_mean <- function(.x, .f, ...) {
  (map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}

map_cbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(cbind)
}

map_rbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(rbind)
}
