# Find the locally optimum resampling weights
# given Lagrange multiplier lambda and data x.
argmin_lagrangian <- function(x, lambda) {
  as.vector(1 / (x %*% lambda))
}

# For the real optimization.
objective_function <- function(x, lambda, target) {
  new_q <- argmin_lagrangian(x, lambda)
  constraint_value <- apply(x * new_q, 2, sum)
  squares <- (constraint_value - target) ^ 2
  sum(squares)
}

# Start the optimization
# at a uniform resampling distribution.
# starting_values <- function(x){
#   y <- rep(nrow(x), nrow(x))
#   fit <- lm(as.formula("y ~ x - 1"))
#   fit$coefficients
# }
