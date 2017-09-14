argmin_lagrangian <- function(x, lambda) {
  as.vector(1 / (x %*% lambda))
}

constraint_value <- function(x, lambda) {
  new_pi <- argmin_lagrangian(x, lambda)
  apply(x * new_pi, 2, sum)
}
