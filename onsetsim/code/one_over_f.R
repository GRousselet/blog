# Edit `one_over_f` function from `primer` package to control variance (Stevens, 2009). 
# Original function is available on [GitHub](https://github.com/HankStevens/primer).
# Copyright Hank Stevens.

one_over_f <- function(gamma = 1, N = 200, outvar = 1){
  N.2 <- N/2
  sine.waves <- matrix(NA, nrow = N, ncol = N.2)
  steps = 2 * pi * (1:N)/N
  phase <- stats::runif(N.2, 0, 2 * pi)
  for (i in 1:N.2) {
    freq <- i
    weight <- 1/(freq^gamma)
    y <- weight * sin(freq * steps + phase[i])
    sine.waves[, i] <- y
  }
  out <- rowSums(sine.waves)
  # force var = outvar and mean = 0
  out <- out - mean(out);
  out <- out / sd(out);
  out <- out * sqrt(outvar);
  return(out)
}