make_D_matrix <- function(N, k=1) {
  D <- matrix(0, N, N)
  j <- 0
  if (N > 1) {
    while(j <= k) {
      rows <- 1:(N-j)
      cols <- rows + j
      indices <- cbind(rows, cols)
      D[indices] <- (-1) ^ (k-j) * choose(k, k-j)
      j <- j + 1
    }
    D <- D[1:(N-k),]
  } else {
    D <- matrix(1)
  }
  return(as(D, 'dgCMatrix'))
}

make_D_matrix_hier <- function(p, N, k=1) {
  D <- make_D_matrix(N, k)
  return(as(as.matrix(Matrix::bdiag(lapply(1:p, function(d){D}))), 'dgCMatrix'))
}

