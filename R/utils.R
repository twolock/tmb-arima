make_TMB_map <- function(par, x=NULL, append_to=NULL) {
  if (!is.null(append_to)) {
    par_names <- x[!(x %in% names(append_to))]
  } else {
    par_names <- x
  }
  tmp_l <- lapply(par_names, function(xx) {
    factor(rep(NA, length(par[[xx]])))
  })
  names(tmp_l) <- par_names
  return(c(append_to, tmp_l))
}

organize_means <- function(sd.obj) {
  if (is.null(sd.obj$jointPrecision)) {
    unique.names <- unique(names(sd_obj$par.fixed))
  } else {
    unique.names <- unique(colnames(sd.obj$jointPrecision))
  }
  all.means <- c()
  for (n in unique.names) {
    par.type <- ifelse(n %in% names(sd.obj$par.fixed), 'par.fixed', 'par.random')
    all.means <- c(all.means, sd.obj[[par.type]][names(sd.obj[[par.type]]) == n])
  }
  return(all.means)
}

scale_gmrf_precision <- function(Q,
                                 A = matrix(1, ncol = ncol(Q)),
                                 eps = sqrt(.Machine$double.eps)) {
  
  nb <- spdep::mat2listw(abs(Q))$neighbours
  comp <- spdep::n.comp.nb(nb)
  
  for (k in seq_len(comp$nc)) {
    idx <- which(comp$comp.id == k)
    Qc <- Q[idx, idx, drop = FALSE]
    
    if (length(idx) == 1) {
      ## set marginal variance for islands = 1
      Qc[1, 1] <- 1
    } else {
      Ac <- A[ , idx, drop = FALSE]
      Qc_eps <- Qc + Matrix::Diagonal(ncol(Qc)) * max(Matrix::diag(Qc)) * eps
      Qc_inv <- qinv(Qc_eps, A = Ac)
      scaling_factor <- exp(mean(log(Matrix::diag(Qc_inv))))
      Qc <- scaling_factor * Qc
    }
    
    Q[idx, idx] <- Qc
  }
  
  Q
}


qinv <- function(Q, A = NULL) {
  Sigma <- Matrix::solve(Q)
  if (is.null(A))
    return(Sigma)
  else {
    A <- matrix(1,1, nrow(Sigma))
    W <- Sigma %*% t(A)
    Sigma_const <- Sigma - W %*% Matrix::solve(A %*% W) %*% Matrix::t(W)
    return(Sigma_const)
  }
}