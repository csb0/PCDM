nhpp = function(cc, w, cyclic = TRUE) {
  n = length(cc)
  if (cyclic == TRUE) {
    m = n * w
    Q = diag(rep(12, m))
    Q[cbind(2:m, 1:(m - 1))] = Q[cbind(1:(m - 1), 2:m)]
                             = Q[m, 1] = Q[1, m] = -8
    Q[cbind(3:m, 1:(m - 2))] = Q[cbind(1:(m - 2), 3:m)] = Q[m - 1, 1]
                             = Q[m, 2] = Q[1, m - 1] = Q[2, m] = 2
  } else {
    m = n * w + 1
    Q = diag(rep(12, m))
    Q[cbind(2:m, 1:(m - 1))] = Q[cbind(1:(m - 1), 2:m)] = -8
    Q[cbind(3:m, 1:(m - 2))] = Q[cbind(1:(m - 2), 3:m)] = 2
    Q[1, 1] = Q[m, m] = 2
    Q[2, 2] = Q[m - 1, m - 1] = 10
    Q[1, 2] = Q[2, 1] = Q[m, m - 1] = Q[m - 1, m] = -4
  }
  Q = Q*w^2
  if (cyclic == TRUE) {
    E = matrix(0, n, m + 1)
    rowpat = c(1, rep(2, w - 1), 1)
    for (i in 1:n) {
      j = (i - 1) * w + 1
      E[i, j:(j + w)] = rowpat
    }
    E = E[ , -(m + 1)]
    E[n, 1] = 1
  } else {
    E = matrix(0, n, m)
    rowpat = c(1, rep(2, w - 1), 1)
    for (i in 1:n) {
      j = (i - 1) * w + 1
      E[i, j:(j + w)] = rowpat
    }
  }
  A = cbind(Q, t(E))
  B = cbind(E, matrix(0, n, n))
  C = rbind(A, B)
  rhs = c(rep(0, m), 2 * w * cc)
  yg = solve(C, rhs)
  return(yg[1:m])
}
