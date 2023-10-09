
checkUni <- function (w, delta, x) {
  x <- drop (x)
  n <- length (x)
  vinv <- solve (smacofVmat (w) + (1 / n)) - (1 / n)
  return (drop (vinv %*% rowSums (w * delta * sign (outer (
    x, x, "-"
  )))))
}

matchMe <- function (x,
                     itmax = 100,
                     eps = 1e-10,
                     verbose = FALSE) {
  m <- length (x)
  y <- sumList (x) / m
  itel <- 1
  fold <- sum (sapply (x, function (z)
    (z - y) ^ 2))
  repeat {
    for (j in 1:m) {
      u <- crossprod (x[[j]], y)
      s <- svd (u)
      r <- tcrossprod (s$u, s$v)
      x[[j]] <- x[[j]] %*% r
    }
    y <- sumList (x) / m
    fnew <- sum (sapply (x, function (z)
      (z - y) ^ 2))
    if (verbose) {
      
    }
    if (((fold - fnew) < eps) || (itel == itmax))
      break
    itel <- itel + 1
    fold <- fnew
  }
  return (x)
}

sumList <- function (x) {
  m <- length (x)
  y <- x[[1]]
  for (j in 2:m) {
    y <- y + x[[j]]
  }
  return (y)
}


smacofLoss <- function (d, w, delta) {
  return (sum (w * (delta - d) ^ 2) / 4)
}

smacofBmat <- function (d, w, delta) {
  dd <- ifelse (d == 0, 0, 1 / d)
  b <- -dd * w * delta
  diag (b) <- -rowSums (b)
  return(b)
}

smacofVmat <- function (w) {
  v <- -w
  diag(v) <- -rowSums(v)
  return (v)
}

smacofGuttman <- function (x, b, vinv) {
  return (vinv %*% b %*% x)
}

columnCenter <- function (x) {
  return (apply (x, 2, function (z)
    z - mean (z)))
}

smacofComplement <- function (y, v) {
  return (sum (v * tcrossprod (y)) / 4)
}

smacofPenalty <-
  function (w,
            delta,
            p = 2,
            lbd = 0,
            zold = columnCenter (diag (nrow (delta))),
            itmax = 10000,
            eps = 1e-10,
            verbose = FALSE) {
    itel <- 1
    n <- nrow (zold)
    vmat <- smacofVmat (w)
    vinv <- solve (vmat + (1 / n)) - (1 / n)
    dold <- as.matrix (dist (zold))
    mold <- sum (w * delta * dold) / sum (w * dold * dold)
    zold <- zold * mold
    dold <- dold * mold
    yold <- zold [, (p + 1):n]
    sold <- smacofLoss (dold, w, delta)
    bold <- smacofBmat (dold, w, delta)
    told <- smacofComplement (yold, vmat)
    uold <- sold + lbd * told
    repeat {
      znew <- smacofGuttman (zold, bold, vinv)
      ynew <- znew [, (p + 1):n] / (1 + lbd)
      znew [, (p + 1):n] <- ynew
      xnew <- znew [, 1:p]
      dnew <- as.matrix (dist (znew))
      bnew <- smacofBmat (dnew, w, delta)
      tnew <- smacofComplement (ynew, vmat)
      snew <- smacofLoss (dnew, w, delta)
      unew <- snew + lbd * tnew
      if (verbose) {
        cat(
          "itel ",
          formatC(itel, width = 4, format = "d"),
          "sold ",
          formatC(
            sold,
            width = 10,
            digits = 6,
            format = "f"
          ),
          "snew ",
          formatC(
            snew,
            width = 10,
            digits = 6,
            format = "f"
          ),
          "told ",
          formatC(
            told,
            width = 10,
            digits = 6,
            format = "f"
          ),
          "tnew ",
          formatC(
            tnew,
            width = 10,
            digits = 6,
            format = "f"
          ),
          "uold ",
          formatC(
            uold,
            width = 10,
            digits = 6,
            format = "f"
          ),
          "unew ",
          formatC(
            unew,
            width = 10,
            digits = 6,
            format = "f"
          ),
          "\n"
        )
      }
      if (((uold - unew) < eps) || (itel == itmax)) {
        break
      }
      itel <- itel + 1
      zold <- znew
      bold <- bnew
      sold <- snew
      told <- tnew
      uold <- unew
    }
    zpri <- znew %*% svd(znew)$v
    xpri <- zpri[, 1:p]
    return (list (
      x = xpri,
      z = zpri,
      b = bnew,
      l = lbd,
      s = snew,
      t = tnew,
      itel = itel
    ))
  }

plotMe2 <- function(hList, labels, s = 1, t = 2) {
  n <- nrow(hList[[1]]$x)
  m <- length (hList)
  par(pty = "s")
  hMatch <- matchMe (lapply (hList, function(r)
    r$x))
  hMat <- matrix (0, 0, 2)
  for (j in 1:m) {
    hMat <- rbind(hMat, hMatch[[j]][, c(s, t)])
  }
  plot(
    hMat,
    xlab = "dim 1",
    ylab = "dim 2",
    col = c(rep("RED", n * (m - 1)), rep("BLUE", n)),
    cex = c(rep(1, n * (m - 1)), rep(2, n))
  )
  for (i in 1:n) {
    hLine <- matrix (0, 0, 2)
    for (j in 1:m) {
      hLine <- rbind (hLine, hMatch[[j]][i, c(s, t)])
    }
    lines(hLine)
  }
  text(hMatch[[m]], labels, cex = .75)
}

plotMe1 <- function(hList, labels) {
  n <- length (hList[[1]]$x)
  m <- length (hList)
  blow <- function (x) {
    n <- length (x)
    return (matrix (c(1:n, x), n, 2))
  }
  hMat <- matrix (0, 0, 2)
  for (j in 1:m) {
    hMat <- rbind(hMat, blow(hList[[j]]$x))
  }
  plot(
    hMat,
    xlab = "index",
    ylab = "x",
    col = c(rep("RED", n * (m - 1)), rep("BLUE", n)),
    cex = c(rep(1, n * (m - 1)), rep(2, n))
  )
  for (i in 1:n) {
    hLine <- matrix (0, 0, 2)
    for (j in 1:m) {
      hLine <- rbind (hLine, blow(hList[[j]]$x)[i,])
      lines(hLine)
    }
  }
  text(blow(hList[[m]]$x), labels, cex = 1.00)
  for (i in 1:n) {
    abline(h = hList[[m]]$x[i])
  }
}




runPenalty <-
  function (w,
            delta,
            lbd,
            p = 2,
            itmax = 10000,
            eps = 1e-10,
            cut = 1e-6,
            write = TRUE,
            verbose = FALSE) {
    m <- length (lbd)
    hList <- as.list (1:m)
    hList[[1]] <-
      smacofPenalty(
        w,
        delta,
        p,
        lbd = lbd[1],
        itmax = itmax,
        eps = eps,
        verbose = verbose
      )
    for (j in 2:m) {
      hList[[j]] <-
        smacofPenalty(
          w,
          delta,
          p,
          zold = hList[[j - 1]]$z,
          lbd = lbd[j],
          itmax = itmax,
          eps = eps,
          verbose = verbose
        )
    }
    mm <- m
    for (i in 1:m) {
      if (write) {
        cat(
          "itel",
          formatC(hList[[i]]$itel, width = 4, format = "d"),
          "lambda",
          formatC(
            hList[[i]]$l,
            width = 10,
            digits = 6,
            format = "f"
          ),
          "stress",
          formatC(
            hList[[i]]$s,
            width = 8,
            digits = 6,
            format = "f"
          ),
          "penalty",
          formatC(
            hList[[i]]$t,
            width = 8,
            digits = 6,
            format = "f"
          ),
          "\n"
        )
      }
      if (hList[[i]]$t < cut) {
        mm <- i
        break
      }
    }
    return(hList[1:mm])
  }

writeSelected <- function (hList, ind) {
  m <- length (hList)
  n <- length (ind)
  mn <- sort (union (union (1:3, ind), m - (2:0)))
  for (i in mn) {
    if (i > m) {
      next
    }
    cat(
      "itel",
      formatC(hList[[i]]$itel, width = 4, format = "d"),
      "lambda",
      formatC(
        hList[[i]]$l,
        width = 10,
        digits = 6,
        format = "f"
      ),
      "stress",
      formatC(
        hList[[i]]$s,
        width = 8,
        digits = 6,
        format = "f"
      ),
      "penalty",
      formatC(
        hList[[i]]$t,
        width = 8,
        digits = 6,
        format = "f"
      ),
      "\n"
    )
  }
}

rhofun <- function (a) {
  rhomax <- 0
  rhoval <- c(0, 0)
  n <- nrow (a)
  for (i in 1:n) {
    z <- a[i, 1] * xbase + a[i, 2] * ybase
    rho <- sum (delta * dist (z))
    if (rho > rhomax) {
      rhomax <- rho
      rhoval <- a[i,]
    }
  }
  return(list(rhomax = rhomax, rhoval = rhoval))
}

rhomax2Plot <- function (inpoints) {
  par(pty = "s")
  s <- (0:500) / 500 * 2 * pi
  x <- sin(s)
  y <- cos(s)
  plot(
    x,
    y,
    type = "l",
    col = "RED",
    xlim = c(-1.25, 1.25),
    ylim = c(-1.25, 1.25),
    lwd = 3
  )
  points(0, 0, cex = 1.2)
  n <- nrow(inpoints)
  outpoints <- matrix(0, n, 2)
  for (i in 1:n) {
    a <- inpoints[i, ]
    if (i == n) {
      b <- inpoints[1, ]
    } else {
      b <- inpoints[i + 1, ]
    }
    d <- a[1] * b[2] - a[2] * b[1]
    outpoints[i, 1] <- (b[2] - a[2]) / d
    outpoints[i, 2] <- (a[1] - b[1]) / d
  }
  lines (
    x = inpoints[, 1],
    y = inpoints[, 2],
    col = "BLUE",
    lwd = 2
  )
  lines (
    x = inpoints[c(1, n), 1],
    y = inpoints[c(1, n), 2],
    col = "BLUE",
    lwd = 2
  )
  lines (
    x = outpoints[, 1],
    y = outpoints[, 2],
    col = "BLUE",
    lwd = 2
  )
  lines (
    x = outpoints[c(1, n), 1],
    y = outpoints[c(1, n), 2],
    col = "BLUE",
    lwd = 2
  )
}


rhomax2Comp <-
  function (inpoints,
            itmax = 100,
            eps = 1e-8,
            verbose = TRUE) {
    itel = 1
    repeat {
      n <- nrow(inpoints)
      outpoints <- matrix(0, n, 2)
      for (i in 1:n) {
        a <- inpoints[i, ]
        if (i == n) {
          b <- inpoints[1, ]
        } else {
          b <- inpoints[i + 1, ]
        }
        d <- a[1] * b[2] - a[2] * b[1]
        outpoints[i, 1] <- (b[2] - a[2]) / d
        outpoints[i, 2] <- (a[1] - b[1]) / d
      }
      infun <- rhofun (inpoints)
      oufun <- rhofun (outpoints)
      inmax <- infun$rhomax
      oumax <- oufun$rhomax
      inval <- infun$rhoval
      ouval <- oufun$rhoval
      for (i in 1:n) {
        outpoints[i, ] <- outpoints[i, ] / sqrt (sum (outpoints[i, ] ^ 2))
      }
      impoints <- matrix (0, 2 * n, 2)
      impoints[seq.int(1, (2 * n) - 1, by = 2),] <- inpoints
      impoints[seq.int(2, 2 * n, by = 2),] <- outpoints
      if (verbose) {
        cat(
          "itel ",
          formatC(itel, width = 6, format = "d"),
          "vertices ",
          formatC(n, width = 6, format = "d"),
          "innermax ",
          formatC(
            inmax,
            digits = 8,
            width = 15,
            format = "f"
          ),
          "outermax ",
          formatC(
            oumax,
            digits = 8,
            width = 15,
            format = "f"
          ),
          "\n"
        )
      }
      if ((itel == itmax) || ((oumax - inmax) < eps)) {
        break
      }
      itel <- itel + 1
      inpoints <- impoints
    }
    return (list (
      itel = itel,
      vertices = n,
      inmax = inmax,
      oumax = oumax,
      inval = inval,
      ouval = ouval
    ))
  }
