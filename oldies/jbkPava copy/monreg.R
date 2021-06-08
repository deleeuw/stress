monreg <-
  function(x,
           w = rep(1, length(x)),
           block = weighted.mean) {
    is.up.satisfied <- function(x, i)
      (i == length(x)) || (x[i] <= x[i + 1])
    is.down.satisfied <- function(x, i)
      (i == 1) || (x[i - 1] <= x[i])
    put.back <- function(n, blocklist, blockvalues) {
      x <- rep(0, n)
      nb <- length(blockvalues)
      for (i in 1:nb) {
        x[blocklist[i, 1]:blocklist[i, 2]] <- blockvalues[i]
      }
      return(x)
    }
    merge.block.up <-
      function(blocklist, blockvalues, x, w, i, block) {
        n <- length(blockvalues)
        nn <- 1:n
        ii <- which(i + 1 != nn)
        blocklist[i,] <- c(blocklist[i, 1], blocklist[i + 1, 2])
        indi <- blocklist[i, 1]:blocklist[i + 1, 2]
        blockvalues[i] <- block(x[indi], w[indi])
        blocklist <- blocklist[ii,]
        if (length(ii) == 1)
          dim(blocklist) <- c(1, 2)
        blockvalues <- blockvalues[ii]
        list(v = blockvalues, l = blocklist)
      }
    nblock <- length (x)
    n <- length (x)
    blocklist <- array(1:n, c(n, 2))
    blockvalues <- x
    active <- 1
    repeat {
      if (!is.up.satisfied(blockvalues, active)) {
        blockmerge <-
          merge.block.up(blocklist, blockvalues, x, w, active, block)
        blockvalues <- blockmerge$v
        blocklist <- blockmerge$l
        nblock <- nblock - 1
        while (!is.down.satisfied(blockvalues, active)) {
          blockmerge <-
            merge.block.up(blocklist, blockvalues, x, w, active - 1, block)
          blockvalues <- blockmerge$v
          blocklist <- blockmerge$l
          nblock <- nblock - 1
          active <- active - 1

        }
      }
      else if (active == nblock)
        break()
      else
        active <- active + 1
    }
    put.back(n, blocklist, blockvalues)
  }

monregP <-
  function(x,
           y,
           w = rep(1, length(x)),
           block = weighted.mean) {
    o <- order(x, y)
    r <- order(o)
    return(pava(y[o], w[o])[r])
  }

monregS <-
  function(x,
           y,
           w = rep(1, length(x)),
           block = weighted.mean) {
    wag <- tapply(w, x, sum)
    yag <- tapply(y, x, mean)
    xag <- tapply(x, x, mean)
    o <- order(xag)
    r <- order(o)
    e <- pava(yag[o], wag[o])[r]
    return(ifelse(outer(x, xag, "=="), 1, 0) %*% e)
  }

monregT <-
  function(x,
           y,
           w = rep(1, length(x)),
           block = weighted.mean) {
    wag <- tapply(w, x, sum)
    yag <- tapply(y, x, mean)
    xag <- tapply(x, x, mean)
    o <- order(xag)
    r <- order(o)
    e <- pava(yag[o], wag[o])[r]
    return(y + ifelse(outer(x, xag, "=="), 1, 0) %*% (e - yag[o]))
  }
