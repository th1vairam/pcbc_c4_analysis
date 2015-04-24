# Modified from http://stackoverflow.com/questions/13088770/how-to-write-linearly-dependent-column-in-a-matrix-in-terms-of-linearly-independ
# Function to find linearly dependednt columns of a matrix
linColumnFinder <- function(mat){
  
  # If the matrix is full rank then we're done
  if(qr(mat)$rank == ncol(mat)){
    print("Matrix is of full rank")
    return(invisible(seq(ncol(mat))))
  }
  
  m <- ncol(mat)
  # cols keeps track of which columns are linearly independent
  cols <- 1
  for(i in seq(2, m)){
    ids <- c(cols, i)
    mymat <- mat[, ids]
    if(qr(mymat)$rank != length(ids)){
      # Regression the column of interest on the previous columns to figure out the relationship
      o <- lm(mat[,i] ~ as.matrix(mat[,cols]) + 0)
      # Construct the output message
      start <- paste0(colnames(mat)[i], " = ")
      # Which coefs are nonzero
      nz <- !(abs(coef(o)) <= .Machine$double.eps^0.5)
      tmp <- colnames(mat)[cols[nz]]
      vals <- paste(coef(o)[nz], tmp, sep = "*", collapse = " + ")
      message <- paste0(start, vals)
      print(message)
    } else {
      # If the matrix subset was of full rank
      # then the newest column in linearly independent
      # so add it to the cols list
      cols <- ids
    }
  }
  return(invisible(cols))
}