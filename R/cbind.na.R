#' bind non equal colunm
#' 
#' @usage cbind.na(..., deparse.level = 1) 
#' @param ...  (generalized) vectors or matrices. 
#' @param deparse.level  integer controlling the construction of labels in 
#' the case of non-matrix-like arguments (for the default method):
#' deparse.level = 0 constructs no labels; the default,
#' deparse.level = 1 or 2 constructs labels from the argument names.
#' 
#' 
#' @return a data frame with merged columns
#'
#' @examples
#' \dontrun{
#' col1 <- c("a","b","c","d")
#' col2 <- c("A", "B", "C")
#' col3 <- cbind.na(col1, col2)
#'}
#'
#'@export

cbind.na <- function(..., deparse.level = 1) 
{
    na <- nargs() - (!missing(deparse.level))    
    deparse.level <- as.integer(deparse.level)
    stopifnot(0 <= deparse.level, deparse.level <= 2)
    argl <- list(...)   
    while (na > 0 && is.null(argl[[na]])) {
        argl <- argl[-na]
        na <- na - 1
    }
    if (na == 0) 
        return(NULL)
    if (na == 1) {         
        if (isS4(..1)) 
            return(cbind2(..1))
        else return(matrix(...))  ##.Internal(cbind(deparse.level, ...)))
    }
    if (deparse.level) {       
        symarg <- as.list(sys.call()[-1L])[1L:na]
        Nms <- function(i) {
            if (is.null(r <- names(symarg[i])) || r == "") {
                if (is.symbol(r <- symarg[[i]]) || deparse.level == 
                        2) 
                    deparse(r)
            }
            else r
        }
    }   
    ## deactivated, otherwise no fill in with two arguments
    if (na == 0) {
        r <- argl[[2]]
        fix.na <- FALSE
    }
    else {
        nrs <- unname(lapply(argl, nrow))
        iV <- sapply(nrs, is.null)
        fix.na <- identical(nrs[(na - 1):na], list(NULL, NULL))
        ## deactivated, otherwise data will be recycled
        #if (fix.na) {
        #    nr <- max(if (all(iV)) sapply(argl, length) else unlist(nrs[!iV]))
        #    argl[[na]] <- cbind(rep(argl[[na]], length.out = nr), 
        #        deparse.level = 0)
        #}       
        if (deparse.level) {
            if (fix.na) 
                fix.na <- !is.null(Nna <- Nms(na))
            if (!is.null(nmi <- names(argl))) 
                iV <- iV & (nmi == "")
            ii <- if (fix.na) 
                2:(na - 1)
            else 2:na
            if (any(iV[ii])) {
                for (i in ii[iV[ii]]) if (!is.null(nmi <- Nms(i))) 
                    names(argl)[i] <- nmi
            }
        }
        
        ## filling with NA's to maximum occuring nrows
        nRow <- as.numeric(sapply(argl, function(x) NROW(x)))
        maxRow <- max(nRow, na.rm = TRUE)  
        argl <- lapply(argl, function(x)  if (is.null(nrow(x))) c(x, rep(NA, maxRow - length(x)))
                       else rbind.na(x, matrix(, maxRow - nrow(x), ncol(x))))
        r <- do.call(cbind, c(argl[-1L], list(deparse.level = deparse.level)))
    }
    d2 <- dim(r)
    r <- cbind2(argl[[1]], r)
    if (deparse.level == 0) 
        return(r)
    ism1 <- !is.null(d1 <- dim(..1)) && length(d1) == 2L
    ism2 <- !is.null(d2) && length(d2) == 2L && !fix.na
    if (ism1 && ism2) 
        return(r)
    Ncol <- function(x) {
        d <- dim(x)
        if (length(d) == 2L) 
            d[2L]
        else as.integer(length(x) > 0L)
    }
    nn1 <- !is.null(N1 <- if ((l1 <- Ncol(..1)) && !ism1) Nms(1))
    nn2 <- !is.null(N2 <- if (na == 2 && Ncol(..2) && !ism2) Nms(2))
    if (nn1 || nn2 || fix.na) {
        if (is.null(colnames(r))) 
            colnames(r) <- rep.int("", ncol(r))
        setN <- function(i, nams) colnames(r)[i] <- if (is.null(nams)) 
            ""
        else nams
        if (nn1) 
            setN(1, N1)
        if (nn2) 
            setN(1 + l1, N2)
        if (fix.na) 
            setN(ncol(r), Nna)
    }
    r
}
