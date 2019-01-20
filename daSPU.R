#' Doubley Adaptive Sum of Powered Score tests (daSPU) test for single trait with multiple weights for each genetic marker
#'
#' It returns p-values and test statistics
#'
#' @param U Score vector for the marker set we are interested in. (N by 1 matrix)
#'
#' @param V Corresponding covariance matrix for the score vector. (N by N matrix)
#'
#' @param weight Multiple weights for each genetic markers. (N by M matrix)
#'
#' @param pow1 SNP or single CpG sites specific power(gamma values) used in daSPU test.
#'
#' @param pow2 Specific power(gamma values) used for different weight.
#'
#' @param n.perm number of permutations.
#'
#' @export
#' @return P-values for SPUMpath tests and aSPUMpath test.
#'
#' @author Chong Wu and Wei Pan


daSPU <- function(U,V,weight, pow1 = c(1:8,Inf),pow2 = c(1:8,Inf),n.perm = 1000) {
    weight = as.matrix(weight)
    if(n.perm < 500000) {
        sim.aSPUw2(U,V,weight,pow1,pow2,n.perm)
    } else {
        sim.aSPUw4(U,V,weight,pow1,pow2,n.perm)
    }
}



# sim.aSPUw1, sim.aSPUw2 are for the situations with the number of permutations less than half million
# sim.aSPUw3, sim.aSPUw4 are for the situations with the number of permutations large than half million
sim.aSPUw2 <- function(U,V,weight, pow1 = c(1:8,Inf),pow2 = c(1:8,Inf),n.perm = 1000){
    
    Ts <- rep(0, length(pow1) * length(pow2))
    Ts.tmp <- matrix(0,dim(weight)[2],length(pow1))
    for(i in 1:dim(weight)[2]) {
        for(j in 1:length(pow1)) {
            if (pow1[j] < Inf) {
                a <- sum( (U * weight[,i])^pow1[j] )
                Ts.tmp[i,j] <- sign(a) * ( (abs(a))^(1/pow1[j]) )
            } else {
                Ts.tmp[i,j] <- max(abs(U * weight[,i]))
            }
        }
    }
    
    for(i in 1:length(pow1)) {
        for(j in 1:length(pow2)) {
            if (pow2[j] < Inf) {
                Ts[(j-1) * length(pow1) + i] <- sum(Ts.tmp[,i]^pow2[j])
            } else {
                Ts[(j-1) * length(pow1) + i] <- max(abs(Ts.tmp[,i]))
            }
        }
    }
    
    eV <- eigen(V)
    eV$values[ eV$values < 0 ] = 0

    CovSsqrt <- t(eV$vectors %*% (t(eV$vectors) * sqrt(eV$values)))
    pow1[pow1==Inf] <- 0 # pass 0 as infitiy
    pow2[pow2==Inf] <- 0 # pass 0 as infitiy
    
    T0sC <- calcT0Wsim2(as.matrix(CovSsqrt),as.matrix(weight), as.vector(pow1), as.vector(pow2), n.perm)
    
    # We declare the variables before using it.
    T0s <- T0sC$T0s
    T0s[1:5,]
    pPerm0 <- rep(NA,length(pow1) * length(pow2))
    P0s <- rep(NA,n.perm)
    minp0 <- rep(NA,n.perm)
    
    T0s.abs <- abs(T0s)
    Ts.abs <- abs(Ts)
    
    for ( j in 1:(length(pow1)* length(pow2))) {
        pPerm0[j] <- sum( Ts.abs[j] <= T0s.abs[,j] ) / n.perm
        P0s <- ( ( n.perm - rank( T0s.abs[,j] ) ) + 1 ) / (n.perm)
        if (j == 1 ) {
            minp0 <- P0s
        } else {
            minp0[which(minp0>P0s)] <- P0s[which(minp0>P0s)]
        }
    }
    
    Paspu <- (sum(minp0 <= min(pPerm0)) + 1) / (n.perm + 1)
    pvs <- c(pPerm0, Paspu)
    
    Ts <- c(Ts, min(pPerm0))
    
    if(min(pow1) == 0) {
        pow1[which(pow1 == 0 )] = Inf
    }
    
    if(min(pow2) == 0) {
        pow2[which(pow2 == 0)] = Inf
    }
    
    nmvec <- NULL;
    for(ii in pow2) {
    	   for(jj in pow1) {
               nmvec <- c(nmvec, paste("SPU(",jj,",",ii,")",sep=""))
           }
    }
    
    nmvec <- c(nmvec, "daSPU")
    names(pvs) <- nmvec
    names(Ts) <- nmvec
    
    return(pvs)
}




sim.aSPUw4 <- function(U,V,weight, pow1 = c(1:8,Inf),pow2 = c(1:8,Inf),n.perm = 1000){
    
    Ts <- rep(0, length(pow1) * length(pow2))
    Ts.tmp <- matrix(0,dim(weight)[2],length(pow1))
    for(i in 1:dim(weight)[2]) {
        for(j in 1:length(pow1)) {
            if (pow1[j] < Inf) {
                a <- sum( (U * weight[,i])^pow1[j] )
                Ts.tmp[i,j] <- sign(a) * ( (abs(a))^(1/pow1[j]) )
            } else {
                Ts.tmp[i,j] <- max(abs(U * weight[,i]))
            }
        }
    }
    
    for(i in 1:length(pow1)) {
        for(j in 1:length(pow2)) {
            if (pow2[j] < Inf) {
                Ts[(j-1) * length(pow1) + i] <- sum(Ts.tmp[,i]^pow2[j])
            } else {
                Ts[(j-1) * length(pow1) + i] <- max(abs(Ts.tmp[,i]))
            }
        }
    }
    
    eV <- eigen(V)
    eV$values[ eV$values < 0 ] = 0
    
    CovSsqrt <- t(eV$vectors %*% (t(eV$vectors) * sqrt(eV$values)))
    pow1[pow1==Inf] <- 0 # pass 0 as infitiy
    pow2[pow2==Inf] <- 0 # pass 0 as infitiy

    len.pow = length(pow1) * length(pow2)
    
    T0s <- big.matrix(n.perm,len.pow,init = 0.0,type = "double")
    
    Ts.abs <- abs(Ts)


    pvalue = daSPU_calcT0(as.matrix(CovSsqrt),as.matrix(weight), as.vector(pow1), as.vector(pow2), as.vector(Ts.abs), T0s@address, n.perm)
    
    pvalue = as.vector(pvalue)
    # We declare the variables before using it.
    
    if(min(pow1) == 0) {
        pow1[which(pow1 == 0 )] = Inf
    }
    
    if(min(pow2) == 0) {
        pow2[which(pow2 == 0)] = Inf
    }
    
    nmvec <- NULL;
    for(ii in pow2) {
    	   for(jj in pow1) {
               nmvec <- c(nmvec, paste("SPU(",jj,",",ii,")",sep=""))
           }
    }
    
    nmvec <- c(nmvec, "daSPU")
    names(pvalue) <- nmvec
    
    pvalue
}

