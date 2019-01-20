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
aSPUO2 <- function(U,V,weight, pow=c(1:8,Inf),n.perm = 1000){
    
    weight = as.matrix(weight)
    
    # remove SNPs corresponding to zero weight
    weight.tmp = abs(weight)
    index = rowSums(weight.tmp) > 0
    U = as.matrix(U)
    U = U[index,]
    V = V[index,]
    V = V[,index]
    weight = weight[index,]
    
    weight = as.matrix(weight)
    
    npow = length(pow)
    nweight = dim(weight)[2]
    Ts <- rep(0,npow * nweight)

    for(i in 1:nweight) {
        for(j in 1:npow) {
            if (pow[j] < Inf) {
                Ts[(i-1) * npow + j] <- sum((U * weight[,i])^pow[j])
            } else {
                Ts[(i-1) * npow + j] <- max(abs(U * weight[,i]))
            }
        }
    }
    
    eV <- eigen(V)
    eV$values[ eV$values < 0 ] = 0
    
    CovSsqrt <- t(eV$vectors %*% (t(eV$vectors) * sqrt(eV$values)))
    pow[pow==Inf] <- 0 # pass 0 as infitiy
    
    T0s.size = npow * nweight
    T0s <- big.matrix(n.perm,T0s.size,type = "double")
    
    Ts.abs <- abs(Ts)
    
    minp0 <- big.matrix(n.perm,nweight,type = "double")
    minp0_sign <- big.matrix(n.perm,nweight,type = "double")

    RES = daSPU_calcT0Wsim4(as.matrix(CovSsqrt), as.matrix(weight),as.matrix(pow),as.matrix(Ts.abs),T0s@address,minp0@address,minp0_sign@address, n.perm)
    
    pPerm0 = RES$pPerm0
    p.cov = RES$cov
    
    p.cov.inv = solve(p.cov)
    
    df = qr(p.cov)$rank
    
    min.aSPU = rep(0,nweight)
    for( i in 1:nweight) {
        minp.tmp = min(pPerm0[(1+(i-1) *npow):(i*npow)])
        min.aSPU[i] = abs(qnorm(minp.tmp/2)) * sign(Ts[1 +(i-1) * npow])
    }
    
    minp = t(as.matrix(min.aSPU)) %*% p.cov.inv %*% as.matrix(min.aSPU)
    pvalue = 1 - pchisq(minp, df = df)
    
    minp.final = daSPU_calc_aSPU(p.cov.inv,minp0@address,n.perm)
    
    Paspu <- (sum(minp.final >= as.numeric(minp)) ) / (n.perm)
    pvs <- c(as.numeric(pPerm0), Paspu)
    
    Ts <- c(Ts, min(pPerm0))
    pow[pow==0] <- Inf
    
    tmp.name = c(paste("SPU(",1,",",pow,")", sep=""))

    for(i in 2:nweight) {
        tmp.name = c(tmp.name, paste("SPU(",i,",",pow,")", sep=""))
    }
    
    names(Ts) <- c(tmp.name,"aSPU")
    names(pvs) <- names(Ts)
    
    return(list(Ts = Ts, pvs = pvs,ptheory = pvalue))
}
