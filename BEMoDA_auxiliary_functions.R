# BEMoDA v1.0 - BiowaivEr aid for Model Dependent-Independent Approach script for in-vitro dissolution profile comparison
# 
# Model Dependent Approach script for in-vitro dissolution profile comparison as proposed by Sathe et al. in 1996
# (Sathe PM, Tsong Y, Shah VP. In-vitro dissolution profile comparison: statistics and analysis, model dependent approach. Pharm Res. 1996 Dec;13(12):1799-803).
# 
# Copyright (C) 2017 Jakub Szlęk, Aleksander Mendyk
# 
# Authors: 
# Jakub Szlęk, Aleksander Mendyk
# 
# Affiliation: 
# Jagiellonian University Medical College,
# Faculty of Pharmacy,
# Department of Pharmaceucial Technology and Biopharmaceutics,
# Medyczna 9 st.,
# 30-688 Kraków
# Poland
# 
# Bugs, issues, please e-mail to maintainer
# Jakub Szlęk: j.szlek@uj.edu.pl
# 
# Copyright (C) 2017 Jakub Szlęk, Aleksander Mendyk
# 
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General 
# Public License as published by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program.
# If not, see <http://www.gnu.org/licenses/>.
# 
# File contains auxiliary functions used by the script:
# BEMoDA_Dep.R
# 


#### Time points extraction from data table

extract.timepoints <- function(dataframe){

    time <- gsub("X","",colnames(dataframe))
    time <- gsub(".min","",time)
    time <- as.numeric(time)

    return(time)
}


#### To find out the points on the borded of critical region (CR)

draw.ellipse <- function(x, y){
#
# Estimate the parameters.
#
X <- as.data.frame(cbind(One=1, x, y, xy=x*y, x2=x^2, y2=y^2))
fit <- lm(One ~ . - 1, X)
beta.hat <- coef(fit)
#
# Plot the estimate, the point, and the original axes.
#
evaluate <- function(x, y, beta) {
  if (missing(y)) {
    y <- x[, 2]; x <- x[, 1]
  }
  as.vector(cbind(x, y, x*y, x^2, y^2) %*% beta - 1)
}
e.x <- diff(range(x)) / 40
e.y <- diff(range(y)) / 40
n.x <- 100
n.y <- 60
u <- seq(min(x)-e.x, max(x)+e.x, length.out=n.x)
v <- seq(min(y)-e.y, max(y)+e.y, length.out=n.y)
z <- matrix(evaluate(as.matrix(expand.grid(u, v)), beta=beta.hat), n.x)

return(list(u=u,v=v,z=z))
}

#### Function to search the par[1] and par[2] via nloptr function for critical region (CR)

cr.ellipse <- function(par,S.cr,k.val,mean.diff,p,nr,level){

    f.val <- qf(level,p,2*nr-p-1)
    
    
# The part of an below equation 
# 
#   left.eq <- S[1,1]*(par[1]-mean.diff[1])^2 + S[1,2]*(par[1]-mean.diff[1])*S[2,1]*(par[2]-mean.diff[2]) + S[2,2]*(par[2]-mean.diff[2])^2    
# 
#  can be substituted by: 
# 
    left.eq <- t(par-mean.diff) %*% solve(S.cr) %*% (par-mean.diff)
   
    right.eq <- f.val / k.val # k.val <- ((nt+nr-p-1)/((nt+nr-2)*p))*((nt*nr)/(nt+nr)) 

    res <- mean((left.eq-right.eq)^2)
    
    if (is.na(res)){
        res<-Inf
    }
    
    
    return(res)
}


#### Function to search the par[1] and par[2] via nloptr function for similarity region (SR)

sr.ellipse <- function(par,S, k.val, p, nr, level){

    alpha <- par[1]
    beta <- par[2]
    f.val <- qf(level,p,2*nr-p-1)
    
# The part of an below equation 
# 
#   left.eq <- S[1,1]*(par[1]-mean.diff[1])^2 + S[1,2]*(par[1]-mean.diff[1])*S[2,1]*(par[2]-mean.diff[2]) + S[2,2]*(par[2]-mean.diff[2])^2    
# 
#  can be substituted by: 
# 
    left.eq <- t(par) %*% solve(S) %*% (par)
   
    right.eq <- (f.val / k.val)

    res <- mean((left.eq-right.eq)^2)
    
    if (is.na(res)){
        res<-Inf
    }
    
    
    return(res)
}


### Weibull model function to optimize

weibull <- function(params, t, X){
    
    alpha <- params[1]
    beta <- params[2]
    pred <- 1-exp(-alpha*t^(beta))
    res <- mean((X-pred)^2)
    
    if (is.na(res)){
        res<-Inf
    }
    
    return(res)
}


### Weibull model with fixed parameters of alpha and beta

weibull.fit <- function(t, alpha.calc,beta.calc){
    
    alpha <- alpha.calc
    beta <- beta.calc
    
    diss <- 1-exp(-alpha*t^(beta))
    
    return(diss)
}

### RMSE function # Function that returns Root Mean Squared Error

RMSE <- function(predicted, observed){
    
    error <- observed - predicted
    res <- sqrt(mean(error^2))
    
    return(res)
}


### Taken after package sparsediscrim, to calculate covariance matrix for multi-sample pooled covariance

cov_pool <- function (x, y){
    x <- as.matrix(x)
    y <- as.factor(y)
    n <- length(y)
    scatter_matrices <- tapply(seq_len(n), y, function(i) {
        (length(i) - 1) * cov(as.matrix(x[i, ]))
    })
    as.matrix(Reduce("+", scatter_matrices)/(n-nlevels(y))) ### MODIFIED PUT (n-nlevels(y)) to reflect the Sathe et al.
}
