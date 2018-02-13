# BEMoDA v1.0 - Bioequivalence Model Dependent-Independent Approach script for in-vitro dissolution profile comparison
# 
# Bioequivalence Model Dependent Approach script for in-vitro dissolution profile comparison as proposed by Sathe et al. in 1996
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


#######################################################
### OPTIONS FOR TEST AND REFERENCE PRODUCTS        ####
#######################################################

# dissolution data for test & reference, format in columns t1 t2 t3 t4 ...
# rows contains ref_1 ref_2 ref_3 etc.
# TAB-delimited file with column names and rownames

# filename of dissolution of test product
filename_test <- c("m_test.csv")

# filename of dissolution of reference product
filename_ref <- c("m_ref.csv")

# filename of dissolution of standard batches - they are needed to establish similarity region (SR)
filenames_std <- list.files(pattern="std.*.*")

###################################################################
### OPTIONS FOR WEIBULL OPTIMIZATION                           ####
###################################################################

# Optimization method for Weibull model fitting, same as for method parameter in optimx()
optim.method <- c("Nelder-Mead") # for optimx ! c("BFGS") or c("Nelder-Mead")

optimx.method <- FALSE
nloptr.method <- TRUE
gensa.method <- FALSE
nls.method <- FALSE

maxit_nloptr.weibull <-10000
max_iter_gensa.weibull <-5000
maxit_NM.weibull <-50000



###################################################################
### OPTIONS FOR nloptr TO SEARCH FOR CR AND SR ELLIPSES        ####
###################################################################

opti_trace<-FALSE #save drive space if FALSE

params.no <- 2 # how many parameters in the equation
starting.params <- rnorm(params.no)/10 # generate starting params for optim functions
lower.boundary <- rep(-100*max(abs(starting.params)),times=params.no) # for alpha and beta
upper.boundary <- rep(100*max(abs(starting.params)),times=params.no) # for alpha and beta

# tolerace limit for optimization
optim_rel_tol <- 1e-20

# no of steps in optimization
maxit_nloptr<-10000


sr.level <- 0.99 # similarity region confidence interval - I think this corresponds to 3SD - as 99% of delta_mu alpha/beta should fall inside this region
cr.level <- 0.9  # critical region confidence interval
ellipse.sr.npts <- 50   # number of points to discover at the similarity region ellipse
ellipse.cr.npts <- 50   # number of points to discover at the critical region ellipse


# Decide to draw ellipse or rectangle as a similarity region
draw.ellipse.SR <- TRUE
draw.rect.SR <- FALSE



####################################################################
###                         VALIDATION MODE                     ####
####################################################################

# If validation.mode = TRUE ONLY for comparison with the article by Sathe et al.
# 
# Files should inculde data from Sathe P M, Tsong Y, Shah V P. Pharm Res, 1996, 13(12): 1799-1803.
# of alpha, beta parameter for reference and test product
validation.mode <- FALSE

# filename of Weibull parameters of test product
val.test <- c("val_test1.csv")

# filename of Weibull parameters of reference product
val.ref <- c("val_ref.csv")

# filename of Weibull parameters of standard batches
val.std <- c("val_s1.csv")


# Additonally check with MM
check.MM <- FALSE # NEED TO PROVIDE "val_test2.csv" file

if(validation.mode==FALSE){
    check.MM <- FALSE
}


################################
#### LOAD REQUIRED PACKAGES ####
################################
require(dplyr)
require(optimx)
require(ggplot2)
require(nloptr)
require(GenSA)
require(MASS)


#######################
#### LOAD THE DATA ####
#######################

# read dissolution data for test and reference batches
mt <- read.csv(file=paste(filename_test),header=TRUE,row.names=1,sep="\t")
mr <- read.csv(file=paste(filename_ref),header=TRUE,row.names=1,sep="\t")

# read dissolution data from standard batches
for(i in 1:length(filenames_std)){

stdname <- paste("std",i,sep="")

assign(stdname, read.csv(file=paste(filenames_std[i]),header=TRUE,row.names=1,sep="\t"))

}

std.no <- length(filenames_std)

#################################################
##### SOURCE FILE WITH AUXILIARY FUNCTIONS ######
#################################################

# Put here:
source("./BEMoDA_auxiliary_functions.R")

##############################
##### CHECK DATA FORMAT ######
#############################

nt <- nrow(mt)
nr <- nrow(mr)
ns <- nrow(std1)

if(nt != nr || ns != nr){
    stop("Number of observations (dissolution profiles) in test and/or reference/standard differ! Please check the data!")
}


p <- ncol(mt)

if(p != ncol(mr)){
    stop("Number of parameters (time points) in test and reference differ! Please check the data!")
}

# WHEN REPORT.TXT DEPRECATED
# 
# ###########################################
# ###         PRINT OUT LOADED DATA      ####
# ###########################################
# 
# cat("REFERENCE DATA:","\n\n")
# mr
# cat("TEST DATA:","\n\n")
# mt

###########################
### TRANSFORM THE DATA ####
###########################



# Prepare data frame for writing Weibull parameters
par.ref.df <- matrix(ncol = 3, nrow = nr)
par.test.df <- matrix(ncol = 3, nrow = nt)


for(i in 1:length(filenames_std)){

par.std.df.name <- paste("par.std.df",i,sep="")

assign(par.std.df.name, setNames(data.frame(matrix(ncol=3, nrow = ns)), c("Formulation","Weibull alpha","Weibull beta")))

}


tmp.x <- c("Formulation","Weibull alpha","Weibull beta")
colnames(par.ref.df) <- tmp.x
colnames(par.test.df) <- tmp.x


# Prepare matrix to save RMSE
diss.rmse.ref <- matrix(ncol = 2, nrow = nr)
diss.rmse.test <- matrix(ncol = 2, nrow = nt)

tmp.x <- c("Formulation","RMSE")
colnames(diss.rmse.ref) <- tmp.x
colnames(diss.rmse.test) <- tmp.x


#####################################################
##### REPORT - SECTION 1 - BEFORE CALCULATIONS ######
#####################################################

sink("report.txt",append=TRUE)

cat("DATE:","\n",date(),"\n")
cat("\n")

cat("SESSION INFO:","\n")
sessionInfo()
cat("\n")
cat("Provided filenames","\n")
cat("filename_test = ", filename_test,"\n")
cat("filename_ref = ", filename_ref,"\n")
cat("filename_std = ", filenames_std,"\n")
cat("\n")
cat("\n")

cat("OPTIONS FOR WEIBULL OPTIMIZATION","\n")
cat("optim.method = ",optim.method,"\n")
cat("optimx.method = ", optimx.method,"\n")
cat("nloptr.method = ", nloptr.method,"\n")
cat("gensa.method = ", gensa.method,"\n")
cat("nls.method = ", nls.method,"\n")
cat("maxit_nloptr.weibull = ", maxit_nloptr.weibull, "\n")
cat("max_iter_gensa.weibull = ", max_iter_gensa.weibull, "\n")
cat("maxit_NM.weibull = ", maxit_NM.weibull, "\n")
cat("\n")
cat("\n")

cat("OPTIONS FOR nloptr TO SEARCH FOR CR AND SR ELLIPSES","\n")
cat("opti_trace = ", opti_trace, "\n")
cat("params.no = ", params.no, "\n")
cat("starting.params = ","\n",starting.params,"\n")
cat("lower.boundary = ", "\n", lower.boundary,"\n")
cat("upper.boundary = ", "\n", upper.boundary,"\n")
cat("optim_rel_tol = ", optim_rel_tol, "\n")
cat("maxit_nloptr = ", maxit_nloptr, "\n")
cat("sr.level = ", sr.level, "\n")
cat("cr.level = ", cr.level, "\n")
cat("ellipse.sr.npts = ", ellipse.sr.npts, "\n")
cat("ellipse.cr.npts = ", ellipse.cr.npts, "\n")
cat("draw.ellipse.SR = ", draw.ellipse.SR, "\n")
cat("draw.rect.SR = ", draw.rect.SR, "\n")
cat("\n")
cat("\n")

cat("VALIDATION OPTIONS","\n")
cat("validation.mode = ", validation.mode, "\n")
cat("val.test = ", val.test, "\n")
cat("val.ref = ", val.ref, "\n")
cat("val.std = ", val.std, "\n")
cat("check.MM = ", check.MM, "\n")
cat("\n")
cat("\n")

if(validation.mode == FALSE){
        cat("LOADED DATA","\n")
        cat("REFERENCE DATA:","\n")
        print(mr)
        cat("\n")
        cat("\n")
        cat("TEST DATA:","\n")
        print(mt)
        cat("\n")
        cat("\n")


        for(i in 1:length(filenames_std)){
        cat("STANDARD NO ",i," DATA:","\n")
        print(get(paste("std",i,sep="")))
        cat("\n")
        cat("\n")

        }
} 

sink()

############################################################
##### END OF REPORT - SECTION 1 - BEFORE CALCULATIONS ######
############################################################



#################################################################
###        WEIBULL MODEL FITTING FOR REF AND TEST            ####
#################################################################

for(i in 1:nt){

# Prepare new data frame
# Prepare data points - acquire them from the table header

dat.r <- extract.timepoints(mr)
dat.t <- extract.timepoints(mt)

dat.r <-rbind(t=dat.r,X=mr[i,])
dat.t <-rbind(t=dat.t,X=mt[i,])



# Make fractions of Q and transform time into hours
dat.r[1,] <- dat.r[1,]/60
dat.r[2,] <- dat.r[2,]/100

dat.t[1,] <- dat.t[1,]/60
dat.t[2,] <- dat.t[2,]/100


# Transpose data frame
dat.r <- as.data.frame(t(dat.r))
dat.t <- as.data.frame(t(dat.t))


# optimx optim method
if(optimx.method == TRUE){

weibull.optimx.r <- optimx(starting.params,weibull,t=dat.r[,1],X=dat.r[,2],method=optim.method, control=list(trace=opti_trace,maxit=maxit_NM))
weibull.optimx.t <- optimx(starting.params,weibull,t=dat.t[,1],X=dat.t[,2],method=optim.method, control=list(trace=opti_trace,maxit=maxit_NM))

# Get the alpha and beta parameters from the model
r.alpha.calc <- weibull.optimx.r$p1
r.beta.calc <- weibull.optimx.r$p2

t.alpha.calc <- weibull.optimx.t$p1
t.beta.calc <- weibull.optimx.t$p2

}


# nloptr optim method

if(nloptr.method == TRUE){

weibull.nloptr.r <- nloptr(x0=starting.params, eval_f=weibull, lb=lower.boundary, ub=upper.boundary,
                            t=dat.r[,1],X=dat.r[,2],
                            opts=list(algorithm="NLOPT_GN_CRS2_LM",xtol_rel=optim_rel_tol,maxeval=maxit_nloptr,print_level=1,local_opts=list(algorithm="NLOPT_LD_MMA",xtol_rel=optim_rel_tol))
                            )

weibull.nloptr.t <- nloptr(x0=starting.params, eval_f=weibull, lb=lower.boundary, ub=upper.boundary,
                            t=dat.t[,1],X=dat.t[,2],
                            opts=list(algorithm="NLOPT_GN_CRS2_LM",xtol_rel=optim_rel_tol,maxeval=maxit_nloptr,print_level=1,local_opts=list(algorithm="NLOPT_LD_MMA",xtol_rel=optim_rel_tol))
                            )
                            
                            
  r.alpha.calc <- weibull.nloptr.r$solution[1]
  r.beta.calc <- weibull.nloptr.r$solution[2]
  
  t.alpha.calc <- weibull.nloptr.t$solution[1]
  t.beta.calc <- weibull.nloptr.t$solution[2]
                            
}

# GenSA optim method

if(gensa.method == TRUE){

  weibull.gensa.r <- GenSA(par=starting.params, fn=weibull,
                            t=dat.r[,1],X=dat.r[,2],lower=lower.boundary,
                            upper=upper.boundary,
                            control=list(smooth=FALSE,maxit=max_iter_gensa,verbose=TRUE,nb.stop.improvement=max_iter_gensa)
                            )

  
  weibull.gensa.t <- GenSA(par=starting.params, fn=weibull,
                            t=dat.t[,1],X=dat.t[,2],lower=lower.boundary,
                            upper=upper.boundary,
                            control=list(smooth=FALSE,maxit=max_iter_gensa,verbose=TRUE,nb.stop.improvement=max_iter_gensa)
                            )

  r.alpha.calc <- weibull.gensa.r$par[1] 
  r.beta.calc <- weibull.gensa.r$par[2]
  
  t.alpha.calc <- weibull.gensa.t$par[1] 
  t.beta.calc <- weibull.gensa.t$par[2]
  
}


# nls optim method
if(nls.method == TRUE){

# Predict nonlinear least square model
fit.r <- nls(X ~ 1-exp(-alpha*t^(beta)),start=list(alpha=1,beta=1),data=as.data.frame(dat.r),control = list(maxiter = 5000,warnOnly=T),trace=TRUE)
fit.t <- nls(X ~ 1-exp(-alpha*t^(beta)),start=list(alpha=1,beta=1),data=as.data.frame(dat.t),control = list(maxiter = 5000,warnOnly=T),trace=TRUE)

# Get the alpha and beta parameters from the model
r.alpha.calc <- coef(fit.r)[1]
r.beta.calc <- coef(fit.r)[2]

t.alpha.calc <- coef(fit.t)[1]
t.beta.calc <- coef(fit.t)[2]

}



# Save parameters in data frame par.ref.df and par.test.df
par.ref.df[i,1] <- i
par.ref.df[i,2] <- r.alpha.calc
par.ref.df[i,3] <- r.beta.calc

par.test.df[i,1] <- i
par.test.df[i,2] <- t.alpha.calc
par.test.df[i,3] <- t.beta.calc


# Prepare new data for time values
new.t.ref <- seq(0,max(dat.r[,1]),0.05)
new.t.test <- seq(0,max(dat.t[,1]),0.05)


# Calculate new Q values with new Weibull model
new.X.ref <- weibull.fit(new.t.ref,r.alpha.calc,r.beta.calc)
new.X.test <- weibull.fit(new.t.test,t.alpha.calc,t.beta.calc)


# Calculate RMSE for optimized Weibull equation

# First calculate predicted dissolution based on Weibull model
diss.ref <- weibull.fit(dat.r[,1],r.alpha.calc,r.beta.calc)
diss.test <- weibull.fit(dat.t[,1],t.alpha.calc,t.beta.calc)

# Calculate % of dissolution
diss.ref <- 100*diss.ref
diss.test <- 100*diss.test

dat.r.diss <- 100*dat.r[,2]
dat.t.diss <- 100*dat.t[,2]

# RMSE
rmse.ref <- RMSE(diss.ref, dat.r.diss)
rmse.test <- RMSE(diss.test, dat.t.diss)

# Save RMSE in matrix
diss.rmse.ref[i,1] <- i
diss.rmse.ref[i,2] <- rmse.ref

diss.rmse.test[i,1] <- i
diss.rmse.test[i,2] <- rmse.test

cat("\n")
cat("=========================","\n")
cat("ref ",i,"\n")
cat("alpha = ", r.alpha.calc,"\n")
cat("beta = ", r.beta.calc,"\n")
cat("=========================","\n")

cat("\n")
cat("=========================","\n")
cat("test ",i,"\n")
cat("alpha = ", t.alpha.calc,"\n")
cat("beta = ", t.beta.calc,"\n")
cat("=========================","\n")



if(validation.mode == FALSE){
# Plot the points and fitted Weibull curve

ggplot(data=dat.r, aes(x=t,y=X),color='blue') + geom_point() +
geom_line(color='red',data=data.frame(t=new.t.ref,X=new.X.ref),aes(x=t,y=X)) +
labs(x="Time [hours]",y="Q [fraction dissolved]") + ggtitle(paste("REFERENCE ",i," Weibull model")) + 
theme(plot.title = element_text(lineheight=.8, face="bold"))
ggsave(paste("REF_",i,"_Weibull_plot.pdf",sep=""))

ggplot(data=dat.t, aes(x=t,y=X),color='blue') + geom_point() +
geom_line(color='red',data=data.frame(t=new.t.test,X=new.X.test),aes(x=t,y=X)) +
labs(x="Time [hours]",y="Q [fraction dissolved]") + ggtitle(paste("TEST ",i," Weibull model")) + 
theme(plot.title = element_text(lineheight=.8, face="bold"))
ggsave(paste("TEST_",i,"_Weibull_plot.pdf",sep=""))
    }
    
}

par.test.df <- as.data.frame(par.test.df[,2:3])
par.ref.df <- as.data.frame(par.ref.df[,2:3])

cat("","\n")
cat("Weibull model parameters for REFERENCE")
print(par.ref.df)
cat("\n")
cat("\n")
cat("Weibull model parameters for TEST")
print(par.test.df)
cat("\n")
cat("\n")


########################################################
###        WEIBULL MODEL FITTING FOR STD            ####
########################################################

for(k in 1:length(filenames_std)){

tmp.std <- get(paste("std",k,sep=""))   

tmp.ns <- nrow(tmp.std)

par.std.df.name <- paste("par.std.df",k,sep="")

assign(paste("tmp.par"),setNames(data.frame(matrix(nrow=tmp.ns,ncol=3)),c("Formulation","Weibull alpha","Weibull beta")))


for(i in 1:ns){

# Prepare new data frame
# Prepare data points - acquire them from the table header

dat.s <- extract.timepoints(tmp.std)

dat.s <-rbind(t=dat.s,X=tmp.std[i,])

# Make fractions of Q and transform time into hours
dat.s[1,] <- dat.s[1,]/60
dat.s[2,] <- dat.s[2,]/100

# Transpose data frame
dat.s <- as.data.frame(t(dat.s))

# optimx optim method
if(optimx.method == TRUE){

weibull.optimx.s <- optimx(starting.params,weibull,t=dat.s[,1],X=dat.s[,2],method=optim.method, control=list(trace=opti_trace,maxit=maxit_NM.weibull))

# Get the alpha and beta parameters from the model
s.alpha.calc <- weibull.optimx.s$p1
s.beta.calc <- weibull.optimx.s$p2

}


# nloptr optim method

if(nloptr.method == TRUE){

weibull.nloptr.s <- nloptr(x0=starting.params, eval_f=weibull, lb=lower.boundary, ub=upper.boundary,
                            t=dat.s[,1],X=dat.s[,2],
                            opts=list(algorithm="NLOPT_GN_CRS2_LM",xtol_rel=optim_rel_tol,maxeval=maxit_nloptr.weibull,print_level=1,local_opts=list(algorithm="NLOPT_LD_MMA",xtol_rel=optim_rel_tol))
                            )
            
  s.alpha.calc <- weibull.nloptr.s$solution[1]
  s.beta.calc <- weibull.nloptr.s$solution[2]
  
}

# GenSA optim method

if(gensa.method == TRUE){

  weibull.gensa.s <- GenSA(par=starting.params, fn=weibull,
                            t=dat.s[,1],X=dat.s[,2],lower=lower.boundary,
                            upper=upper.boundary,
                            control=list(smooth=FALSE,maxit=max_iter_gensa.weibull,verbose=TRUE,nb.stop.improvement=max_iter_gensa.weibull)
                            )

  s.alpha.calc <- weibull.gensa.s$par[1] 
  s.beta.calc <- weibull.gensa.s$par[2]
  
}


# nls optim method
if(nls.method == TRUE){

# Predict nonlinear least square model
fit.s <- nls(X ~ 1-exp(-alpha*t^(beta)),start=list(alpha=1,beta=1),data=as.data.frame(dat.s),control = list(maxiter = 5000,warnOnly=T),trace=TRUE)

# Get the alpha and beta parameters from the model
s.alpha.calc <- coef(fit.s)[1]
s.beta.calc <- coef(fit.s)[2]

}

tmp.par[i,1] <- i
tmp.par[i,2] <- s.alpha.calc
tmp.par[i,3] <- s.beta.calc




# Prepare new data for time values
new.t.s <- seq(0,max(dat.s[,1]),0.05)

# Calculate new Q values with new Weibull model
new.X.s <- weibull.fit(new.t.s,s.alpha.calc,s.beta.calc)

# Calculate RMSE for optimized Weibull equation

# First calculate predicted dissolution based on Weibull model
diss.std <- weibull.fit(dat.s[,1],s.alpha.calc,s.beta.calc)

# Calculate % of dissolution
diss.std <- 100*diss.std

dat.s.diss <- 100*dat.s[,2]

# RMSE
rmse.std <- RMSE(diss.std, dat.s.diss)

print("RMSE")
print(rmse.std)
# 
# # Save RMSE in matrix
# diss.rmse.std[i,1] <- i
# diss.rmse.std[i,2] <- rmse.std

cat("\n")
cat("=========================","\n")
cat("ref ",i,"\n")
cat("alpha = ", s.alpha.calc,"\n")
cat("beta = ", s.beta.calc,"\n")
cat("=========================","\n")


if(validation.mode == FALSE){
# Plot the points and fitted Weibull curve

# pdf(paste(i,"_mt_plot.pdf"))

ggplot(data=dat.s, aes(x=t,y=X),color='blue') + geom_point() +
geom_line(color='red',data=data.frame(t=new.t.s,X=new.X.s),aes(x=t,y=X)) +
labs(x="Time [hours]",y="Q [fraction dissolved]") + ggtitle(paste("STANDARD ",i," Weibull model")) + 
theme(plot.title = element_text(lineheight=.8, face="bold"))
ggsave(paste("STD_",k,"_",i,"_Weibull_plot.pdf",sep=""))

    }
}

assign(par.std.df.name, tmp.par[,2:3])

cat("","\n")
cat("Weibull model parameters for STD_",k)
print(tmp.par)
cat("\n")

}

# Calculate log-normal matrix for std batches

for(i in 1:length(filenames_std)){

par.std.df.name <- paste("par.std.df",i,sep="")

assign(par.std.df.name, log(get(par.std.df.name)))

}



# For the future - correct the code and insert this transformation into a list much earlier in the code
# I would be much more efficient

std.list <- vector(mode = "list", length = length(filenames_std))

for (i in 1:length(filenames_std)) {
    std.list[[i]] <- get(paste("par.std.df",i,sep=""))
}

all.std.par <- bind_rows(std.list,.id="id")



### S.pooled - variance/covariance matrix of standard batches (std) 

S.sr <- cov_pool(all.std.par[,2:3],all.std.par[,1])


if(validation.mode == TRUE){

    # read parameters data for test and reference batches
        mt <- read.csv(file=paste(val.test),header=TRUE,row.names=1,sep="\t")
        mr <- read.csv(file=paste(val.ref),header=TRUE,row.names=1,sep="\t")
        mv <- read.csv(file=paste(val.std),header=TRUE,row.names=1,sep="\t")
        
        S.sr <- var(mv)

}

# # S.pooled for SR
# n.std1 <- nrow(par.std.df1)
# n.std2 <- nrow(par.std.df2)
# n.std3 <- nrow(par.std.df3)

# S.pooled - variance/covariance matrix of standard batches (std) 
# S.sr <- cov_pool(all.std.par[,2:3],all.std.par[,1])
# S.sr2 <- ((n.std1 - 1) * cov(par.std.df1) + (n.std2 - 1) * cov(par.std.df2) + (n.std3 - 1) * cov(par.std.df3)) / (n.std1  + n.std2 + n.std3 - 3)
# S.sr == S.sr2
#               Weibull alpha Weibull beta
# Weibull alpha          TRUE         TRUE
# Weibull beta           TRUE         TRUE

#####################################################################





#####################################################
##### REPORT - SECTION 2 - WEIBULL PARAMETERS  ######
#####################################################

sink("report.txt",append=TRUE)

cat("DATE:","\n",date(),"\n")
cat("\n")

cat()

if(validation.mode == FALSE){
        cat("WEIBULL PARAMETERS","\n")
        cat("REFERENCE PRODUCT:","\n")
        print(par.ref.df)
        cat("\n")
        cat("\n")
        cat("TEST PRODUCT:","\n")
        print(par.test.df)
        cat("\n")
        cat("\n")
        cat("STANDARD LOG-NORMAL WEIBULL PARAMETERS:","\n")
        print(std.list)
        cat("\n")
        cat("\n")

        
} else {

        cat("LOADED VALIDATION DATA","\n")
        cat("VALIDATION REFERENCE DATA:","\n")
        print(mr)
        cat("\n")
        cat("\n")
        cat("VALIDATION TEST DATA:","\n")
        print(mt)
        cat("\n")
        cat("\n")
        cat("VALIDATION STANDARD DATA:","\n")
        print(mv)
        cat("\n")
        cat("\n")

}

cat("\n")

sink()

############################################################
##### END OF REPORT - SECTION 2 - WEIBULL PARAMETERS  ######
############################################################


# As an option? rectangle SR?

# Univariate SD.pooled
# In In Vitro-In Vivo Correlations, Eds.: David B. Young,John G. Devane,Jackie Butler
# In Vitro Dissolution Profile Comparison And IVIVR
# 
# # matrix of mean variances, intra-lot variances

sd.intra.lot.alpha <- var(all.std.par[,1])
sd.intra.lot.beta <- var(all.std.par[,2])

sd.log.std.alpha <- sapply(std.list, function(x) var(x[,1]))
sd.log.std.beta <- sapply(std.list, function(x) var(x[,2]))

sqr.sd.log.std.alpha <- sd.log.std.alpha^2
sqr.sd.log.std.beta <- sd.log.std.beta^2

SD.pooled.alpha <- sqrt((sum(sqr.sd.log.std.alpha)/length(std.list))+sd.intra.lot.alpha)
SD.pooled.beta <- sqrt((sum(sqr.sd.log.std.beta)/length(std.list))+sd.intra.lot.beta)

if(validation.mode == TRUE){

    SD.pooled.alpha <- var(mv[,1]) 
    SD.pooled.beta <- var(mv[,2])
}


cat("SD.pooled.alpha","\n")
SD.pooled.alpha
cat("SD.pooled.beta","\n")
SD.pooled.beta

SR1 <- data.frame(x1=SD.pooled.alpha, x2=-SD.pooled.alpha, y1=SD.pooled.beta, y2=-SD.pooled.beta)
SR2 <- data.frame(x1=2*SD.pooled.alpha, x2=-2*SD.pooled.alpha, y1=2*SD.pooled.beta, y2=-2*SD.pooled.beta)
SR3 <- data.frame(x1=3*SD.pooled.alpha, x2=-3*SD.pooled.alpha, y1=3*SD.pooled.beta, y2=-3*SD.pooled.beta)

# About log-transformation
# After Zhang et al. An introduction to the approaches used by DDSolver
# (4) Because the MSD is calculated under the assumption that the model parameters are multivariate normally distributed, it should be determined whether the model parameters need to be log-transformed before the comparison is performed, which depends mainly on the statistical distribution properties of the compared parameters. Take the widely used Weibull model, for example: some studies have concluded similarity between two sets of parameters based on natural-logarithm-transformed data (31), while some other studies have used untransformed data (32). To ensure rational application of model-dependent approaches, this issue obviously needs to be explicitly clarified.

log.mt <- log(par.test.df)
log.mr <- log(par.ref.df)

# Check once again number of parameters and tablets
nt <- nrow(log.mt)
nr <- nrow(log.mr)
p <- ncol(log.mt)

if(validation.mode == TRUE){
    
    mt <- read.csv(file=paste(val.test),header=TRUE,row.names=1,sep="\t")
    mr <- read.csv(file=paste(val.ref),header=TRUE,row.names=1,sep="\t")
    
    log.mt <- log(mt)
    log.mr <- log(mr)
    
    nt <- nrow(mt)
    nr <- nrow(mr)
    p <- ncol(mt)


}


# Scaling factor
k.val <- ((nt+nr-p-1)/((nt+nr-2)*p))*((nt*nr)/(nt+nr)) 


# F-distribution value of critical region
F.cr <- qf(0.90,p,(2*nt-p-1))

# F-distribution value of similarity region 
F.sr <- qf(0.99,p,(2*nt-p-1))

# Covariance matrix for S_pooled 
S.cr <- (cov(log.mt)+cov(log.mr))/2

# Means of paramters for test and reference dissolution profiles
mean.log.mt <- colMeans(log.mt)
mean.log.mr <- colMeans(log.mr)
mean.log.diff <- mean.log.mt - mean.log.mr

# Difference of parameters
# log.diff <- log.mt - log.mr
log.diff <- log.mt - log.mr

# Mean of parameters difference
# mean.log.sr <- colMeans(rbind(colMeans(log.std1), colMeans(log.std2),colMeans(log.std3)))

###############################################################################################                            
#####     CALCULATE MAHALANOBIS DISTANCE AND HOTTELING T^2 FOR REFERENCE AND TEST          ####
###############################################################################################


# Mahalanobis distance, 'M' Distance
M.dist.cr <- sqrt(t(mean.log.diff) %*% solve(S.cr) %*% mean.log.diff)

# Hotelling's T^2, as Scaled 'M' Distance as in Sathe et. al - Hotelling's two-sample t-squared statistic
# https://en.wikipedia.org/wiki/Hotelling%27s_T-squared_distribution
# http://sites.stat.psu.edu/~ajw13/stat505/fa06/11_2sampHotel/01_2sampHotel.html
# 
M.dist.cr.scaled <- (nt*nr)/(nt+nr)*M.dist.cr^2



######################################################################################################################                            
#####     SEARCH FOR POINTS ON THE BORDER OF CRITICAL REGION (OPTIMIZE FUNCTION BY NLOPTR) AND SIMILARITY REGION  ####
######################################################################################################################




ellipse.sr <- matrix(nrow=ellipse.sr.npts,ncol=2,data=NA)
ellipse.cr <- matrix(nrow=ellipse.cr.npts,ncol=2,data=NA)




for (i in 1:ellipse.sr.npts){
fit <- nloptr(x0=starting.params, eval_f=sr.ellipse, lb=lower.boundary, ub=upper.boundary,
                            S=S.sr, k.val=k.val, p=p, nr=nr, level=sr.level, 
                            opts=list(algorithm="NLOPT_GN_CRS2_LM",xtol_rel=optim_rel_tol,maxeval=maxit_nloptr,print_level=1,local_opts=list(algorithm="NLOPT_LD_MMA",xtol_rel=optim_rel_tol))
                            )
ellipse.sr[i,] <- fit$solution
                            
                            
}
                            

for (i in 1:ellipse.cr.npts){
fit <- nloptr(x0=starting.params, eval_f=cr.ellipse, lb=lower.boundary, ub=upper.boundary,
                            S.cr=S.cr, mean.diff=mean.log.diff, k.val=k.val, p=p, nr=nr, level=cr.level, 
                            opts=list(algorithm="NLOPT_GN_CRS2_LM",xtol_rel=optim_rel_tol,maxeval=maxit_nloptr,print_level=1,local_opts=list(algorithm="NLOPT_LD_MMA",xtol_rel=optim_rel_tol))
                            )
ellipse.cr[i,] <- fit$solution
                            
                            
}
                            


# To check if MM-ref will fall outside SR                            

if(check.MM == TRUE){

# read dissolution data for test and reference batches
filename_MM <- c("val_test2.csv")
MMt <- read.csv(file=filename_MM,header=TRUE,row.names=1,sep="\t")


log.MMt <- log(MMt)
S.MMt <- (cov(log.MMt)+cov(log.mr))/2

log.diff.MMt <- log.MMt - log.mr

mean.log.diff.MMt <- colMeans(log.diff.MMt)


#################################################################################################                            
#####     CALCULATE MAHALANOBIS DISTANCE AND HOTTELING T^2 FOR REF AND MM AFTER SATHE et. al ####
#################################################################################################

# Mahalanobis distance, 'M' Distance
M.dist.MMt <- sqrt(t(mean.log.diff.MMt) %*% solve(S.MMt) %*% mean.log.diff.MMt)

# Hotelling's T^2, as Scaled 'M' Distance as in Sathe et. al - Hotelling's two-sample t-squared statistic
# https://en.wikipedia.org/wiki/Hotelling%27s_T-squared_distribution
# http://sites.stat.psu.edu/~ajw13/stat505/fa06/11_2sampHotel/01_2sampHotel.html
# 
M.dist.MMt.scaled <- nt*nr/(nt+nr)*M.dist.MMt^2


ellipse.MM <- matrix(nrow=ellipse.cr.npts,ncol=2,data=NA)
                   
for (i in 1:ellipse.cr.npts){
fit <- nloptr(x0=starting.params, eval_f=cr.ellipse, lb=lower.boundary, ub=upper.boundary,
                            S.cr=S.MMt, mean.diff=mean.log.diff.MMt, k.val=k.val, p=p, nr=nr, level=cr.level, 
                            opts=list(algorithm="NLOPT_GN_CRS2_LM",xtol_rel=optim_rel_tol,maxeval=maxit_nloptr,print_level=1,local_opts=list(algorithm="NLOPT_LD_MMA",xtol_rel=optim_rel_tol))
                            )
ellipse.MM[i,] <- fit$solution
                            
                            
    }
}                            
                            
                            
                            
                            
###################################################################                            
#####     FIT ELLIPSE BASED ON THE POINTS GENERETED BY NLOPTR  ####
###################################################################



ellipse1 <- draw.ellipse(ellipse.cr[,1], ellipse.cr[,2])
ellipse2 <- draw.ellipse(ellipse.sr[,1], ellipse.sr[,2])
ellipse1.pts <- contourLines(ellipse1$u, ellipse1$v, ellipse1$z, levels=0)
ellipse2.pts <- contourLines(ellipse2$u, ellipse2$v, ellipse2$z, levels=0)



xlim <- c(min=min(ellipse.cr[,1],ellipse.sr[,1]), max=max(ellipse.cr[,1],ellipse.sr[,1]))
ylim <- c(min=min(ellipse.cr[,2],ellipse.sr[,2]), max=max(ellipse.cr[,2],ellipse.sr[,2]))


if(check.MM == TRUE){
    ellipse3 <- draw.ellipse(ellipse.MM[,1], ellipse.MM[,2])
    ellipse3.pts <- contourLines(ellipse3$u, ellipse3$v, ellipse3$z, levels=0)
    
    xlim <- c(min=min(ellipse.cr[,1],ellipse.sr[,1],ellipse.MM[,1]), max=max(ellipse.cr[,1],ellipse.sr[,1],ellipse.MM[,1]))
    ylim <- c(min=min(ellipse.cr[,2],ellipse.sr[,2],ellipse.MM[,2]), max=max(ellipse.cr[,2],ellipse.sr[,2],ellipse.MM[,2]))
    
    
}


# Plot and save resulting graph!

png("output.png", width=1800, height=1800, res=300)

contour(ellipse1$u, ellipse1$v, ellipse1$z, levels=0, lwd=2, xlab="Ln(alpha)", ylab="Ln(beta)", asp=1,col="red", xlim=c(xlim[1]*1.10,xlim[2]*1.10),ylim=c(ylim[1]*1.10,ylim[2]*1.10), label="mm-REF", labcex=0.5)


if(draw.ellipse.SR == TRUE){
contour(ellipse2$u, ellipse2$v, ellipse2$z, levels=0, lwd=2, xlab="Ln(alpha)", ylab="Ln(beta)", asp=1,col="blue",add=TRUE, label="Similarity\n region\n (SR)", labcex=0.5)
}

if(draw.rect.SR == TRUE){
rect(xleft=SR1$x2, ybottom=SR1$y2, xright=SR1$x1, ytop=SR1$y1, border = "black", lty=1) # 1 SD
rect(xleft=SR2$x2, ybottom=SR2$y2, xright=SR2$x1, ytop=SR2$y1, border = "black", lty=2) # 2 SD
rect(xleft=SR3$x2, ybottom=SR3$y2, xright=SR3$x1, ytop=SR3$y1, border = "black", lty=3) # 3 SD

}


if(check.MM == TRUE){
    contour(ellipse3$u, ellipse3$v, ellipse3$z, levels=0, lwd=2, xlab="Ln(alpha)", ylab="Ln(beta)", asp=1,col="black",add=TRUE,label="MM-REF",labcex=0.5)
}


dev.off()


#####################################################
##### REPORT - SECTION 3 - MAHALANOBIS         ######
#####################################################


sink("report.txt",append=TRUE)

cat("MAHALANOBIS CALCULATIONS:","\n")

if(validation.mode == FALSE){
        cat("LOG-NORMAL WEIBULL PARAMETERS","\n")
        cat("REFERENCE PRODUCT:","\n")
        print(log.mr)
        cat("\n")
        cat("\n")
        cat("TEST PRODUCT:","\n")
        print(log.mt)
        cat("\n")
        cat("\n")
        
        for(i in 1:length(filenames_std)){
        cat("LOG-NORMAL STANDARD NO ",i," WEIBULL PARAMETERS:","\n")
        print(get(paste("par.std.df",i,sep="")))
        cat("\n")
        cat("\n")
        }
        
        cat("NUMBER OF ROWS IN REFERENCE DATA = ",nr,"\n")
        cat("NUMBER OF ROWS IN TEST DATA = ",nt,"\n")
        cat("NUMBER OF PARAMETERS = ",p,"\n")
        
        cat("SCALING FACTOR (K) = ", k.val,"\n")
        cat("F-DIST VALUE FOR CRITICAL REGION, 2N-P-1 (0.9) = ", F.cr,"\n")
        cat("F-DIST VALUE FOR SIMILARITY REGION, 2N-P-1 (0.99) = ", F.sr,"\n")
        cat("\n")
        cat("VARIANCE/COVARIANCE MATRIX FOR CRITICAL REGION CALCULATIONS: ","\n")
        print(S.cr)
        cat("\n")
        cat("VARIANCE/COVARIANCE MATRIX FOR SIMILARITY REGION CALCULATIONS: ","\n")
        print(S.sr)
        cat("\n")
        cat("MEAN LOG-NORMAL VALUES FOR REFERENCE: ","\n")
        print(mean.log.mr)
        cat("\n")
        cat("MEAN LOG-NORMAL VALUES FOR TEST: ","\n")
        print(mean.log.mt)
        cat("\n")
        cat("MEAN LOG-NORMAL DIFFERENCE (TEST-REF): ","\n")
        print(mean.log.diff)
        cat("\n")
        cat("\n")
        
        cat("MAHALANOBIS DISTANCE (D^2) FOR TEST-REFERENCE = ",M.dist.cr, "\n")
        cat("\n")
        cat("\n")
        cat("SCALED MAHALANOBIS DISTANCE (HOTTELING'S T^2) FOR TEST-REFERENCE = ",M.dist.cr.scaled, "\n")
        cat("\n")
        cat("\n")
        cat("CRITICAL REGION MIN-MAX","\n")
        cat("LOG-NORMAL ALPHA","\n")
        cat("(",min(ellipse1.pts[[1]]$x),", ", max(ellipse1.pts[[1]]$x) ,")","\n")
        cat("\n")
        cat("LOG-NORMAL BETA","\n")
        cat("(",min(ellipse1.pts[[1]]$y),", ", max(ellipse1.pts[[1]]$y) ,")","\n")
        cat("\n")
        cat("SIMILARITY REGION MIN-MAX","\n")
        cat("LOG-NORMAL ALPHA","\n")
        cat("(",min(ellipse2.pts[[1]]$x),", ", max(ellipse2.pts[[1]]$x) ,")","\n")
        cat("\n")
        cat("LOG-NORMAL BETA","\n")
        cat("(",min(ellipse2.pts[[1]]$y),", ", max(ellipse2.pts[[1]]$y) ,")","\n")
        
        
} else {

        cat("LOG-NORMAL WEIBULL VALIDATION PARAMETERS","\n")
        cat("VALIDATION REFERENCE:","\n")
        print(log.mr)
        cat("\n")
        cat("\n")
        cat("VALIDATION TEST (mm):","\n")
        print(log.mt)
        cat("\n")
        cat("\n")
        if(check.MM == TRUE){
        cat("VALIDATION TEST (MM):","\n")
        print(log.MMt)
        cat("\n")
        cat("\n")

        }
        
        cat("NUMBER OF ROWS IN REFERENCE DATA = ",nr,"\n")
        cat("NUMBER OF ROWS IN TEST DATA = ",nt,"\n")
        cat("NUMBER OF PARAMETERS = ",p,"\n")
        cat("SCALING FACTOR (K) = ", k.val,"\n")
        cat("F-DIST VALUE FOR CRITICAL REGION, 2N-P-1 (0.9) = ", F.cr,"\n")
        cat("F-DIST VALUE FOR SIMILARITY REGION, 2N-P-1 (0.99) = ", F.sr,"\n")
        cat("\n")
        cat("VARIANCE/COVARIANCE MATRIX FOR CRITICAL REGION CALCULATIONS (mm-REF): ","\n")
        print(S.cr)
        cat("\n")
        if(check.MM == TRUE){
        cat("\n")
        cat("VARIANCE/COVARIANCE MATRIX FOR CRITICAL REGION CALCULATIONS (MM-REF): ","\n")
        print(S.MMt)
        cat("\n")
        }
        
        cat("\n")
        cat("VARIANCE/COVARIANCE MATRIX FOR SIMILARITY REGION CALCULATIONS: ","\n")
        print(S.sr)
        cat("\n")
        cat("\n")
        cat("MEAN LOG-NORMAL VALUES FOR REFERENCE: ","\n")
        print(mean.log.mr)
        cat("\n")
        cat("\n")
        cat("MEAN LOG-NORMAL VALUES FOR TEST: ","\n")
        print(mean.log.mt)
        cat("\n")
        cat("\n")
        cat("MEAN LOG-NORMAL DIFFERENCE (mm-REF): ","\n")
        print(mean.log.diff)
        cat("\n")
        cat("\n")
        if(check.MM == TRUE){
        cat("MEAN LOG-NORMAL DIFFERENCE (MM-REF): ","\n")
        print(mean.log.diff.MMt)
        cat("\n")
        cat("\n")
        }
        
        
        cat("MAHALANOBIS DISTANCE (D^2) FOR TEST-REFERENCE = ",M.dist.cr, "\n")
        cat("\n")

        cat("SCALED MAHALANOBIS DISTANCE (HOTTELING'S T^2) FOR TEST-REFERENCE = ",M.dist.cr.scaled, "\n")
        cat("\n")
        
        if(check.MM == TRUE){
        cat("\n")
        cat("MAHALANOBIS DISTANCE (D^2) FOR TEST (MM)-REFERENCE = ",M.dist.MMt, "\n")
        cat("\n")
        cat("SCALED MAHALANOBIS DISTANCE (HOTTELING'S T^2) FOR TEST (MM)-REFERENCE = ",M.dist.MMt.scaled, "\n")
        cat("\n")

        }
        
        
        cat("CRITICAL REGION MIN-MAX","\n")
        cat("LOG-NORMAL ALPHA","\n")
        cat("(",min(ellipse1.pts[[1]]$x),", ", max(ellipse1.pts[[1]]$x) ,")","\n")
        cat("\n")

        cat("LOG-NORMAL BETA","\n")
        cat("(",min(ellipse1.pts[[1]]$y),", ", max(ellipse1.pts[[1]]$y) ,")","\n")
        cat("\n")

        
        cat("SIMILARITY REGION MIN-MAX","\n")
        cat("LOG-NORMAL ALPHA","\n")
        cat("(",min(ellipse2.pts[[1]]$x),", ", max(ellipse2.pts[[1]]$x) ,")","\n")
        cat("\n")

        cat("LOG-NORMAL BETA","\n")
        cat("(",min(ellipse2.pts[[1]]$y),", ", max(ellipse2.pts[[1]]$y) ,")","\n")
        cat("\n")
        
        if(check.MM == TRUE){
        
        cat("CRITICAL REGION (MM) MIN-MAX","\n")
        cat("LOG-NORMAL ALPHA","\n")
        cat("(",min(ellipse3.pts[[1]]$x),", ", max(ellipse3.pts[[1]]$x) ,")","\n")
        cat("\n")

        cat("LOG-NORMAL BETA","\n")
        cat("(",min(ellipse3.pts[[1]]$y),", ", max(ellipse3.pts[[1]]$y) ,")","\n")
        cat("\n")
        
        }

}

cat("\n")

sink()


############################################################
##### END OF REPORT - SECTION 3 - MAHALANOBIS         ######
############################################################
