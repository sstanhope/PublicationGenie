# y ~ the dependent variable
# D ~ searchable dichotomous covariates
# C ~ searchable continuous covariates
# Z ~ nonsearchable and not tested covariates
# S ~ the covariates of interest
# its ~ number of search iterations (1000)
# family ~ the regression model family ('binomial','gamma','gaussian','poisson')
# min.obs ~ a lower bound on the number of observations post-inclusion/exclusion criteria
#	    search.
# min.fn ~ an option on the fucntion used to compare p-values for the covariates of interest
#	   0 ~ max(p1,p2,..,pn)
#	   1 ~ sum((pi)^2)
# stop.crit ~ the number of lookback iterations used to determine early convergence
# trace ~ should search results be displayed (0=no, 1=yeah)
#
# The algorithm works by specifying a code for each covariate:
# 'I' - include
# 'E0','E1' - dichotomous covariates - exclude 0,1 observations
# 'T01','T05' - continous covariates - exclude 1,5% tail observations
# 'X' - all covariates - exclude the covariate from the model
# and doing a (1+1) genetic algorithm search to minimize the p-value.
#
# Return list:
# p0 ~ the minimized p-value
# s.D0 ~ the code for the D covariates
# s.C0 ~ the code for the C covariates.
# y0 ~ the dependent variables after observation inclusion / exclusion.
# D0 ~ the D matrix after variable exclusion and observation inclusion / exclusion.
# C0 ~ the C matrix after variable exclusion and observation inclusion / exclusion.
# Z0 ~ the nonsearchable covariates after observation inclusion / exclusion. 
# S0 ~ the covariates of interest after observation inclusion / exclusion.
# mod0 ~ the lm() object from the regression y0 ~ D0 + C0 + Z0 + S0
# p.trace ~ the trace of the objective function


publication.genie <- function(y=NULL,D=NULL,C=NULL,Z=NULL,S=NULL,
		              its=1000,family='gaussian',min.obs=0,min.fn=0,
			      stop.crit=100,trace=0)  {

if (is.null(y) | is.null(S))  {
	ret <- list(p0=NULL,s.D0=NULL,s.C0=NULL,y0=NULL,D0=NULL,C0=NULL,Z0=NULL,S0=NULL,mod0=NULL,
		    p.lookback=NULL)
	return(ret)
	}

p.lookback <- rep(NA,its+1)

if (!is.null(D))  D <- as.matrix(D)
if (!is.null(C))  C <- as.matrix(C)
if (!is.null(Z))  Z <- as.matrix(Z)
if (!is.null(S))  S <- as.matrix(S)

# define inclusion/exclusion cutpoints for continuous variables based on quartiles

if (!is.null(C))  {
	cut.01 <- rep(0,ncol(C))
	cut.05 <- rep(0,ncol(C))
	cut.95 <- rep(0,ncol(C))
	cut.99 <- rep(0,ncol(C))
	for (i in 1:ncol(C))  {
		cut.01[i] <- sort(C[,i])[0.01*nrow(C)]
		cut.05[i] <- sort(C[,i])[0.05*nrow(C)]
		cut.95[i] <- sort(C[,i])[0.95*nrow(C)]
		cut.99[i] <- sort(C[,i])[0.99*nrow(C)]
		}
	}

# iteration 1 - include all searchable covariates

if (!is.null(D))  s.D <- rep('I',ncol(D))
if (!is.null(C))  s.C <- rep('I',ncol(C))

# for dichotomous variables included in the models (cols), remove observations (rows)
# depending on 0/1 exlusions

ind.row.D <- rep(TRUE,NROW(y))
if (!is.null(D))  {
	ind.col.D <- s.D=='I' | s.D=='E0' | s.D=='E1'
	for (i in 1:ncol(D))  {	
		omit <- rep(FALSE,nrow(D))
		if (s.D[i]=='E0')  {
			omit <- D[,i] == 0
			}
		if (s.D[i]=='E1')  {
			omit <- D[,i] == 1
			}
		ind.row.D[omit] <- FALSE
		}
	}

# for continous variables included in the models (cols), remove observations (rows)
# depending on 1,5% tail exlusions

ind.row.C <- rep(TRUE,NROW(y))
if (!is.null(C))  {
	ind.col.C <- s.C=='I' | s.C=='T01' | s.C=='T05'
	for (i in 1:ncol(C))  {	
		omit <- rep(FALSE,nrow(C))
		if (s.C[i]=='T01')  {
			omit <- (C[,i] < cut.01[i]) | (C[,i] > cut.99[i])
			}
		if (s.C[i]=='T05')  {
			omit <- (C[,i] < cut.05[i]) | (C[,i] > cut.95[i])
			}
		ind.row.C[omit] <- FALSE
		}
	}

# included observations (rows) are those that are included based on all criteria

ind.row <- ind.row.D & ind.row.C

# specify submatrices based on included covariates and observations, run the 
# regression

y1 <- y[ind.row]
D1 <- D
if (!is.null(D))  D1 <- D[ind.row,ind.col.D]
C1 <- C
if (!is.null(C))  C1 <- C[ind.row,ind.col.C]
Z1 <- Z
if (!is.null(Z))  Z1 <- Z[ind.row,]
S1 <- S[ind.row,]
M1 <- cbind(D1,C1,Z1,S1)
if (family=='binomial')  mod <- glm(y1~1+M1,family='binomial')
if (family=='gamma')  mod <- glm(y1~1+M1,family='gamma')
if (family=='gaussian') mod <- lm(y1~1+M1)
if (family=='poisson')  mod <- glm(y1~1+M1,family='poisson')

# get the p-values for the targeted covariates. save it, the combinations
# of covariates and inclusion / exclusions factors that produced it, the resulting
# submatrices, and the regression model. 

s.D0 <- NULL
if (!is.null(D))  s.D0 <- s.D
s.C0 <- NULL
if (!is.null(C))  s.C0 <- s.C
y0 <- y1
D0 <- D1
C0 <- C1
Z0 <- Z1
S0 <- S1
N.S <- nrow(summary(mod)$coef)-(NCOL(S0)-1)
p0 <- summary(mod)$coef[N.S:nrow(summary(mod)$coef),4]
mod0 <- mod

# trace.

p.lookback[1] <- (1-min.fn)*max(p0)+min.fn*sum(p0^2)
if (trace==1)  {
	print(paste('Iteration 0'))
	print(s.D0)
	print(s.C0)
	print(p0)
	print(p.lookback[1])
	}

# Now, use (1+1) genetic algorithm search to investigate combinations of 
# covariates and inclusion / exclusion criteria.

for (it in 1:its)  {

	# adjust one of the dichotomous covariate model inclusion or 
	# observations inclusions / exclusions at random.

	s.D <- NULL
	if (!is.null(s.D0))  {
		i <- sample(1:NROW(s.D0),1)
		r <- sample(c('I','E0','E1','X'),1)
		s.D <- s.D0
		s.D[i] <- r
		}

	# adjust one of the continous covariate model inclusion or 
	# observations inclusions / exclusions at random.

	s.C <- NULL
	if (!is.null(s.C0))  {
		i <- sample(1:NROW(s.C0),1)
		r <- sample(c('I','T01','T05','X'),1)
		s.C <- s.C0
		s.C[i] <- r
		}
	
	# for dichotomous variables included in the model (cols), 
	# remove observations (rows) depending on 0/1 exclusions

	ind.row.D <- rep(TRUE,NROW(y))
	if (!is.null(s.D))  {
		ind.col.D <- s.D=='I' | s.D=='E0' | s.D=='E1'
		for (i in 1:ncol(D))  {	
			omit <- rep(FALSE,nrow(D))
			if (s.D[i]=='E0')  {
				omit <- D[,i] == 0
				}
			if (s.D[i]=='E1')  {
				omit <- D[,i] == 1
				}
			ind.row.D[omit] <- FALSE
			}
		}

	# for dichotomous variables included in the model (cols), 
	# remove observations (rows) depending on 1,5% tail exclusions

	ind.row.C <- rep(TRUE,NROW(y))
	if (!is.null(s.C))  {
		ind.col.C <- s.C=='I' | s.C=='T01' | s.C=='T05'
		for (i in 1:ncol(C))  {	
			omit <- rep(FALSE,nrow(C))
			if (s.C[i]=='T01')  {
				omit <- (C[,i] < cut.01[i]) | (C[,i] > cut.99[i])
				}
			if (s.C[i]=='T05')  {
				omit <- (C[,i] < cut.05[i]) | (C[,i] > cut.95[i])
				}
			ind.row.C[omit] <- FALSE
			}
		}

	# included observations (rows) are those that are included based on all criteria

	ind.row <- ind.row.D & ind.row.C

	# specify submatrices based on included covariates and observations, run the 
	# regression. If the regression fails, just set the p-value to 1.

	y1 <- y[ind.row]
	D1 <- D
	if (!is.null(D))  D1 <- D[ind.row,ind.col.D]
	C1 <- C
	if (!is.null(C))  C1 <- C[ind.row,ind.col.C]
	Z1 <- Z
	if (!is.null(Z))  Z1 <- Z[ind.row,]	
	S1 <- S[ind.row,]
	M1 <- cbind(D1,C1,Z1,S1)
	if (family=='binomial')  mod <- glm(y1~1+M1,family='binomial')
	if (family=='gamma')  mod <- glm(y1~1+M1,family='gamma')
	if (family=='gaussian') mod <- lm(y1~1+M1)
	if (family=='poisson')  mod <- glm(y1~1+M1,family='poisson')

	#print(s.D)
	#print(s.C)
	#print(dim(D1))
	#print(dim(C1))
	#print(dim(Z1))
	#print(dim(M1))
	#print(dim(S1))
	#print(N.S)
	#print(summary(mod)$coef)
	#print(dim(summary(mod)$coef))

	N.S <- nrow(summary(mod)$coef)-(NCOL(S1)-1)
	p <- rep(NCOL(S1),1)
	if (NROW(y1) > min.obs)  p <- summary(mod)$coef[N.S:nrow(summary(mod)$coef),4]
	p[is.nan(p)] <- 1
	
	# If the new objective function is less than the old, save it and the 
	# related information
	

	if ((1-min.fn)*max(p)+min.fn*sum(p^2) < (1-min.fn)*max(p0)+min.fn*sum(p0^2))  {
		p0 <- p
		s.D0 <- s.D
		s.C0 <- s.C
		y0 <- y1
		D0 <- D1
		C0 <- C1
		Z0 <- Z1
		S0 <- S1
		mod0 <- mod
		}

	# trace.

	p.lookback[it+1] <- (1-min.fn)*max(p0)+min.fn*sum(p0^2)
	if (trace==1)  {
		print(paste('Iteration',it))
		#print(s.D)
		#print(s.C)
		#print(p)
		print(s.D0)
		print(s.C0)
		print(p0)
		print(p.lookback[it+1])
		}

	if (it > stop.crit)  {
		if (p.lookback[it]==p.lookback[it-stop.crit]) break
		} 
	}

# specify the return list

ret <- list(p0=p0,s.D0=s.D0,s.C0=s.C0,y0=y0,D0=D0,C0=C0,Z0=Z0,S0=S0,mod0=mod0,
	    p.lookback=p.lookback)

ret
}
