## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
knitr::opts_chunk$set(fig.width = 6, fig.height=6)
# knitr::opts_chunk$set(tidy.opts=list(width.cutoff=80), tidy=TRUE)

## ----setup--------------------------------------------------------------------
library(imprecise101)

## -----------------------------------------------------------------------------
op <- idm(nj=1, s=2, N=6, k=4)
c(op$p.lower, op$p.upper)

## -----------------------------------------------------------------------------
op <- idm(nj=1, s=1, N=6, k=4)
round(c(op$p.lower, op$p.upper),3)

## -----------------------------------------------------------------------------
op <- idm(nj=1, s=0, N=6, k=4)
round(c(op$p.lower, op$p.upper),3)

## -----------------------------------------------------------------------------
r1 <- c(idm(nj=1, s=1, N=6, k=4)$p, 
        idm(nj=1, s=2, N=6, k=4)$p, 
        idm(nj=1, s=4/2, N=6, k=4)$p, 
        idm(nj=1, s=4, N=6, k=4)$p) # Omega 1

r2 <- c(idm(nj=1, s=1, N=6, k=2)$p, 
        idm(nj=1, s=2, N=6, k=2)$p, 
        idm(nj=1, s=2/2, N=6, k=2)$p, 
        idm(nj=1, s=2, N=6, k=2)$p) # Omega 2 

r3 <- c(idm(nj=1, s=1, N=6, k=3, cA=2)$p, 
        idm(nj=1, s=2, N=6, k=3, cA=2)$p, 
        idm(nj=1, s=3/2, N=6, k=3, cA=2)$p, 
        idm(nj=1, s=3, N=6, k=3, cA=2)$p) # Omega 3

tb1 <- rbind(r1, r2, r3)
rownames(tb1) <- c("Omega1", "Omega2", "Omega3")
colnames(tb1) <- c("s=1", "s=2", "s=k/2", "s=k")
round(tb1,3)

## -----------------------------------------------------------------------------
mat <- cbind(
  unlist(pbetabinom(M=6, x=1, s=1, N=6, y=0)),
  unlist(pbetabinom(M=6, x=1, s=1, N=6, y=1)),
  unlist(pbetabinom(M=6, x=1, s=1, N=6, y=2)),
  unlist(pbetabinom(M=6, x=1, s=1, N=6, y=3)),
  unlist(pbetabinom(M=6, x=1, s=1, N=6, y=4)),
  unlist(pbetabinom(M=6, x=1, s=1, N=6, y=5)),
  unlist(pbetabinom(M=6, x=1, s=1, N=6, y=6))
)
colnames(mat) <- c("y=0", "y=1", "y=2", "y=3", "y=4", "y=5", "y=6")
round(mat, 3)

## -----------------------------------------------------------------------------
op <- idm(nj=1, s=2, N=6, k=4)
round(c(op$p.upper, op$p.lower, op$s.upper, op$s.lower),3)

op <- idm(nj=1, s=1, N=6, k=4)
round(c(op$p.upper, op$p.lower, op$s.upper, op$s.lower),3)

## -----------------------------------------------------------------------------
round(hpd(alpha=3, beta=5, p=0.95),4) # s=2
round(hpd(alpha=3, beta=5, p=0.90),4) # s=2
round(hpd(alpha=3, beta=5, p=0.50),4) # s=2 (required for message of failure)
round(hpd(alpha=2, beta=5, p=0.95),4) # s=1
round(hpd(alpha=2, beta=5, p=0.90),4) # s=1
round(hpd(alpha=2, beta=5, p=0.50),4) # s=1

## -----------------------------------------------------------------------------
x <- pscl::betaHPD(alpha=2, beta=6, p=0.95, plot=FALSE)
round(x,4)

## -----------------------------------------------------------------------------
fn <- function(x) choose(7,1)*(1-x)^6
integrate(f=fn, lower=1/2, upper=1)$value

fn <- function(x) dbeta(x, 3, 5)
integrate(f=fn, lower=1/2, upper=1)$value 

## ----out.width="45%"----------------------------------------------------------
x <- seq(-0.99, 0.99, 0.02)
ymax <- ymin <- numeric(length(x))
for(i in 1:length(x)) ymin[i] <- dbetadif(x=x[i], a1=9,b1=2,a2=8,b2=4)
for(i in 1:length(x)) ymax[i] <- dbetadif(x=x[i], a1=11,b1=0.01,a2=6,b2=6)

plot(x=x, y=cumsum(ymin)/sum(ymin), type="l", ylab="F(z)", xlab="z",
     main=expression(paste("Fig 1. Posterior upper and lower CDFs for ", psi, 
                           "=", theta[e]-theta[c])))
points(x=x, y=cumsum(ymax)/sum(ymax), type="l")

## ----out.width="45%", echo=FALSE----------------------------------------------
tc <- seq(0,1,0.1)
s <- 2

op <- ibm(n=10, m=6, showplot=TRUE, xlab1="z", main1=expression(paste("Fig 1. Posterior ", theta[c], " based on ", s==2)))

op <- ibm(n=9, m=9, showplot=TRUE, xlab1="z", main1=expression(paste("Fig 2. Posterior ", theta[e], " based on ", s==2)))

## ----out.width="45%", echo=FALSE----------------------------------------------
x <- seq(-0.99, 0.99, 0.02)
ymax <- ymin <- numeric(length(x))
for(i in 1:length(x)) ymin[i] <- dbetadif(x=x[i], a1=9,b1=2,a2=8,b2=4)
for(i in 1:length(x)) ymax[i] <- dbetadif(x=x[i], a1=11,b1=0.01,a2=6,b2=6)

plot(x=x, y=1-cumsum(ymin)/sum(ymin), type="l", ylab="G(z)", xlab="z", main=expression(paste("Fig 3. ", psi, "=", theta[e]-theta[c], " Walley 1996 p.496" )))
points(x=x, y=1-cumsum(ymax)/sum(ymax), type="l")

dmin <- 1-cumsum(ymin/sum(ymin))
dmax <- 1-cumsum(ymax/sum(ymax))

