## Generate and fit a "small-mix" data set a la Banerjee et al.
# install.packages('CircStats')
# library(CircStats)
# Generate 100 observations from a von Mises distribution.
# with mean direction 0 and concentration 3.
# data.vm <- rvm(100, 0, 3)
# Plot data set. All points do not fit on plot.
# circ.plot(data.vm, stack=TRUE, bins=150)
# Shrink the plot so that all points fit.
# circ.plot(data.vm, stack=TRUE, bins=150, shrink=1.5)

install.packages('circular'); 
library(circular);
x <- rvonmises(20, circular(0), 10)
y <- runif(20, 0.5, 1)
plot(x, shrink=2)
lines(x, y)

dev.set(which = dev.next())
# generate mixed von Mises
data.vm <- rmixedvonmises(n=100, mu1 = circular(0), mu2=circular(pi/3),
kappa1=3, kappa2 = 3, prop = 0.4)
data.cart <- cbind(cos(data.vm),sin(data.vm))
dev.new()
plot(data.vm, stack=TRUE, bins=150)
dev.new()
plot(data.cart[,1], data.cart[,2], xlim = c(-1,1), ylim = c(-1,1), asp = 1)
lines(data.cart[,1], data.cart[,2], join = FALSE, nosort = FALSE, offset=1, shrink=1,
plot.info = NULL, zero = NULL, rotation = NULL)

library(movMF)
mu <- rbind(c(-0.251, -0.968),
c(0.399, 0.917))
kappa <- c(4, 4)
theta <- kappa * mu
theta
alpha <- c(0.48, 0.52)
## Generate a sample of size n = 50 from the von Mises-Fisher mixture
## with the above parameters.
set.seed(123)
x <- rmovMF(50, theta, alpha)
## Fit a von Mises-Fisher mixture with the "right" number of components,
## using 10 EM runs.
y2 <- movMF(x, 2, nruns = 10)
## Inspect the fitted parameters:
y2
## Compare the fitted classes to the true ones:
table(True = attr(x, "z"), Fitted = predict(y2))
## To use a common kappa:
y2cv <- movMF(x, 2, nruns = 10, kappa = list(common = TRUE))
## To use a common kappa fixed to the true value of 4:
y2cf <- movMF(x, 2, nruns = 10, kappa = 4)
## Comparing solutions via BIC:
sapply(list(y2, y2cf, y2cv), BIC)

library(circular)
# Generate 100 observations from a von Mises distribution.
# with mean direction 0 and concentration 3.
data.vm <- rvonmises(n=100, mu=circular(0), kappa=3)
# Plot data set. All points do not fit on plot.
plot(data.vm, stack=TRUE, bins=150)
# Shrink the plot so that all points fit.
plot(data.vm, stack=TRUE, bins=150, shrink=1.5)
# Recentering the figure in a different place

# fit this using movMF package - bummer, only >= 2D. Craps
# fit this using skmeans
#y3 <- skmeans(data.cart,2)
data.cart <- cbind(cos(data.vm),sin(data.vm))
y3 <- movMF(data.cart,2, nruns = 100)
y3cv <- movMF(data.cart, 2, nruns = 100, kappa = list(common = TRUE))
y3cf <- movMF(data.cart, 2, nruns = 10, kappa = 3)

x <- cbind(x1 = 3, x2 = c(4:1, 2:5))
dimnames(x)[[1]] <- letters[1:8]
apply(x, 2, mean, trim = .2)
col.sums <- apply(x, 2, sum)
row.sums <- apply(x, 1, sum)
rbind(cbind(x, Rtot = row.sums), Ctot = c(col.sums, sum(col.sums)))

## plot circular density
dev.new()
set.seed(1234)
x <- rvonmises(n=100, mu=circular(pi), kappa=2)
res25x <- density(x, bw=25)
plot(res25x, points.plot=TRUE, xlim=c(-1.5,1))
plot(x,stack=TRUE, bins=150, shrink=1.5,xlim=c(-1.5,1))
res50x <- density(x, bw=25, adjust=2)
lines(res50x, col=2) # higher bandwidth supposed to make the kde smoother

dev.new()
resp25x <- plot(res25x, points.plot=F, xlim=c(-1.3, 1.3), ylim=c(-1.5,1.5),bins = 150,
template="none", main="Plotting density estimate for two data set", shrink = 1)
par(new=T)
plot(x,stack=TRUE, bins = 150, shrink = 1.5)

y <- rvonmises(n=100, mu=circular(pi/2), kappa=2,
control.circular=list(template="none"))
res25y <- density(y, bw=25)
lines(res25y, points.plot=TRUE, plot.info=resp25x, col=2, points.col=2)

dev.new()
plot(res25x, plot.type="line", points.plot=TRUE, xlim=c(-1, 1.3), ylim=c(-1.5,1.2),
template="geographics", main="Plotting density estimate for two data set")
lines(res25y, plot.type="line", points.plot=TRUE, col=2, points.col=2)

# plot empirical distribution
data1 <- rvonmises(n=10, mu=circular(0), kappa=3)
data2 <- rvonmises(n=10, mu=circular(0), kappa=1)
dev.new()
plot.edf(data1, xlab="Data", ylab="EDF", main="Plots of Two EDF's")
lines.edf(data2, lty=2, col=2)


x <- rvonmises(n=50, mu=circular(0), kappa=5)
mle.vonmises(x) # estimation of mu and kappa
mle.vonmises(x, mu=circular(0)) # estimation of kappa only

set.seed(1234)
x <- cbind(rnorm(20), rnorm(20))
y <- coord2rad(x)

####################################################################
####################################################################
setwd('/Users/lun5/Research/github/HE-tumor-object-segmentation/')
library(circular)
#install.packages('movMF')
library(movMF)
mu <- rbind(c(-0.251, -0.968),
c(0.399, 0.917))
kappa <- c(4, 4)
theta <- kappa * mu
theta
alpha <- c(0.48, 0.52)
## Generate a sample of size n = 50 from the von Mises-Fisher mixture
## with the above parameters.
set.seed(123)
x <- rmovMF(50, theta, alpha)
## Fit a von Mises-Fisher mixture with the "right" number of components,
## using 10 EM runs.
y2 <- movMF(x, 2, nruns = 10)
## Inspect the fitted parameters:
y2
## Compare the fitted classes to the true ones:
table(True = attr(x, "z"), Fitted = predict(y2))
## To use a common kappa:
y2cv <- movMF(x, 2, nruns = 10, kappa = list(common = TRUE))
## To use a common kappa fixed to the true value of 4:
y2cf <- movMF(x, 2, nruns = 10, kappa = 4)
## Comparing solutions via BIC:
sapply(list(y2, y2cf, y2cv), BIC)

# generate mixed von Mises
mu.rad <- coord2rad(mu)
data.vm <- rmixedvonmises(n=500, mu1 = circular(mu.rad[1]), 
mu2 = circular(mu.rad[2]),kappa1= kappa[1], kappa2 = kappa[2], prop = alpha[1])
data.cart <- cbind(cos(data.vm),sin(data.vm))
# fit this using skmeans
#y3 <- skmeans(data.cart,2)
# fit this using movMF package
num.comps <-3
y3 <- movMF(data.cart,num.comps, nruns = 10)
y3cv <- movMF(data.cart, num.comps, nruns = 10, kappa = list(common = TRUE))
y3cf <- movMF(data.cart, num.comps, nruns = 10, kappa = 4)
sapply(list(y3, y3cf, y3cv), BIC)

kappa.fitted <- cbind(norm(y3$theta[1,],'2'),norm(y3$theta[2,],'2'),
norm(y3$theta[3,],'2'));kappa.fitted
mu.fitted = rbind(y3$theta[1,]/kappa.fitted[1], y3$theta[2,]/kappa.fitted[2],
 y3$theta[3,]/kappa.fitted[3]); mu.fitted

# plot 2 distribution together
dev.new()
plot(data.vm,stack=TRUE, bins=500, shrink=1.5,xlim=c(-2.5,2.5))

set.seed(13)
x1 <- rvonmises(n=100, mu=circular(coord2rad(mu.fitted[1,1],mu.fitted[1,2])), 
kappa=kappa.fitted[1])
res25x1 <- density(x1, bw=25)
lines(res25x1, points.plot=F, col = 4)

#par(new = T)
x2 <- rvonmises(n=100, mu=circular(coord2rad(mu.fitted[2,1],mu.fitted[2,2])), 
kappa=kappa.fitted[2])
res25x2 <- density(x2, bw=25)
lines(res25x2, points.plot=F, col = 2)

x3 <- rvonmises(n=100, mu=circular(coord2rad(mu.fitted[3,1],mu.fitted[3,2])), 
kappa=kappa.fitted[3])
res25x3 <- density(x2, bw=25)
lines(res25x3, points.plot=F, col = 3)

###
data.vm <- read.csv('theta.csv',sep =',', header=F)
data.vm <- t(as.circular(data.vm))
