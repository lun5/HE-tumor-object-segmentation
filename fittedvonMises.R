####################################################################
####################################################################
library(circular)
#install.packages('movMF')
library(movMF)
setwd('/Users/lun5/Research/github/HE-tumor-object-segmentation/results')
#setwd('C:/Users/luong_nguyen/Documents/GitHub/HE-tumor-object-segmentation/results')
data.vm <- read.csv('gland3_snip_opp.csv',sep =',', header=F);dim(data.vm)
# data.vm <- read.csv('gland3_snip_lch.csv',sep =',', header=F);dim(data.vm)
# data.vm <- read.csv('gland3_snip_oppNat.csv',sep =',', header=F);dim(data.vm)
# data.vm <- read.csv('tp10-611_gland1snip_theta_oppCol.csv',sep =',', header=F);dim(data.vm)
data.vm <- read.csv('tp10-867-1gland21_snip.csv',sep =',', header=F);dim(data.vm)

data.vm <- t(as.circular(data.vm))

dist.matrix <- dist.circular(data.vm, method = "chord", diag = FALSE, upper = FALSE)

data.cart <- cbind(cos(data.vm),sin(data.vm))
# fit this using skmeans
#y3 <- skmeans(data.cart,2)
# fit this using movMF package

num.comps <- 3
y3 <- movMF(data.cart,num.comps, nruns = 10)
y3cv <- movMF(data.cart, num.comps, nruns = 10, kappa = list(common = TRUE)) # same kappa
y3cf <- movMF(data.cart, num.comps, nruns = 10, kappa = 4) # given some kappa
sapply(list(y3, y3cf, y3cv), BIC) # compare the three fits

# fitted values
kappa.fitted <- apply(y3$theta, 1, function(x) norm(x,'2')); kappa.fitted
mu.fitted <- t(apply(y3$theta, 1, function(x) x/norm(x,'2'))); mu.fitted
angles.fitted <- coord2rad(mu.fitted); angles.fitted
# plot multiple distributions together
dev.new()
plot.range <- 1.5
circ.hist <- plot(data.vm,stack=TRUE, bins=150, shrink=1.5,
xlim=c(-plot.range,plot.range),ylim=c(-plot.range,plot.range),se=5e-3,col='gray')

set.seed(13)
x1 <- rvonmises(n=100, mu=circular(coord2rad(mu.fitted[1,1],mu.fitted[1,2])), 
kappa=kappa.fitted[1])
res25x1 <- density(x1, bw=25)
lines(res25x1, points.plot=F, col = 'pink')

#par(new = T)
set.seed(135)
x2 <- rvonmises(n=100, mu=circular(coord2rad(mu.fitted[2,1],mu.fitted[2,2])), 
kappa=kappa.fitted[2])
res25x2 <- density(x2, bw=25)
lines(res25x2, points.plot=F, col = 'purple')

set.seed(1113)
x3 <- rvonmises(n=100, mu=circular(coord2rad(mu.fitted[3,1],mu.fitted[3,2])), 
kappa=kappa.fitted[3])
res25x3 <- density(x3, bw=25)
lines(res25x3, points.plot=F, col = 'pink')

dev.copy(tiff,paste(getwd(),'/results/theta_vonMises.tiff',sep = ''))
dev.off()

####################################################################
####################################################################
# fit kde for spherical data
angle.pairs <- read.csv('gland3_snip_samplepairs.csv',sep =',', header=F);
dim(angle.pairs)
angle.pairs <- as.circular(angle.pairs)
density(angle.pairs,bw=25)

data(coope)
res <- lsfit.circle(x=x.coope, y=y.coope)
res
plot(res)
par(mfcol=c(1,2))
plot(res$angles)
hist(res$radius)
plot(circular(0), type="n", xlim=c(-5.2, 5.2), ylim=c(-5.2, 5.2),
xlab="The Radius of the circle \n is measured from the base line of the axes.")
lines(x=res$angles, y=res$radius, join=TRUE, type="b")
ff <- function(x) sqrt((res$coefficients[1]*cos(x))^2+(res$coefficients[1]*sin(x))^2)
curve.circular(ff, add=TRUE, join=TRUE, nosort=FALSE, col=2)
windrose(x=res$angles, y=res$radius)
