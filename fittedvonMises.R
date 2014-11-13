####################################################################
####################################################################
library(circular)
#install.packages('movMF')
library(movMF)
#setwd('/Users/lun5/Research/github/HE-tumor-object-segmentation/results')
setwd('C:/Users/luong_nguyen/Documents/GitHub/HE-tumor-object-segmentation/results')
data.vm <- read.csv('gland3_snip_opp.csv',sep =',', header=F);dim(data.vm)
data.vm <- read.csv('gland3_snip_lch.csv',sep =',', header=F);dim(data.vm)
data.vm <- read.csv('gland3_snip_oppNat.csv',sep =',', header=F);dim(data.vm)

data.vm <- t(as.circular(data.vm))

dist.matrix <- dist.circular(data.vm, method = "chord", diag = FALSE, upper = FALSE)

data.cart <- cbind(cos(data.vm),sin(data.vm))
# fit this using skmeans
#y3 <- skmeans(data.cart,2)
# fit this using movMF package

num.comps <- 2
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
plot(data.vm,stack=TRUE, bins=150, shrink=1.5,xlim=c(-1.5,1.5),
ylim=c(-1.5,1.5),se=5e-3,col='gray')

set.seed(13)
x1 <- rvonmises(n=100, mu=circular(coord2rad(mu.fitted[1,1],mu.fitted[1,2])), 
kappa=kappa.fitted[1])
res25x1 <- density(x1, bw=25)
lines(res25x1, points.plot=F, col = 'purple')

#par(new = T)
set.seed(135)
x2 <- rvonmises(n=100, mu=circular(coord2rad(mu.fitted[2,1],mu.fitted[2,2])), 
kappa=kappa.fitted[2])
res25x2 <- density(x2, bw=25)
lines(res25x2, points.plot=F, col = 'pink')

set.seed(1113)
x3 <- rvonmises(n=100, mu=circular(coord2rad(mu.fitted[3,1],mu.fitted[3,2])), 
kappa=kappa.fitted[3])
res25x3 <- density(x2, bw=25)
lines(res25x3, points.plot=F, col = 1)

dev.copy(tiff,paste(getwd(),'/results/theta_vonMises.tiff',sep = ''))
dev.off()
###
