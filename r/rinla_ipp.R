#############################
## Point process with R-INLA
############################
# Based on Krainski et al. (2019)
# https://becarioprecario.bitbucket.io/spde-gitbook/ch-lcox.html

# Set of functions to ease the use of INLA
source("r/functions/spde-book-files/R/spde-book-functions.R")
library(INLA)
library(rgeos)
library(spatstat)

#------------
## Simulation
#------------
## Spatial domain
win <- owin(c(0, 3), c(0, 3))
npix <- 300
spatstat.options(npixel = npix)

## Parameters value
beta0 <- 5
# exp(beta0) * diff(range(win$x)) * diff(range(win$y)) # expected value of points
sigma2x <- 0.5
range <- 3
nu <- 1

## Simulate Log-Gaussian Cox Process
library(RandomFields)
set.seed(1)
lg.s <- rLGCP('matern', beta0, var = sigma2x,
              scale = range / sqrt(8), nu = nu, win = win)
xy <- cbind(lg.s$x, lg.s$y)[, 2:1] # coord of the observed events
# (n <- nrow(xy)) # number of points simulated

Lam <- attr(lg.s, 'Lambda')
rf.s <- log(Lam$v)
summary(as.vector(rf.s))

# plot
book.plot.field(list(x = Lam$yrow, y = Lam$xcol, z = rf.s), 
                xlim=c(0,3), ylim=c(0,3))
points(xy, pch = 19) 

#-------------------------------------
## Fitting the model with NO covariate
#-------------------------------------
## Build the mesh
loc.d <- 3 * cbind(c(0, 1, 1, 0, 0), c(0, 0, 1, 1, 0))
mesh <- inla.mesh.2d(loc.domain = loc.d, offset = c(0.3, 1), 
                     max.edge = c(0.3, 0.7), cutoff = 0.05)
nv <- mesh$n

# plot mesh
par(mar = c(0, 0, 0, 0))
plot(mesh, asp = 1, main = '')
points(xy, col = 4, pch = 19)
lines(loc.d, col = 3)

spde <- inla.spde2.pcmatern(mesh = mesh,
                            # PC-prior on range: P(practic.range < 0.05) = 0.01
                            prior.range = c(0.05, 0.01),
                            # PC-prior on sigma: P(sigma > 1) = 0.01
                            prior.sigma = c(1, 0.01)) 

## Get weights for estimating the point process (integrating concstant - Simpson et al. (2016))
dmesh <- book.mesh.dual(mesh)
domain.polys <- Polygons(list(Polygon(loc.d)), '0')
domainSP <- SpatialPolygons(list(domain.polys))

w <- sapply(1:length(dmesh), function(i) {
  if (gIntersects(dmesh[i, ], domainSP))
    return(gArea(gIntersection(dmesh[i, ], domainSP)))
  else return(0)
})

# Plot integration points
par(mar = c(2, 2, 1, 1), mgp = 2:0)
plot(mesh$loc, asp = 1, col = (w == 0) + 1, pch = 19, xlab = '', ylab = '') 
plot(dmesh, add = TRUE)
lines(loc.d, col = 3)

## shape data
n <- nrow(xy)
y.pp <- rep(0:1, c(nv, n))
e.pp <- c(w, rep(0, n)) 
imat <- Diagonal(nv, rep(1, nv))
lmat <- inla.spde.make.A(mesh, xy)
A.pp <- rbind(imat, lmat)

stk.pp <- inla.stack(
  data = list(y = y.pp, e = e.pp), 
  A = list(1, A.pp),
  effects = list(list(b0 = rep(1, nv + n)), list(i = 1:nv)),
  tag = 'pp')


## Fit model
pp.res <- inla(y ~ 0 + b0 + f(i, model = spde), 
               family = 'poisson', data = inla.stack.data(stk.pp), 
               control.predictor = list(A = inla.stack.A(stk.pp)), 
               E = inla.stack.data(stk.pp)$e)

## Plot parameters posterior
par(mfrow = c(1, 3), mar = c(3, 3, 1, 0.3), mgp = c(2, 1, 0)) 
plot(pp.res$marginals.fixed[[1]], type = 'l', xlab = expression(beta[0]),
     ylab = 'Density')
abline(v = beta0, col = 2)
plot(pp.res$marginals.hyperpar[[2]], type = 'l', xlab = expression(sigma),
     ylab = 'Density', xlim = c(0,2))
abline(v = sqrt(sigma2x), col = 2)
plot(pp.res$marginals.hyperpar[[1]], type = 'l', xlab = 'Nominal range',
     ylab = 'Density', xlim = c(0, 8))
abline(v = range, col = 2)


#----------------------------------
## Fitting the model with covariate
#----------------------------------
# simulate LGCP as a function of a covariate and latent field
# Use expanded range
x0 <- seq(min(mesh$loc[, 1]), max(mesh$loc[, 1]), length = npix)
y0 <- seq(min(mesh$loc[, 2]), max(mesh$loc[, 2]), length = npix)
gridcov <- outer(x0, y0, function(x,y) cos(x) - sin(y - 2))

set.seed(1)
beta1 <- 0.5
lg.s.c <- rLGCP('matern', im(beta0 + beta1 * gridcov, xcol = x0,
                             yrow = y0), var = sigma2x, scale = range / sqrt(8), 
                nu = 1, win = win)
xy.c <- cbind(lg.s.c$x, lg.s.c$y)[, 2:1]
n.c <- nrow(xy.c)

par(mfrow = c(1, 2), mar = c(2, 1, 0.5, 4.5), mgp = c(1, 0.5, 0))
book.plot.field(list(x = x0, y = y0, z = gridcov), xlim = c(0, 3),
                ylim = c(0, 3))
title(main = "Covariate")
book.plot.field(list(x = x0, y = y0, z = log(attr(lg.s.c, 'Lambda')$v)),
                xlim = c(0, 3), ylim = c(0, 3))  
points(xy.c, pch = 19)
title(main = "Latent field and points")

covariate.im <- im(gridcov, x0, y0)
covariate <- interp.im(covariate.im, 
                       x = c(mesh$loc[, 1], xy.c[, 1]),
                       y = c(mesh$loc[, 2], xy.c[, 2]))

y.pp.c <- rep(0:1, c(nv, n.c))
e.pp.c <- c(w, rep(0, n.c))

lmat.c <- inla.spde.make.A(mesh, xy.c)

A.pp.c <- rbind(imat, lmat.c)

stk.pp.c <- inla.stack(
  data = list(y = y.pp.c, e = e.pp.c), 
  A = list(1, A.pp.c), 
  effects = list(list(b0 = 1, covariate = covariate), 
                 list(i = 1:nv)),
  tag = 'pp.c')

## Fit model
pp.c.res <- inla(y ~ 0 + b0 + covariate + f(i, model = spde),
                 family = 'poisson', data = inla.stack.data(stk.pp.c), 
                 control.predictor = list(A = inla.stack.A(stk.pp.c)), 
                 E = inla.stack.data(stk.pp.c)$e)

## Parameters estimates
par(mfrow = c(2, 2), mar = c(3, 3, 1, 0.3), mgp = c(2, 1, 0)) 
plot(pp.c.res$marginals.fix[[1]], type = 'l', ylab = 'Density', 
     xlab = expression(beta[0]))
abline(v = beta0, col = 2)
plot(pp.c.res$marginals.fix[[2]], type = 'l', ylab = 'Density', 
     xlab = expression(beta[1]))
abline(v = beta1, col = 2)
plot(pp.c.res$marginals.hyperpar[[2]], type = 'l', ylab = 'Density', 
     xlab = expression(sigma), xlim = c(0, 2))
abline(v = sqrt(sigma2x), col = 2)
plot(pp.c.res$marginals.hyperpar[[1]], type = 'l', ylab = 'Density',
     xlab = "Spatial range", xlim = c(0, 10))
abline(v = range, col = 2)

