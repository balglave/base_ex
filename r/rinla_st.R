#####################################
## Spatio-temporal models with R-INLA
#####################################
# B. Alglave based on Krainski et al. (2019)
# https://becarioprecario.bitbucket.io/spde-gitbook/ch-spacetime.html

## Packages and functions
library(INLA)
library(splancs)

source("r/functions/spde-book-files/R/spde-book-functions.R")

## Data
data(PRborder) # study region
data(PRprec) # location of the points
coords <- as.matrix(PRprec[sample(1:nrow(PRprec)), 1:2])

## Space time domain
k <- 12 # time steps
# domain and mesh
prdomain <- inla.nonconvex.hull(as.matrix(PRprec[, 1:2]),
                                convex = -0.03, concave = -0.05,
                                resolution = c(100, 100))
prmesh1 <- inla.mesh.2d(boundary = prdomain, 
                        max.edge = c(0.7, 0.7), cutoff = 0.35,
                        offset = c(-0.05, -0.05))

#------------------------------------------------------------------
#------------------- Simulation step ------------------------------
#------------------------------------------------------------------
## space-time latent field parameterization
params <- c(variance = 1, kappa = 1) 
rho <- 0.7
set.seed(1)
# k independent spatial components
x.k <- book.rspde(coords, range = sqrt(8) / params[2], # function for simulating spatial field
                  sigma = sqrt(params[1]), n = k, mesh = prmesh1,
                  return.attributes = TRUE)

# add the autoregressive part
x <- x.k
for (j in 2:k) x[, j] <- rho * x[, j - 1] + sqrt(1 - rho^2) * x.k[, j]

# plot simulation
x11()
par(mfrow = c(5, 3), mar = c(0, 0, 0.7, 0))

# Values for scaling
x.min <- min(as.vector(x))
x.max <- max(as.vector(x))
x.range <- x.max - x.min

c100 <- book.color.c(101)
for (j in 1:k) {
  cols <- c100[1 + round(100 * (x[, j] - x.min)) / x.range ]
  plot(coords, col = cols, axes = FALSE, asp = 1, pch = 19, cex = 0.5,
       main = paste0("Time: ", j))
}

# add a categorial variable
n <- nrow(coords)
set.seed(2)
ccov <- factor(sample(LETTERS[1:3], n * k, replace = TRUE))
beta <- -1:1 # coefficient

# compute response variable
sd.y <- 0.1 # observation variance
y <- beta[unclass(ccov)] + x + rnorm(n * k, 0, sd.y) # covariate effect + space-time random effect + white noise with variance sd.y

# create the data dataframe
isel <- sample(1:(n * k), n * k / 2) # select datapoints to be included in model fitting
dat <- data.frame(y = as.vector(y), w = ccov, 
                  time = rep(1:k, each = n), 
                  xcoo = rep(coords[, 1], k), 
                  ycoo = rep(coords[, 2], k))[isel, ]

spde <- inla.spde2.pcmatern(mesh = prmesh1, 
                            prior.range = c(0.5, 0.01), # P(range < 0.05) = 0.01
                            prior.sigma = c(1, 0.01)) # P(sigma > 1) = 0.01

#----------------------------------------------------------
##-------------- Data stack preparation -------------------
#----------------------------------------------------------
## SPDE objects
spde <- inla.spde2.pcmatern(mesh = prmesh1, 
                            prior.range = c(0.5, 0.01), # P(range < 0.05) = 0.01
                            prior.sigma = c(1, 0.01)) # P(sigma > 1) = 0.01

iset <- inla.spde.make.index('i', n.spde = spde$n.spde, # index set for time component (specific to space-time analysis)
                             n.group = k)

A <- inla.spde.make.A(mesh = prmesh1, # matrix design (relate mesh to datapoints)
                      loc = cbind(dat$xcoo, dat$ycoo), group = dat$time) 

sdat <- inla.stack(
  data = list(y = dat$y), 
  A = list(A, 1), 
  effects = list(iset, w = dat$w), # index set and categorical variable
  tag = 'stdata') 

#----------------------------------------------------------
##-------------- Fitting the model ------------------------
#----------------------------------------------------------
## Configuration
# PC prior
h.spec <- list(rho = list(prior = 'pc.cor1', param = c(0, 0.9)))

# Model formula
formulae <- y ~ 0 + w + f(i, model = spde, group = i.group, 
                          control.group = list(model = 'ar1', hyper = h.spec)) 

# PC prior on the autoreg. param.
prec.prior <- list(prior = 'pc.prec', param = c(1, 0.01))

## Model fitting
res <- inla(formulae,  data = inla.stack.data(sdat), 
            control.predictor = list(compute = TRUE,
                                     A = inla.stack.A(sdat)), 
            control.family = list(hyper = list(prec = prec.prior)), 
            control.fixed = list(expand.factor.strategy = 'inla'))

## Model parameters - estimated vs. simulated values
x11()
par(mfrow = c(2, 2), mar = c(3, 3, 1, 0.1), mgp = 2:0)
for (j in 1:4) {
  plot(res$marginals.hyper[[j]], type = 'l', 
       xlab = names(res$marginals.hyper)[j], ylab = 'Density')
  abline(v = c(1 / sd.y^2, sqrt(8) / params[1], params[2]^0.5, rho)[j],
         col = 2)
}

## Plot latent field values
stepsize <- 4 * 1 / 111
nxy <- round(c(diff(range(coords[, 1])), 
               diff(range(coords[, 2]))) / stepsize)
projgrid <- inla.mesh.projector( # projection grid
  prmesh1, xlim = range(coords[, 1]), 
  ylim = range(coords[, 2]), dims = nxy)

xmean <- list() # compute predictions
for (j in 1:k){
  xmean[[j]] <- inla.mesh.project(
    projgrid, res$summary.random$i$mean[iset$i.group == j])
}

xy.in <- inout(projgrid$lattice$loc, 
               cbind(PRborder[, 1], PRborder[, 2]))

x11()
par(mfrow = c(4,3), mar = c(0, 0, 0, 0))
for (j in 1:k) {
  xmean[[j]][!xy.in] <- NA
  book.plot.field(list(x = projgrid$x, y = projgrid$y, z = xmean[[j]]),
                  zlim = round(range(unlist(xmean), na.rm = TRUE), 1) )
}

#----------------------------------------------------------
##--------------------- Validation ------------------------
#----------------------------------------------------------
## Make validation dataset
vdat <- data.frame(r = as.vector(y), w = ccov,
                   t = rep(1:k, each = n), x = rep(coords[, 1], k),
                   y = rep(coords[, 2], k))
vdat <- vdat[-isel, ]

## Make projection matrix and data stack for validation data
Aval <- inla.spde.make.A(prmesh1, 
                         loc = cbind(vdat$x, vdat$y), group = vdat$t) 
stval <- inla.stack(
  data = list(y = NA), # NA: no data, only enable predictions
  A = list(Aval, 1), 
  effects = list(iset, w = vdat$w),
  tag = 'stval')  

## join fit and validation stack 
stfull <- inla.stack(sdat, stval)

## Re-run simulation
vres <- inla(formulae,  data = inla.stack.data(stfull), 
             control.predictor = list(compute = TRUE,
                                      A = inla.stack.A(stfull)), 
             control.family = list(hyper = list(prec = prec.prior)), 
             control.fixed = list(expand.factor.strategy = 'inla'), 
             control.mode = list(theta = res$mode$theta, restart = FALSE))

## Plot observed values vs. predicted
ival <- inla.stack.index(stfull, 'stval')$data 
par(mfrow = c(1, 1), mar = c(3, 3, 0.5, 0.5), mgp = c(1.75, 0.5, 0))
plot(vres$summary.fitted.values$mean[ival], vdat$r, asp = 1,
     xlab = 'Posterior mean', ylab = 'Observed') 
abline(0:1, col = gray(0.7)) 
