#####################################################################
# Data: NY8_utm18 TCE NY_nb.gal NY8cities Sy_GeoDa1.GAL Sy_GeoDa4.GAL
#####################################################################

# Set up: clear all and set work directory 
rm(list = ls())
# setwd("")

library(rgdal)
library(spdep)

NY8 <- readOGR(".", "NY8_utm18")
TCE <- readOGR(".", "TCE")
cities <- readOGR(".", "NY8cities")

##########################################################
# Let's look at how the data is stored

# How is the data stored?
class(NY8)

# It's a SpatialPolygonsDataFrame - a data.frame on polygon
getClass("SpatialPolygonsDataFrame")

# Spatial polygons contains polygons and Spatial characteristics
getClass("SpatialPolygons")

# a polygon is a sequence of closed lines;
# points coordinates where the first point equals the last
getClass("Polygon")

# labpt - label point, centroid of polygon 
# a line is an ordered list of coordinates 
getClass("Line")

# finally.. Spatial is the mother class of all Spatial classes
# used in the sp package
getClass("Spatial")

##########################################################
# What does the data contain?
summary(NY8)

##########################################################
# What about the cities data?
summary(cities)

##########################################################
# Let's make some plots
par(mfrow=c(1,2))
plot(NY8, border="grey60", axes=TRUE)
text(coordinates(cities), labels=as.character(cities$names), font=2, cex=0.9)
text(bbox(NY8)[1,1], bbox(NY8)[2,2], labels="a)", cex=0.8)
plot(NY8, border="grey60", axes=TRUE)
points(TCE, pch=1, cex=0.7)
points(TCE, pch=3, cex=0.7)
text(coordinates(TCE), labels=as.character(TCE$name), cex=0.7,
 font=1, pos=c(4,1,4,1,4,4,4,2,3,4,2), offset=0.3)
text(bbox(NY8)[1,1], bbox(NY8)[2,2], labels="b)", cex=0.8)

# let's plot one of the features - percent age > 65
spplot(NY8, c("PCTAGE65P"))
spplot(NY8, c("PCTAGE65P"), col="transparent")

# Let's make a different plot, with a new color palette

library("RColorBrewer")

#color palette creator function
rds <- colorRampPalette(brewer.pal(8, "RdBu"))
#get a range for the values
tr_at <- seq(min(NY8$PCTAGE65P), max(NY8$PCTAGE65P), length.out=20)
#create a color interpolating function taking the required
#number of shades as argument
tr_rds <- rds(20)
#parameters
# at - at which values colors change
# col.regions - specify fill colors 
tr_pl <- spplot(NY8, c("PCTAGE65P"), at=tr_at, col="transparent",
                col.regions=tr_rds, main=list(label="Age>65", cex=0.8))
plot(tr_pl)

# reads a GAL lattice file into a neighbors list
NY_nb <- read.gal("NY_nb.gal", region.id=row.names(NY8))

summary(NY_nb) #which states are neighbors?

par(mfrow=c(1,1))
plot(NY8, border="grey60", axes=TRUE)
plot(NY_nb, coordinates(NY8), pch=19, cex=0.6, add=TRUE)

# vignette
vignette("nb", package = "spdep")

library(spdep)

Syracuse <- NY8[NY8$AREANAME == "Syracuse city",]
Sy0_nb <- subset(NY_nb, NY8$AREANAME == "Syracuse city")
summary(Sy0_nb)

coords <- coordinates(Syracuse)
IDs <- row.names(Syracuse)
Sy8_nb <- knn2nb(knearneigh(coords, k=1), row.names=IDs)
Sy9_nb <- knn2nb(knearneigh(coords, k=2), row.names=IDs)
Sy10_nb <- knn2nb(knearneigh(coords, k=4), row.names=IDs)
dsts <- unlist(nbdists(Sy8_nb, coords))
Sy11_nb <- dnearneigh(coords, d1=0, d2=0.75*max(dsts), row.names=IDs)

Sy0_lw_W <- nb2listw(Sy0_nb)
Sy0_lw_W

names(Sy0_lw_W)
names(attributes(Sy0_lw_W))

1/rev(range(card(Sy0_lw_W$neighbours)))
summary(unlist(Sy0_lw_W$weights))
summary(sapply(Sy0_lw_W$weights, sum))

dsts <- nbdists(Sy0_nb, coordinates(Syracuse))
idw <- lapply(dsts, function(x) 1/(x/1000))
Sy0_lw_B <- nb2listw(Sy0_nb, style="B")
summary(unlist(Sy0_lw_B$weights))
summary(sapply(Sy0_lw_B$weights, sum))

Sy0_lw_idwB <- nb2listw(Sy0_nb, glist=idw, style="B")
summary(unlist(Sy0_lw_idwB$weights))
summary(sapply(Sy0_lw_idwB$weights, sum))

library(RColorBrewer)

# Three spatial weights representations for Syracuse
pal <- brewer.pal(9, "Reds")
oopar <- par(mfrow=c(1,3), mar=c(1,1,3,1)+0.1)
z <- t(listw2mat(Sy0_lw_W))
brks <- c(0,0.1,0.143,0.167,0.2,0.5,1)
nbr3 <- length(brks)-3
image(1:63, 1:63, z[,ncol(z):1], breaks=brks, col=pal[c(1,(9-nbr3):9)],
      main="W style", axes=FALSE)
box()
z <- t(listw2mat(Sy0_lw_B))
image(1:63, 1:63, z[,ncol(z):1], col=pal[c(1,9)], main="B style", axes=FALSE)
box()
z <- t(listw2mat(Sy0_lw_idwB))
brks <- c(0,0.35,0.73,0.93,1.2,2.6)
nbr3 <- length(brks)-3
image(1:63, 1:63, z[,ncol(z):1], breaks=brks, col=pal[c(1,(9-nbr3):9)],
      main="IDW B style", axes=FALSE)
box()
par(oopar)

try(Sy0_lw_D1 <- nb2listw(Sy11_nb, style="B"))

Sy0_lw_D1 <- nb2listw(Sy11_nb, style="B", zero.policy=TRUE)
print(Sy0_lw_D1, zero.policy=TRUE)

Sy14_nb <- read.gal("Sy_GeoDa1.GAL")
isTRUE(all.equal(Sy0_nb, Sy14_nb, check.attributes=FALSE))

Sy16_nb <- read.gwt2nb("Sy_GeoDa4.GWT")
isTRUE(all.equal(Sy10_nb, Sy16_nb, check.attributes=FALSE))

# Using Weights to Simulate Spatial Autocorrelation

set.seed(987654)
n <- length(Sy0_nb)
uncorr_x <- rnorm(n)
rho <- 0.5
autocorr_x <- invIrW(Sy0_lw_W, rho) %*% uncorr_x

# Simulating spatial autocorrelation: spatial lag plots, 
# showing a locally weighted smoother line
oopar <- par(mfrow=c(1,2), mar=c(4,4,3,2)+0.1)
plot(uncorr_x, lag(Sy0_lw_W, uncorr_x), xlab="", cex.lab=0.8,
 ylab="spatial lag", main="Uncorrelated random variable", cex.main=0.8)
lines(lowess(uncorr_x, lag(Sy0_lw_W, uncorr_x)), lty=2, lwd=2)
plot(autocorr_x, lag(Sy0_lw_W, autocorr_x),
 xlab="", ylab="",
 main="Autocorrelated random variable", cex.main=0.8, cex.lab=0.8)
lines(lowess(autocorr_x, lag(Sy0_lw_W, autocorr_x)), lty=2, lwd=2)
par(oopar)

moran.test(uncorr_x, listw=Sy0_lw_W)
moran.test(autocorr_x, listw=Sy0_lw_W)
moran.test(autocorr_x, listw=nb2listw(Sy9_nb, style="W"))

et <- coords[,1] - min(coords[,1])
trend_x <- uncorr_x + 0.00025 * et
moran_t <- moran.test(trend_x, listw=Sy0_lw_W)
moran_t1 <- lm.morantest(lm(trend_x ~ et), listw=Sy0_lw_W)
# K <- moran(NY8$Cases, listw=nb2listw(NY_nb, style="B"), n=length(NY8$Cases), S0=Szero(nb2listw(NY_nb, style="B")))$K

moran.test(NY8$Cases, listw=nb2listw(NY_nb))
lw_B <- nb2listw(NY_nb, style="B")

moran.test(NY8$Cases, listw=lw_B)
moran.test(NY8$Cases, listw=lw_B, randomisation=FALSE)

lm.morantest(lm(Cases ~ 1, NY8), listw=lw_B)
lm.morantest.sad(lm(Cases ~ 1, NY8), listw=lw_B)
lm.morantest.exact(lm(Cases ~ 1, NY8), listw=lw_B)

# Monte Carlo test
set.seed(1234)
bperm <- moran.mc(NY8$Cases, listw=lw_B, nsim=999)
bperm

r <- sum(NY8$Cases)/sum(NY8$POP8)
rni <- r*NY8$POP8
CR <- function(var, mle) rpois(length(var), lambda=mle)
MoranI.pboot <- function(var, i, listw, n, S0, ...) {
  return(moran(x=var, listw=listw, n=n, S0=S0)$I)
}

set.seed(1234)

library(boot)

boot2 <- boot(NY8$Cases, statistic=MoranI.pboot, R=999, sim="parametric",
  ran.gen=CR, listw=lw_B, n=length(NY8$Cases), S0=Szero(lw_B), mle=rni)

pnorm((boot2$t0 - mean(boot2$t))/sd(boot2$t[,1]), lower.tail=FALSE)

oopar <- par(mfrow=c(1,2))
xlim <- range(c(bperm$res, boot2$t[,1]))
hist(bperm$res[-length(bperm$res)], main="Permutation bootstrap", xlab=expression(I[std]), xlim=xlim, density=15, angle=45, ylim=c(0,260))
abline(v=bperm$statistic, lty=2)
hist(boot2$t, col=rgb(0.4,0.4,0.4), main="Parametric bootstrap", xlab=expression(I[CR]), xlim=xlim, ylim=c(0,260))
hist(bperm$res[-length(bperm$res)], density=15, angle=45, add=TRUE)
abline(v=boot2$t0, lty=2)
par(oopar)

rni <- fitted(glm(Cases ~ 1 + offset(log(POP8)), data=NY8, family="poisson"))

set.seed(1234)
EBImoran.mc(n=NY8$Cases, x=NY8$POP8, listw=nb2listw(NY_nb, style="B"), nsim=999)

cor8 <- sp.correlogram(neighbours=NY_nb, var=NY8$Cases, order=8, method="I", style="C")

library(pgirmess)
corD <- correlog(coordinates(NY8), NY8$Cases, method="Moran")

oopar <- par(mfrow=c(1,2))
plot(cor8, main="Contiguity lag orders")
plot(corD, main="Distance bands")
par(oopar)

print(cor8, p.adj.method="holm") 

oopar <- par(mfrow=c(1,2))
msp <- moran.plot(NY8$Cases, listw=nb2listw(NY_nb, style="C"), quiet=TRUE)
title("Moran scatterplot")
infl <- ifelse(packageVersion("spdep") > "1.1.4", msp$is_inf, apply(msp$is.inf, 1, any))
x <- NY8$Cases
lhx <- cut(x, breaks=c(min(x), mean(x), max(x)), labels=c("L", "H"), include.lowest=TRUE)
wx <- lag(nb2listw(NY_nb, style="C"), NY8$Cases)
lhwx <- cut(wx, breaks=c(min(wx), mean(wx), max(wx)), labels=c("L", "H"), include.lowest=TRUE)
lhlh <- interaction(lhx, lhwx, infl, drop=TRUE)
cols <- rep(1, length(lhlh))
cols[lhlh == "H.L.TRUE"] <- 2
cols[lhlh == "L.H.TRUE"] <- 3
cols[lhlh == "H.H.TRUE"] <- 4
plot(NY8, col=brewer.pal(4, "Accent")[cols])
legend("topright", legend=c("None", "HL", "LH", "HH"), fill=grey.colors(4, 0.95, 0.55, 2.2), bty="n", cex=0.8, y.intersp=0.8)
title("Tracts with influence")
par(oopar)

lm1 <- localmoran(NY8$Cases, listw=nb2listw(NY_nb, style="C"))
lm2 <- as.data.frame(localmoran.sad(lm(Cases ~ 1, NY8), nb=NY_nb, style="C"))
lm3 <- as.data.frame(localmoran.exact(lm(Cases ~ 1, NY8), nb=NY_nb, style="C"))

r <- sum(NY8$Cases)/sum(NY8$POP8)
rni <- r*NY8$POP8
lw <- nb2listw(NY_nb, style="C")
sdCR <- (NY8$Cases - rni)/sqrt(rni)
wsdCR <- lag(lw, sdCR)
I_CR <- sdCR * wsdCR

library(RColorBrewer)
gry <- c(rev(brewer.pal(8, "Reds")[1:6]), brewer.pal(6, "Blues"))
NY8$Standard <- lm1[,1]
NY8$"Constant_risk" <- I_CR
nms <- match(c("Standard", "Constant_risk"), names(NY8))
spplot(NY8, c("Standard", "Constant_risk"), at=c(-2.5,-1.4,-0.6,-0.2,0,0.2,0.6,4,7), col.regions=colorRampPalette(gry)(8))

set.seed(1234)
nsim <- 999
N <- length(rni)
sims <- matrix(0, ncol=nsim, nrow=N)
for (i in 1:nsim) {
  y <- rpois(N, lambda=rni)
  sdCRi <- (y - rni)/sqrt(rni)
  wsdCRi <- lag(lw, sdCRi)
  sims[,i] <- sdCRi * wsdCRi 
}
xrank <- apply(cbind(I_CR, sims), 1, function(x) rank(x)[1])
diff <- nsim - xrank
diff <- ifelse(diff > 0, diff, 0)
pval <- punif((diff + 1)/(nsim + 1))

NY8$Normal <- lm2[,3]
NY8$Randomisation <- lm1[,5]
NY8$Saddlepoint <- lm2[,5]
NY8$Exact <- lm3[,5]
NY8$Constant_risk <- pval
gry <- c(rev(brewer.pal(6, "Reds")), brewer.pal(6, "Blues"))
spplot(NY8, c("Normal", "Randomisation", "Saddlepoint", "Exact", "Constant_risk"), at=c(0,0.01,0.05,0.1,0.9,0.95,0.99,1), col.regions=colorRampPalette(gry)(7))
spplot(NY8, c("Normal", "Exact", "Constant_risk"), xlim=c(405200, 432200), ylim=c(4652700, 4672000), at=c(0,0.01,0.05,0.1,0.9,0.95,0.99,1), col.regions=colorRampPalette(gry)(7))

nylm <- lm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8)
summary(nylm)
NY8$lmresid <- residuals(nylm)

NYlistw <- nb2listw(NY_nb, style = "B")
lm.morantest(nylm, NYlistw)

NYlistwW <- nb2listw(NY_nb, style = "W")
aple(residuals(nylm), listw=NYlistwW)
spautolm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, listw=NYlistwW)$lambda

nysar <- spautolm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, listw=NYlistw)
summary(nysar)

nylam1 <- c(nysar$lambda)
nylam2 <- c(LR1.Spautolm(nysar)$p.value)

NY8$sar_trend <- nysar$fit$signal_trend
NY8$sar_stochastic <- nysar$fit$signal_stochastic
rds <- colorRampPalette(brewer.pal(8, "RdBu"))
tr_at <- seq(-1, 1.3, length.out=21)
tr_rds <- rds(sum(tr_at >= 0)*2)[-(1:(sum(tr_at >= 0)-sum(tr_at < 0)))]
tr_pl <- spplot(NY8, c("sar_trend"), at=tr_at, col="transparent", col.regions=tr_rds, main=list(label="Trend", cex=0.8))
st_at <- seq(-0.16, 0.39, length.out=21)
st_rds <- rds(sum(st_at >= 0)*2)[-(1:(sum(st_at >= 0)-sum(st_at < 0)))]
st_pl <- spplot(NY8, c("sar_stochastic"), at=st_at, col="transparent", col.regions=st_rds, main=list(label="Stochastic", cex=0.8))
plot(tr_pl, split=c(1,1,2,1), more=TRUE)
plot(st_pl, split=c(2,1,2,1), more=FALSE)

nylmw <- lm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, weights=POP8)
summary(nylmw)
NY8$lmwresid <- residuals(nylmw)

library(RColorBrewer)
gry <- c(rev(brewer.pal(6, "Reds")[1:4]), colorRampPalette(brewer.pal(5, "Blues"))(9))
TCEpts <- list("sp.points", TCE, pch=16, col="grey5")
spplot(NY8, c("lmresid", "lmwresid"), sp.layout=list(TCEpts), col.regions=gry, col="transparent", lwd=0.5, at=seq(-2,4.5,0.5))

lm.morantest(nylmw, NYlistw)
nysarw <- spautolm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME , data=NY8, listw=NYlistw, weights=POP8)
summary(nysarw)

NY8$sarw_trend <- nysarw$fit$signal_trend
NY8$sarw_stochastic <- nysarw$fit$signal_stochastic
tr_pl <- spplot(NY8, c("sarw_trend"), at=tr_at, col="transparent", col.regions=tr_rds, main=list(label="Trend", cex=0.8))
st_pl <- spplot(NY8, c("sarw_stochastic"), at=st_at, col="transparent", col.regions=st_rds, main=list(label="Stochastic", cex=0.8))
plot(tr_pl, split=c(1,1,2,1), more=TRUE)
plot(st_pl, split=c(2,1,2,1), more=FALSE)

nycar <- spautolm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME , data=NY8, family="CAR",listw=NYlistw)
summary(nycar)

nylam1 <- c(nycar$lambda)
nycarw <- spautolm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, family="CAR",listw=NYlistw, weights=POP8)
summary(nycarw)

nysarwM <- spautolm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, family="SAR",listw=NYlistw, weights=POP8, method="Matrix")
summary(nysarwM)

1/range(eigenw(NYlistw))
nysar_ll <- spautolm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, family="SAR",listw=NYlistw, llprof=100)
nysarw_ll <- spautolm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, family="SAR",listw=NYlistw, weights=POP8, llprof=100)

ylim <- range(c(nysarw_ll$llprof$ll, nysar_ll$llprof$ll), na.rm=TRUE)
plot(nysarw_ll$llprof$lambda, nysarw_ll$llprof$ll, type="l", xlab=expression(lambda), ylab="log likelihood", ylim=ylim, lwd=2)
abline(v=nysarw_ll$lambda)
abline(h=nysarw_ll$LL)
lines(nysar_ll$llprof$lambda, nysar_ll$llprof$ll, lty=2, lwd=2)
abline(v=nysar_ll$lambda, lty=2)
abline(h=nysar_ll$LL, lty=2)
legend("bottom", legend=c("weighted SAR", "SAR"), lty=c(1,2), lwd=2, bty="n")

nysmaw <- spautolm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, family="SMA",listw=NYlistw, weights=POP8)
summary(nysmaw)