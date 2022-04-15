if (!require("pacman")) install.packages("pacman")
pacman::p_load(httr,EcoHydRology,curl,data.table,multisensi)
# KEEP OPEN AS YOU WILL BE WALKING THROUGH IT FOR LAB	
vignette("multisensi-vignette")
#
# Let’s get started as normal. 
setwd("~")
objects()
rm(list=objects())
#
dir.create("~/Week11")
setwd("~/Week11/")
list.files(all.files = T)
objects()   # Should be empty

PET_fromTemp <- function (Jday, Tmax_C, Tmin_C, lat_radians, 
                          AvgT = (Tmax_C + Tmin_C)/2, albedo = 0.18, 
                          TerrestEmiss = 0.97, aspect = 0, slope = 0, 
                          forest = 0, PTconstant=1.26, 
                          AEparams=list(vp=NULL, opt="linear"))

{
  cloudiness <- EstCloudiness(Tmax_C, Tmin_C)
  DailyRad <- NetRad(lat_radians, Jday, Tmax_C, Tmin_C, albedo, forest, slope, aspect, AvgT, cloudiness, TerrestEmiss, AvgT, AEparams=AEparams)
  potentialET <- PTpet(DailyRad, AvgT, PTconstant)
  potentialET[which(potentialET < 0)] <- 0
  potentialET[which(Tmax_C == -999 | Tmin_C == -999)] <- (-999)
  return(potentialET)
}


T <- seq(from = 5, to = 365, by = 5)
PET_fromTemp2 <- function(X, t = T) {
  out <- matrix(nrow = nrow(X), ncol = length(t), NA)
  for (i in 1:nrow(X)) {
    out[i, ] <- PET_fromTemp(lat_radians=X$lat[i], 
                                    Tmax_C=X$Tx[i], 
                                    Tmin_C=(X$Tx[i]-X$Trange[i]), 
                                    albedo=X$albedo[i],
                                    slope=X$slope[i],
                                    aspect=X$aspect[i],
                                    Jday=t)
  }
  out <- as.data.frame(out)
  names(out) <- paste("t", t, sep = "")
  return(out)
}

n <- 10
set.seed(1234)
X <- data.frame(lat = runif(n, min = 0, max =pi/3 ), 
                Tx = runif(n, min = 1,max = 40),
                Trange = runif(n, min = 1,max = 10),
                albedo = runif(n, min = 0,max = 1),
                slope = runif(n, min = 0,max = 0.2),
                aspect = runif(n, min = 0, max =pi*2))
Y <- PET_fromTemp2(X)
par(cex.axis = 0.7, cex.lab = 0.8)
plot(T, Y[1, ], type = "l", xlab = "Time", ylab = "PET (m)",
     ylim = c(0, max(Y)))
for (i in 2:n) {
  lines(T, Y[i, ], type = "l", col = i)
}

X <- expand.grid(Tx = c(5,15,25), 
                 Trange = c(2,9,16), 
                 slope = c(0.1,0.2,0.3),
                 aspect = c(0.1,.5,1.0),
                 albedo = c(0.1,0.5,1.0),
                 lat=c(0.1,.77,1.1))

Y <- PET_fromTemp2(X) ## this part can be performed outside R if necessary
PET.seq <- multisensi(design=X, model=Y, reduction=NULL, center=FALSE)

plot(PET.seq, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)
title(xlab = "Time in half-decades.")
plot(PET.seq, normalized = FALSE, color = terrain.colors, gsi.plot = FALSE)
title(xlab = "Time in half-decades.")

#putting above 3 lines in one code
PET.seq <- multisensi(model=PET_fromTemp2, reduction=NULL, center=FALSE,
                        design.args = list(Tx = c(5,15,25), 
                                           Trange = c(2,9,16), 
                                           slope = c(0.1,0.2,0.3),
                                           aspect = c(0.1,.5,1.0),
                                           albedo = c(0.1,0.5,1.0),
                                           lat=c(0.1,.77,1.1)))
#PCA
PET.pca <- multisensi(design=X, model=Y, reduction=basis.ACP, scale=FALSE)
summary(PET.pca,digits=2)

plot(PET.pca, graph = 1)
plot(PET.pca, graph = 2)
plot(PET.pca, graph = 3)

#Sobel method
library(sensitivity)
m <- 10000
Xb <- data.frame(lat = runif(m, min = 0, max =pi/3 ), 
                 Tx = runif(m, min = 1,max = 40),
                 Trange = runif(m, min = 1,max = 10),
                 albedo = runif(m, min = 0,max = 1),
                 slope = runif(m, min = 0,max = 0.2),
                 aspect = runif(m, min = 0, max =pi*2))  
PET.seq.sobol <- multisensi(design = sobol2007, model = PET_fromTemp2,
                              reduction = NULL, analysis = analysis.sensitivity, center = TRUE,
                              design.args = list(X1 = Xb[1:(m/2), ], 
                                                 X2 = Xb[(1 + m/2):m, ], 
                                                 nboot = 100),
                              analysis.args = list(keep.outputs = FALSE))
## [*] Design
## [*] Response simulation
## [*] Analysis + Sensitivity Indices
print(PET.seq.sobol, digits = 2)
plot(PET.seq.sobol, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)
title(xlab = "Time in half-decades")

#Fast99
PET.seq.fast <- multisensi(design = fast99, model = PET_fromTemp2,
                             center = FALSE, reduction = NULL, analysis = analysis.sensitivity,
                             design.args=list( factors=c("lat","Tx","Trange"
                                                         ,"albedo","slope","aspect"),
                                               n=1000, q = "qunif",
                                               q.arg = list(
                                                 list(min = 0, max =pi/4), 
                                                 list(min = 1,max = 40),
                                                 list(min = 1,max = 10),
                                                 list(min = 0,max = 1),
                                                 list(min = 0,max = 0.2),
                                                 list(min = 0, max =pi*2))),
                             
                             analysis.args=list(keep.outputs=FALSE))
print(PET.seq.fast,digits=2)
dev.off()
plot(PET.seq.fast, normalized = FALSE, color = terrain.colors, gsi.plot = TRUE)
title(xlab = "Time in half-decades")

