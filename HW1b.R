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

NetRad(lat, Jday, Tx, Tn, albedo = 0.18, forest = 0, slope = 0, 
       aspect = 0, airtemp = (Tn+Tx)/2, cloudiness = "Estimate", 
       surfemissivity = 0.97, surftemp = (Tn+Tx)/2, units = "kJm2d", 
       AEparams=list(vp=NULL, opt="linear"))

T <- seq(from = 5, to = 365, by = 5)
NetRad2 <- function(X, t = T) {
  out <- matrix(nrow = nrow(X), ncol = length(t), NA)
  for (i in 1:nrow(X)) {
    out[i, ] <- NetRad(lat=X$lat[i], 
                             Tx=X$Tx[i], 
                             Tn=(X$Tx[i]-X$Trange[i]), 
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
Y <- NetRad2(X)
par(cex.axis = 0.7, cex.lab = 0.8)
plot(T, Y[1, ], type = "l", xlab = "Time", ylab = "Net radiation (KJ m2/d)",
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

Y <- NetRad2(X) ## this part can be performed outside R if necessary
NetRad.seq <- multisensi(design=X, model=Y, reduction=NULL, center=FALSE)

plot(NetRad.seq, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)
title(xlab = "Time in half-decades.")
plot(NetRad.seq, normalized = FALSE, color = terrain.colors, gsi.plot = FALSE)
title(xlab = "Time in half-decades.")

#putting above 3 lines in one code
NetRad.seq <- multisensi(model=NetRad2, reduction=NULL, center=FALSE,
                      design.args = list(Tx = c(5,15,25), 
                                         Trange = c(2,9,16), 
                                         slope = c(0.1,0.2,0.3),
                                         aspect = c(0.1,.5,1.0),
                                         albedo = c(0.1,0.5,1.0),
                                         lat=c(0.1,.77,1.1)))
#PCA
NetRad.pca <- multisensi(design=X, model=Y, reduction=basis.ACP, scale=FALSE)
summary(NetRad.pca,digits=2)

plot(NetRad.pca, graph = 1)
plot(NetRad.pca, graph = 2)
plot(NetRad.pca, graph = 3)

#Sobel method
library(sensitivity)
m <- 10000
Xb <- data.frame(lat = runif(m, min = 0, max =pi/3 ), 
                 Tx = runif(m, min = 1,max = 40),
                 Trange = runif(m, min = 1,max = 10),
                 albedo = runif(m, min = 0,max = 1),
                 slope = runif(m, min = 0,max = 0.2),
                 aspect = runif(m, min = 0, max =pi*2))  
NetRad.seq.sobol <- multisensi(design = sobol2007, model = NetRad2,
                            reduction = NULL, analysis = analysis.sensitivity, center = TRUE,
                            design.args = list(X1 = Xb[1:(m/2), ], 
                                               X2 = Xb[(1 + m/2):m, ], 
                                               nboot = 100),
                            analysis.args = list(keep.outputs = FALSE))
## [*] Design
## [*] Response simulation
## [*] Analysis + Sensitivity Indices
print(NetRad.seq.sobol, digits = 2)
plot(NetRad.seq.sobol, normalized = TRUE, color = terrain.colors, gsi.plot = FALSE)
title(xlab = "Time in half-decades")

#Fast99
NetRad.seq.fast <- multisensi(design = fast99, model = NetRad2,
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
print(NetRad.seq.fast,digits=2)
plot(NetRad.seq.fast, normalized = FALSE, color = terrain.colors, gsi.plot = FALSE)
title(xlab = "Time in half-decades")

