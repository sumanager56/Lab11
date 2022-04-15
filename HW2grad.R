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

SoilStorage(S_avg, field_capacity, soil_water_content, porosity)

T <- seq(from = 50, to = 500, by = 25)
SStorage2 <- function(X, t = T) {
  out <- matrix(nrow = nrow(X), ncol = length(t), NA)
  for (i in 1:nrow(X)) {
    out[i, ] <- EcoHydRology::SoilStorage(S_avg=t, 
                                    field_capacity=X$fc[i], 
                                    soil_water_content=X$swc[i], 
                                    porosity=X$nfrac[i])
  }
  out <- as.data.frame(out)
  names(out) <- paste("t", t, sep = "")
  return(out)
}

n <- 10
set.seed(1234)
X <- data.frame(#Savg = runif(n, min = 0, max = 20), 
                fc = runif(n, min = 0.25,max = 0.4),
                swc = runif(n, min = 0.25,max = 0.4),
                nfrac = runif(n, min = 0.1,max = 0.45)
                )
Y <- SStorage2(X)
par(cex.axis = 0.7, cex.lab = 0.8)
plot(T, Y[1, ], type = "l", xlab = "S_avg (mm)", ylab = "S value (mm)",
     ylim = c(0, max(Y)))
for (i in 2:n) {
  lines(T, Y[i, ], type = "l", col = i)
}

 
#Fast99
SStorage.seq.fast <- multisensi(design = fast99, model = SStorage2,
                             center = FALSE, reduction = NULL, analysis = analysis.sensitivity,
                             design.args=list( factors=c("fc","swc"
                                                         ,"nfrac"),
                                               n=1000, q = "qunif",
                                               q.arg = list(
                                                 list(min = 0.25,max = 0.4),
                                                 list(min = 0.25,max = 0.4),
                                                 list(min = 0.1,max = 0.45)),
                             analysis.args=list(keep.outputs=FALSE)))

print(SStorage.seq.fast,digits=2)
plot(SStorage.seq.fast, normalized = TRUE, color = terrain.colors, gsi.plot = TRUE)
title(xlab = "Average estimated S value (mm)")
