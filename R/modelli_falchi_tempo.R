rm(list = ls())
load("audio_segments.RData")

library(seewave)
library(tuneR)
library(refund)
library(dplyr)

oscillo(audio_segments[[1]][[1]]$segments[[6]])


# Define a common time grid (normalized time from 0 to 1)
# NOTE: augment n_points 
n_points <- 100  # resolution of interpolation
t_common <- seq(0, 1, length.out = n_points)

# Function to time-normalize each waveform
normalize_time <- function(wave) {
  w <- normalize(wave, unit = "1") # normalize amplitude to [-1,1]
  y <- w@left  
  dur <- duration(w)
  t_original <- seq(0, 1, length.out = length(y))  # map to [0,1]
  approx(t_original, y, xout = t_common)$y  # interpolate to common time grid
}

# Apply to all waveforms
normalized_matrix <- sapply(audio_segments[[1]][[1]]$segments, normalize_time)
dim(normalized_matrix)

for (sample_id in 1:length(audio_segments)){
  if(length(audio_segments[[sample_id]]) > 0){
    
    for(i in 1:length(audio_segments[[sample_id]])){
      audio_segments[[sample_id]][[i]]$norm_matrix = t(sapply(audio_segments[[sample_id]][[i]]$segments,
                                                      normalize_time))
      }
    }
}


# TO DO: check normalization

plot(audio_segments[[1]][[1]]$norm_matrix[2,], type = "l")


# Modello ------------------------------------------------------

temp_matr_1 <- audio_segments[[1]][[1]]$norm_matrix

# concatenated super matrix
temp_matr_Y = temp_matr_1[-1,]
temp_matr_X = temp_matr_1[-NROW(temp_matr_1),]

# initialize
for(i in 2:length(audio_segments[[1]])){
  temp_matr = audio_segments[[1]][[i]]$norm_matrix
  temp_matr_Y = rbind(temp_matr_Y, temp_matr[-1,])
  temp_matr_X = rbind(temp_matr_X, temp_matr[-NROW(temp_matr),])
}

# populate
for (sample_id in 2:length(audio_segments)){
  
  print("sample_id")
  print(sample_id)
  if(length(audio_segments[[sample_id]]) > 0){
    
    for(i in 1:length(audio_segments[[sample_id]])){
      if(length(audio_segments[[sample_id]][[i]]) > 1){
        print("i")
        print(i)
        temp_matr = audio_segments[[sample_id]][[i]]$norm_matrix
        print(dim(temp_matr))
        temp_matr_Y = rbind(temp_matr_Y, temp_matr[-1,])
        temp_matr_X = rbind(temp_matr_X, temp_matr[-NROW(temp_matr),])
      }
      
    }
  }
}

# just to have the correct data structure
# to subscribe
data1 <- pffrSim(scenario = "ff", n = NROW(temp_matr_Y),
                 nxgrid = n_points, nygrid = n_points)

data1$Y = temp_matr_Y
data1$X1 = temp_matr_X

dim(temp_matr_Y)


# REML
fit = pffr(Y ~ ff(X1, xind = t_common,,
                  basistype = "te",
                  splinepars = list(bs="ps",
                                    m = list(c(2, 1),c(2,1)),
                                    k = 10)),
           yind = t_common,
           data = data1,
           bs.yindex = list(bs = "ps", k = 10, m = c(2, 1)),
           bs.int =  list(bs="ps", k = 10, m = c(2, 1)),
           method = "REML")

plot(fit, scale = 0, scheme = 3, select = 1)
plot(fit, scale = 0, scheme = 3, select = 2)
plot(fit, scale = 0, scheme = 3, select = 2, shade = T)

fit_summary = summary(fit)
fit$coefficients

# Bande di confidenza simultanee ------------------
residuals <- residuals(fit, type = "response")
dim(residuals)

plot(seq(0,1, length = 100),
     residuals[1,],
     pch = 16,
     ylim = c(-1,1),
     xlab = "Time", ylab = "Residual")

for(i in 2:NROW(residuals)){
  points(seq(0,1, length = 100), residuals[i,],
         pch = 16)
}

res_mean <- colMeans(residuals)
res_sd <- apply(residuals, 2, sd)

plot(seq(0,1, length = 100), res_mean, type = "l", ylim = c(-1,1),
     xlab = "Time", ylab = "Mean Residual", main = "Mean ± 2SD Residuals")
lines(seq(0,1, length = 100), res_mean + 2 * res_sd, col = "red", lty = 2)
lines(seq(0,1, length = 100), res_mean - 2 * res_sd, col = "red", lty = 2)
abline(h = 0, col = "gray")


fitted_vals <- fitted(fit)

plot(fitted_vals)

fit_coef <- coef(fit)

str(fit_coef$smterms$`ff(X1,t_common,`)
str(fit_coef$smterms$`Intercept(t_common`)

plot(fit_coef$smterms$`ff(X1,t_common,t_common`)

refund:::coef.pffr

fit_coef$smterms$`ff(X1,t_common,t_common`$coef$L.X1

pffr.check(fit)









