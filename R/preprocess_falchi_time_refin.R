library(tuneR)
library(seewave)
library(av)

rm(list = ls())


load(file = "data/audio_segments.RData")

# only first time
# audio_segments_refin = as.list(audio_segments)

load(file = "data/audio_segments_refin.RData")



# Actual procedure ----------------------------------------------

#start with 
# i = 1
# j = 1
# k = 1

i = 3
j = 1
k = 1

i = i + 1
j = 1 # reset j
k = 1 # reset k

j = j + 1
k = 1 # reset k

k = k + 1


oscillo(audio_segments_refin[[i]][[j]]$segments[[k]], fastdisp = T)

# original duration
TLIM = c(0, seewave::duration(audio_segments_refin[[i]][[j]]$segments[[k]]))

TLIM = c(0.02,0.04)

cutted_wave = cutw(audio_segments_refin[[i]][[j]]$segments[[k]],
                   from = TLIM[1], to = TLIM[2],
                   output = "Wave")


oscillo(cutted_wave, fastdisp = T)

audio_segments_refin[[i]][[j]]$segments[[k]] = cutted_wave

cat("i: ", i, "\n")
cat("j: ", j, "\n")
cat("k: ", k, "\n")


save(audio_segments_refin,
     file = "data/audio_segments_refin.RData")






