rm(list = ls())

library(tuneR)
library(seewave)
library(av)

input_folder <- "falchi_audio"
output_folder <- "falchi_audio_wav"

mp3_files <- list.files(input_folder, pattern = "\\.mp3$", full.names = TRUE)


# Loop and convert
# if FALSE is assumed to already have made the conversione to wav
BOOL_SAVE = F

if (BOOL_SAVE){
  i = 0
  for (file in mp3_files) {
    # Create new filename in output folder
    wav_filename <- paste0(tools::file_path_sans_ext(basename(file)), ".wav")
    output_path <- file.path(output_folder, wav_filename)
    
    # Convert MP3 to WAV
    av_audio_convert(file, output = output_path)
    
    cat("Converted:", output_path, "\n")
    i = i + 1
    print(i)
  }
}


# Plot ----------------------
wav_files <- list.files(output_folder, pattern = "\\.wav$", full.names = TRUE)
# 
# for (file in wav_files) {
#   wave <- readWave(file)
#   oscillo(wave, fastdisp = T)  # Plot waveform
#   title(main = paste("Amplitude Plot of", basename(file)))  # Add title
# }


# Save images

output_oscillo_folder = "falchi_oscillo"

# If FALSE assume plot have already been saved
BOOL_PLOT = F

if(BOOL_PLOT){
  for (file in wav_files) {
    wave <- readWave(file)
    
    # Normalize the signal (optional but good for consistency)
    wave <- normalize(wave, unit = "16")
    
    # Calculate RMS
    rms_value <- rms(wave@left)
    rms_label <- paste("RMS:", round(rms_value, 4))
    
    # Output PNG file
    png_filename <- paste0(tools::file_path_sans_ext(basename(file)), "_amplitude.png")
    png_path <- file.path(output_oscillo_folder, png_filename)
    
    # Open PNG graphics device
    png(png_path, width = 1000, height = 400)
    
    # Plot waveform
    oscillo(wave)
    
    # Add title and RMS label
    title(main = paste("Amplitude Plot of", basename(file)))
    text(x = 0, y = max(wave@left), labels = rms_label, pos = 4, col = "blue", cex = 1.2)
    
    # Close PNG device
    dev.off()
    
    cat("Saved:", png_path, "\n")
  }
}


# try filtering out most noisy audios -------------

# save only the selected samples -------------------------
id_chosen <- read.table("falchi_clear.txt", header = F)
id_chosen$V1

# save selected audio in folder
falchi_selected_wav_folder = "falchi_selected_audio_wav"
falchi_selected_oscillo_folder = "falchi_selected_oscillo"

# save selected audio images in folder

# Convert IDs to strings
id_strings <- as.character(id_chosen$V1)

# Filter files that contain any of the IDs
matching_files <- wav_files[
  sapply(wav_files, function(f) any(grepl(paste0("\\b", id_strings, "\\b", collapse = "|"), f)))
]

file.copy(from = matching_files,
          to = file.path(falchi_selected_wav_folder, basename(matching_files)))

BOOL_PLOT = F

if(BOOL_PLOT){
  for (file in matching_files) {
    wave <- readWave(file)
    
    # Normalize the signal (optional but good for consistency)
    wave <- normalize(wave, unit = "16")
    
    # Calculate RMS
    rms_value <- rms(wave@left)
    rms_label <- paste("RMS:", round(rms_value, 4))
    
    # Output PNG file
    png_filename <- paste0(tools::file_path_sans_ext(basename(file)), "_amplitude.png")
    png_path <- file.path(falchi_selected_oscillo_folder, png_filename)
    
    # Open PNG graphics device
    png(png_path, width = 1000, height = 400)
    
    # Plot waveform
    oscillo(wave, fastdisp = T)
    
    # Add title and RMS label
    title(main = paste("Amplitude Plot of", basename(file)))
    text(x = 0, y = max(wave@left), labels = rms_label, pos = 4, col = "blue", cex = 1.2)
    
    # Close PNG device
    dev.off()
    
    cat("Saved:", png_path, "\n")
  }
}


# Blocks Preprocess -----------------------------------

# select one simple audio
selected_wav = list.files(falchi_selected_wav_folder, pattern = "\\.wav$", full.names = TRUE)

load("audio_segments.RData")


# test ------------------------------
test_wave = readWave(selected_wav[[9]])
oscillo(test_wave, fastdisp = T)
timer(test_wave, threshold = 5, envt = "abs",
      ssmooth = 60, dmin = 0.03,
      tlim = c(0, 2.5))

# Actual procedure ----------------------------------------------

# start with i = 1
i = length(audio_segments) + 1


wave <- readWave(selected_wav[i])

# Normalize the signal (optional but good for consistency)
wave <- normalize(wave, unit = "16")

TLIM = c(0,duration(wave))

oscillo(wave, fastdisp = T)

TLIM = c(0,3.5)

cutted_wave = cutw(wave, from = TLIM[1], to = TLIM[2],
                   output = "Wave")

timer_result <- timer(cutted_wave, threshold = 10, envt = "abs",
                      ssmooth = 60, dmin = 0.03)

# Cut

segments <- lapply(1:length(timer_result$s.start), function(index) {
  cutw(cutted_wave,
    from = timer_result$s.start[index],
       to = timer_result$s.end[index],
       output = "Wave")
})

audio_segments[[i]][[j]] = list("file_name" = selected_wav[i],
                             "segments" = segments,
                           "block_index" = j)
# increment block index
j = j + 1



# increment index
i = i + 1
audio_segments[[i]] = list()
j = 1



print(i)

save(audio_segments,
     file = "audio_segments.RData")























