library(av)
library(tuneR)
library(seewave)
library(soundgen)
library(fda)
library(warbleR)
library(suwo)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(terra)
library(tidyverse)

load("data_1.RData")
# save(
#   falchi,
#   falchi_meanspec_amps,
#   falchi_meanspec_freqs,
#   gufi,
#   gufi_meanspec_amps,
#   gufi_meanspec_freqs,
#   gabbiani,
#   gabbiani_meanspec_amps,
#   gabbiani_meanspec_freqs,
#   file = "data_1.RData"
# )

dwnl = F # Scelgo se scaricare o no i file
# ╭──────╮
# │Falchi│
# ╰──────╯
falchi_dir_name = "animali_audio/falchi_audio/old"
if (!dir.exists(falchi_dir_name)) dir.create(falchi_dir_name, recursive = T)
falchi = query_xc(
  'Falco tinnunculus q:">C" area:europe len:"<31"',
  download = dwnl,
  path = falchi_dir_name,
  parallel = 4
)
falchi = falchi |>
  # Tolgo gli audio con più specie nello stesso momento
  filter(if_all(contains("Other_species"), is.na)) |>
  # Calcolo minuti e secondi per poi calcolare la durata degli audio
  separate(Length, into = c("Minutes", "Seconds"), sep = ":", convert = T) |>
  mutate(
    Scientific_name = paste(Genus, Specific_epithet, sep = "-"),
    # Sistemo le coordinate mancanti
    Latitude = case_match(
      Locality,
      "Holmnäs, Umeå, Västerbottens län" ~ 63.8007,
      "Mungialdea (near  Iturribaltzaga), Bizkaia, Euskadi" ~ 43.3833,
      .default = Latitude
    ),
    Longitude = case_match(
      Locality,
      "Holmnäs, Umeå, Västerbottens län" ~ 19.85004,
      "Mungialdea (near  Iturribaltzaga), Bizkaia, Euskadi" ~ -2.863,
      .default = Longitude
    ),
    # Calcolo durata audio in secondi
    Duration = Minutes * 60 + Seconds,
    English_name = as.factor(English_name)
  ) |>
  # Tolgo colonne inutili
  select(
    Recording_ID,
    English_name,
    Latitude,
    Longitude,
    Vocalization_type,
    Audio_file,
    Quality,
    Scientific_name,
    Duration
  )

# Elimino audio non selezionati e converto quelli rimasti in un formato comune
remove_useless_files = function(dir_name, dataset) {
  all_files = list.files(dir_name, pattern = "\\.mp3$", full.names = TRUE)
  file_ids = sapply(
    all_files,
    function(f)
      as.numeric(gsub("\\D+", "", tools::file_path_sans_ext(basename(f))))
  )
  file.remove(all_files[!file_ids %in% dataset$Recording_ID])
}
convert_files = function(filenames, old_dir, new_dir) {
  for (f in filenames) {
    print(f)
    system(paste0(
      "ffmpeg -i ",
      file.path(old_dir, f),
      " -ar 44100 -ac 1 -ab 256k ",
      file.path(new_dir, f)
    ))
  }
}
remove_and_convert_audio = function(dir_name, dataset) {
  remove_useless_files(dir_name, dataset)
  filenames = basename(list.files(
    dir_name,
    pattern = "\\.mp3$",
    full.names = TRUE
  ))
  new_dir_name = gsub("/old", "", dir_name)
  convert_files(
    filenames,
    dir_name,
    new_dir_name
  )
}
if (dwnl) {
  remove_and_convert_audio(falchi_dir_name, falchi)
}

# Mappa degli audio
map_europe <- ne_countries(scale = "medium", returnclass = "sf")
ggplot(data = map_europe) +
  geom_sf(fill = "gray90", color = "white") +
  geom_point(
    data = falchi,
    aes(x = Longitude, y = Latitude),
    size = 3,
    colour = "red"
  ) +
  theme_minimal() +
  coord_sf(ylim = c(25, 73), xlim = c(-30, 45)) +
  labs(title = "Mappa audio Gheppio Comune (Falco tinnunculus)")

# Assegno ogni audio ad una certa zona climatica
kg_map = rast("koppen_geiger_0p00833333.tif")
crop(kg_map, ext(c(-30, 45, 25, 73))) |> plot()
points(falchi$Longitude, falchi$Latitude, pch = 17)
falchi_points = st_as_sf(
  falchi,
  coords = c("Longitude", "Latitude"),
  crs = 4326
)
climate_codes = terra::extract(
  kg_map,
  vect(falchi_points),
  search_radius = 5000
)
crop(kg_map, ext(c(-30, 45, 25, 73))) |> plot()
plot(st_geometry(falchi_points[is.na(climate_codes[, 2]), ]), add = T, pch = 17)
kg_zones = read.csv("climate_classifications_hex.csv")
climate_codes = left_join(
  x = climate_codes,
  y = kg_zones,
  by = join_by(koppen_geiger_0p00833333 == ID)
)
table(climate_codes$Description)
falchi = cbind(falchi, climate_codes)[, -c(10, 12:13)]
ggplot(data = map_europe) +
  geom_sf(fill = "gray90", color = "white") +
  geom_point(
    data = falchi,
    aes(x = Longitude, y = Latitude, colour = Code),
    size = 3
  ) +
  scale_color_manual(values = setNames(falchi$ColorHex, falchi$Code)) +
  theme_minimal() +
  coord_sf(ylim = c(25, 73), xlim = c(-30, 45)) +
  labs(title = "Mappa clima Gheppio Comune (Falco tinnunculus)")
falchi = falchi |>
  mutate(
    Climate_zone = case_match(
      Code,
      c("BSh", "BSk", "BWh", "Csa", "Csb") ~ "Hot",
      c("Cfa", "Cfb") ~ "Temperate",
      c("Dfa", "Dfb", "Dfc") ~ "Cold",
    )
  ) |>
  mutate(Climate_zone = as.factor(Climate_zone))

# Carico gli audio
falchi_audio = lapply(
  file.path(
    gsub("/old", "", falchi_dir_name),
    paste0(falchi$Scientific_name, "-", falchi$Recording_ID, ".mp3")
  ),
  readMP3
)
# i = 7
# spectro(falchi_audio[[i]])
# play(falchi_audio[[i]])
# par(mfrow = c(1, 1))
# oscillo(falchi_audio[[i]])
# meanspec(falchi_audio[[i]], wl = 1024)

# Filtro gli audio dai 2.5 kHz ai 7.5 kHz
falchi_audio_filt = list()
for (i in 1:nrow(falchi)) {
  af = ffilter(
    falchi_audio[[i]],
    from = 2500,
    to = 7500,
    bandpass = TRUE,
    output = "Wave",
    ovlp = 50,
    rescale = T
  )
  falchi_audio_filt[[i]] = af
}
# i = 7
# spectro(falchi_audio_filt[[i]])
# play(falchi_audio_filt[[i]])
# par(mfrow = c(1, 1))
# oscillo(falchi_audio_filt[[i]])
# meanspec(falchi_audio_filt[[i]], wl = 1024)

# Calcolo gli spettri di frequenza medi
falchi_meanspec_matrix = meanspec(
  falchi_audio_filt[[1]],
  wl = 1024,
  plot = F,
  norm = T
)
for (i in 2:nrow(falchi)) {
  falchi_meanspec_matrix = cbind(
    falchi_meanspec_matrix,
    meanspec(falchi_audio_filt[[i]], wl = 1024, plot = F, norm = T)[, 2]
  )
}
# Seleziono solo le frequenze non filtrate
falchi_meanspec_matrix = falchi_meanspec_matrix[
  which(
    (2.3 <= falchi_meanspec_matrix[, 1]) &
      (falchi_meanspec_matrix[, 1] <= 7.7)
  ),
]
falchi_meanspec_freqs = falchi_meanspec_matrix[, 1]
falchi_meanspec_amps = falchi_meanspec_matrix[, -1]
# Grafico spettri normalizzati per ogni uccello
par(mfrow = c(1, 1))
plot(
  falchi_meanspec_freqs,
  falchi_meanspec_amps[, 1],
  ylim = c(0, max(falchi_meanspec_amps)),
  type = "l"
)
for (i in 2:ncol(falchi_meanspec_amps)) {
  lines(
    falchi_meanspec_freqs,
    falchi_meanspec_amps[, i],
    col = falchi$Climate_zone[i]
  )
}

# ╭────╮
# │Gufi│
# ╰────╯
gufi_dir_name = "animali_audio/gufi_audio/old"
if (!dir.exists(gufi_dir_name)) dir.create(gufi_dir_name, recursive = T)
gufi = query_xc(
  'Strix aluco q:">C" area:europe len:"<31"',
  download = dwnl,
  path = gufi_dir_name,
  parallel = 4
)
gufi = gufi |>
  # Tolgo gli audio con più specie nello stesso momento
  filter(if_all(contains("Other_species"), is.na)) |>
  # Calcolo minuti e secondi per poi calcolare la durata degli audio
  separate(Length, into = c("Minutes", "Seconds"), sep = ":", convert = T) |>
  mutate(
    Scientific_name = paste(Genus, Specific_epithet, sep = "-"),
    Duration = Minutes * 60 + Seconds,
    English_name = as.factor(English_name)
  ) |>
  # Tolgo colonne inutili
  select(
    Recording_ID,
    English_name,
    Latitude,
    Longitude,
    Vocalization_type,
    Audio_file,
    Quality,
    Scientific_name,
    Duration
  ) |>
  filter(!is.na(Latitude))

# Elimino audio non selezionati e converto quelli rimasti in un formato comune
if (dwnl) {
  remove_and_convert_audio(gufi_dir_name, gufi)
}

# Mappa degli audio
ggplot(data = map_europe) +
  geom_sf(fill = "gray90", color = "white") +
  geom_point(
    data = gufi,
    aes(x = Longitude, y = Latitude),
    size = 3,
    color = "red"
  ) +
  theme_minimal() +
  coord_sf(ylim = c(25, 73), xlim = c(-30, 45)) +
  labs(title = "Mappa audio Allocco Comune (Strix aluco)")

# Assegno ogni audio ad una certa zona climatica
crop(kg_map, ext(c(-30, 45, 25, 73))) |> plot()
points(gufi$Longitude, gufi$Latitude, pch = 17)
gufi_points = st_as_sf(
  gufi,
  coords = c("Longitude", "Latitude"),
  crs = 4326
)
climate_codes = terra::extract(
  kg_map,
  vect(gufi_points),
  search_radius = 5000
)
crop(kg_map, ext(c(-30, 45, 25, 73))) |> plot()
plot(st_geometry(gufi_points[is.na(climate_codes[, 2]), ]), add = T, pch = 17)
climate_codes = left_join(
  x = climate_codes,
  y = kg_zones,
  by = join_by(koppen_geiger_0p00833333 == ID)
)
table(climate_codes$Description)
gufi = cbind(gufi, climate_codes)[, -c(10, 12:13)]
gufi = gufi[gufi$Code != "Dsb", ]
ggplot(data = map_europe) +
  geom_sf(fill = "gray90", color = "white") +
  geom_point(
    data = gufi,
    aes(x = Longitude, y = Latitude, colour = Code),
    size = 3
  ) +
  scale_color_manual(values = setNames(gufi$ColorHex, gufi$Code)) +
  theme_minimal() +
  coord_sf(ylim = c(25, 73), xlim = c(-30, 45)) +
  labs(title = "Mappa clima Allocco Comune (Strix aluco)")
gufi = gufi |>
  mutate(
    Climate_zone = case_match(
      Code,
      c("BSk", "Csa", "Csb") ~ "Hot",
      c("Cfa", "Cfb") ~ "Temperate",
      c("Dfb", "Dfc") ~ "Cold",
    )
  ) |>
  mutate(Climate_zone = as.factor(Climate_zone))
# Carico gli audio
gufi_audio = lapply(
  file.path(
    gsub("/old", "", gufi_dir_name),
    paste0(gufi$Scientific_name, "-", gufi$Recording_ID, ".mp3")
  ),
  readMP3
)
# i = 7
# spectro(gufi_audio[[i]])
# play(gufi_audio[[i]])
# par(mfrow = c(1, 1))
# oscillo(gufi_audio[[i]])
# meanspec(gufi_audio[[i]], wl = 1024)

# Filtro gli audio dai 0 Hz ai 5 kHz
gufi_audio_filt = list()
for (i in 1:nrow(gufi)) {
  af = ffilter(
    gufi_audio[[i]],
    from = 0,
    to = 5000,
    bandpass = TRUE,
    output = "Wave",
    ovlp = 50,
    rescale = T
  )
  gufi_audio_filt[[i]] = af
}
# i = 7
# spectro(gufi_audio_filt[[i]])
# play(gufi_audio_filt[[i]])
# par(mfrow = c(1, 1))
# oscillo(gufi_audio_filt[[i]])
# meanspec(gufi_audio_filt[[i]], wl = 1024)

# Calcolo gli spettri di frequenza medi
gufi_meanspec_matrix = meanspec(
  gufi_audio_filt[[1]],
  wl = 1024,
  plot = F,
  norm = T
)
for (i in 2:nrow(gufi)) {
  gufi_meanspec_matrix = cbind(
    gufi_meanspec_matrix,
    meanspec(gufi_audio_filt[[i]], wl = 1024, plot = F, norm = T)[, 2]
  )
}
# Seleziono solo le frequenze non filtrate
gufi_meanspec_matrix = gufi_meanspec_matrix[
  which(gufi_meanspec_matrix[, 1] <= 5.2),
]
gufi_meanspec_freqs = gufi_meanspec_matrix[, 1]
gufi_meanspec_amps = gufi_meanspec_matrix[, -1]
# Grafico spettri normalizzati per ogni uccello
par(mfrow = c(1, 1))
plot(
  gufi_meanspec_freqs,
  gufi_meanspec_amps[, 1],
  ylim = c(0, max(gufi_meanspec_amps)),
  type = "l"
)
for (i in 2:ncol(gufi_meanspec_amps)) {
  lines(
    gufi_meanspec_freqs,
    gufi_meanspec_amps[, i],
    col = gufi$Climate_zone[i]
  )
}

# ╭────────╮
# │Gabbiani│
# ╰────────╯
gabbiani_dir_name = "animali_audio/gabbiani_audio/old"
if (!dir.exists(gabbiani_dir_name)) dir.create(gabbiani_dir_name, recursive = T)
gabbiani = query_xc(
  'Larus michahellis q:">C" area:europe len:"<31"',
  download = dwnl,
  path = gabbiani_dir_name,
  parallel = 4
)
gabbiani = gabbiani |>
  # Tolgo gli audio con più specie nello stesso momento
  filter(if_all(contains("Other_species"), is.na)) |>
  # Calcolo minuti e secondi per poi calcolare la durata degli audio
  separate(Length, into = c("Minutes", "Seconds"), sep = ":", convert = T) |>
  mutate(
    Scientific_name = paste(Genus, Specific_epithet, sep = "-"),
    # Calcolo durata audio in secondi
    Duration = Minutes * 60 + Seconds,
    English_name = as.factor(English_name)
  ) |>
  # Tolgo colonne inutili
  select(
    Recording_ID,
    English_name,
    Latitude,
    Longitude,
    Vocalization_type,
    Audio_file,
    Quality,
    Scientific_name,
    Duration
  ) |>
  filter(Latitude < 50, Longitude > -20, Longitude < 25)

# Elimino audio non selezionati e converto quelli rimasti in un formato comune
if (dwnl) {
  remove_and_convert_audio(gabbiani_dir_name, gabbiani)
}

# Mappa degli audio
ggplot(data = map_europe) +
  geom_sf(fill = "gray90", color = "white") +
  geom_point(
    data = gabbiani,
    aes(x = Longitude, y = Latitude),
    size = 3,
    color = "red"
  ) +
  theme_minimal() +
  coord_sf(ylim = c(25, 73), xlim = c(-30, 45)) +
  labs(title = "Mappa audio Gabbiano Reale Zampegialle (Larus michahellis)")

# Assegno ogni audio ad una certa zona climatica
crop(kg_map, ext(c(-30, 45, 25, 73))) |> plot()
points(gabbiani$Longitude, gabbiani$Latitude, pch = 17)
gabbiani_points = st_as_sf(
  gabbiani,
  coords = c("Longitude", "Latitude"),
  crs = 4326
)
climate_codes = terra::extract(
  kg_map,
  vect(gabbiani_points),
  search_radius = 5000
)
crop(kg_map, ext(c(-30, 45, 25, 73))) |> plot()
plot(st_geometry(gufi_points[is.na(climate_codes[, 2]), ]), add = T, pch = 17)
climate_codes = left_join(
  x = climate_codes,
  y = kg_zones,
  by = join_by(koppen_geiger_0p00833333 == ID)
)
table(climate_codes$Description)
gabbiani = cbind(gabbiani, climate_codes)[, -c(10, 12:13)]
ggplot(data = map_europe) +
  geom_sf(fill = "gray90", color = "white") +
  geom_point(
    data = gabbiani,
    aes(x = Longitude, y = Latitude, colour = Code),
    size = 3
  ) +
  geom_vline(xintercept = 16) +
  geom_hline(yintercept = 45.95) +
  scale_color_manual(values = setNames(gabbiani$ColorHex, gabbiani$Code)) +
  theme_minimal() +
  coord_sf(ylim = c(25, 73), xlim = c(-30, 45)) +
  labs(title = "Mappa clima Gabbiano Reale Zampegialle (Larus michahellis)")
gabbiani = gabbiani |>
  mutate(
    Cluster = case_when(
      Longitude > -20 & Longitude < -10 & Latitude > 25 & Latitude < 30 ~
        "Canary",
      Longitude > -10 & Longitude < -4.7 & Latitude > 35 & Latitude < 45 ~
        "Atlantic1",
      Longitude > -4.7 & Longitude < -0.1 & Latitude > 43 & Latitude < 48 ~
        "Atlantic2",
      Longitude > -0.1 & Longitude < 16 & Latitude > 45.95 & Latitude < 50 ~
        "Centre",
      .default = "Mediterrean"
    )
  ) |>
  mutate(
    Cluster = ifelse(
      Cluster == "Atlantic1" | Cluster == "Atlantic2",
      "Atlantic",
      Cluster
    )
  ) |>
  mutate(Cluster = as.factor(Cluster))
ggplot(data = map_europe) +
  geom_sf(fill = "gray90", color = "white") +
  geom_point(
    data = gabbiani,
    aes(x = Longitude, y = Latitude, colour = Cluster),
    size = 3
  ) +
  theme_minimal() +
  coord_sf(ylim = c(25, 73), xlim = c(-30, 45)) +
  labs(title = "Mappa zone Gabbiano Reale Zampegialle (Larus michahellis)")
table(gabbiani$Cluster)

# Carico gli audio
gabbiani_audio = lapply(
  file.path(
    gsub("/old", "", gabbiani_dir_name),
    paste0(gabbiani$Scientific_name, "-", gabbiani$Recording_ID, ".mp3")
  ),
  readMP3
)
# i = 5
# spectro(gabbiani_audio[[i]])
# play(gabbiani_audio[[i]])
# par(mfrow = c(1, 1))
# oscillo(gabbiani_audio[[i]])
# meanspec(gabbiani_audio[[i]], wl = 1024)

# Filtro gli audio dai 0 Hz ai 5 kHz
gabbiani_audio_filt = list()
for (i in 1:nrow(gabbiani)) {
  af = ffilter(
    gabbiani_audio[[i]],
    from = 0,
    to = 5000,
    bandpass = TRUE,
    output = "Wave",
    ovlp = 50,
    rescale = T
  )
  gabbiani_audio_filt[[i]] = af
}
# i = 7
# spectro(gabbiani_audio_filt[[i]])
# play(gabbiani_audio_filt[[i]])
# par(mfrow = c(1, 1))
# oscillo(gabbiani_audio_filt[[i]])
# meanspec(gabbiani_audio_filt[[i]], wl = 1024)

# Calcolo gli spettri di frequenza medi
gabbiani_meanspec_matrix = meanspec(
  gabbiani_audio_filt[[1]],
  wl = 1024,
  plot = F,
  norm = T
)
for (i in 2:nrow(gabbiani)) {
  gabbiani_meanspec_matrix = cbind(
    gabbiani_meanspec_matrix,
    meanspec(gabbiani_audio_filt[[i]], wl = 1024, plot = F, norm = T)[, 2]
  )
}
# Seleziono solo le frequenze non filtrate
gabbiani_meanspec_matrix = gabbiani_meanspec_matrix[
  which(gabbiani_meanspec_matrix[, 1] <= 5.2),
]
gabbiani_meanspec_freqs = gabbiani_meanspec_matrix[, 1]
gabbiani_meanspec_amps = gabbiani_meanspec_matrix[, -1]
# Grafico spettri normalizzati per ogni uccello
par(mfrow = c(1, 1))
plot(
  gabbiani_meanspec_freqs,
  gabbiani_meanspec_amps[, 1],
  ylim = c(0, max(gabbiani_meanspec_amps)),
  type = "l"
)
for (i in 2:ncol(gabbiani_meanspec_amps)) {
  lines(
    gabbiani_meanspec_freqs,
    gabbiani_meanspec_amps[, i],
    col = gabbiani$Cluster[i]
  )
}
