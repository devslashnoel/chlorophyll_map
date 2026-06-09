# Ocean map: OVL-style Sentinel-3A/3B OLCI Chl-a + sample site fluorescence
# Data product: CMEMS Global Ocean Colour, OLCI 300m, Level 3 Multi-Year daily
# Dataset: cmems_obs-oc_glo_bgc-plankton_my_l3-olci-300m_P1D (S3A + S3B merged)

library(ncdf4)
library(ggplot2)
library(ggnewscale)
library(ggrepel)
library(ggspatial)
library(rnaturalearth)
library(rnaturalearthhires)
library(sf)
library(dplyr)
library(scales)

# ── Sample sites (from metadata email, Jan 2025) ──────────────────────────────
sites <- data.frame(
  station = c("S319A-001", "S319A-002", "S319A-003",
              "S319A-004", "S319A-005", "S319A-006"),
  site    = c("Mangawhai (N)", "Leigh (D)",   "Te Arai (D)",
              "Te Arai (N)",   "Leigh (N)",   "Mangawhai (D)"),
  lat     = c(-36.078,  -36.243,  -36.1615, -36.2274, -36.1696, -36.075),
  lon     = c(174.6683, 174.8191, 174.702,  174.8229, 174.7099, 174.6478),
  fluor   = c(5.7, 6.78, 12.02, 8.34, 8.99, 8.42)
)

# ── Bounding box matching the OVL image ───────────────────────────────────────
lat_wide <- c(-37.5, -35.4)
lon_wide <- c(173.8, 176.5)

# ── 300m OLCI Chl-a from Copernicus Marine Service (Multi-Year reprocessed) ────
# Dataset: cmems_obs-oc_glo_bgc-plankton_my_l3-olci-300m_P1D
# Sentinel-3A and 3B are merged into this single daily composite (2016–present).
# 7-day window used to maximise cloud-free coverage.
# Data download is delegated to fetch_cmems.sh (copernicusmarine CLI wrapper).
# Requires: copernicusmarine Python package + free CMEMS credentials
#   Install: pip install copernicusmarine
#   Authenticate once: copernicusmarine login
#   Register: https://data.marine.copernicus.eu/

fetch_cmems_mean <- function(dataset_id, var, time_start, time_end,
                             lon_range, lat_range) {
  tmp <- tempfile(fileext = ".nc")

  # Resolve fetch_cmems.sh relative to this R script; fall back to working directory
  args       <- commandArgs(trailingOnly = FALSE)
  file_args  <- grep("--file=", args, value = TRUE)
  script_dir <- if (length(file_args))
                  dirname(normalizePath(sub("--file=", "", file_args[1]), mustWork = FALSE))
                else
                  "."
  script <- file.path(script_dir, "fetch_cmems.sh")

  cmd <- paste(
    shQuote(script),
    shQuote(dataset_id),
    shQuote(var),
    shQuote(time_start),
    shQuote(time_end),
    lat_range[1], lat_range[2],
    lon_range[1], lon_range[2],
    shQuote(tmp)
  )
  exit_code <- system(cmd, wait = TRUE)
  if (exit_code != 0) {
    stop(
      "fetch_cmems.sh failed (exit code ", exit_code, ").\n",
      "Check credentials (run: copernicusmarine login) and that the dataset ID is correct.\n",
      "Command: ", cmd
    )
  }

  nc      <- nc_open(tmp)
  lon_v   <- ncvar_get(nc, "longitude")
  lat_v   <- ncvar_get(nc, "latitude")
  chl_arr <- ncvar_get(nc, var)          # dims in R: [lon × lat × time]
  nc_close(nc)
  file.remove(tmp)

  chl_mean <- apply(chl_arr, c(1, 2), mean, na.rm = TRUE)
  expand.grid(longitude = lon_v, latitude = lat_v) |>
    mutate(chlor_a = as.vector(chl_mean)) |>
    filter(!is.na(chlor_a) & is.finite(chlor_a))
}

chla_df <- fetch_cmems_mean(
  dataset_id  = "cmems_obs-oc_glo_bgc-plankton_my_l3-olci-300m_P1D",
  var         = "CHL",
  time_start  = "2025-01-01T00:00:00",
  time_end    = "2025-01-07T23:59:59",
  lon_range   = lon_wide,
  lat_range   = lat_wide
)

# ── NZ high-resolution coastline ──────────────────────────────────────────────
nz_coast <- ne_countries(country = "New Zealand", scale = "large", returnclass = "sf")

# ── OVL-style rainbow colour palette ──────────────────────────────────────────
ovl_colors <- c("#6600CC", "#0000FF", "#0099FF", "#00CCCC",
                "#00FF99", "#99FF00", "#FFFF00", "#FFAA00",
                "#FF5500", "#FF0000", "#CC0000")

# ── Plot ──────────────────────────────────────────────────────────────────────
ggplot() +

  # Satellite Chl-a background
  geom_raster(
    data = chla_df,
    aes(x = longitude, y = latitude, fill = log10(chlor_a))
  ) +
  scale_fill_gradientn(
    colours  = ovl_colors,
    limits   = c(-1.5, 1.5),
    oob      = squish,
    na.value = NA,
    name     = "Chl-a\n(log₁₀ mg/m³)"
  ) +

  # Grey land
  geom_sf(data = nz_coast, fill = "grey50", colour = "grey30",
          linewidth = 0.3, inherit.aes = FALSE) +

  # Sample sites coloured by in-situ fluorescence
  new_scale_fill() +
  geom_point(
    data  = sites,
    aes(x = lon, y = lat, fill = fluor),
    size = 5, shape = 21, colour = "white", stroke = 1.5
  ) +
  scale_fill_gradient(
    low  = "#ffffb2", high = "#bd0026",
    name = "In-situ\nFluorescence"
  ) +
  geom_label_repel(
    data        = sites,
    aes(x = lon, y = lat, label = paste0(site, " (", fluor, ")")),
    size        = 2.8, fontface = "bold",
    fill        = "white", colour = "black",
    box.padding = 0.4, label.size = 0, alpha = 0.9
  ) +

  annotation_scale(
    location   = "bl",
    width_hint = 0.25,
    text_col   = "white",
    line_col   = "white",
    bar_cols   = c("white", "grey30")
  ) +
  annotation_north_arrow(
    location = "bl",
    pad_y    = unit(0.6, "cm"),
    style    = north_arrow_fancy_orienteering(
      fill     = c("white", "grey30"),
      line_col = "grey30",
      text_col = "white"
    )
  ) +
  coord_sf(xlim = lon_wide, ylim = lat_wide, expand = FALSE) +
  labs(
    title    = "Hauraki Gulf – Sentinel-3 OLCI Chlorophyll-a (300m, S3A+S3B)",
    subtitle = "CMEMS 7-day mean, 1–7 January 2025  |  Sample sites coloured by in-situ fluorescence",
    x        = "Longitude",
    y        = "Latitude",
    caption  = "Satellite: Copernicus Marine Service / Sentinel-3A & 3B OLCI, 300m Multi-Year reprocessed"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold"),
    legend.position  = "right",
    panel.background = element_rect(fill = "grey85", colour = NA),
    panel.grid       = element_line(colour = "grey70", linewidth = 0.2)
  )

ggsave("Hauraki_Gulf_Fluorescence_Map.png", width = 12, height = 10)

# ── Zoomed-in map: focus on the 6 sample sites ───────────────────────────────
lat_zoom <- c(-36.45, -35.95)
lon_zoom <- c(174.50, 175.05)

ggplot() +
  geom_raster(
    data = chla_df,
    aes(x = longitude, y = latitude, fill = log10(chlor_a))
  ) +
  scale_fill_gradientn(
    colours  = ovl_colors,
    limits   = c(-1.5, 1.5),
    oob      = squish,
    na.value = NA,
    name     = "Chl-a\n(log₁₀ mg/m³)"
  ) +
  geom_sf(data = nz_coast, fill = "grey50", colour = "grey30",
          linewidth = 0.3, inherit.aes = FALSE) +
  new_scale_fill() +
  geom_point(
    data  = sites,
    aes(x = lon, y = lat, fill = fluor),
    size = 6, shape = 21, colour = "white", stroke = 1.8
  ) +
  scale_fill_gradient(
    low  = "#ffffb2", high = "#bd0026",
    name = "In-situ\nFluorescence"
  ) +
  geom_label_repel(
    data        = sites,
    aes(x = lon, y = lat, label = paste0(site, "\n(", fluor, ")")),
    size        = 3, fontface = "bold",
    fill        = "white", colour = "black",
    box.padding = 0.5, label.size = 0.2, alpha = 0.92,
    min.segment.length = 0
  ) +
  annotation_scale(
    location   = "br", width_hint = 0.25,
    text_col   = "white", line_col  = "white",
    bar_cols   = c("white", "grey30")
  ) +
  annotation_north_arrow(
    location = "br", pad_y = unit(0.7, "cm"),
    style    = north_arrow_fancy_orienteering(
      fill = c("white", "grey30"), line_col = "grey30", text_col = "white"
    )
  ) +
  coord_sf(xlim = lon_zoom, ylim = lat_zoom, expand = FALSE) +
  labs(
    title    = "Hauraki Gulf – Sentinel-3 OLCI Chlorophyll-a (300m, S3A+S3B)",
    subtitle = "CMEMS 7-day mean, 1–7 January 2025  |  Sample sites coloured by in-situ fluorescence",
    x = "Longitude", y = "Latitude",
    caption  = "Satellite: Copernicus Marine Service / Sentinel-3A & 3B OLCI, 300m Multi-Year reprocessed"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold"),
    legend.position  = "right",
    panel.background = element_rect(fill = "grey85", colour = NA),
    panel.grid       = element_line(colour = "grey70", linewidth = 0.2)
  )

ggsave("Hauraki_Gulf_Zoomed_Map.png", width = 10, height = 9)
