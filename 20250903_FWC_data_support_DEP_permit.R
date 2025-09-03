
##20250903 Code to support DEP permit application for CKOA
# ----------------------------
# Packages
# ----------------------------
library(readxl)
library(dplyr)
library(ggplot2)
library(sf)   ##

# ----------------------------
# Inputs
# ----------------------------
file_path     <- "C:/Users/flori/Documents/SWCA/Levy/Levy_awarded/Levy Permits/FWC Data Gap/density shapefile data.xlsx"
sheet_name    <- "Summarized Shapefile Data"  # exact tab name
dist_thresh_m <- 100                          # <-- distance threshold (meters), edit as needed

# ----------------------------
# Read & clean
# ----------------------------
df_all <- read_excel(file_path, sheet = sheet_name)
names(df_all) <- trimws(names(df_all))

# Coerce Date if present (handles Excel weird dates if those creep in)
if ("Date" %in% names(df_all)) {
  df_all$Date <- if (is.numeric(df_all$Date)) {
    as.Date(df_all$Date, origin = "1899-12-30")
  } else {
    as.Date(df_all$Date)
  }
}

# Keep only columns we need 
keep_cols <- c("FixedLocationID","Latitude","Longitude","Site","Section","Subsection",
               "SHA","Date","N_Quadrats","OY_m_mean","OY_m_SD","Legal_m_mean",
               "Legal_m_SD","Habitat")
df_subset <- df_all[, intersect(keep_cols, names(df_all)), drop = FALSE]


###########


###########

# ----------------------------
# Define CKOA sites (decimal degrees)
# ----------------------------
ckoa_sites <- data.frame(
  Name      = c("Thompson_Gap", "Sheephead_Hole", "Crooked_Gap_Cowpen", "Herman_Allen_Gap"),
  Latitude  = c(29.15763, 29.15948, 29.15756, 29.15567),
  Longitude = c(-83.00696, -83.00220, -82.99719, -82.99000),
  stringsAsFactors = FALSE
)

# ----------------------------
# Subset data: CK and WAC
# ----------------------------
CK  <- subset(df_subset, Site == "SS" & Section == "S" &
                !is.na(Latitude) & !is.na(Longitude))
WAC <- subset(df_subset, Site == "WC" & Section == "N" &
                !is.na(Latitude) & !is.na(Longitude))

# Combine for distance calculations; tag group
pts <- bind_rows(
  transform(CK,  Group = "CK (SS, S)"),
  transform(WAC, Group = "WAC (WC, N)")
)

# ----------------------------
# Project to meters (UTM 17N) → simple Euclidean distances
# ----------------------------
pts_sf  <- st_as_sf(pts,  coords = c("Longitude","Latitude"), crs = 4326) |> st_transform(32617)
ckoa_sf <- st_as_sf(ckoa_sites, coords = c("Longitude","Latitude"), crs = 4326) |> st_transform(32617)

# Distances from each sample point to each CKOA site
D <- st_distance(pts_sf, ckoa_sf)  # matrix, units = meters
min_d       <- apply(D, 1, min)         # nearest distance per sample
nearest_idx <- apply(D, 1, which.min)   # which CKOA site is nearest

# add points with nearest CKOA + distance (m)
pts$nearest_ckoa  <- ckoa_sites$Name[nearest_idx]
pts$dist_m        <- as.numeric(min_d)
pts$within_thresh <- pts$dist_m <= dist_thresh_m

# Split back to CK/WAC with distance 
CK_d  <- subset(pts, Group == "CK (SS, S)")
WAC_d <- subset(pts, Group == "WAC (WC, N)")

# ----------------------------
# reports
# ----------------------------
cat(sprintf("CK:  %d of %d points within %.0f m of a CKOA site.\n",
            sum(CK_d$within_thresh), nrow(CK_d), dist_thresh_m))
cat(sprintf("WAC: %d of %d points within %.0f m of a CKOA site.\n",
            sum(WAC_d$within_thresh), nrow(WAC_d), dist_thresh_m))

if (any(CK_d$within_thresh)) {
  cat("\nCK points within threshold:\n")
  print(CK_d[CK_d$within_thresh,
             c("FixedLocationID","Latitude","Longitude","Habitat","nearest_ckoa","dist_m")],
        row.names = FALSE)
}
if (any(WAC_d$within_thresh)) {
  cat("\nWAC points within threshold:\n")
  print(WAC_d[WAC_d$within_thresh,
              c("FixedLocationID","Latitude","Longitude","Habitat","nearest_ckoa","dist_m")],
        row.names = FALSE)
}

# Nearest sampled point to EACH CKOA site
nn_idx <- st_nearest_feature(ckoa_sf, pts_sf)
closest_tbl <- data.frame(
  CKOA_site       = ckoa_sites$Name,
  FixedLocationID = pts$FixedLocationID[nn_idx],
  Site            = pts$Site[nn_idx],
  Section         = pts$Section[nn_idx],
  Subsection      = pts$Subsection[nn_idx],
  Habitat         = pts$Habitat[nn_idx],
  dist_m          = as.numeric(st_distance(ckoa_sf, pts_sf[nn_idx, ], by_element = TRUE)),
  stringsAsFactors = FALSE
)
cat("\nClosest sampled point to each CKOA site:\n")
print(closest_tbl, row.names = FALSE)

# ----------------------------
# How far are CKOA sites from each other?
# ----------------------------
D_units <- st_distance(ckoa_sf)               # units: m
D_mat   <- round(as.matrix(D_units), 1)
rownames(D_mat) <- colnames(D_mat) <- ckoa_sites$Name
cat("\nPairwise distances among CKOA sites (m):\n")
print(D_mat)

# nearest neighbor among CKOA sites
D_num <- as.matrix(D_units); diag(D_num) <- Inf
nn_idx_ckoa <- apply(D_num, 1, which.min)
nn_tbl <- data.frame(
  site         = ckoa_sites$Name,
  nearest_site = ckoa_sites$Name[nn_idx_ckoa],
  dist_m       = round(mapply(function(i,j) as.numeric(D_units[i,j]),
                              seq_len(nrow(ckoa_sites)), nn_idx_ckoa), 1),
  stringsAsFactors = FALSE
)
cat("\nNearest neighbor for each CKOA site (m):\n")
print(nn_tbl, row.names = FALSE) 

# ----------------------------
# Plots with CKOA stars; longitude (y) reversed
# ----------------------------
p_ck <- ggplot(CK, aes(x = Latitude, y = Longitude, shape = Habitat)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_point(data = ckoa_sites, aes(x = Latitude, y = Longitude),
             shape = 8, size = 4, inherit.aes = FALSE) +
  geom_text(data = ckoa_sites, aes(x = Latitude, y = Longitude, label = Name),
            nudge_x = 0.0006, nudge_y = -0.0006, size = 3, hjust = 0, inherit.aes = FALSE) +
  coord_equal() +
  scale_y_reverse() +  # put larger negative longitudes at the top
  labs(title = "CK (Site = SS, Section = S)",
       subtitle = "Stars = CKOA sites; axes in decimal degrees",
       x = "Latitude (decimal degrees)", y = "Longitude (decimal degrees)",
       shape = "Habitat", caption = paste0("N = ", nrow(CK))) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right", plot.title.position = "plot",
        legend.title = element_text(face = "bold"))

p_wac <- ggplot(WAC, aes(x = Latitude, y = Longitude, shape = Habitat)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_point(data = ckoa_sites, aes(x = Latitude, y = Longitude),
             shape = 8, size = 4, inherit.aes = FALSE) +
  geom_text(data = ckoa_sites, aes(x = Latitude, y = Longitude, label = Name),
            nudge_x = 0.0006, nudge_y = -0.0006, size = 3, hjust = 0, inherit.aes = FALSE) +
  coord_equal() +
  scale_y_reverse() +
  labs(title = "WAC (Site = WC, Section = N)",
       subtitle = "Stars = CKOA sites; axes in decimal degrees",
       x = "Latitude (decimal degrees)", y = "Longitude (decimal degrees)",
       shape = "Habitat", caption = paste0("N = ", nrow(WAC))) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right", plot.title.position = "plot",
        legend.title = element_text(face = "bold"))

print(p_ck)
print(p_wac)

# save

# Ensure the folder exists
out_dir <- "C:/Users/flori/Documents/SWCA/Levy/Levy_awarded/Levy Permits/FWC Data Gap"

# Save the plots
ggsave(filename = file.path(out_dir, "CK_points_with_CKOA.jpeg"),
       plot = p_ck, width = 6, height = 5, units = "in", dpi = 300)

ggsave(filename = file.path(out_dir, "WAC_points_with_CKOA.jpeg"),
       plot = p_wac, width = 6, height = 5, units = "in", dpi = 300)


ggsave("CK_points_with_CKOA.jpeg",  p_ck,  width = 6, height = 5, dpi = 300)
ggsave("WAC_points_with_CKOA.jpeg", p_wac, width = 6, height = 5, dpi = 300)


###



# Combine CK and WAC
wac_ck <- bind_rows(
  CK  %>% mutate(Group = "CK"),
  WAC %>% mutate(Group = "WAC")
) %>%
  filter(!is.na(Date)) %>%
  mutate(
    Year  = as.integer(format(Date, "%Y")),
    Month = factor(format(Date, "%b"), levels = month.abb, ordered = TRUE)
  )

# Summary: unique sites per Year × Month for each group
summary_sites <- wac_ck %>%
  #group_by(Group, Year, Month) %>%
  group_by(Group) %>%
  summarise(
    n_sites   = n_distinct(FixedLocationID),
    n_records = dplyr::n(),              # (optional) number of rows/samples
    .groups = "drop"
  ) %>%
  arrange(Group, Year, Month)

# View per group
CK_summary  <- filter(summary_sites, Group == "CK")
WAC_summary <- filter(summary_sites, Group == "WAC")

print(CK_summary)
print(WAC_summary)


    
    
    