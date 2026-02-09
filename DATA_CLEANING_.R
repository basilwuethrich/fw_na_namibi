###############################################################################
# 01_clean_data.R  —  DATA CLEANING (SWARM + DADA2)  [GitHub-ready]
#
# Input:
#   - data/raw/swarm/*.csv
#   - data/raw/dada2/*.csv
#   - data/raw/sample_coordinates.csv
#
# Output:
#   - data/processed/sample_meta.rds
#   - data/processed/swarm_clean_long.rds
#   - data/processed/dada2_clean_long.rds
#
# Config:
#   - config/config.yml
###############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(yaml)
})

# ----------------------------
# 00) Helpers
# ----------------------------
assert_files_exist <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) {
    stop("Missing file(s):\n - ", paste(missing, collapse = "\n - "), call. = FALSE)
  }
}

assert_has_cols <- function(df, cols, df_name = "data") {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) {
    stop(df_name, " is missing column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  }
}

safe_first <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) NA_character_ else x[[1]]
}

# ----------------------------
# 01) Load config
# ----------------------------
cfg <- yaml::read_yaml(here::here("config", "config.yml"))

SWARM_DIR <- here::here(cfg$paths$swarm_dir)
DADA2_DIR <- here::here(cfg$paths$dada2_dir)
COORDS_CSV <- here::here(cfg$paths$coords_csv)
OUT_DIR <- here::here(cfg$paths$out_dir)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

selected_samples <- unlist(cfg$selected_samples)

# removal lists
rm_sw_species <- unlist(cfg$remove$swarm$species)
rm_sw_genus   <- unlist(cfg$remove$swarm$genus)
rm_sw_family  <- unlist(cfg$remove$swarm$family)
rm_sw_class_l <- tolower(str_trim(unlist(cfg$remove$swarm$classes)))

rm_da_species <- unlist(cfg$remove$dada2$species)
rm_da_genus   <- unlist(cfg$remove$dada2$genus)
rm_da_family  <- unlist(cfg$remove$dada2$family)
rm_da_drop_sci <- unlist(cfg$remove$dada2$drop_scientific_names)

keep_taxon_ranks <- tolower(unlist(cfg$filters$keep_taxon_ranks))

# input files
swarm_files <- here::here(cfg$paths$swarm_dir, unlist(cfg$inputs$swarm_files))
dada2_files <- here::here(cfg$paths$dada2_dir, unlist(cfg$inputs$dada2_files))

assert_files_exist(c(COORDS_CSV, swarm_files, dada2_files))

# ----------------------------
# 02) Sample metadata (coords)
# ----------------------------
sample_meta <- readr::read_csv(COORDS_CSV, show_col_types = FALSE)

# expected columns in your coords file
assert_has_cols(sample_meta, c("eventID", "decimalLatitude", "decimalLongitude"), "sample_coordinates")

sample_meta <- sample_meta %>%
  filter(eventID %in% selected_samples) %>%
  rename(sample_id = eventID) %>%
  arrange(decimalLongitude, decimalLatitude) %>%    # West → Ost
  mutate(site_nr = row_number())

if (nrow(sample_meta) == 0) stop("sample_meta is empty after filtering selected_samples.", call. = FALSE)

# ----------------------------
# 03) SWARM import + cleaning
# ----------------------------
swarm_raw <- swarm_files %>%
  set_names() %>%
  purrr::map_dfr(~ readr::read_csv(.x, show_col_types = FALSE))

# Required columns (adapt if your raw schema differs)
assert_has_cols(
  swarm_raw,
  c("sample_name","marker","definition","count_reads",
    "scientific_name_ncbi_corrected","rank_ncbi_corrected",
    "class_name","genus_name_corrected","family_name_corrected"),
  "swarm_raw"
)

swarm_clean_long <- swarm_raw %>%
  mutate(
    sample_name = stringr::str_remove(sample_name, "_\\d+$"),
    taxon       = stringr::str_trim(scientific_name_ncbi_corrected),
    taxon_rank  = tolower(stringr::str_trim(rank_ncbi_corrected)),
    genus_from_taxon = stringr::word(taxon, 1)
  ) %>%
  filter(sample_name %in% selected_samples) %>%
  left_join(sample_meta, by = c("sample_name" = "sample_id")) %>%
  filter(
    !(taxon_rank == "species" & taxon %in% rm_sw_species),
    !(taxon_rank == "genus"   & taxon %in% rm_sw_genus),
    !(taxon_rank == "species" & genus_from_taxon %in% rm_sw_genus),
    !(taxon_rank == "family"  & taxon %in% rm_sw_family)
  ) %>%
  transmute(
    pipeline   = "SWARM",
    sample_id  = sample_name,
    marker     = as.character(marker),
    feature_id = as.character(definition),
    count      = as.numeric(count_reads),
    taxon      = taxon,
    taxon_rank = taxon_rank,
    class      = as.character(class_name),

    species = ifelse(taxon_rank == "species", taxon, NA_character_),
    genus   = ifelse(taxon_rank == "genus",   taxon, as.character(genus_name_corrected)),
    family  = ifelse(taxon_rank == "family",  taxon, as.character(family_name_corrected)),

    decimalLatitude,
    decimalLongitude
  ) %>%
  group_by(sample_id, marker, feature_id) %>%
  summarise(
    pipeline   = first(pipeline),
    count      = sum(count, na.rm = TRUE),
    taxon      = first(taxon),
    taxon_rank = first(taxon_rank),
    class      = first(class),
    species    = first(species),
    genus      = first(genus),
    family     = first(family),
    decimalLatitude  = first(decimalLatitude),
    decimalLongitude = first(decimalLongitude),
    .groups = "drop"
  ) %>%
  mutate(class_l = tolower(str_trim(class))) %>%
  filter(is.na(class_l) | !(class_l %in% rm_sw_class_l)) %>%
  select(-class_l)

# ----------------------------
# 04) DADA2 import + cleaning
# ----------------------------
dada2_raw <- dada2_files %>%
  set_names() %>%
  purrr::map_dfr(~ readr::read_csv(.x, show_col_types = FALSE))

assert_has_cols(
  dada2_raw,
  c("eventID","scientificName","TaxonRank","organismQuantity","pcr_primer_name_forward",
    "class","order","family"),
  "dada2_raw"
)

dada2_filt <- dada2_raw %>%
  mutate(
    TaxonRank_l = tolower(TaxonRank),
    sci = str_trim(scientificName),
    genus_from_sci = stringr::word(sci, 1)
  ) %>%
  filter(TaxonRank_l %in% keep_taxon_ranks) %>%
  filter(eventID %in% selected_samples) %>%
  filter(
    !(TaxonRank_l == "species" & sci %in% rm_da_species),
    !(TaxonRank_l == "genus"   & sci %in% rm_da_genus),
    !(TaxonRank_l == "species" & genus_from_sci %in% rm_da_genus),
    !(TaxonRank_l == "family"  & sci %in% rm_da_family),
    !(sci %in% rm_da_drop_sci)
  ) %>%
  left_join(sample_meta, by = c("eventID" = "sample_id")) %>%
  select(-TaxonRank_l, -genus_from_sci, -sci)

dada2_clean_long <- dada2_filt %>%
  mutate(marker = as.character(pcr_primer_name_forward)) %>%
  group_by(eventID, marker, TaxonRank, scientificName) %>%
  summarise(
    count = sum(as.numeric(organismQuantity), na.rm = TRUE),
    class  = safe_first(class),
    order  = safe_first(order),
    family = safe_first(family),
    genus  = if ("genus" %in% names(dada2_filt)) safe_first(genus) else NA_character_,
    decimalLatitude  = first(decimalLatitude),
    decimalLongitude = first(decimalLongitude),
    .groups = "drop"
  ) %>%
  filter(count > 0) %>%
  transmute(
    pipeline   = "DADA2",
    sample_id  = eventID,
    marker     = marker,
    feature_id = as.character(scientificName),
    count      = count,
    taxon      = as.character(scientificName),
    taxon_rank = as.character(TaxonRank),
    class, order,
    species = ifelse(tolower(taxon_rank) == "species", taxon, NA_character_),
    genus   = ifelse(tolower(taxon_rank) == "genus",   taxon, genus),
    family  = ifelse(tolower(taxon_rank) == "family",  taxon, family),
    decimalLatitude, decimalLongitude
  )

# ----------------------------
# 05) QC summary
# ----------------------------
message("SWARM: Samples=", n_distinct(swarm_clean_long$sample_id),
        " | Features=", n_distinct(swarm_clean_long$feature_id),
        " | Rows=", nrow(swarm_clean_long),
        " | Missing coords rows=", sum(is.na(swarm_clean_long$decimalLongitude)))

message("DADA2: Samples=", n_distinct(dada2_clean_long$sample_id),
        " | Features=", n_distinct(dada2_clean_long$feature_id),
        " | Rows=", nrow(dada2_clean_long),
        " | Missing coords rows=", sum(is.na(dada2_clean_long$decimalLongitude)))

missing_swarm <- setdiff(selected_samples, unique(swarm_clean_long$sample_id))
missing_dada2 <- setdiff(selected_samples, unique(dada2_clean_long$sample_id))
if (length(missing_swarm) > 0) message("WARN SWARM missing samples: ", paste(missing_swarm, collapse = ", "))
if (length(missing_dada2) > 0) message("WARN DADA2 missing samples: ", paste(missing_dada2, collapse = ", "))

# ----------------------------
# 06) Save outputs
# ----------------------------
saveRDS(sample_meta,       file.path(OUT_DIR, "sample_meta.rds"))
saveRDS(swarm_clean_long,  file.path(OUT_DIR, "swarm_clean_long.rds"))
saveRDS(dada2_clean_long,  file.path(OUT_DIR, "dada2_clean_long.rds"))

message("Saved to: ", OUT_DIR)



