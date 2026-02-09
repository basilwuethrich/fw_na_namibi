###############################################################################
# 01_clean_data.R  —  DATA CLEANING (SWARM + DADA2)
# Ziel: pro Pipeline ein "clean long format" Objekt erzeugen, das später
#       für PA-Matrix + PCoA wiederverwendet wird.
#
# Output:
#   - swarm_clean_long   (sample_id, feature_id, count, marker, coords + taxon)
#   - dada2_clean_long   (sample_id, feature_id, count, primer + coords + taxonrank)
#   - sample_meta        (sample_id, decimalLatitude, decimalLongitude, ...)
#
# Optional speichern:
#   saveRDS(swarm_clean_long, "data/swarm_clean_long.rds")
#   saveRDS(dada2_clean_long, "data/dada2_clean_long.rds")
###############################################################################
suppressPackageStartupMessages({
  library(tidyverse)  # includes readr, dplyr, tidyr, tibble, stringr
})

dada2_clean_long %>%
  filter(class == "Polychaeta") %>%
  count(taxon, taxon_rank, marker, sort = TRUE)

dada2_clean_long %>%
  filter(class == "Polychaeta") %>%
  summarise(min=min(count), median=median(count), max=max(count), n=n())
dada2_clean_long %>%
  filter(class == "Polychaeta") %>%
  count(sample_id, sort = TRUE)

###############################################################################
# 00) CONFIG (bitte hier alles sammeln, was du später evtl. änderst)
###############################################################################

# --- Pfade (anpassen) ---
path_swarm <- "C:/Users/basil/OneDrive/Desktop/Master_Thesis_2.0/DATA/Basil_Namibia_SWARM"
path_dada2 <- "C:/Users/basil/OneDrive/Desktop/Master_Thesis_2.0/DATA/taxadata"
path_coords <- "C:/Users/basil/OneDrive/Desktop/Master_Thesis_2.0/DATA/sample_coordinates.csv"


# --- Samples (deine 20) ---
selected_samples <- c(
  "SPY221856","SPY221857","SPY221858","SPY221859","SPY221860",
  "SPY221861","SPY221862","SPY221863","SPY221864","SPY221865",
  "SPY221866","SPY221867","SPY221868","SPY221869","SPY221872",
  "SPY221873","SPY221870","SPY221871","SPY221653","SPY221654"
)

# --- Entfernen (Afrika-Fremde / Kontamination etc.) ---
species_to_remove_swarm <- c(
  "Chloebia gouldiae","Iranocichla hormuzensis",
  "Plegadis chihi","Pycnonotus taivanus","Pycnonotus xanthorrhous",
  "Trachurus trachurus","Turdus dissimilis","Homo sapiens",
  "Heteromormyrus longilateralis","Hyperopisus bebe"
)
genus_to_remove_swarm <- c(
  "Acanthorhynchus","Homo","Plagioscion","Systomus",
  "Bos","Eudocimus","Trachurus"
)
family_to_remove_swarm <- c(
  "Loricariidae","Ampharetidae"
)

# --- DADA2 removal lists ---
species_to_remove_dada2 <- c(
  "Chloebia gouldiae","Eupera ferruginea","Iranocichla hormuzensis",
  "Plegadis chihi","Pycnonotus taivanus","Pycnonotus xanthorrhous",
  "Trachurus trachurus","Turdus dissimilis",
  "Creagrutus melanzonus","Butastur liventer","Anhinga melanogaster",
  "Cairina moschata","Equus asinus","Columba livia"
)

genus_to_remove_dada2 <- c(
  "Acanthorhynchus","Homo","Plagioscion","Systomus",
  "Gallus","Gorilla","Pan","Trachurus","Cairina","Aphaniops","Amphisamytha"
)

family_to_remove_dada2 <- c(
  "Loricariidae","Carangidae","Fundulidae","Atherinopsidae",
  "Serrasalmidae","Stevardiidae"
)
# --- Klassen entfernen
classes_to_remove_swarm <- c("Flabellinia","Polychaeta","Scyphozoa")
classes_to_remove_swarm_l <- tolower(classes_to_remove_swarm)
  

# --- DADA2: welche TaxonRanks bleiben? ---
keep_taxon_ranks <- c("species", "genus", "family")

###############################################################################
# 01) SAMPLE METADATA (Koordinaten)  —  EINMAL laden, für beide nutzen
###############################################################################

sample_meta <- read_csv(path_coords, show_col_types = FALSE) %>%
  filter(eventID %in% selected_samples) %>%
  rename(sample_id = eventID)

sample_meta <- sample_meta %>%
  arrange(decimalLongitude, decimalLatitude) %>%  # West → Ost
  mutate(site_nr = row_number())

# Quick QC
stopifnot(nrow(sample_meta) > 0)

###############################################################################
# 02) SWARM  —  import + cleaning → swarm_clean_long
###############################################################################

# 02.1 Import (4 Marker) + Merge
swarm_files <- c(
  file.path(path_swarm, "swarm_Euka3_clean.csv"),
  file.path(path_swarm, "swarm_Unio_clean.csv"),
  file.path(path_swarm, "swarm_V05_clean.csv"),
  file.path(path_swarm, "swarm_teleo_clean.csv")

)

swarm_raw <- swarm_files %>%
  set_names() %>%
  map_dfr(~ read_csv(.x, show_col_types = FALSE))

 # 02.2 Cleaning / Harmonisierung
# - sample_name suffix "_1" etc entfernen
# - auf selected samples filtern
# - Koordinaten joinen
# - species_to_remove filtern
# - in "long canonical" bringen: sample_id, feature_id, count, marker, taxonomy...
swarm_clean_long <- swarm_raw %>%
  mutate(
    sample_name = str_remove(sample_name, "_\\d+$"),
    taxon      = str_trim(scientific_name_ncbi_corrected),
    taxon_rank = tolower(str_trim(rank_ncbi_corrected)),
    genus_from_taxon = stringr::word(taxon, 1)
  ) %>%
  filter(sample_name %in% selected_samples) %>%
  left_join(sample_meta, by = c("sample_name" = "sample_id")) %>%
  filter(
    # 1) Species entfernen (nur wenn taxon_rank == "species")
    !(taxon_rank == "species" & taxon %in% species_to_remove_swarm),
    
    # 2) Genus entfernen:
    #    a) wenn taxon_rank == "genus" und taxon selbst in Liste
    !(taxon_rank == "genus"   & taxon %in% genus_to_remove_swarm),
    #    b) wenn taxon_rank == "species" und Genus (aus Binomialname) in Liste
    !(taxon_rank == "species" & genus_from_taxon %in% genus_to_remove_swarm),
    
    # 3) Family entfernen (nur wenn taxon_rank == "family")
    !(taxon_rank == "family"  & taxon %in% family_to_remove_swarm)
  ) %>%
  transmute(
    pipeline   = "SWARM",
    sample_id  = sample_name,
    marker     = marker,
    feature_id = definition,
    count      = count_reads,
    taxon      = taxon,
    taxon_rank = taxon_rank,
    class      = class_name,
    
    species = ifelse(taxon_rank == "species", taxon, NA_character_),
    genus   = ifelse(taxon_rank == "genus",   taxon, genus_name_corrected),
    family  = ifelse(taxon_rank == "family",  taxon, family_name_corrected),
    
    decimalLatitude,
    decimalLongitude
  )

 # 02.3 variable count: richtig zusammenfassen
swarm_clean_long <- swarm_clean_long %>%
  group_by(sample_id, marker, feature_id) %>%
  summarise(
    pipeline   = first(pipeline),
    count      = sum(count, na.rm = TRUE),
    taxon      = first(taxon),
    taxon_rank = first(taxon_rank),
    
    class  = first(class),
    species = first(species),
    genus   = first(genus),
    family  = first(family),
    
    decimalLatitude  = first(decimalLatitude),
    decimalLongitude = first(decimalLongitude),
    .groups = "drop"
  )
swarm_clean_long <- swarm_clean_long %>%
  mutate(class_l = tolower(str_trim(class))) %>%
  filter(
    is.na(class_l) | !(class_l %in% classes_to_remove_swarm_l)
  ) %>%
  select(-class_l)

message("SWARM classes (top):")
print(sort(table(swarm_clean_long$class), decreasing = TRUE)[1:20])

# SWARM QC
message("SWARM: Samples=", n_distinct(swarm_clean_long$sample_id),
        " | Features=", n_distinct(swarm_clean_long$feature_id),
        " | Rows=", nrow(swarm_clean_long))
message("SWARM: Missing coords rows=", sum(is.na(swarm_clean_long$decimalLongitude)))

###############################################################################
# 03) DADA2  —  import + cleaning → dada2_clean_long
###############################################################################
# 03.1 Import (5 CSVs) + Merge
dada2_files <- c(
  file.path(path_dada2, "gbif_fw_na_namibia_2022_euka3.csv"),
  file.path(path_dada2, "gbif_fw_na_namiibia_2022_teleo.csv"),
  file.path(path_dada2, "gbif_fw_na_namiibia_2022_unio.csv"),
  file.path(path_dada2, "gbif_fw_na_namiibia_2022_v05.csv"),
  file.path(path_dada2, "gbif_fw_na_namiibia_2022_vene.csv")
)

dada2_raw <- dada2_files %>%
  set_names() %>%
  map_dfr(~ read_csv(.x, show_col_types = FALSE))

# 03.2 Filter TaxonRank + selected samples
dada2_filt <- dada2_raw %>%
  mutate(
    TaxonRank_l = tolower(TaxonRank),
    sci = str_trim(scientificName),
    genus_from_sci = word(sci, 1)  # für Species-Zeilen "Genus species"
  ) %>%
  filter(TaxonRank_l %in% keep_taxon_ranks) %>%
  filter(eventID %in% selected_samples) %>%
  filter(
    # 1) remove species (nur wenn Rank == species)
    !(TaxonRank_l == "species" & sci %in% species_to_remove_dada2),
    
    # 2) remove genus (wenn Rank == genus ODER Species gehört zu Genus)
    !(TaxonRank_l == "genus"   & sci %in% genus_to_remove_dada2),
    !(TaxonRank_l == "species" & genus_from_sci %in% genus_to_remove_dada2),
    
    # 3) remove family (wenn Rank == family)
    !(TaxonRank_l == "family" & sci %in% family_to_remove_dada2)
  ) %>%
  filter(sci != "Hippopotamus") %>%  # optional: Genus entfernen, wenn du willst
  left_join(sample_meta, by = c("eventID" = "sample_id")) %>%
  mutate(
    decimalLatitude  = decimalLatitude.x,
    decimalLongitude = decimalLongitude.x
  ) %>%
  select(-decimalLatitude.x, -decimalLongitude.x, -TaxonRank_l, -genus_from_sci, -sci)


# 03.3 "Merge duplicates" sauber: pro sample_id + marker + scientificName + TaxonRank summieren
# (Marker bleibt drin, damit nicht Primer gemischt werden)

safe_first <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) NA_character_ else x[[1]]
}

dada2_clean_long <- dada2_filt %>%
  mutate(marker = pcr_primer_name_forward) %>%
  group_by(eventID, marker, TaxonRank, scientificName) %>%
  summarise(
    count = sum(organismQuantity, na.rm = TRUE),
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
    feature_id = scientificName,
    count      = count,
    taxon      = scientificName,
    taxon_rank = TaxonRank,
    class, order,
    species = ifelse(tolower(taxon_rank) == "species", taxon, NA_character_),
    genus   = ifelse(tolower(taxon_rank) == "genus",   taxon, genus),
    family  = ifelse(tolower(taxon_rank) == "family",  taxon, family),
    decimalLatitude, decimalLongitude
  )




# DADA2 QC
message("DADA2: Samples=", n_distinct(dada2_clean_long$sample_id),
        " | Features=", n_distinct(dada2_clean_long$feature_id),
        " | Rows=", nrow(dada2_clean_long))
message("DADA2: Missing coords rows=", sum(is.na(dada2_clean_long$decimalLongitude)))

###############################################################################
# 04) OPTIONAL: Konsistenz-Checks (empfohlen)
###############################################################################

# Check: alle selected samples vorhanden?
missing_swarm <- setdiff(selected_samples, unique(swarm_clean_long$sample_id))
missing_dada2 <- setdiff(selected_samples, unique(dada2_clean_long$sample_id))

if (length(missing_swarm) > 0) message("WARN SWARM missing samples: ", paste(missing_swarm, collapse = ", "))
if (length(missing_dada2) > 0) message("WARN DADA2 missing samples: ", paste(missing_dada2, collapse = ", "))

###############################################################################
# 05) OPTIONAL: Speichern (empfohlen statt CSV)
###############################################################################
# dir.create("data", showWarnings = FALSE)
saveRDS(sample_meta,       "C:/Users/basil/OneDrive/Desktop/Master_Thesis_2.0/DATA/sample_meta.rds")
saveRDS(swarm_clean_long,  "C:/Users/basil/OneDrive/Desktop/Master_Thesis_2.0/DATA/swarm_clean_long.rds")
saveRDS(dada2_clean_long,  "C:/Users/basil/OneDrive/Desktop/Master_Thesis_2.0/DATA/dada2_clean_long.rds")
###############################################################################

unique(dada2_clean_long$species)
unique(dada2_clean_long$genus)
unique(swarm_clean_long$species)
unique(dada2_raw$scientificName)

list(
  species = sort(unique(na.omit(swarm_clean_long$species))),
  genus   = sort(unique(na.omit(swarm_clean_long$genus))),
  family  = sort(unique(na.omit(swarm_clean_long$family)))
)

print(unique(swarm_clean_long$class))

str(swarm_clean_long)
str(dada2_clean_long)

sort(unique(swarm_clean_long$marker))
sort(unique(dada2_clean_long$marker))
str(swarm_raw)
str(dada2_raw)

library(dplyr)

primer_lookup_dada2 <- dada2_raw %>%
  distinct(
    marker = pcr_primer_name_forward,
    primer_F_name = pcr_primer_name_forward,
    primer_F_seq  = pcr_primer_forward,
    primer_R_name = pcr_primer_name_reverse,
    primer_R_seq  = pcr_primer_reverse,
    lib_layout
  ) %>%
  arrange(marker)

primer_lookup_dada2

swarm_raw %>% filter(marker=="Euka3") %>% summarise(n=n(), min=min(length_sequence), med=median(length_sequence), max=max(length_sequence))
dada2_raw %>% filter(pcr_primer_name_forward=="Pera02F") %>% summarise(n=n(), min=min(nchar(DNA_sequence)), med=median(nchar(DNA_sequence)), max=max(nchar(DNA_sequence)))

library(dplyr)

pera_raw <- dada2_raw %>%
  filter(pcr_primer_name_forward == "Pera02F")

pera_raw %>%
  summarise(
    n_rows = n(),
    reads_sum = sum(organismQuantity, na.rm = TRUE),
    n_unique_scientificName = n_distinct(scientificName),
    n_unique_class = n_distinct(class),
    n_unique_phylum = n_distinct(phylum),
    n_unique_kingdom = n_distinct(kingdom)
  )

pera_raw %>%
  count(TaxonRank, sort = TRUE)

pera_raw %>%
  mutate(class = if_else(is.na(class) | class=="", "NA", class)) %>%
  group_by(class) %>%
  summarise(
    reads_sum = sum(organismQuantity, na.rm = TRUE),
    n_taxa = n_distinct(scientificName)
  ) %>%
  arrange(desc(reads_sum)) %>%
  head(25)

str(swarm_clean_long)

library(dplyr)
library(tidyr)

# 1) Filter auf Bivalvia
biv <- swarm_clean_long %>%
  filter(class == "Bivalvia")

# 2) Basis-Zahlen
summary_biv <- biv %>%
  summarise(
    n_rows = n(),
    n_samples = n_distinct(sample_id),
    n_markers = n_distinct(marker),
    n_otus = n_distinct(feature_id),
    
    # wie viele OTUs haben überhaupt eine ID auf dem jeweiligen Rang?
    otus_with_family = n_distinct(feature_id[!is.na(family)]),
    otus_with_genus  = n_distinct(feature_id[!is.na(genus)]),
    otus_with_species= n_distinct(feature_id[!is.na(species)]),
    
    # wie viele verschiedene Namen (Taxa) gibt es je Rang?
    n_family_names  = n_distinct(family, na.rm = TRUE),
    n_genus_names   = n_distinct(genus, na.rm = TRUE),
    n_species_names = n_distinct(species, na.rm = TRUE)
  )

summary_biv

