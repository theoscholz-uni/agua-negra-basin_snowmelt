#================================================================================================================================================================================================
###AGUA NEGRA BECKEN - VORWORT###
#================================================================================================================================================================================================
##DIESER CODE IST IM RAHMEN EINER BACHELORARBEIT AN DER UNIVERSITÄT BONN IM FACHBEREICH GEOGRAPHIE ENTSTANDEN##

##DIESER IST ALS FORTSETZUNG DER ARBEIT VON LUTZ KLEIN, DOKTORAND DER AG HYDROLOGIE, ZU VERSTEHEN!##

##BENÖTIGTE PACKAGES SOWIE DATEN KÖNNEN IN DER README-DATEI DES GITHUB-REPOSITORIES AUFGEFUNDEN WERDEN##

##INHALT:##
#ABSCHNITT A: HYDROLOGISCHE ABGRENZUNG DES AGUA-NEGRA-TEILEINZUGSGEBIETS
#ABSCHNITT B: DOWNLOAD UND PROZESSIERUNG DER KLIMATOLOGISCHEN DATEN DER ERA5-REANALYSE
#ABSCHNITT C: ERSTELLEN EINER BINÄREN SCHNEEMASKE UND ABLEITUNG EINER FRAKTIONALEN SCHNEEBEDECKUNG
#ABSCHNITT D: ABLEITUNG DER POTENZIELLEN DIREKTSTRAHLUNG AUS SONNENGEOMETRISCHEN GRUNDGLEICHUNGEN 
#ABSCHNITT E: MONTE-CARLO-BASIERTE KALIBRIERUG DES MODELLS GEGEN DIE BEOBACHTETE FSCA-ZEITREIHE
#ABSCHNITT F: EXTERNE PLAUSIBILISIERUNG DER ERGEBNISSE MIT HILFE VON SATELLITENBASIERTEN SCHNEEBEDECKUNGSPRODUKTEN (VIIRS)

##ALLE PFADE MÜSSEN MANUELL ANGEPASST WERDEN!##

##KI-NUTZUNG: DIE ERSTELLUNG DES CODES WURDE DURCH DIE TOOLS GITHUB COPILOT, CODEX UND CHATGPT UNTERSTÜTZT!##

##VERWENDETE UND WEITERFÜHRENDE LITERATUR UND QUELLEN SIND IN DER README-DATEI DES GITHUB-REPOSITORIES AUFGEFÜHRT!## 

###WICHTIG: DIE TECHNISCHEN VORAUSSETZUNGEN AN DIESES SKRIPT SIND BEI KOMPLETTEN DURCHLAUF HOCH!##
##DESWEGEN EMPFOHLEN: 32 GB ARBEITSSPEICHER UND SOLIDE CPU##
##IST DIES NICHT GEGEBEN: EMPFEHLENSWERT DIE BERECHNUNGS-ZEITRÄUME ZU MINIMIEREN!##
#================================================================================================================================================================================================
###ABSCHNITT A: HYDROLOGISCHE ABGRENZUNG DES AGUA-NEGRA-TEILEINZUGSGEBIETS###
#================================================================================================================================================================================================
##ANMERKUNG: DAS DEM WIRD HYDROLOGISCH KONDITONIERT, IN SEINE TEILBECKEN ABGEGRENZT UND DER FINALE BEZUGSRAUM ANHAND DER KOORDINATEN DER SCHNEEKAMERA DEFINIERT##

#A.0 Kamerastandort als sf- und SpatVector-Objekt erstellen
VALLEY_observer_coords <- data.frame(lon = -69.800472, lat = -30.178385)
VALLEY_location_sf <- st_as_sf(VALLEY_observer_coords, coords = c("lon", "lat"), crs = 4326)
VALLEY_location_vect <- vect(VALLEY_location_sf)

#A.1 Pfade für Input- und Output-Daten
DEM_AguaNegra_path <- "PFAD/ZUM/DEM/DEM.tif" #PFAD ANPASSEN!
basin_output_dir   <- "PFAD/FÜR/BASIN/OUTPUT/DIRECTORY" #PFAD ANPASSEN!

#A.2 Einlesen des DEMs
DEM_AguaNegra <- rast(DEM_AguaNegra_path)

#A.3 Senkenbehandlung, Fließrichtung, Fließakkumulation und Teilbeckenberechnung
DEM_cellsize_mean <- mean(res(DEM_AguaNegra))
breach_distance_cells <- max(5L, as.integer(round(200 / DEM_cellsize_mean)))

#A.4 Pfade für die Zwischenergebnisse der hydrologischen Abgrenzung
DEM_filled_path      <- file.path(basin_output_dir, "DEM_filled.tif")
flow_direction_path  <- file.path(basin_output_dir, "DEM_FlowD.tif")
flow_accum_path      <- file.path(basin_output_dir, "DEM_FAC_cells.tif")
streams_path         <- file.path(basin_output_dir, "DEM_streams.tif")
subbasins_path       <- file.path(basin_output_dir, "DEM_subbasins.tif")

#A.5 Senkenbehandlung mit Least-Cost-Breaching
wbt_breach_depressions_least_cost(
  dem    = DEM_AguaNegra_path,
  output = DEM_filled_path,
  dist   = breach_distance_cells)

#A.6 Berechnung der Fließrichtung
wbt_d8_pointer(
  dem    = DEM_filled_path,
  output = flow_direction_path)

#A.7 Berechnung der Fließakkumulation in Zellen
wbt_d8_flow_accumulation(
  input    = DEM_filled_path,
  output   = flow_accum_path,
  out_type = "cells")

#A.8 Extraktion der Fließwege mit einem Schwellenwert von 2 km² Einzugsgebiet
stream_threshold_cells <- round((2.0 * 1e6) / DEM_cellsize_mean)

#A.9 Extraktion der Fließwege
wbt_extract_streams(
  flow_accum = flow_accum_path,
  output     = streams_path,
  threshold  = stream_threshold_cells)

#A.10 Berechnung der Teilbecken
wbt_subbasins(
  d8_pntr = flow_direction_path,
  streams = streams_path,
  output  = subbasins_path)

#A.11 Teilbecken an der Position des VALLEY-Beobachtungspunktes extrahieren
AguaNegra_subbasin_raster <- rast(subbasins_path)
VALLEY_location_DEM_crs <- project(VALLEY_location_vect, crs(AguaNegra_subbasin_raster))

VALLEY_subbasin_extract <- terra::extract(AguaNegra_subbasin_raster, VALLEY_location_DEM_crs)
VALLEY_subbasin_id      <- VALLEY_subbasin_extract[1, 2]

VALLEY_basin_mask <- AguaNegra_subbasin_raster == VALLEY_subbasin_id
VALLEY_basin_mask <- ifel(VALLEY_basin_mask, 1, NA)

plot(VALLEY_basin_mask)
points(VALLEY_location_DEM_crs, col = "blue", pch = 16)

#A.12 Umwandlung der Maske in ein Polygon
VALLEY_basin_poly <- as.polygons(
  VALLEY_basin_mask,
  dissolve = TRUE,
  na.rm = TRUE)

#A.13 Export des Einzugsgebietspolygons als Shapefile
#OPTIONALER SCHRITT!
writeVector(
  VALLEY_basin_poly,
  file.path(basin_output_dir, "VALLEY_basin.shp"),
  overwrite = TRUE)

#A.13 Zuschneiden und Maskieren des DEMs auf das VALLEY-Becken
DEM_VALLEY <- crop(DEM_AguaNegra, VALLEY_basin_poly)
DEM_VALLEY <- mask(DEM_VALLEY, VALLEY_basin_poly)
DEM_VALLEY <- trim(DEM_VALLEY)

#0.14 Export des zugeschnittenen und maskierten DEMs als GeoTIFF
#OPTIONALER SCHRITT!
writeRaster(
  DEM_VALLEY,
  file.path(basin_output_dir, "DEM_VALLEY.tif"),
  overwrite = TRUE)

#MAIN OUTPUT: VALLEY_basin_poly als SpatVector-Objekt, das das hydrologisch abgegrenzte Einzugsgebiet des VALLEY-Beobachtungspunktes repräsentiert

#================================================================================================================================================================================================
###ABSCHNITT B: PROZESSIERUNG DER KLIMATOLOGISCHEN DATEN DER ERA5-REANALYSE###
#================================================================================================================================================================================================
#================================================================================================================================================================================================
###SCHRITT B.1: VORBEREITUNGEN###
#================================================================================================================================================================================================
##ANMERKUNG: FÜR DEN DOWNLOAD DER ERA5-DATEN MUSS ZUERST DAS EINZUGSGEBIET ABGEGRENZT UND DIE BOUNDING BOX IN WGS 84 BESTIMMT WERDEN!##
#DIE BOUNDING BOX KANN AUCH FÜR DEN DOWNLOAD DER VIIRS SCHNEEBEDECKUNG IN ABSCHNITT F GENUTZT WERDEN#

#B.1.1 Einzugsgebiet als Vektor bereitstellen (für die spätere Verarbeitung der ERA5-Daten)
VALLEY_basin_vect <- VALLEY_basin_poly 

#B.1.2 Bounding Box des bereits abgegrenzten Einzugsgebietes in WGS 84 bestimmen
VALLEY_basin_wgs84 <- st_transform(st_as_sf(VALLEY_basin_poly), 4326)
AN_basin_bbox <- st_bbox(VALLEY_basin_wgs84)

AN_basin_bbox_cds <- c(
  North = unname(AN_basin_bbox["ymax"]),
  West  = unname(AN_basin_bbox["xmin"]),
  South = unname(AN_basin_bbox["ymin"]),
  East  = unname(AN_basin_bbox["xmax"]))

##WICHTIG: DAS SKRIPT ERWARTET DIE DATEN IN FORM VON .NC-FILES VORLIEGEND IN ZWEI ORDNERN:##
#1. FÜR DEN GESAMTNIEDERSCHLAG (TOTAL PRECIPITATION)
#2. FÜR DIE LUFTTEMPERATUR (TEMPERATURE 2m)

##DIESE MÜSSEN MONATLICH ALS NETCDF-DATEIEN HERUNTERGELADEN UND IN DEN JEWEILIGEN ORDNERN ABGELEGT WERDEN!##

#VERFÜGBAR NACH: MUÑOZ SABATER, J. (2019): ERA5-Land hourly data from 1950 to present. Copernicus Climate Change Service (C3S) Climate Data Store (CDS). [Online-Publikation]. [Zugriff am: 11.04.2026]. Verfügbar unter: DOI: 10.24381/cds.e2161bac.

#================================================================================================================================================================================================
###SCHRITT B.2: AUFRUFEN DER DATEN UND VORBEREITUNG DES EINZUGSGEBIETS-LAYERS###
#================================================================================================================================================================================================

#B.2.1 Directories der heruntergeladenen ERA5-Daten festlegen
ERA5_precip_dir <- "ORNDER/PFAD/ZU/ERA5 - PRECIPITATION" #PFAD ANPASSEN!
ERA5_temp_dir   <- "ORDNER/PFAD/ZU/ERA5 - TEMPERATURE" #PFAD ANPASSEN!

#B.2.2 Listen aller .nc-Dateien im jeweiligen Ordner
ERA5_precip_files <- list.files(ERA5_precip_dir, pattern = "\\.nc$", full.names = TRUE)
ERA5_temp_files   <- list.files(ERA5_temp_dir, pattern = "\\.nc$", full.names = TRUE)

#B.2.3 Referenzlayer der ERA5-Daten einlesen und Einzugsgebiet in dessen CRS projizieren
ERA5_reference_raster <- rast(ERA5_temp_files[1])
crs(ERA5_reference_raster) # MUSS WGS 84 sein!

VALLEY_basin_ERA5_vect <- project(VALLEY_basin_vect, crs(ERA5_reference_raster))
crs(VALLEY_basin_ERA5_vect) # MUSS WGS 84 sein!

#================================================================================================================================================================================================
###SCHRITT B.3: FUNKTIONEN ZUR PROZESSIERUNG DER ERA5-DATEN###
#================================================================================================================================================================================================

#B.3.1 Zeitinformationen aus den ERA5-Layernamen extrahieren und in POSIXct umwandeln
ERA5_extract_layer_time <- function(layer_names) {
  valid_time_seconds <- as.numeric(sub(".*valid_time=", "", layer_names))
  as.POSIXct(valid_time_seconds, origin = "1970-01-01", tz = "UTC")}

#B.3.2 Tageswerte aus akkumulierten ERA5-Variablen ableiten (NIEDERSCHLAG)
ERA5_daily_from_accumulated <- function(monthly_raster, variable_pattern, basin_vect, output_col) {
  variable_raster <- subset(monthly_raster, grepl(variable_pattern, names(monthly_raster)))
  layer_time <- ERA5_extract_layer_time(names(variable_raster))
  
  variable_raster <- crop(variable_raster, basin_vect)
  basin_mean_values <- terra::extract(variable_raster, basin_vect, fun = mean, na.rm = TRUE)[1, -1]
  
  daily_df <- data.frame(
    time = layer_time,
    value_m = as.numeric(basin_mean_values))
  
  daily_df$date <- as.Date(daily_df$time)
  daily_df <- aggregate(value_m ~ date, data = daily_df, max, na.rm = TRUE)
  daily_df[[output_col]] <- daily_df$value_m * 1000
  daily_df$value_m <- NULL
  
  daily_df}

#B.3.3 Tageswerte aus instantanen ERA5-Variablen ableiten (TEMPERATUR)
ERA5_daily_from_instantaneous <- function(monthly_raster, variable_pattern, basin_vect, output_col, multiplier = 1, daily_fun = mean) {
  variable_raster <- subset(monthly_raster, grepl(variable_pattern, names(monthly_raster)))
  layer_time <- ERA5_extract_layer_time(names(variable_raster))
  
  variable_raster <- crop(variable_raster, basin_vect)
  basin_mean_values <- terra::extract(variable_raster, basin_vect, fun = mean, na.rm = TRUE)[1, -1]
  
  daily_df <- data.frame(
    time = layer_time,
    value = as.numeric(basin_mean_values) * multiplier)
  
  daily_df$date <- as.Date(daily_df$time)
  daily_df <- aggregate(value ~ date, data = daily_df, daily_fun, na.rm = TRUE)
  names(daily_df)[2] <- output_col
  
  daily_df}


#B.3.4 Monatliche ERA5-Niederschlagsdatei verarbeiten
ERA5_process_precip_file <- function(nc_file, basin_vect) {
  monthly_raster <- rast(nc_file)
  
  daily_total_precip <- ERA5_daily_from_accumulated(
    monthly_raster   = monthly_raster,
    variable_pattern = "^tp_",
    basin_vect       = basin_vect,
    output_col       = "tp_mm")
  
  daily_total_precip}

#B.3.5 Monatliche ERA5-Temperaturdatei verarbeiten
ERA5_process_temp_file <- function(nc_file, basin_vect) {
  monthly_raster <- rast(nc_file)
  
  daily_temperature_k <- ERA5_daily_from_instantaneous(
    monthly_raster   = monthly_raster,
    variable_pattern = "^t2m_",
    basin_vect       = basin_vect,
    output_col       = "t2m_k",
    multiplier       = 1,
    daily_fun        = mean)
  
  daily_temperature_k$t2m_c <- as.numeric(daily_temperature_k$t2m_k) - 273.15
  daily_temperature_k$t2m_k <- NULL
  
  daily_temperature_k}

#================================================================================================================================================================================================
###SCHRITT B.4: UMSETZUNG DER PROZESSIERUNG UND FERTIGSTELLEN DER ENDGÜLTIGEN DATAFRAMES###
#================================================================================================================================================================================================

#B.4.1 Alle monatlichen Niederschlagsdateien verarbeiten
ERA5_precip_monthly_list <- lapply(
  ERA5_precip_files,
  ERA5_process_precip_file,
  basin_vect = VALLEY_basin_ERA5_vect)

#B.4.2 Dataframe der Niederschlagsdaten erstellen
ERA5_precip_df <- do.call(rbind, ERA5_precip_monthly_list)
ERA5_precip_df <- ERA5_precip_df[order(ERA5_precip_df$date), ]

#B.4.3 Alle monatlichen Temperaturdateien verarbeiten
ERA5_temp_monthly_list <- lapply(
  ERA5_temp_files,
  ERA5_process_temp_file,
  basin_vect = VALLEY_basin_ERA5_vect)

#4.4 Dataframe der Temperaturdaten erstellen
ERA5_temp_df <- do.call(rbind, ERA5_temp_monthly_list)
ERA5_temp_df <- ERA5_temp_df[order(ERA5_temp_df$date), ]

#4.5 Niederschlags- und Temperaturdaten zusammenführen
ERA5_all_df <- merge(ERA5_precip_df, ERA5_temp_df, by = "date", all = TRUE)

#4.6 Spalten in die richtigen Datentypen umwandeln
ERA5_all_df$date  <- as.Date(ERA5_all_df$date)
ERA5_all_df$tp_mm <- as.numeric(ERA5_all_df$tp_mm)
ERA5_all_df$t2m_c <- as.numeric(ERA5_all_df$t2m_c)

#MAIN OUTPUT: ERA5_all_df mit den Spalten date, tp_mm und t2m_c##

#================================================================================================================================================================================================
###ABSCHNITT C: ABLEITUNG DER FSCA AUS AUFNAHMEN DER SCHNEEKAMERA### 
#================================================================================================================================================================================================
#================================================================================================================================================================================================
###SCHRITT C.0: GEOTIFFS EINMALIG AUF AOI CROPPEN###
#================================================================================================================================================================================================
##ANMERKUNG: IN DIESEM SCHRITT GEHT ES UM DAS CROPPING DER GROSSEN GEOTIFFS AUF DIE AOI DES SICHTFELDES DER KAMERA!##
##DA NACH ABSPRACHE LEDIGLICH DIE BEREITS GECROPPTEN DATEIEN ZUR VERFÜGUNG STEHEN, KANN DIESER SCHRITT ÜBERSPRUNGEN WERDEN!##
##FALLS DOCH GETESTET WERDEN SOLL, BITTE ROHDATEN ANFRGAFEN UND SAFETY-SCHALTER AUF TRUE SETZEN!##

#C.0.1 Directories der unbearbeiteten und bereits gecroppten Blue-Kanal-Geotiffs
VALLEY_blue_input_dir   <- "PFAD/ZU/UNBEARBEITETEN/BLUE/KANAL/GEOTIFFS" #PFAD ANPASSEN! (ORDNER)
VALLEY_blue_cropped_dir <- "PFAD/ZU/GECROPPTEN/BLUE/KANAL/GEOTIFFS" #PFAD ANPASSEN! (ORDNER)

#C.0.2 Wert definieren, der Pixel außerhalb des Sichtfeldes repräsentiert
VALLEY_viewshed_invalid_value <- 0

#C.0.3 AOI des sichtbaren Kamerabereichs einlesen (SHAPEFILE WURDE SEPERAT IN QGIS ERSTELLT)
VALLEY_aoi_path <- "PFAD/ZUR/SHAPEFILE/DES/SiCHTFELDES/DER/KAMERA/AOI_mask_VALLEY_vis.shp" #PFAD ANPASSEN!
VALLEY_aoi_sf   <- st_read(VALLEY_aoi_path, quiet = TRUE)
VALLEY_aoi_vect <- vect(VALLEY_aoi_sf)

#C.0.4 Safety-Schalter: Schritt nicht versehentlich erneut ausführen (ZUM AUSFÜHREN = TRUE)
VALLEY_crop_blue_geotiffs_once <- FALSE

if (VALLEY_crop_blue_geotiffs_once) {
  
  #C.0.5 Liste der unbearbeiteten Blue-Kanal-Geotiffs
  VALLEY_blue_files_raw <- list.files(
    VALLEY_blue_input_dir,
    pattern = "\\.tif$",
    full.names = TRUE)
  
  for (VALLEY_blue_file in VALLEY_blue_files_raw) {
    
    VALLEY_blue_raster <- rast(VALLEY_blue_file)[[1]]
    
    #AOI in das CRS des jeweiligen Geotiffs transformieren
    VALLEY_aoi_same_crs_sf   <- st_transform(VALLEY_aoi_sf, crs(VALLEY_blue_raster))
    VALLEY_aoi_same_crs_vect <- vect(VALLEY_aoi_same_crs_sf)
    
    #Raster auf die AOI zuschneiden und maskieren
    VALLEY_blue_raster_cropped <- crop(VALLEY_blue_raster, VALLEY_aoi_same_crs_vect)
    VALLEY_blue_raster_cropped <- mask(VALLEY_blue_raster_cropped, VALLEY_aoi_same_crs_vect)
    
    #Pixel außerhalb des Sichtfeldes als NA setzen
    VALLEY_blue_raster_cropped[VALLEY_blue_raster_cropped == VALLEY_viewshed_invalid_value] <- NA
    
    #Gecropptes GeoTIFF speichern
    VALLEY_cropped_output_path <- file.path(
      VALLEY_blue_cropped_dir,
      basename(VALLEY_blue_file))
    
    writeRaster(
      VALLEY_blue_raster_cropped,
      VALLEY_cropped_output_path,
      overwrite = TRUE,
      gdal = c("COMPRESS=LZW", "TILED=YES"))}}


#================================================================================================================================================================================================
###SCHRITT C.1: EINLESEN DER GECROPPTEN BLUE-KANAL-GEOTIFFS UND DEFINITION DER OUTPUT-PFADE###
#================================================================================================================================================================================================

#C.1.1 Directory der bereits gecroppten Blue-Kanal-Geotiffs festlegen (erneut falls SCHRITT C.0 nicht durchgeführt wird, ansonsten bereits definiert)
VALLEY_blue_cropped_dir <- "PFAD/ZU/GECROPPTEN/BLUE/KANAL/GEOTIFFS" #PFAD ANPASSEN! (ORDNER)

#C.1.2 Output-Pfade für Schneemasken und fSCA-Zeitreihe festlegen, OPTIONAL: Direktes Erstellen der Ornder
VALLEY_snowmask_output_dir <- "OUTPUT/FÜR/SCHNEEMASKE/FESTLEGEN/ORDNER"
VALLEY_fSCA_output_csv     <- "OUTPUT/FÜR/FSCA/ZEITREIHE/fSCA_VALLEY_results.csv"

dir.create(VALLEY_snowmask_output_dir, recursive = TRUE, showWarnings = FALSE) #OPTIONAL!
dir.create(dirname(VALLEY_fSCA_output_csv), recursive = TRUE, showWarnings = FALSE) #OPTIONAL!

#C.1.3 Liste der bereits gecroppten Blue-Kanal-Geotiffs erstellen
VALLEY_blue_files_cropped <- list.files(
  VALLEY_blue_cropped_dir,
  pattern = "\\.tif$",
  full.names = TRUE)

#================================================================================================================================================================================================
###SCHRITT C.2: HILFSFUNKTIONEN ZUR DATEI- UND HISTOGRAMMVERARBEITUNG###
#================================================================================================================================================================================================
##ANMERKUNG: DIE AUFNAHMEN WURDEN MANUELL NACH FOLGENDER STRUKTUR BENANNT: ID_VALLEY_blue_DDMMYYYY.tif (BEISPIEL: 1_VALLEY_blue_18022023.tif)##

#C.2.1 Bild-ID aus dem Dateinamen extrahieren 
extract_image_id_from_filename <- function(file_path) {
  file_name <- basename(file_path)
  leading_digits <- regexpr("^\\d+", file_name)
  
  if (leading_digits[1] == -1) return(NA_integer_)
  
  as.integer(regmatches(file_name, leading_digits))}

#C.2.2 Aufnahmedatum aus dem Dateinamen extrahieren
extract_image_date_from_filename <- function(file_path) {
  file_name <- basename(file_path)
  date_match <- regexpr("\\d{8}", file_name)
  
  if (date_match[1] == -1) return(as.Date(NA))
  
  date_string <- regmatches(file_name, date_match)
  as.Date(date_string, format = "%d%m%Y")}

#C.2.3 5-Punkt-Glättung des DN-Histogramms im Blaukanal (nach SALVATORI et al. 2011)
smooth_blue_histogram_5point <- function(hist_counts) {
  as.numeric(stats::filter(hist_counts, rep(1 / 5, 5), sides = 2))}

#C.2.4 Schwellenwert nach SALVATORI et al. (2011) bestimmen
#Erster lokaler Minimalwert jenseits von DN >= 127, andernfalls Fallback auf 127
find_salvatori_blue_threshold <- function(blue_values, min_dn = 127L) {
  
  blue_histogram <- hist(
    blue_values,
    breaks = seq(-0.5, 255.5, by = 1),
    plot = FALSE)
  
  hist_counts_raw <- blue_histogram$counts
  hist_counts_smoothed <- smooth_blue_histogram_5point(hist_counts_raw)
  hist_counts_smoothed[is.na(hist_counts_smoothed)] <- hist_counts_raw[is.na(hist_counts_smoothed)]
  
  dn_values <- 0:255
  candidate_indices <- which(dn_values >= min_dn)
  
  selected_threshold <- NA_integer_
  
  for (hist_index in candidate_indices) {
    if (hist_index <= 1 || hist_index >= length(hist_counts_smoothed)) next
    
    if (hist_counts_smoothed[hist_index] < hist_counts_smoothed[hist_index - 1] &&
        hist_counts_smoothed[hist_index] < hist_counts_smoothed[hist_index + 1]) {
      selected_threshold <- dn_values[hist_index]
      break}}
  
  if (is.na(selected_threshold)) {
    selected_threshold <- min_dn}
  
  list(
    threshold        = as.integer(selected_threshold),
    hist_counts_raw  = hist_counts_raw,
    hist_counts_smoothed = hist_counts_smoothed,
    hist_bins        = dn_values)}

#C.2.5 Binäre Schneemaske nach Schwellenwert aus dem Blaukanal ableiten
classify_snow_from_blue_channel <- function(blue_raster, threshold_blue) {
  ifel(is.na(blue_raster), NA, ifel(blue_raster >= threshold_blue, 1, 0))}

#================================================================================================================================================================================================
###SCHRITT C.3: SICHERHEITSCHECK DES WERTEBEREICHS DER BLUE-KANAL-DATEN###
#================================================================================================================================================================================================

#C.3.1 Testweise ein gecropptes Blue-Kanal-Geotiff einlesen
VALLEY_blue_test_raster <- rast(VALLEY_blue_files_cropped[1])[[1]]

#C.3.2 Werte extrahieren und auf gültige numerische Werte beschränken
VALLEY_blue_test_values <- values(VALLEY_blue_test_raster, mat = FALSE)
VALLEY_blue_test_values <- VALLEY_blue_test_values[is.finite(VALLEY_blue_test_values)]

#C.3.3 Wertebereich prüfen
range(VALLEY_blue_test_values, na.rm = TRUE)

#================================================================================================================================================================================================
###SCHRITT C.4: SCHWELLENWERTBESTIMMUNG, BINÄRE SCHNEEKLASSIFIKATION UND BERECHNUNG DER FSCA###
#================================================================================================================================================================================================

#C.4.1 Liste zur Speicherung aller bildspezifischen Ergebnisse vorbereiten
VALLEY_fSCA_results <- vector("list", length(VALLEY_blue_files_cropped))

#C.4.2 Progress-Bar für die Verarbeitung aller gecroppten Blue-Kanal-Geotiffs
VALLEY_fSCA_progress <- txtProgressBar(
  min = 0,
  max = length(VALLEY_blue_files_cropped),
  style = 3)

#C.4.3 Alle gecroppten Blue-Kanal-Geotiffs durchlaufen
for (file_index in seq_along(VALLEY_blue_files_cropped)) {
  setTxtProgressBar(VALLEY_fSCA_progress, file_index)
  
  VALLEY_blue_file <- VALLEY_blue_files_cropped[file_index]
  VALLEY_blue_raster <- rast(VALLEY_blue_file)[[1]]
  
  VALLEY_blue_values <- values(VALLEY_blue_raster, mat = FALSE)
  VALLEY_blue_values <- as.numeric(VALLEY_blue_values)
  VALLEY_blue_values <- VALLEY_blue_values[is.finite(VALLEY_blue_values)]
  
  #Falls keine gültigen Pixel vorliegen, wird ein leerer Ergebnisdatensatz erzeugt
  if (length(VALLEY_blue_values) == 0) {
    VALLEY_fSCA_results[[file_index]] <- data.frame(
      img_id         = extract_image_id_from_filename(VALLEY_blue_file),
      date           = extract_image_date_from_filename(VALLEY_blue_file),
      file           = basename(VALLEY_blue_file),
      threshold_blue = NA_integer_,
      fSCA           = NA_real_,
      n_valid        = 0L,
      qc_flag        = "no_valid_pixels")
    next}
  
  #Schwellenwert nach SALVATORI et al. (2011) auf Basis des Blaukanals bestimmen
  VALLEY_threshold_result <- find_salvatori_blue_threshold(
    blue_values = VALLEY_blue_values,
    min_dn = 127L)
  
  VALLEY_blue_threshold <- VALLEY_threshold_result$threshold
  
  #Binäre Schneemaske aus dem Blaukanal ableiten
  VALLEY_snowmask_raster <- classify_snow_from_blue_channel(
    blue_raster = VALLEY_blue_raster,
    threshold_blue = VALLEY_blue_threshold)
  
  #fSCA als Anteil der als Schnee klassifizierten gültigen Pixel berechnen
  VALLEY_fSCA_value <- global(
    VALLEY_snowmask_raster,
    "mean",
    na.rm = TRUE)[1, 1]
  
  #Binäre Schneemaske als GeoTIFF speichern
  VALLEY_snowmask_output_path <- file.path(
    VALLEY_snowmask_output_dir,
    paste0(
      tools::file_path_sans_ext(basename(VALLEY_blue_file)),
      "_snowmask.tif"))
  
  writeRaster(
    VALLEY_snowmask_raster,
    VALLEY_snowmask_output_path,
    overwrite = TRUE,
    gdal = c("COMPRESS=LZW", "TILED=YES")) #KOMPRESSION UND TILING FÜR EFFIZIENTE SPEICHERUNG
  
  #Ergebniszeile für das aktuelle Bild speichern
  VALLEY_fSCA_results[[file_index]] <- data.frame(
    img_id         = extract_image_id_from_filename(VALLEY_blue_file),
    date           = extract_image_date_from_filename(VALLEY_blue_file),
    file           = basename(VALLEY_blue_file),
    threshold_blue = VALLEY_blue_threshold,
    fSCA           = as.numeric(VALLEY_fSCA_value),
    n_valid        = length(VALLEY_blue_values),
    qc_flag        = "ok")}

close(VALLEY_fSCA_progress)

#MAIN OUTPUT: VALLEY_fSCA_results als Liste mit einem Dataframe pro Bild, der die Bild-ID, das Aufnahmedatum, den Dateinamen, den ermittelten Schwellenwert im Blaukanal, den berechneten fSCA-Wert, die Anzahl der gültigen Pixel und eine QC-Flag enthält#

#================================================================================================================================================================================================
###SCHRITT C.5: ZUSAMMENFÜHREN, SPEICHERN UND FORMATIEREN DER FSCA-ZEITREIHE###
#================================================================================================================================================================================================

#C.5.1 Alle bildspezifischen Ergebnisse in einem Dataframe zusammenführen und sortieren (Nach Date und Bild-ID)
VALLEY_fSCA_df <- bind_rows(VALLEY_fSCA_results) %>%
  arrange(date, img_id)

#C.5.2 fSCA-Zeitreihe als .CSV exportieren (OPTIONAL!)
write.csv(VALLEY_fSCA_df, VALLEY_fSCA_output_csv, row.names = FALSE)

#C.5.3 fSCA-Werte zusätzlich mit sechs Nachkommastellen formatieren (OPTIONAL, da fSCA-Werte bereits numerisch vorliegen, aber für die Darstellung in formatiert werden können)
VALLEY_fSCA_df$fSCA_decimal <- formatC(
  VALLEY_fSCA_df$fSCA,
  format = "f",
  digits = 6)

#C.5.4 Kurze Kontrolle der Ergebnisse
print(head(VALLEY_fSCA_df, 10))
summary(VALLEY_fSCA_df$fSCA)
summary(VALLEY_fSCA_df$threshold_blue)

#================================================================================================================================================================================================
###SCHRITT C.6: EINBRINGEN DER MANUELLEN QC-FLAGS BEZÜGLICH SCHATTEN, NEBEL UND KONTRAST###
#================================================================================================================================================================================================
##ANMERKUNG: DIE MANUELL DEFINIERTEN QC-FLAGS BASIEREN AUF EINER SORGFÄLTIGEN DURCHSICHTUNG ALLER AUFNAHMEN UND EINEM ABGLEICH DIESER MIT DER FSCA!##

#6.1 Manuell definierte QC-Tage für Schatten, Nebel und Kontrast
VALLEY_shadow_drop_dates <- as.Date(c(
  "2023-02-23", "2023-03-02", "2023-03-18",
  "2023-05-27", "2023-06-18", "2023-07-18", "2023-08-28", "2023-09-01", "2023-10-17",
  "2023-10-21", "2023-10-26", "2023-11-12", "2023-11-14",
  "2024-01-24", "2024-01-27", "2024-01-28", "2024-01-29", "2024-01-30",
  "2024-02-05", "2024-02-06", "2024-02-08", "2024-02-09"))

VALLEY_shadow_keep_dates <- as.Date(c(
  "2023-02-21", "2023-03-26", "2023-04-06", "2023-04-16", "2023-04-20", "2023-04-22",
  "2023-04-29", "2023-05-21", "2023-06-04", "2023-06-09", "2023-06-15", "2023-07-07",
  "2023-07-08", "2023-07-17", "2023-07-24", "2023-07-31", "2023-08-25", "2023-09-17",
  "2023-09-23", "2023-10-02", "2023-10-15", "2023-10-23", "2023-10-25", "2023-10-27",
  "2023-11-02", "2023-11-16", "2024-01-25", "2024-01-26", "2024-02-01", "2024-02-03",
  "2024-02-13"))

VALLEY_fog_drop_dates <- as.Date(c(
  "2023-06-10", "2023-08-07", "2023-08-17", "2023-09-08", "2023-10-29",
  "2023-11-01", "2023-11-10", "2023-12-14", "2023-12-26", "2023-12-30"))

VALLEY_fog_keep_dates <- as.Date(c(
  "2023-05-22", "2023-07-19", "2023-08-23",
  "2023-08-26", "2023-09-02", "2023-12-28"))

VALLEY_contrast_drop_dates <- as.Date("2023-03-28")

#C.6.2 Hilfsfunktion zur Erstellung der manuellen QC-Flags
create_manual_qc_flag_df <- function(flag_dates, flag_type, keep_status) {
  data.frame(
    date = flag_dates,
    qc_type = flag_type,
    qc_keep = keep_status,
    qc_label = paste0(
      "QC_", flag_type, "_",
      ifelse(keep_status == "yes", "keep", "drop")))}

#C.6.3 Alle manuellen QC-Flags in einem Dataframe zusammenführen
VALLEY_manual_qc_flags <- bind_rows(
  create_manual_qc_flag_df(VALLEY_shadow_drop_dates, "shadow", "no"),
  create_manual_qc_flag_df(VALLEY_shadow_keep_dates, "shadow", "yes"),
  create_manual_qc_flag_df(VALLEY_fog_drop_dates, "fog", "no"),
  create_manual_qc_flag_df(VALLEY_fog_keep_dates, "fog", "yes"),
  create_manual_qc_flag_df(VALLEY_contrast_drop_dates, "contrast", "no")) %>%
  arrange(date, qc_type)

#C.6.4 Manuelle QC-Flags an die rohe fSCA-Zeitreihe anhängen
VALLEY_fSCA_df <- VALLEY_fSCA_df %>%
  left_join(VALLEY_manual_qc_flags, by = "date") %>%
  mutate(
    qc_flag_manual = ifelse(is.na(qc_label), "none", qc_label),
    qc_type = ifelse(is.na(qc_type), "none", qc_type),
    qc_keep = ifelse(is.na(qc_keep), "yes", qc_keep)) %>%
  select(-qc_label)

#C.6.5 Rohwerte und QC-gefilterte Werte in einem separaten Dataframe bereitstellen
VALLEY_fSCA_df_qc <- VALLEY_fSCA_df %>%
  mutate(
    fSCA_raw = fSCA,
    fSCA_qc_filtered = ifelse(qc_keep == "yes", fSCA, NA_real_))

#C.6.6 Anzahl der verfügbaren Roh- und QC-gefilterten Werte bestimmen
VALLEY_fSCA_n_raw <- sum(!is.na(VALLEY_fSCA_df_qc$fSCA_raw))
VALLEY_fSCA_n_qc  <- sum(!is.na(VALLEY_fSCA_df_qc$fSCA_qc_filtered))

#C.6.7 Vergleichs-Dataframe für die Darstellung der ungefilterten und QC-gefilterten fSCA erstellen
VALLEY_fSCA_comparison_df <- bind_rows(
  VALLEY_fSCA_df_qc %>%
    transmute(
      date = date,
      fSCA = fSCA_raw,
      panel = "raw"),
  VALLEY_fSCA_df_qc %>%
    filter(!is.na(fSCA_qc_filtered)) %>%
    transmute(
      date = date,
      fSCA = fSCA_qc_filtered,
      panel = "qc")) %>%
  group_by(panel) %>%
  mutate(n_panel = sum(!is.na(fSCA))) %>%
  ungroup() %>%
  mutate(
    panel_label = case_when(
      panel == "raw" ~ paste0("Ungefilterte fSCA | n = ", n_panel),
      panel == "qc"  ~ paste0("QC-gefilterte fSCA | n = ", n_panel)),
    panel_label = factor(
      panel_label,
      levels = c(
        paste0("Ungefilterte fSCA | n = ", VALLEY_fSCA_n_raw),
        paste0("QC-gefilterte fSCA | n = ", VALLEY_fSCA_n_qc))))

#C.6.8 Vergleich der ungefilterten und QC-gefilterten fSCA-Zeitreihe visualisieren
ggplot(VALLEY_fSCA_comparison_df, aes(x = date, y = fSCA)) +
  geom_line(linewidth = 0.5, color = "black", na.rm = FALSE) +
  facet_wrap(~ panel_label, ncol = 1, scales = "fixed") +
  labs(
    x = "Datum",
    y = "fraktionale Schneebedeckung (-)"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 11, face = "bold"),
    panel.grid.minor = element_blank())

#C.6.9 Vollständige Datumsskala für den Untersuchungszeitraum anlegen
VALLEY_fSCA_full_dates <- data.frame(
  date = seq(as.Date("2023-02-18"), as.Date("2024-02-25"), by = "day"))

#C.6.10 Tage identifizieren, die nach der QC-Filterung ausgeschlossen wurden
VALLEY_qc_excluded_days <- VALLEY_fSCA_df_qc %>%
  filter(is.na(fSCA_qc_filtered)) %>%
  mutate(
    fSCA_raw_decimal = formatC(fSCA_raw, format = "f", digits = 6)) %>%
  select(
    date, file, qc_flag_manual, qc_type, qc_keep,
    fSCA_raw, fSCA_raw_decimal, fSCA_qc_filtered) %>%
  arrange(date)

#C.6.11 Tage ohne verfügbares Bild im Untersuchungszeitraum identifizieren
VALLEY_missing_days <- VALLEY_fSCA_full_dates %>%
  anti_join(
    VALLEY_fSCA_df_qc %>% select(date),
    by = "date") %>%
  arrange(date)

print(VALLEY_qc_excluded_days)
print(VALLEY_missing_days)

#MAIN OUTPUT: VALLEY_fSCA_df_qc mit den Spalten date, file, qc_flag_manual, qc_type, qc_keep, fSCA_raw, fSCA_qc_filtered und den QC-Flags für Schatten, Nebel und Kontrast#

#================================================================================================================================================================================================
###SCHRITT C.7: AUSWERTUNG DER BILDSPEZIFISCHEN SCHWELLENWERTE###
#================================================================================================================================================================================================

#C.7.1 Absolute Häufigkeit aller tatsächlich vorkommenden Schwellenwerte
VALLEY_threshold_frequency <- VALLEY_fSCA_df_qc %>%
  filter(!is.na(threshold_blue)) %>%
  count(threshold_blue, name = "n_days") %>%
  arrange(threshold_blue)

print(VALLEY_threshold_frequency)

#C.7.2 Histogramm der qualitätsgeprüften Schwellenwerte mit Quartilen
VALLEY_threshold_histogram_df <- VALLEY_fSCA_df_qc %>%
  filter(qc_keep == "yes", !is.na(threshold_blue))

VALLEY_threshold_xmin <- floor(min(VALLEY_threshold_histogram_df$threshold_blue, na.rm = TRUE) / 5) * 5
VALLEY_threshold_xmax <- ceiling(max(VALLEY_threshold_histogram_df$threshold_blue, na.rm = TRUE) / 5) * 5

VALLEY_threshold_quartiles <- quantile(
  VALLEY_threshold_histogram_df$threshold_blue,
  probs = c(0.25, 0.50, 0.75),
  na.rm = TRUE)

VALLEY_threshold_hist_ymax <- hist(
  VALLEY_threshold_histogram_df$threshold_blue,
  breaks = seq(VALLEY_threshold_xmin, VALLEY_threshold_xmax + 5, by = 5),
  plot = FALSE)$counts |> max(na.rm = TRUE)

VALLEY_threshold_quartiles_df <- data.frame(
  x = as.numeric(VALLEY_threshold_quartiles),
  x_label = c(
    as.numeric(VALLEY_threshold_quartiles[1]) - 1.5,
    as.numeric(VALLEY_threshold_quartiles[2]) + 1.2,
    as.numeric(VALLEY_threshold_quartiles[3]) + 1.2),
  hjust = c(1, 0, 0),
  label = c(
    paste0("Q1 = ", round(VALLEY_threshold_quartiles[1], 0)),
    paste0("Median = ", round(VALLEY_threshold_quartiles[2], 0)),
    paste0("Q3 = ", round(VALLEY_threshold_quartiles[3], 0))),
  y = rep(VALLEY_threshold_hist_ymax * 1.15, 3))

VALLEY_threshold_q1_q3_df <- VALLEY_threshold_quartiles_df[c(1, 3), ]
VALLEY_threshold_median_df <- VALLEY_threshold_quartiles_df[2, ]

#C.7.3 Plot des Histogramms der Schwellenwerte mit markierten Quartilen und Median
ggplot(VALLEY_threshold_histogram_df, aes(x = threshold_blue)) +
  geom_histogram(
    binwidth = 5,
    boundary = VALLEY_threshold_xmin,
    color = "black",
    fill = "grey80",
    linewidth = 0.3) +
  geom_vline(
    data = VALLEY_threshold_q1_q3_df,
    aes(xintercept = x),
    linetype = "dashed",
    linewidth = 1.0,
    color = "red") +
  geom_vline(
    data = VALLEY_threshold_median_df,
    aes(xintercept = x),
    linetype = "solid",
    linewidth = 1.0,
    color = "red") +
  geom_text(
    data = VALLEY_threshold_q1_q3_df,
    aes(x = x_label, y = y, label = label, hjust = hjust),
    vjust = 1,
    size = 3.6) +
  geom_text(
    data = VALLEY_threshold_median_df,
    aes(x = x_label, y = y, label = label, hjust = hjust),
    vjust = 1,
    size = 3.6,
    fontface = "bold") +
  scale_x_continuous(
    breaks = seq(VALLEY_threshold_xmin, VALLEY_threshold_xmax, by = 10),
    limits = c(VALLEY_threshold_xmin, VALLEY_threshold_xmax)) +
  expand_limits(y = VALLEY_threshold_hist_ymax * 1.08) +
  labs(
    x = "Bildspezifischer Schwellenwert im Blaukanal (DN)",
    y = "Anzahl der Tage") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    panel.grid.minor = element_blank())

#C.7.4 Prüfung der Anzahl der unterschiedlichen Schwellenwerte und der häufigsten Schwellenwerte
VALLEY_threshold_n_distinct <- VALLEY_fSCA_df_qc %>%
  filter(!is.na(threshold_blue)) %>%
  summarise(n_distinct_thresholds = n_distinct(threshold_blue))

print(VALLEY_threshold_n_distinct)

VALLEY_threshold_mode_table <- VALLEY_threshold_frequency %>%
  arrange(desc(n_days), threshold_blue)

print(head(VALLEY_threshold_mode_table, 10))

#C.7.5 Prüfung des Fallback-Werts 127 DN
VALLEY_threshold_127_n <- sum(VALLEY_fSCA_df_qc$threshold_blue == 127, na.rm = TRUE)

VALLEY_threshold_127_days <- VALLEY_fSCA_df_qc %>%
  filter(threshold_blue == 127) %>%
  select(date, file, threshold_blue, fSCA, qc_flag_manual, qc_keep) %>%
  arrange(date)

VALLEY_threshold_127_n
print(VALLEY_threshold_127_days)

#C.7.6 Höchste und niedrigste Schwellenwerte
VALLEY_threshold_highest_days <- VALLEY_fSCA_df_qc %>%
  arrange(desc(threshold_blue)) %>%
  select(date, file, threshold_blue, fSCA, qc_flag_manual, qc_keep)

print(head(VALLEY_threshold_highest_days, 15))

VALLEY_threshold_lowest_days <- VALLEY_fSCA_df_qc %>%
  arrange(threshold_blue) %>%
  select(date, file, threshold_blue, fSCA, qc_flag_manual, qc_keep)

print(head(VALLEY_threshold_lowest_days, 15))

#================================================================================================================================================================================================
###SCHRITT C.8: BEISPIELPLOT DES GEGLÄTTETEN HISTOGRAMMS FÜR EINEN TAG###
#================================================================================================================================================================================================

#C.8.1 Beispielhafter Tag für die Darstellung des geglätteten Blaukanal-Histogramms (Kann problemlos variiert werden!)
VALLEY_example_date <- as.Date("2023-06-4")

VALLEY_example_file <- VALLEY_blue_files_cropped[
  extract_image_date_from_filename(VALLEY_blue_files_cropped) == VALLEY_example_date][1]

VALLEY_example_blue_raster <- rast(VALLEY_example_file)[[1]]

VALLEY_example_blue_values <- values(VALLEY_example_blue_raster, mat = FALSE)
VALLEY_example_blue_values <- as.numeric(VALLEY_example_blue_values)
VALLEY_example_blue_values <- VALLEY_example_blue_values[is.finite(VALLEY_example_blue_values)]

#C.8.2 Schwellenwertbestimmung und Histogrammauswertung für den Beispieltag
VALLEY_example_threshold_result <- find_salvatori_blue_threshold(
  blue_values = VALLEY_example_blue_values,
  min_dn = 127L)

VALLEY_example_histogram_df <- data.frame(
  dn = VALLEY_example_threshold_result$hist_bins,
  counts_raw = VALLEY_example_threshold_result$hist_counts_raw,
  counts_smoothed = VALLEY_example_threshold_result$hist_counts_smoothed)

VALLEY_example_threshold <- VALLEY_example_threshold_result$threshold

#C.8.3 Plot des geglätteten Blaukanal-Histogramms mit markiertem Schwellenwert und Fallback-Wert
ggplot(VALLEY_example_histogram_df, aes(x = dn, y = counts_smoothed)) +
  geom_col(width = 1, fill = "grey70", color = "grey40") +
  geom_vline(
    xintercept = VALLEY_example_threshold,
    linewidth = 0.6,
    linetype = "solid") +
  geom_vline(
    xintercept = 127,
    linewidth = 0.6,
    linetype = "dashed") +
  annotate(
    "text",
    x = VALLEY_example_threshold,
    y = max(VALLEY_example_histogram_df$counts_smoothed, na.rm = TRUE),
    label = paste0("Gewählter Schwellenwert = ", VALLEY_example_threshold),
    hjust = -0.05,
    vjust = -0.3,
    size = 3) +
  annotate(
    "text",
    x = 127,
    y = max(VALLEY_example_histogram_df$counts_smoothed, na.rm = TRUE) * 0.92,
    label = "Minimaler Schwellenwert = 127",
    hjust = -0.05,
    vjust = -0.3,
    size = 3) +
  labs(
    x = "Blau-Kanal DN-Werte",
    y = "Anzahl der Pixel") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank())


#================================================================================================================================================================================================
###SCHRITT C.9: FESTLEGUNG DER FINALEN FSCA-BEOBACHTUNGSREIHE###
#================================================================================================================================================================================================

#C.9.1 Finales QC-gefiltertes fSCA-Observationsfeld festlegen 
VALLEY_fSCA_df_qc <- VALLEY_fSCA_df_qc %>%
  mutate(
    fSCA_obs = fSCA_qc_filtered)

summary(VALLEY_fSCA_df_qc$fSCA_obs)
sum(!is.na(VALLEY_fSCA_df_qc$fSCA_obs))

#================================================================================================================================================================================================
###ABSCHNITT D: ABLEITUNG DER POTENZIELLEN DIREKTSTRAHLUNG####
#================================================================================================================================================================================================
#================================================================================================================================================================================================
###SCHRITT D.1: VORBEREITUNGEN (PFADE + ZEITRAUM + PARAMETER)###
#================================================================================================================================================================================================

#D.1.1 DEM, AOI und Output-Verzeichnis definieren
VALLEY_Ipot_dem_base <- DEM_VALLEY
VALLEY_Ipot_aoi_vect <- VALLEY_basin_poly
VALLEY_Ipot_aoi_vect <- project(VALLEY_Ipot_aoi_vect, crs(VALLEY_Ipot_dem_base))
VALLEY_Ipot_aoi_sf <- st_as_sf(VALLEY_Ipot_aoi_vect)

#D.1.2 Outpit-Directory für Raster der potenziellen Direktstrahlung (OPTIONAL! Schreiben der Files kann nachfolgend deaktiviert werden!)
VALLEY_Ipot_output_dir <- "C:/Users/theos/Desktop/Bachelorarbeit_AguaNegra/Snow Cameras/images/VALLEY/10_VALLEY_radiation_daily"
dir.create(VALLEY_Ipot_output_dir, recursive = TRUE, showWarnings = FALSE)

#D.1.2 Untersuchungszeitraum auf Tagesbasis festlegen
VALLEY_Ipot_date_start <- as.Date("2023-02-18")
VALLEY_Ipot_date_end   <- as.Date("2024-02-25")
VALLEY_Ipot_dates_seq  <- seq(VALLEY_Ipot_date_start, VALLEY_Ipot_date_end, by = "day")

#D.1.3 Stündliche Zeitschritte festlegen
VALLEY_Ipot_hours_seq <- 0:23

#D.1.4 Solarkonstante und psi alpha festlegen
VALLEY_Ipot_solar_constant <- 1368
VALLEY_Ipot_psi_alpha_default  <- 0.75

#D.1.5 Zeitzone festlegen
VALLEY_Ipot_tz <- "America/Argentina/San_Juan" 

make_local_time <- function(date, hour, tz = VALLEY_Ipot_tz) {
  as.POSIXct(paste(date, sprintf("%02d:00:00", hour)), tz = tz)}

#================================================================================================================================================================================================
###SCHRITT D.2: DEM LADEN, AOI MASKIEREN, SLOPE/ASPECT BERECHNEN (BOGENMAß)###
#================================================================================================================================================================================================

#D.2.1 DEM auf die AOI zuschneiden und maskieren
VALLEY_Ipot_dem <- crop(VALLEY_Ipot_dem_base, VALLEY_Ipot_aoi_vect) |> mask(VALLEY_Ipot_aoi_vect)

#D.2.2 Hangneigung und Exposition in Grad berechnen
VALLEY_Ipot_slope_deg  <- terrain(VALLEY_Ipot_dem, v = "slope", unit = "degrees")
VALLEY_Ipot_aspect_deg <- terrain(VALLEY_Ipot_dem, v = "aspect", unit = "degrees")

#D.2.3 Hilfsfunktion zur Umrechnung von Grad in Radiant
deg2rad <- function(x) x * pi / 180

#D.2.4 Hangneigung und Exposition in Radiant umrechnen
VALLEY_Ipot_slope_rad  <- deg2rad(VALLEY_Ipot_slope_deg)
VALLEY_Ipot_aspect_rad <- deg2rad(VALLEY_Ipot_aspect_deg)

#D.2.5 Trigonometrische Raster für die spätere Berechnung vorberechnen
VALLEY_Ipot_cos_slope  <- cos(VALLEY_Ipot_slope_rad)
VALLEY_Ipot_sin_slope  <- sin(VALLEY_Ipot_slope_rad)
VALLEY_Ipot_cos_aspect <- cos(VALLEY_Ipot_aspect_rad)
VALLEY_Ipot_sin_aspect <- sin(VALLEY_Ipot_aspect_rad)

#================================================================================================================================================================================================
###SCHRITT D.3: REPRÄSENTATIVE LAT/LON AUS AOI-CENTROID (WGS84)###
#================================================================================================================================================================================================
##ANMERKUNG: FESTLEGEN EINES ZENTRALEN REFERENZPUNKTES!##

#D.3.1 AOI-Zentrum bestimmen
VALLEY_Ipot_aoi_centroid <- st_centroid(st_union(VALLEY_Ipot_aoi_sf))

#D.3.2 AOI-Zentrum nach WGS 84 transformieren
VALLEY_Ipot_aoi_centroid_wgs84 <- st_transform(VALLEY_Ipot_aoi_centroid, 4326)

#D.3.3 Geographische Länge und Breite des AOI-Zentrums extrahieren
VALLEY_Ipot_centroid_coords <- st_coordinates(VALLEY_Ipot_aoi_centroid_wgs84)

VALLEY_Ipot_lon_deg <- as.numeric(VALLEY_Ipot_centroid_coords[1])
VALLEY_Ipot_lat_deg <- as.numeric(VALLEY_Ipot_centroid_coords[2])

#================================================================================================================================================================================================
###SCHRITT D.4: FUNKTIONEN ZUR SONNENGEOMETRIE UND LUFTMASSE###
#================================================================================================================================================================================================
##ANMERKUNG: BERECHNUNG AUF BASIS SONNENGEOMETRISCHER NÄHERUNGS-FUNKTIONEN NACH IQBAL (1983) und SZOKOLAY (2007)##

#D.4.1 Solare Deklination und Equation of Time berechnen (in Abhängigkeit von Tag des Jahres und Tagwinkel gamma)
compute_solar_declination_eot <- function(day_of_year) {
  
  #Tagwinkel nach 365-Tage-Näherung
  day_angle_gamma <- 2 * pi * (day_of_year - 1) / 365
  
  #Solare Deklination in Radiant
  solar_declination_rad <- 0.006918 -
    0.399912 * cos(day_angle_gamma) + 0.070257 * sin(day_angle_gamma) -
    0.006758 * cos(2 * day_angle_gamma) + 0.000907 * sin(2 * day_angle_gamma) -
    0.002697 * cos(3 * day_angle_gamma) + 0.00148  * sin(3 * day_angle_gamma)
  
  #Equation of Time in Minuten
  equation_of_time_min <- 229.18 * (
    0.000075 +
      0.001868 * cos(day_angle_gamma) -
      0.032077 * sin(day_angle_gamma) -
      0.014615 * cos(2 * day_angle_gamma) -
      0.040849 * sin(2 * day_angle_gamma))
  
  list(
    solar_declination_rad = solar_declination_rad,
    equation_of_time_min  = equation_of_time_min)}

#4.2 Sonnenhöhe und solaren Azimut berechnen
#Azimut-Konvention: 0 = Norden, pi/2 = Osten, pi = Süden, 3pi/2 = Westen
#Stundenwinkel: 0 = Sonnenhöchststand, positiv am Nachmittag, negativ am Vormittag

compute_solar_position <- function(date, hour, latitude_deg, longitude_deg, tz_offset_hours) {
  
  day_of_year <- yday(date)
  
  solar_geometry <- compute_solar_declination_eot(day_of_year)
  solar_declination_rad <- solar_geometry$solar_declination_rad
  equation_of_time_min  <- solar_geometry$equation_of_time_min
  
  latitude_rad <- deg2rad(latitude_deg)
  
  #Lokaler Standardmeridian der Zeitzone
  local_standard_meridian_deg <- 15 * tz_offset_hours
  
  #Umrechnung von lokaler Uhrzeit in True Solar Time
  clock_time_min <- hour * 60
  solar_time_correction_min <- equation_of_time_min + 4 * (longitude_deg - local_standard_meridian_deg)
  true_solar_time_min <- clock_time_min + solar_time_correction_min
  
  #Stundenwinkel in Radiant
  hour_angle_rad <- deg2rad((true_solar_time_min / 4) - 180)
  
  #Sonnenhöhe in Radiant
  sin_solar_altitude <- sin(latitude_rad) * sin(solar_declination_rad) +
    cos(latitude_rad) * cos(solar_declination_rad) * cos(hour_angle_rad)
  
  solar_altitude_rad <- asin(pmin(pmax(sin_solar_altitude, -1), 1))
  
  #Solarer Azimut in Radiant
  cos_solar_altitude <- cos(solar_altitude_rad)
  
  if (cos_solar_altitude < 1e-6) {
    solar_azimuth_rad <- NA_real_ } else {
      sin_solar_azimuth <- -cos(solar_declination_rad) * sin(hour_angle_rad) / cos_solar_altitude
      cos_solar_azimuth <- (sin(solar_declination_rad) - sin(latitude_rad) * sin(solar_altitude_rad)) / (cos(latitude_rad) * cos_solar_altitude)
      
      solar_azimuth_raw <- atan2(sin_solar_azimuth, cos_solar_azimuth)
      solar_azimuth_rad <- ifelse(solar_azimuth_raw < 0, solar_azimuth_raw + 2 * pi, solar_azimuth_raw)}
  
  list(
    solar_altitude_rad = solar_altitude_rad,
    solar_azimuth_rad  = solar_azimuth_rad)}

#D.4.3 Luftmasse nach Kasten und Young (1989) berechnen
compute_air_mass_kasten_young <- function(solar_altitude_rad) {
  
  solar_altitude_deg <- solar_altitude_rad * 180 / pi
  solar_zenith_deg <- 90 - solar_altitude_deg
  
  air_mass <- 1 / (
    cos(deg2rad(solar_zenith_deg)) +
      0.50572 * (96.07995 - solar_zenith_deg)^(-1.6364))
  air_mass}

#================================================================================================================================================================================================
###SCHRITT D.5a: TÄGLICHE POTENZIELLE DIREKTSTRAHLUNG (MJ/m²/Tag) BERECHNEN###
#================================================================================================================================================================================================
##ANMERKUNG: DURCH DEN OPERATOR "write_tif" KANN EIN EXPORT ALS GEOTIFF STATTFINDEN. FALLS NICHT ERWÜNSCHT = FALSE SETZEN!##

#D.5a.1 Zeitzonen-Offset relativ zu UTC für ein gegebenes Datum bestimmen
get_tz_offset_hours <- function(date, tz = VALLEY_Ipot_tz) {
  reference_time <- as.POSIXct(paste(date, "12:00:00"), tz = tz)
  as.integer(as.numeric(format(reference_time, "%z")) / 100)}

#D.5a.2 Tägliche potenzielle Direktstrahlung für einen einzelnen Tag berechnen
VALLEY_compute_daily_Ipot <- function(
    target_date,
    dem_raster,
    cos_slope, sin_slope,
    cos_aspect, sin_aspect,
    latitude_deg, longitude_deg,
    hours_seq = 0:23,
    solar_constant = 1368,
    psi_alpha = 0.75,
    output_dir,
    timezone_string = "America/Argentina/San_Juan",
    write_tif = FALSE,
    gdal_options = c("COMPRESS=LZW", "TILED=YES")) {
  
  #Tag des Jahres und Exzentrizitätskorrektur der Erde-Sonne-Distanz
  day_of_year <- lubridate::yday(target_date)
  earth_sun_distance_factor <- 1 + 0.033 * cos(2 * pi * day_of_year / 365)
  
  #Zeitzonen-Offset für den jeweiligen Tag bestimmen
  tz_offset_hours <- get_tz_offset_hours(target_date, timezone_string)
  
  #Tagesraster der potenziellen Direktstrahlung initialisieren
  daily_Ipot_raster <- rast(dem_raster)
  values(daily_Ipot_raster) <- 0
  
  #Stündliche Integration über alle Stunden des Tages
  for (hour_value in hours_seq) {
    
    solar_position <- compute_solar_position(
      date = target_date,
      hour = hour_value,
      latitude_deg = latitude_deg,
      longitude_deg = longitude_deg,
      tz_offset_hours = tz_offset_hours)
    
    solar_altitude_rad <- solar_position$solar_altitude_rad
    solar_azimuth_rad  <- solar_position$solar_azimuth_rad
    
    #Nachtstunden oder ungültige Sonnenpositionen überspringen
    if (!is.finite(solar_altitude_rad) || solar_altitude_rad <= 0 || !is.finite(solar_azimuth_rad)) next
    
    sin_solar_altitude <- sin(solar_altitude_rad)
    cos_solar_altitude <- cos(solar_altitude_rad)
    
    #cos(Azimut - Aspect) effizient berechnen
    cos_solar_azimuth <- cos(solar_azimuth_rad)
    sin_solar_azimuth <- sin(solar_azimuth_rad)
    cos_azimuth_minus_aspect <- cos_solar_azimuth * cos_aspect + sin_solar_azimuth * sin_aspect
    
    #Einfallswinkel auf geneigter Fläche
    cos_incidence_angle <- sin_solar_altitude * cos_slope + cos_solar_altitude * sin_slope * cos_azimuth_minus_aspect
    cos_incidence_angle <- clamp(cos_incidence_angle, lower = 0, upper = 1, values = TRUE)
    
    #Luftmasse und atmosphärische Transmissivität
    air_mass <- compute_air_mass_kasten_young(solar_altitude_rad)
    atmospheric_transmissivity <- psi_alpha^air_mass
    
    #Potenzielle Direktstrahlung in W/m² und Umrechnung in MJ/m² pro Stunde
    direct_radiation_Wm2 <- (solar_constant * earth_sun_distance_factor * atmospheric_transmissivity) * cos_incidence_angle
    direct_radiation_MJm2_hour <- direct_radiation_Wm2 * (3600 / 1e6)
    
    #Stündlichen Beitrag zum Tagesraster addieren
    daily_Ipot_raster <- daily_Ipot_raster + direct_radiation_MJm2_hour}
  
  #Tägliches GeoTIFF schreiben
  daily_output_file <- file.path(
    output_dir,
    paste0("I_pot_MJm2_day_", format(target_date, "%Y-%m-%d"), ".tif"))
  
  if (isTRUE(write_tif)) {
    writeRaster(
      daily_Ipot_raster,
      daily_output_file,
      overwrite = TRUE,
      gdal = gdal_options)}
  
  #AOI-Mittelwert der täglichen potenziellen Direktstrahlung berechnen
  daily_Ipot_mean <- terra::global(daily_Ipot_raster, "mean", na.rm = TRUE)[1, 1]
  
  list(
    date = target_date,
    tz_offset_hours = tz_offset_hours,
    output_file = daily_output_file,
    I_pot_raster = daily_Ipot_raster,
    I_pot_MJm2_day_mean = as.numeric(daily_Ipot_mean))}


#================================================================================================================================================================================================
###SCHRITT D.5b: LOOP ÜBER ALLE TAGE - AOI-BASIERTE I_pot-ERGEBNISSE BERECHNEN UND ALS DAILY TABLE SPEICHERN###
#================================================================================================================================================================================================
##ANMERKUNG: DURCH DEN OPERATOR "write_tif" KANN EIN EXPORT ALS GEOTIFF STATTFINDEN. FALLS NICHT ERWÜNSCHT = FALSE SETZEN##

#D.5b.1 Ergebnisliste vorbereiten
VALLEY_Ipot_results <- vector("list", length(VALLEY_Ipot_dates_seq))

#D.5b.2 Progress-Bar für die tägliche Berechnung
VALLEY_Ipot_progress <- txtProgressBar(
  min = 0,
  max = length(VALLEY_Ipot_dates_seq),
  style = 3)

#D:5b.3 Alle Tage des Untersuchungszeitraums durchlaufen
for (date_index in seq_along(VALLEY_Ipot_dates_seq)) {
  setTxtProgressBar(VALLEY_Ipot_progress, date_index)
  
  current_date <- VALLEY_Ipot_dates_seq[date_index]
  
  VALLEY_Ipot_results[[date_index]] <- VALLEY_compute_daily_Ipot(
    target_date = current_date,
    dem_raster = VALLEY_Ipot_dem,
    cos_slope = VALLEY_Ipot_cos_slope,
    sin_slope = VALLEY_Ipot_sin_slope,
    cos_aspect = VALLEY_Ipot_cos_aspect,
    sin_aspect = VALLEY_Ipot_sin_aspect,
    latitude_deg = VALLEY_Ipot_lat_deg,
    longitude_deg = VALLEY_Ipot_lon_deg,
    hours_seq = VALLEY_Ipot_hours_seq,
    solar_constant = VALLEY_Ipot_solar_constant,
    psi_alpha = VALLEY_Ipot_psi_alpha_default,
    output_dir = VALLEY_Ipot_output_dir,
    timezone_string = VALLEY_Ipot_tz,
    write_tif = FALSE)}

close(VALLEY_Ipot_progress)

#D.5b.4 AOI-basierte Tagesreihe der potenziellen Direktstrahlung erstellen
VALLEY_Ipot_daily_summary <- data.frame(
  date = as.Date(VALLEY_Ipot_dates_seq),
  I_pot_MJm2_day_mean = vapply(
    VALLEY_Ipot_results,
    function(result_entry) result_entry$I_pot_MJm2_day_mean,
    numeric(1))) %>%
  mutate(date = as.Date(date)) %>%
  arrange(date)

#5b.5 Kurze Kontrolle der Tagesreihe
print(head(VALLEY_Ipot_daily_summary))
summary(VALLEY_Ipot_daily_summary$I_pot_MJm2_day_mean)

#================================================================================================================================================================================================
###SCHRITT D.6: BASIS-QC DER AOI-BASIERTEN I_pot-ZEITREIHE###
#================================================================================================================================================================================================

#D.6.1 Struktur, Duplikate und Missing Values prüfen
VALLEY_Ipot_basic_qc <- VALLEY_Ipot_daily_summary %>%
  summarise(
    n_rows = n(),
    start_date = min(date, na.rm = TRUE),
    end_date = max(date, na.rm = TRUE),
    n_missing_Ipot = sum(!is.finite(I_pot_MJm2_day_mean)))

VALLEY_Ipot_duplicate_dates <- VALLEY_Ipot_daily_summary %>%
  count(date) %>%
  filter(n > 1)

VALLEY_Ipot_missing_days <- setdiff(VALLEY_Ipot_dates_seq, VALLEY_Ipot_daily_summary$date)

print(VALLEY_Ipot_basic_qc)
print(VALLEY_Ipot_duplicate_dates)
print(VALLEY_Ipot_missing_days)

#D.6.2 Wertebereich der AOI-basierten I_pot-Zeitreihe prüfen
VALLEY_Ipot_value_range <- VALLEY_Ipot_daily_summary %>%
  summarise(
    min_Ipot_MJm2_day = min(I_pot_MJm2_day_mean, na.rm = TRUE),
    max_Ipot_MJm2_day = max(I_pot_MJm2_day_mean, na.rm = TRUE))

print(VALLEY_Ipot_value_range)

#================================================================================================================================================================================================
###SCHRITT D.7a: SAISONALE ZUSAMMENFASSUNG DER AOI-BASIERTEN I_pot-ZEITREIHE###
#================================================================================================================================================================================================

#D.7a.1 Zeitinformationen ergänzen und einfache Jahreszeitenklassen zuweisen
VALLEY_Ipot_seasonal_df <- VALLEY_Ipot_daily_summary %>%
  mutate(
    month = lubridate::month(date),
    season_simple = case_when(
      month %in% c(12, 1, 2) ~ "summer",
      month %in% c(3, 4, 5)  ~ "autumn",
      month %in% c(6, 7, 8)  ~ "winter",
      TRUE                   ~ "spring"))

#D.7a.2 Übersicht über Mittelwert, Streuung und Spannweite der gesamten Zeitreihe
VALLEY_Ipot_summary <- VALLEY_Ipot_seasonal_df %>%
  summarise(
    n_days        = n(),
    I_pot_mean_MJ = mean(I_pot_MJm2_day_mean, na.rm = TRUE),
    I_pot_sd_MJ   = sd(I_pot_MJm2_day_mean, na.rm = TRUE),
    I_pot_min_MJ  = min(I_pot_MJm2_day_mean, na.rm = TRUE),
    I_pot_max_MJ  = max(I_pot_MJm2_day_mean, na.rm = TRUE))

print(VALLEY_Ipot_summary)

#D.7a.3 Saisonale Mittelwerte der AOI-basierten I_pot-Zeitreihe
VALLEY_Ipot_seasonal_summary <- VALLEY_Ipot_seasonal_df %>%
  group_by(season_simple) %>%
  summarise(
    I_pot_mean_MJ = mean(I_pot_MJm2_day_mean, na.rm = TRUE),
    .groups = "drop") %>%
  arrange(factor(season_simple, levels = c("summer", "autumn", "winter", "spring")))

print(VALLEY_Ipot_seasonal_summary)

#D.7a.4 Jahreszeitlichen Kontrast zwischen Winter und Sommer bestimmen
VALLEY_Ipot_winter_summer_contrast <- VALLEY_Ipot_seasonal_summary %>%
  pivot_wider(
    names_from = season_simple,
    values_from = I_pot_mean_MJ) %>%
  mutate(
    winter_to_summer_ratio = winter / summer)

print(VALLEY_Ipot_winter_summer_contrast)


#================================================================================================================================================================================================
###SCHRITT D.7b: PLOT DES SAISONALEN VERLAUFS DER POTENZIELLEN DIREKTSTRAHLUNG###
#================================================================================================================================================================================================

#D.7b.1 Plot-Dataframe vorbereiten
VALLEY_Ipot_plot_df <- VALLEY_Ipot_daily_summary %>%
  transmute(
    date = as.Date(date),
    I_pot_MJm2_day_mean = as.numeric(I_pot_MJm2_day_mean))

#D.7b.2 Gemeinsames y-Limit bestimmen
VALLEY_Ipot_ymax <- ceiling(max(VALLEY_Ipot_plot_df$I_pot_MJm2_day_mean, na.rm = TRUE) / 5) * 5

#D.7b.3 Plot des saisonalen Verlaufs
ggplot(VALLEY_Ipot_plot_df, aes(x = date, y = I_pot_MJm2_day_mean)) +
  geom_line(linewidth = 0.5, color = "black", na.rm = TRUE) +
  scale_y_continuous(
    limits = c(0, VALLEY_Ipot_ymax),
    breaks = seq(0, VALLEY_Ipot_ymax, by = 5)) +
  labs(
    x = "Datum",
    y = expression("Potenzielle Direktstrahlung" ~ (MJ~m^-2~d^-1))) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12))

#================================================================================================================================================================================================
###SCHRITT D.8: SENSITIVITÄTSCHECK FÜR psi alpha###
#================================================================================================================================================================================================

#D.8.1 Zu testende psi_alpha-Werte und Sensitivitätszeitraum festlegen (Festgelegt nach HOCK (1999))
VALLEY_Ipot_psi_alpha_grid <- c(0.75, 0.65, 0.85)

VALLEY_Ipot_sensitivity_dates <- seq(
  min(VALLEY_Ipot_daily_summary$date, na.rm = TRUE),
  max(VALLEY_Ipot_daily_summary$date, na.rm = TRUE),
  by = "day")

#D.8.2 Hilfsfunktion: AOI-basierte I_pot-Zeitreihe für einen psi_alpha-Wert berechnen
compute_Ipot_sensitivity_series <- function(dates_seq, psi_alpha_value) {
  
  sensitivity_results <- vector("list", length(dates_seq))
  
  for (date_index in seq_along(dates_seq)) {
    current_date <- as.Date(dates_seq[date_index])
    
    daily_Ipot_result <- compute_daily_Ipot(
      target_date = current_date,
      dem_raster = VALLEY_Ipot_dem,
      cos_slope = VALLEY_Ipot_cos_slope,
      sin_slope = VALLEY_Ipot_sin_slope,
      cos_aspect = VALLEY_Ipot_cos_aspect,
      sin_aspect = VALLEY_Ipot_sin_aspect,
      latitude_deg = VALLEY_Ipot_lat_deg,
      longitude_deg = VALLEY_Ipot_lon_deg,
      hours_seq = VALLEY_Ipot_hours_seq,
      solar_constant = VALLEY_Ipot_solar_constant,
      psi_alpha = psi_alpha_value,
      output_dir = VALLEY_Ipot_output_dir,
      timezone_string = VALLEY_Ipot_tz,
      write_tif = FALSE)
    
    sensitivity_results[[date_index]] <- data.frame(
      date = current_date,
      psi_alpha = psi_alpha_value,
      I_pot_MJm2_day_mean = daily_Ipot_result$I_pot_MJm2_day_mean)}
  
  bind_rows(sensitivity_results) %>%
    mutate(
      season_simple = case_when(
        lubridate::month(date) %in% c(12, 1, 2) ~ "summer",
        lubridate::month(date) %in% c(3, 4, 5)  ~ "autumn",
        lubridate::month(date) %in% c(6, 7, 8)  ~ "winter",
        TRUE                                    ~ "spring")) %>%
    arrange(psi_alpha, date)}

#D.8.3 psi_alpha-Grid durchrechnen
VALLEY_Ipot_psi_alpha_sensitivity_results <- vector("list", length(VALLEY_Ipot_psi_alpha_grid))

for (psi_alpha_index in seq_along(VALLEY_Ipot_psi_alpha_grid)) {
  current_psi_alpha <- VALLEY_Ipot_psi_alpha_grid[psi_alpha_index]
  
  VALLEY_Ipot_psi_alpha_sensitivity_results[[psi_alpha_index]] <- compute_Ipot_sensitivity_series(
    dates_seq = VALLEY_Ipot_sensitivity_dates,
    psi_alpha_value = current_psi_alpha)}

VALLEY_Ipot_psi_alpha_sensitivity <- bind_rows(VALLEY_Ipot_psi_alpha_sensitivity_results)

#D.8.4 Mittelwerte der potenziellen Direktstrahlung pro psi_alpha und Jahreszeit
VALLEY_Ipot_psi_alpha_scaling_summary <- VALLEY_Ipot_psi_alpha_sensitivity %>%
  group_by(psi_alpha, season_simple) %>%
  summarise(
    I_pot_mean_MJ = mean(I_pot_MJm2_day_mean, na.rm = TRUE),
    .groups = "drop") %>%
  arrange(psi_alpha, factor(season_simple, levels = c("summer", "autumn", "winter", "spring")))

print(VALLEY_Ipot_psi_alpha_scaling_summary)

#D.8.5 Relative Änderung gegenüber psi_alpha = 0.75
VALLEY_Ipot_psi_alpha_reference <- VALLEY_Ipot_psi_alpha_sensitivity %>%
  filter(psi_alpha == 0.75) %>%
  select(date, I_pot_reference = I_pot_MJm2_day_mean)

VALLEY_Ipot_psi_alpha_relative_change <- VALLEY_Ipot_psi_alpha_sensitivity %>%
  left_join(VALLEY_Ipot_psi_alpha_reference, by = "date") %>%
  mutate(
    ratio_to_reference = I_pot_MJm2_day_mean / I_pot_reference) %>%
  group_by(psi_alpha, season_simple) %>%
  summarise(
    ratio_to_reference = mean(ratio_to_reference, na.rm = TRUE),
    rel_change_pct = (mean(ratio_to_reference, na.rm = TRUE) - 1) * 100,
    .groups = "drop") %>%
  arrange(psi_alpha, factor(season_simple, levels = c("summer", "autumn", "winter", "spring")))

print(VALLEY_Ipot_psi_alpha_relative_change)

#D.8.6 Kompakte Gesamtauswertung über alle Jahreszeiten
VALLEY_Ipot_psi_alpha_relative_change_compact <- VALLEY_Ipot_psi_alpha_sensitivity %>%
  left_join(VALLEY_Ipot_psi_alpha_reference, by = "date") %>%
  mutate(
    ratio_to_reference = I_pot_MJm2_day_mean / I_pot_reference) %>%
  group_by(psi_alpha) %>%
  summarise(
    rel_change_mean = mean(ratio_to_reference, na.rm = TRUE),
    rel_change_pct  = (rel_change_mean - 1) * 100,
    .groups = "drop") %>%
  arrange(psi_alpha)

print(VALLEY_Ipot_psi_alpha_relative_change_compact)

#================================================================================================================================================================================================
###ABSCHNITT E: MONTE-CARLO-BASIERTE KALIBRIERUNG DER EMPIRISCHEN PARAMETER UND MODELLIERUNG DER FSCA###
#================================================================================================================================================================================================
#================================================================================================================================================================================================
###SCHRITT E.1: VORBEREITUNG (FIXWERTE, PARAMETERGRENZEN, LHS SETTINGS)###
#================================================================================================================================================================================================

#E.1.1 Reproduzierbarkeit sicherstellen (wird zusätzlich zum LHS-spezifischen Seed gesetzt, um die Reproduzierbarkeit aller Schritte zu gewährleisten)
set.seed(1)

#E.1.2 Kalibrierfenster festlegen
VALLEY_calibration_start <- as.Date("2023-02-18")
VALLEY_calibration_end   <- as.Date("2024-02-25")

#E.1.3 Fixwerte für die Niederschlagspartitionierung
VALLEY_calibration_T_snow_fixed <- -1
VALLEY_calibration_T_rain_fixed <-  2

#E.1.4 Startwert für den vereinfachten fSCA-Beobachtungsoperator
VALLEY_calibration_eps_fsca_start <- 35

#E.1.5 LHS- und Monte-Carlo-Einstellungen
VALLEY_calibration_N_lhs <- 5000
VALLEY_calibration_seed  <- 1 

#E.1.6 Parametergrenzen des Screening-Ranges festlegen
VALLEY_calibration_parameter_mins <- c(kP = 0.6, Tcrit = -1.5, Cdd = 2.0, kR = 0.00)

VALLEY_calibration_parameter_maxs <- c(kP = 1.6,Tcrit = 2.5,Cdd = 8.0,kR = 2.50)

#================================================================================================================================================================================================
###SCHRITT E.2: FINALEN KALIBRIER-DATAFRAME AUS FSCA, KLIMA UND STRAHLUNG AUFBAUEN###
#================================================================================================================================================================================================

#E.2.1 fSCA-Beobachtungsreihe für das Kalibrierfenster vorbereiten
VALLEY_calibration_fSCA_observations <- VALLEY_fSCA_df_qc %>%
  transmute(
    date = as.Date(date),
    fSCA_obs = as.numeric(fSCA_obs),
    qc_flag_manual,
    qc_keep) %>%
  filter(date >= VALLEY_calibration_start, date <= VALLEY_calibration_end)

#E.2.2 Klimadaten für das Kalibrierfenster vorbereiten
VALLEY_calibration_climate_input <- data.frame(
  date = as.Date(ERA5_all_df$date),
  T_mean_C = as.numeric(ERA5_all_df$t2m_c),
  P_total_mm = as.numeric(ERA5_all_df$tp_mm)) %>%
  filter(date >= VALLEY_calibration_start, date <= VALLEY_calibration_end)

#E.2.3 Strahlungsdaten für das Kalibrierfenster vorbereiten
VALLEY_calibration_radiation_input <- VALLEY_Ipot_daily_summary %>%
  transmute(
    date = as.Date(date),
    I_pot_MJm2_day = as.numeric(I_pot_MJm2_day_mean)) %>%
  filter(date >= VALLEY_calibration_start, date <= VALLEY_calibration_end)

#E.2.4 Finalen Kalibrier-Dataframe aus Klima, Strahlung und fSCA zusammenführen
VALLEY_calibration_df <- VALLEY_calibration_climate_input %>%
  left_join(VALLEY_calibration_radiation_input, by = "date") %>%
  left_join(VALLEY_calibration_fSCA_observations, by = "date") %>%
  arrange(date)

#E.2.5 Zusatzspalten für die spätere Modellierung ergänzen
VALLEY_calibration_df <- VALLEY_calibration_df %>%
  mutate(
    year = lubridate::year(date),
    month = lubridate::month(date),
    doy = lubridate::yday(date))

#E.2.6 Quick check
print(head(VALLEY_calibration_df, 10))
summary(VALLEY_calibration_df)

#================================================================================================================================================================================================
###SCHRITT E.3: MODELLFUNKTION FÜR DIE AOI-BASIERTEN FSCA-SIMULATION DEFINIEREN###
#================================================================================================================================================================================================

#E.3.1 lineare Niederschlagspartitionierung nach HARPOLD et al. (2017)
VALLEY_calibration_partition_precip <- function(T_mean_C, P_total_mm, T_snow, T_rain) {
  snow_fraction <- ifelse(
    is.na(T_mean_C) | is.na(P_total_mm), NA_real_,
    ifelse(
      T_mean_C <= T_snow, 1,
      ifelse(
        T_mean_C >= T_rain, 0,
        (T_rain - T_mean_C) / (T_rain - T_snow))))
  
  snow_mm <- P_total_mm * snow_fraction
  rain_mm <- P_total_mm - snow_mm
  
  list(
    snow_fraction = snow_fraction,
    snow_mm = snow_mm,
    rain_mm = rain_mm)}

##AMMERKUNG: IN ANLEHNUNG AN HOCK (1999), ANDERSON (2006), HARPOLD et al. (2017), BEVEN & BINLEY (1992) und ESSERY & POMEROY (2004)##

#E.3.2 Hauptfunktion: Modellierung einer simulierte fSCA Zeitreihe auf Basis der Kalibrier-Inputs und der empirischen Parameter
VALLEY_calibration_run_fsca_model <- function(
    df,
    params,
    T_snow = VALLEY_calibration_T_snow_fixed,
    T_rain = VALLEY_calibration_T_rain_fixed,
    eps_fsca = VALLEY_calibration_eps_fsca_start) {
  
  kP    <- as.numeric(params["kP"])
  Tcrit <- as.numeric(params["Tcrit"])
  Cdd   <- as.numeric(params["Cdd"])
  kR    <- as.numeric(params["kR"])
  
  n_days <- nrow(df)
  
  simulation_df <- df %>%
    dplyr::mutate(
      snow_fraction = NA_real_,
      snow_mm = NA_real_,
      rain_mm = NA_real_,
      snow_input_mm = NA_real_,
      melt_potential_mm = NA_real_,
      melt_actual_mm = NA_real_,
      SWE_mm = NA_real_,
      fSCA_mod = NA_real_)
  
  SWE_previous <- 0
  
  for (day_index in seq_len(n_days)) {
    
    current_temperature <- simulation_df$T_mean_C[day_index]
    current_precipitation <- simulation_df$P_total_mm[day_index]
    current_radiation <- simulation_df$I_pot_MJm2_day[day_index]
    
    #Niederschlagspartitionierung
    partition_result <- VALLEY_calibration_partition_precip(
      T_mean_C = current_temperature,
      P_total_mm = current_precipitation,
      T_snow = T_snow,
      T_rain = T_rain)
    
    current_snow_mm <- partition_result$snow_mm
    current_rain_mm <- partition_result$rain_mm
    current_snow_fraction <- partition_result$snow_fraction
    
    #Schneeinput mit skalierendem Niederschlagsfaktor
    current_snow_input <- ifelse(is.na(current_snow_mm), NA_real_, kP * current_snow_mm)
    
    #Potenzielle Schmelze: Temperaturindex + Strahlung
    current_melt_potential <- ifelse(
      is.na(current_temperature) | is.na(current_radiation),
      NA_real_,
      pmax(0, Cdd * (current_temperature - Tcrit) + kR * current_radiation))
    
    #Speicher vor Schmelze
    SWE_before_melt <- ifelse(
      is.na(current_snow_input),
      SWE_previous,
      SWE_previous + current_snow_input)
    
    #Tatsächliche Schmelze kann Speicher nicht überschreiten
    current_melt_actual <- ifelse(
      is.na(current_melt_potential) | is.na(SWE_before_melt),
      NA_real_,
      min(SWE_before_melt, current_melt_potential))
    
    #Aktualisierter Schneespeicher
    current_SWE <- ifelse(
      is.na(SWE_before_melt) | is.na(current_melt_actual),
      NA_real_,
      pmax(0, SWE_before_melt - current_melt_actual))
    
    #Vereinfachter Beobachtungsoperator (nach ESSERY & POMEROY 2004)
    current_fSCA_mod <- ifelse(
      is.na(current_SWE),
      NA_real_,
      current_SWE / (current_SWE + eps_fsca))
    
    # Ergebnisse speichern
    simulation_df$snow_fraction[day_index]   <- current_snow_fraction
    simulation_df$snow_mm[day_index]         <- current_snow_mm
    simulation_df$rain_mm[day_index]         <- current_rain_mm
    simulation_df$snow_input_mm[day_index]   <- current_snow_input
    simulation_df$melt_potential_mm[day_index] <- current_melt_potential
    simulation_df$melt_actual_mm[day_index]  <- current_melt_actual
    simulation_df$SWE_mm[day_index]          <- current_SWE
    simulation_df$fSCA_mod[day_index]        <- current_fSCA_mod
    
    # State update
    SWE_previous <- ifelse(is.na(current_SWE), SWE_previous, current_SWE)}
  
  simulation_df}

#================================================================================================================================================================================================
###SCHRITT E.4: SCORE-FUNKTIONEN FÜR DIE KALIBRIERUNG DEFINIEREN###
#================================================================================================================================================================================================
#ANMERKUNG: nach KNOBEN et al. (2019), SETIYOWATI et al. (2025) und HAUGHTON et al. (2016)#

#E.4.1 RMSE (Root Mean Square Error)
VALLEY_calibration_rmse <- function(obs, sim) {
  valid_values <- is.finite(obs) & is.finite(sim)
  if (sum(valid_values) == 0) return(NA_real_)
  sqrt(mean((sim[valid_values] - obs[valid_values])^2))}

#4.2 MBE (Mean Bias Error)
VALLEY_calibration_mbe <- function(obs, sim) {
  valid_values <- is.finite(obs) & is.finite(sim)
  if (sum(valid_values) == 0) return(NA_real_)
  mean(sim[valid_values] - obs[valid_values])}

#4.3 Pearson-Korrelationskoeffizient r
VALLEY_calibration_pearson_r <- function(obs, sim) {
  valid_values <- is.finite(obs) & is.finite(sim)
  if (sum(valid_values) < 2) return(NA_real_)
  stats::cor(obs[valid_values], sim[valid_values], method = "pearson")}

#4.4 KGE (Kling-Gupta Efficiency)
VALLEY_calibration_kge <- function(obs, sim) {
  valid_values <- is.finite(obs) & is.finite(sim)
  
  if (sum(valid_values) < 2) return(NA_real_)
  
  obs_valid <- obs[valid_values]
  sim_valid <- sim[valid_values]
  
  r_value <- stats::cor(obs_valid, sim_valid, method = "pearson")
  
  sd_obs <- stats::sd(obs_valid)
  sd_sim <- stats::sd(sim_valid)
  mean_obs <- mean(obs_valid)
  mean_sim <- mean(sim_valid)
  
  if (!is.finite(sd_obs) || !is.finite(sd_sim) || !is.finite(mean_obs) || !is.finite(mean_sim)) {
    return(NA_real_)}
  
  if (sd_obs == 0 || mean_obs == 0) {
    return(NA_real_)}
  
  alpha <- sd_sim / sd_obs
  beta  <- mean_sim / mean_obs
  
  1 - sqrt((r_value - 1)^2 + (alpha - 1)^2 + (beta - 1)^2)}

#4.5 Zentrale Score-Funktion für einen Modelllauf
VALLEY_calibration_score_run <- function(model_run_df) {
  
  obs <- model_run_df$fSCA_obs
  sim <- model_run_df$fSCA_mod
  
  data.frame(
    n_obs = sum(is.finite(obs) & is.finite(sim)),
    RMSE  = VALLEY_calibration_rmse(obs, sim),
    MBE   = VALLEY_calibration_mbe(obs, sim),
    r     = VALLEY_calibration_pearson_r(obs, sim),
    KGE   = VALLEY_calibration_kge(obs, sim))}

#================================================================================================================================================================================================
###SCHRITT E.5: KALIBRIERUNG (LHS/MONTE-CARLO): PARAMETER-SAMPLE AUFBAUEN###
#================================================================================================================================================================================================

#E.5.1 LHS-ähnliches Parameter-Sampling definieren (in Anlehnung an SHIELDS & ZHANG 2016)
VALLEY_calibration_generate_lhs_sample <- function(N, parameter_mins, parameter_maxs) {
  parameter_names <- names(parameter_mins)
  
  lhs_unit_sample <- sapply(parameter_names, function(parameter_name) {
    unit_values <- (seq_len(N) - runif(N)) / N
    sample(unit_values, size = N, replace = FALSE)})
  
  lhs_sample_df <- as.data.frame(lhs_unit_sample)
  
  for (parameter_name in parameter_names) {
    lhs_sample_df[[parameter_name]] <- parameter_mins[parameter_name] +
      lhs_sample_df[[parameter_name]] * (parameter_maxs[parameter_name] - parameter_mins[parameter_name])}
  
  lhs_sample_df}

#E.5.2 LHS-Sample ziehen
set.seed(VALLEY_calibration_seed)

VALLEY_calibration_parameter_sample <- VALLEY_calibration_generate_lhs_sample(
  N = VALLEY_calibration_N_lhs,
  parameter_mins = VALLEY_calibration_parameter_mins,
  parameter_maxs = VALLEY_calibration_parameter_maxs)

#E.5.3 Run-ID ergänzen und Reihenfolge festlegen
VALLEY_calibration_parameter_sample <- VALLEY_calibration_parameter_sample %>%
  mutate(run_id = row_number()) %>%
  select(run_id, kP, Tcrit, Cdd, kR)

#================================================================================================================================================================================================
###SCHRITT E.6: KALIBRIERUNGSLOOP ÜBER ALLE LHS-LÄUFE###
#================================================================================================================================================================================================

#E.6.1 Hauptscore definieren
VALLEY_calibration_main_score <- "KGE"

#E.6.2 Leere Ergebnisliste anlegen
VALLEY_calibration_results <- vector("list", length = nrow(VALLEY_calibration_parameter_sample))

#E.6.3 Progress-Bar für die Kalibrierung
VALLEY_calibration_progress <- txtProgressBar(
  min = 0,
  max = nrow(VALLEY_calibration_parameter_sample),
  style = 3)

#E.6.4 Loop über alle LHS-Läufe
for (run_index in seq_len(nrow(VALLEY_calibration_parameter_sample))) {
  setTxtProgressBar(VALLEY_calibration_progress, run_index)
  
  current_parameter_row <- VALLEY_calibration_parameter_sample[run_index, ]
  
  current_parameter_vector <- c(
    kP    = current_parameter_row$kP,
    Tcrit = current_parameter_row$Tcrit,
    Cdd   = current_parameter_row$Cdd,
    kR    = current_parameter_row$kR)
  
  current_model_run <- VALLEY_calibration_run_fsca_model(
    df = VALLEY_calibration_df,
    params = current_parameter_vector,
    T_snow = VALLEY_calibration_T_snow_fixed,
    T_rain = VALLEY_calibration_T_rain_fixed,
    eps_fsca = VALLEY_calibration_eps_fsca_start)
  
  current_scores <- VALLEY_calibration_score_run(current_model_run)
  
  VALLEY_calibration_results[[run_index]] <- data.frame(
    run_id = current_parameter_row$run_id,
    kP     = current_parameter_row$kP,
    Tcrit  = current_parameter_row$Tcrit,
    Cdd    = current_parameter_row$Cdd,
    kR     = current_parameter_row$kR,
    n_obs  = current_scores$n_obs,
    RMSE   = current_scores$RMSE,
    MBE    = current_scores$MBE,
    r      = current_scores$r,
    KGE    = current_scores$KGE)}

close(VALLEY_calibration_progress)

#E.6.5 Ergebnistabelle zusammenführen und nach Hauptscore sortieren
VALLEY_calibration_results_df <- dplyr::bind_rows(VALLEY_calibration_results) %>%
  mutate(main_score = .data[[VALLEY_calibration_main_score]])

if (VALLEY_calibration_main_score %in% c("KGE", "r")) {
  VALLEY_calibration_results_df <- VALLEY_calibration_results_df %>%
    arrange(desc(main_score))} else {
  VALLEY_calibration_results_df <- VALLEY_calibration_results_df %>%
    arrange(main_score)}

#E.6.5 Quick check
print(head(VALLEY_calibration_results_df, 10))
summary(VALLEY_calibration_results_df)

#================================================================================================================================================================================================
###SCHRITT E.7: BESTEN LAUF EXTRAHIEREN UND ERNEUT DURCHRECHNEN###
#================================================================================================================================================================================================

#E.7.1 Besten Parametersatz auswählen
VALLEY_calibration_best_run_summary <- VALLEY_calibration_results_df %>%
  slice(1)

print(VALLEY_calibration_best_run_summary)

#E.7.2 Parametervektor für den besten Lauf aufbauen
VALLEY_calibration_best_params <- c(
  kP    = VALLEY_calibration_best_run_summary$kP,
  Tcrit = VALLEY_calibration_best_run_summary$Tcrit,
  Cdd   = VALLEY_calibration_best_run_summary$Cdd,
  kR    = VALLEY_calibration_best_run_summary$kR)

print(VALLEY_calibration_best_params)

#E.7.3 Best-Run-Zeitreihe erneut simulieren
VALLEY_calibration_best_run_df <- VALLEY_calibration_run_fsca_model(
  df = VALLEY_calibration_df,
  params = VALLEY_calibration_best_params,
  T_snow = VALLEY_calibration_T_snow_fixed,
  T_rain = VALLEY_calibration_T_rain_fixed,
  eps_fsca = VALLEY_calibration_eps_fsca_start)

#E.7.4 Begleitende Scores des besten Laufs prüfen
VALLEY_calibration_best_run_scores <- VALLEY_calibration_score_run(
  VALLEY_calibration_best_run_df)

print(VALLEY_calibration_best_run_scores)

#================================================================================================================================================================================================
###SCHRITT 8: BEST-RUN-ZEITREIHE ALS DIREKTER VERGLEICH VON BEOBACHTUNG UND SIMULATION###
#================================================================================================================================================================================================

#E.8.1 Panel-Label mit Hauptscore des Best-Runs erzeugen
VALLEY_calibration_panel_label <- paste0(
  "Agua Negra Becken | ",
  VALLEY_calibration_main_score,
  ": ",
  round(VALLEY_calibration_best_run_summary$main_score, 2))

#E.8.2 Gemeinsamen Plot-Dataframe erzeugen
VALLEY_calibration_best_run_plot_df <- bind_rows(
  VALLEY_calibration_best_run_df %>%
    filter(is.finite(fSCA_obs)) %>%
    transmute(
      date = as.Date(date),
      fSCA = as.numeric(fSCA_obs),
      series = "Observed"),
  
  VALLEY_calibration_best_run_df %>%
    transmute(
      date = as.Date(date),
      fSCA = as.numeric(fSCA_mod),
      series = "Simulated")) %>%
  mutate(
    panel = VALLEY_calibration_panel_label,
    series = factor(series, levels = c("Observed", "Simulated")))

#E.8.3 Gemeinsames y-Limit bestimmen
VALLEY_calibration_fSCA_ymax <- max(
  VALLEY_calibration_best_run_plot_df$fSCA,
  na.rm = TRUE)

#E.8.4 Plot des Best-Runs: Beobachtung vs. Simulation
ggplot(VALLEY_calibration_best_run_plot_df, aes(x = date, y = fSCA, color = series)) +
  geom_line(linewidth = 0.9, na.rm = TRUE) +
  facet_wrap(~ panel, ncol = 1, scales = "fixed") +
  scale_color_manual(
    values = c(
      "Observed" = "blue",
      "Simulated"  = "red"),
    name = NULL) +
  scale_y_continuous(
    limits = c(0, VALLEY_calibration_fSCA_ymax)) +
  labs(
    x = "Datum",
    y = "fraktionale Schneebedeckung (-)") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 11, lineheight = 0.9),
    strip.background = element_rect(fill = "grey85", colour = "black"),
    panel.grid.minor = element_blank(),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12),
    legend.text  = element_text(size = 12))

#================================================================================================================================================================================================
###SCHRITT E.9: SENSITIVITÄTSPRÜFUNG FÜR eps_fsca###
#================================================================================================================================================================================================

#E.9.1 Grid für eps_fsca festlegen
VALLEY_calibration_eps_fsca_grid <- c(10, 25, 45, 60)

#E.9.2 Leere Listen für Ergebnisse anlegen
VALLEY_calibration_eps_sensitivity_results <- vector("list", length(VALLEY_calibration_eps_fsca_grid))
VALLEY_calibration_eps_sensitivity_runs    <- vector("list", length(VALLEY_calibration_eps_fsca_grid))

#E.9.3 Loop über eps_fsca: nur Best-Run neu simulieren
for (eps_index in seq_along(VALLEY_calibration_eps_fsca_grid)) {
  
  current_eps_fsca <- VALLEY_calibration_eps_fsca_grid[eps_index]
  
  current_eps_run <- VALLEY_calibration_run_fsca_model(
    df = VALLEY_calibration_df,
    params = VALLEY_calibration_best_params,
    T_snow = VALLEY_calibration_T_snow_fixed,
    T_rain = VALLEY_calibration_T_rain_fixed,
    eps_fsca = current_eps_fsca)
  
  current_eps_scores <- VALLEY_calibration_score_run(current_eps_run)
  
  VALLEY_calibration_eps_sensitivity_runs[[eps_index]] <- current_eps_run
  
  VALLEY_calibration_eps_sensitivity_results[[eps_index]] <- data.frame(
    eps_fsca = current_eps_fsca,
    kP       = unname(VALLEY_calibration_best_params["kP"]),
    Tcrit    = unname(VALLEY_calibration_best_params["Tcrit"]),
    Cdd      = unname(VALLEY_calibration_best_params["Cdd"]),
    kR       = unname(VALLEY_calibration_best_params["kR"]),
    n_obs    = current_eps_scores$n_obs,
    RMSE     = current_eps_scores$RMSE,
    MBE      = current_eps_scores$MBE,
    r        = current_eps_scores$r,
    KGE      = current_eps_scores$KGE)}

#E.9.4 Ergebnisse zusammenführen und sortieren
VALLEY_calibration_eps_sensitivity_df <- dplyr::bind_rows(VALLEY_calibration_eps_sensitivity_results) %>%
  arrange(eps_fsca)

print(VALLEY_calibration_eps_sensitivity_df)

#E.9.5 Plot-Dataframe erzeugen: Best-Run-Referenz + Best-Run-Simulationen für alle eps_fsca
VALLEY_calibration_eps_fsca_default <- VALLEY_calibration_eps_fsca_start
VALLEY_calibration_best_series_label <- paste0("Best-Run | eps_fSCA = ", VALLEY_calibration_eps_fsca_default)

VALLEY_calibration_eps_sensitivity_plot_df <- dplyr::bind_rows(
  VALLEY_calibration_best_run_df %>%
    dplyr::transmute(
      date   = as.Date(date),
      fSCA   = as.numeric(fSCA_mod),
      series = VALLEY_calibration_best_series_label),
  dplyr::bind_rows(lapply(seq_along(VALLEY_calibration_eps_sensitivity_runs), function(eps_index) {
    VALLEY_calibration_eps_sensitivity_runs[[eps_index]] %>%
      dplyr::transmute(
        date   = as.Date(date),
        fSCA   = as.numeric(fSCA_mod),
        series = paste0("eps_fSCA = ", VALLEY_calibration_eps_fsca_grid[eps_index]))}))) %>%
  dplyr::mutate(
    series = factor(
      series,
      levels = c(
        VALLEY_calibration_best_series_label,
        paste0("eps_fSCA = ", VALLEY_calibration_eps_fsca_grid))))

#E.9.6 Plot: Best-Run unter variierendem eps_fsca
ggplot(VALLEY_calibration_eps_sensitivity_plot_df, aes(x = date, y = fSCA, color = series)) +
  geom_line(linewidth = 0.7, na.rm = TRUE) +
  scale_color_manual(
    values = c(
      setNames("black", VALLEY_calibration_best_series_label),
      "eps_fSCA = 10" = "blue",
      "eps_fSCA = 25" = "darkgreen",
      "eps_fSCA = 45" = "orange",
      "eps_fSCA = 60" = "red"),
    labels = c(
      bquote("Best-Run | " * epsilon[fSCA] == .(VALLEY_calibration_eps_fsca_default)),
      expression(epsilon[fSCA] == 10),
      expression(epsilon[fSCA] == 25),
      expression(epsilon[fSCA] == 45),
      expression(epsilon[fSCA] == 60)),
    name = NULL) +
  labs(
    x = "Datum",
    y = "fraktionale Schneebedeckung (-)") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12),
    legend.text  = element_text(size = 14))

#================================================================================================================================================================================================
###SCHRITT E.10: ROBUSTHEIT DER KALIBRIERUNG - TOP-ENSEMBLE (ÄQUIFINALITÄT)###
#================================================================================================================================================================================================

#E.10.1 Top-Ensemble definieren: beste 30 Läufe
VALLEY_calibration_top_n <- min(30, nrow(VALLEY_calibration_results_df))
VALLEY_calibration_top_runs <- VALLEY_calibration_results_df %>%
  slice(seq_len(VALLEY_calibration_top_n))

#E.10.2 Hilfsfunktion für 5-, 50- und 95-%-Quantile
VALLEY_calibration_quantile_summary <- function(x) {
  as.numeric(stats::quantile(x, probs = c(0.05, 0.50, 0.95), na.rm = TRUE))}

#E.10.3 Parameter-Quantile des Top-Ensembles
VALLEY_calibration_top_parameter_quantiles <- data.frame(
  top_n = VALLEY_calibration_top_n,
  
  kP_q05 = VALLEY_calibration_quantile_summary(VALLEY_calibration_top_runs$kP)[1],
  kP_q50 = VALLEY_calibration_quantile_summary(VALLEY_calibration_top_runs$kP)[2],
  kP_q95 = VALLEY_calibration_quantile_summary(VALLEY_calibration_top_runs$kP)[3],
  
  Tcrit_q05 = VALLEY_calibration_quantile_summary(VALLEY_calibration_top_runs$Tcrit)[1],
  Tcrit_q50 = VALLEY_calibration_quantile_summary(VALLEY_calibration_top_runs$Tcrit)[2],
  Tcrit_q95 = VALLEY_calibration_quantile_summary(VALLEY_calibration_top_runs$Tcrit)[3],
  
  Cdd_q05 = VALLEY_calibration_quantile_summary(VALLEY_calibration_top_runs$Cdd)[1],
  Cdd_q50 = VALLEY_calibration_quantile_summary(VALLEY_calibration_top_runs$Cdd)[2],
  Cdd_q95 = VALLEY_calibration_quantile_summary(VALLEY_calibration_top_runs$Cdd)[3],
  
  kR_q05 = VALLEY_calibration_quantile_summary(VALLEY_calibration_top_runs$kR)[1],
  kR_q50 = VALLEY_calibration_quantile_summary(VALLEY_calibration_top_runs$kR)[2],
  kR_q95 = VALLEY_calibration_quantile_summary(VALLEY_calibration_top_runs$kR)[3])

print(VALLEY_calibration_top_parameter_quantiles)

#================================================================================================================================================================================================
###SCHRITT E.11: TOP-VARIANTEN IM KONTEXT DES BEST-RUNS###
#================================================================================================================================================================================================

#E.11.1 Spannweite der Scores und Parameter innerhalb des Top-Ensembles
VALLEY_calibration_top_run_ranges <- data.frame(
  top_n = nrow(VALLEY_calibration_top_runs),
  
  main_score_min = min(VALLEY_calibration_top_runs$main_score, na.rm = TRUE),
  main_score_max = max(VALLEY_calibration_top_runs$main_score, na.rm = TRUE),
  
  RMSE_min = min(VALLEY_calibration_top_runs$RMSE, na.rm = TRUE),
  RMSE_max = max(VALLEY_calibration_top_runs$RMSE, na.rm = TRUE),
  
  MBE_min = min(VALLEY_calibration_top_runs$MBE, na.rm = TRUE),
  MBE_max = max(VALLEY_calibration_top_runs$MBE, na.rm = TRUE),
  
  r_min = min(VALLEY_calibration_top_runs$r, na.rm = TRUE),
  r_max = max(VALLEY_calibration_top_runs$r, na.rm = TRUE),
  
  KGE_min = min(VALLEY_calibration_top_runs$KGE, na.rm = TRUE),
  KGE_max = max(VALLEY_calibration_top_runs$KGE, na.rm = TRUE),
  
  kP_min = min(VALLEY_calibration_top_runs$kP, na.rm = TRUE),
  kP_max = max(VALLEY_calibration_top_runs$kP, na.rm = TRUE),
  
  Tcrit_min = min(VALLEY_calibration_top_runs$Tcrit, na.rm = TRUE),
  Tcrit_max = max(VALLEY_calibration_top_runs$Tcrit, na.rm = TRUE),
  
  Cdd_min = min(VALLEY_calibration_top_runs$Cdd, na.rm = TRUE),
  Cdd_max = max(VALLEY_calibration_top_runs$Cdd, na.rm = TRUE),
  
  kR_min = min(VALLEY_calibration_top_runs$kR, na.rm = TRUE),
  kR_max = max(VALLEY_calibration_top_runs$kR, na.rm = TRUE))

print(VALLEY_calibration_top_run_ranges)

#E.11.2 Median der Scores innerhalb des Top-Ensembles
VALLEY_calibration_top_score_median <- data.frame(
  run_type = "top_median",
  main_score = median(VALLEY_calibration_top_runs$main_score, na.rm = TRUE),
  RMSE = median(VALLEY_calibration_top_runs$RMSE, na.rm = TRUE),
  MBE = median(VALLEY_calibration_top_runs$MBE, na.rm = TRUE),
  r = median(VALLEY_calibration_top_runs$r, na.rm = TRUE),
  KGE = median(VALLEY_calibration_top_runs$KGE, na.rm = TRUE))

#E.11.3 Best-Run und Top-Median direkt gegenüberstellen
VALLEY_calibration_best_vs_top_median <- data.frame(
  run_type = c("best_run", "top_median"),
  main_score = c(
    round(VALLEY_calibration_best_run_summary$main_score, 4),
    round(VALLEY_calibration_top_score_median$main_score, 4)),
  RMSE = c(
    round(VALLEY_calibration_best_run_summary$RMSE, 4),
    round(VALLEY_calibration_top_score_median$RMSE, 4)),
  MBE = c(
    round(VALLEY_calibration_best_run_summary$MBE, 4),
    round(VALLEY_calibration_top_score_median$MBE, 4)),
  r = c(
    round(VALLEY_calibration_best_run_summary$r, 4),
    round(VALLEY_calibration_top_score_median$r, 4)),
  KGE = c(
    round(VALLEY_calibration_best_run_summary$KGE, 4),
    round(VALLEY_calibration_top_score_median$KGE, 4)))

print(VALLEY_calibration_best_vs_top_median)

#================================================================================================================================================================================================
###SCHRITT E.12: TEST DES STRAHLUNGSTERMS (kR = 0) - EINFLUSS DER POTENZIELLEN DIREKTSTRAHLUNG###
#================================================================================================================================================================================================

#E.12.1 Parametersatz des Best-Runs ohne Strahlungsterm ableiten
VALLEY_calibration_best_params_no_radiation <- VALLEY_calibration_best_params
VALLEY_calibration_best_params_no_radiation["kR"] <- 0

print(VALLEY_calibration_best_params_no_radiation)

#E.12.2 Best-Run ohne Strahlungsterm erneut simulieren
VALLEY_calibration_run_without_radiation <- VALLEY_calibration_run_fsca_model(
  df = VALLEY_calibration_df,
  params = VALLEY_calibration_best_params_no_radiation,
  T_snow = VALLEY_calibration_T_snow_fixed,
  T_rain = VALLEY_calibration_T_rain_fixed,
  eps_fsca = VALLEY_calibration_eps_fsca_start)

#E.12.3 Gütemaße des Best-Runs und der Variante ohne Strahlung vergleichen
VALLEY_calibration_scores_best_run <- VALLEY_calibration_score_run(
  VALLEY_calibration_best_run_df) %>%
  dplyr::mutate(run_type = "best_run")

VALLEY_calibration_scores_without_radiation <- VALLEY_calibration_score_run(
  VALLEY_calibration_run_without_radiation) %>%
  dplyr::mutate(run_type = "kR_0")

VALLEY_calibration_radiation_comparison <- dplyr::bind_rows(
  VALLEY_calibration_scores_best_run,
  VALLEY_calibration_scores_without_radiation) %>%
  dplyr::select(run_type, n_obs, RMSE, MBE, r, KGE)

print(VALLEY_calibration_radiation_comparison)

#E.12.4 Differenzen zwischen Best-Run und kR = 0 berechnen
VALLEY_calibration_radiation_differences <- VALLEY_calibration_radiation_comparison %>%
  tidyr::pivot_wider(
    names_from = run_type,
    values_from = c(RMSE, MBE, r, KGE)) %>%
  dplyr::mutate(
    d_RMSE = RMSE_kR_0 - RMSE_best_run,
    d_MBE  = MBE_kR_0  - MBE_best_run,
    d_r    = r_kR_0    - r_best_run,
    d_KGE  = KGE_kR_0  - KGE_best_run)

print(VALLEY_calibration_radiation_differences)

#E.12.5 Plot-Dataframe für den direkten Vergleich beider Simulationen erzeugen
VALLEY_calibration_radiation_panel_label <- "Agua Negra Becken | Einfluss des Strahlungsterms (kR = 0)"

VALLEY_calibration_radiation_plot_df <- dplyr::bind_rows(
  VALLEY_calibration_best_run_df %>%
    dplyr::transmute(
      date = as.Date(date),
      fSCA = as.numeric(fSCA_mod),
      series = "Simulierter Best-Run"),
  VALLEY_calibration_run_without_radiation %>%
    dplyr::transmute(
      date = as.Date(date),
      fSCA = as.numeric(fSCA_mod),
      series = "Simulierter Best-Run kR = 0")) %>%
  dplyr::mutate(
    panel = VALLEY_calibration_radiation_panel_label,
    series = factor(
      series,
      levels = c("Simulierter Best-Run", "Simulierter Best-Run kR = 0")))

VALLEY_calibration_radiation_ymax <- max(
  VALLEY_calibration_radiation_plot_df$fSCA,
  na.rm = TRUE)

#E.12.6 Plot: Best-Run vs. Best-Run ohne Strahlungsterm
ggplot(VALLEY_calibration_radiation_plot_df, aes(x = date, y = fSCA, color = series)) +
  geom_line(linewidth = 0.9, na.rm = TRUE) +
  facet_wrap(~ panel, ncol = 1, scales = "fixed") +
  scale_color_manual(
    values = c(
      "Simulierter Best-Run" = "red",
      "Simulierter Best-Run kR = 0" = "black"),
    name = NULL) +
  scale_y_continuous(
    limits = c(0, VALLEY_calibration_radiation_ymax)) +
    labs(
    x = "Datum",
    y = "fraktionale Schneebedeckung (-)") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 11, lineheight = 0.9),
    strip.background = element_rect(fill = "grey85", colour = "black"),
    panel.grid.minor = element_blank(),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12),
    legend.text  = element_text(size = 12))


#================================================================================================================================================================================================
###ABSCHNITT F: EXTERNE PLAUSIBILISIERUNG DER SAISONALEN SCHNEEDYNAMIK MIT VIIRS VNP10D1F###
#================================================================================================================================================================================================
#================================================================================================================================================================================================
###SCHRITT F.1: INPUT-PFADE UND ZIELPOLYGON FESTLEGEN###
#================================================================================================================================================================================================
##ANMERKUNG: FÜR DIESEN ABSCHNITT IST DER DOWNLOAD VON VON VIIRS DATEN (VNP10D1F, Version 2) ZUR SCHNEEBEDECKUNG NOTWENDIG!##

#VERFÜRGBAR NACH: RIGGS, G. & HALL, D. K. (2025): VIIRS/NPP Snow Cover Daily L3 Global 375m CMG.(VNP10D1F, Version 2). [Online-Publikation]. [Zugriff am: 11.04.2026]. Verfügbar unter: DOI: 10.5067/F3VW9TYIK26A. 

#WICHTIG: ERWARTET WERDEN TÄGLICHE DATEN (01.05.2023 bis 30.09.2023) ALS HS-DATEIEN ZUSAMMENGEFASST IN EINEM ORDNER!#

#F.1.1 Input- und Output-Pfade festlegen
VALLEY_VIIRS_input_dir <- "PFAD/ZU/VIIRS/HS/DATEIEN/ORDNER" #MUSS ANGEPASST WERDEN!
VALLEY_VIIRS_output_csv <- "OUTPUT/PFAD/VALLEY_VIIRS_daily_snow_cover.csv" #MUSS ANGEPASST WERDEN!

#F.1.2 Liste aller VIIRS-HDF5-Dateien erstellen
VALLEY_VIIRS_file_list <- list.files(
  VALLEY_VIIRS_input_dir,
  pattern = "\\.h5$",
  full.names = TRUE)

#F.1.3 Zielpolygon für die externe Plausibilisierung festlegen
VALLEY_VIIRS_target_vect <- project(VALLEY_basin_poly, "EPSG:4326")

#================================================================================================================================================================================================
###SCHRITT F.2: HILFSFUNKTIONEN ZUR VERARBEITUNG DER VIIRS-DATEN###
#================================================================================================================================================================================================

#F.2.1 Datum aus dem VIIRS-Dateinamen extrahieren
VIIRS_extract_date_from_filename <- function(file_path) {
  date_string <- stringr::str_extract(basename(file_path), "A\\d{7}")
  as.Date(substr(date_string, 2, 8), format = "%Y%j")}

#F.2.2 Schneebedeckungs-Layer aus der HDF5-Datei auswählen
VIIRS_get_snow_cover_layer <- function(h5_file, pattern = "NDSI_Snow_Cover$") {
  subdatasets <- terra::sds(h5_file)
  matching_layer <- grep(pattern, names(subdatasets))
  
  if (length(matching_layer) == 0) return(NULL)
  
  subdatasets[[matching_layer[1]]]}

#F.2.3 Einzelne VIIRS-Datei auf Tageswert für das Zielpolygon verarbeiten
VIIRS_process_daily_file <- function(h5_file, target_vect) {
  
  VIIRS_snow_raster <- VIIRS_get_snow_cover_layer(h5_file, "NDSI_Snow_Cover$")
  
  if (is.null(VIIRS_snow_raster)) {
    return(data.frame(
      date = VIIRS_extract_date_from_filename(h5_file),
      file = basename(h5_file),
      viirs_snow = NA_real_,
      n_valid = 0L))}
  
  target_vect_same_crs <- project(target_vect, crs(VIIRS_snow_raster))
  
  VIIRS_snow_raster <- crop(VIIRS_snow_raster, target_vect_same_crs)
  VIIRS_snow_raster <- mask(VIIRS_snow_raster, target_vect_same_crs)
  
  #Nur gültige Werte zwischen 0 und 100 behalten; Qualitätsflags >100 werden als NA gesetzt
  VIIRS_snow_raster <- ifel(VIIRS_snow_raster >= 0 & VIIRS_snow_raster <= 100, VIIRS_snow_raster, NA)
  
  VIIRS_values <- values(VIIRS_snow_raster, mat = FALSE)
  VIIRS_values <- as.numeric(VIIRS_values)
  VIIRS_values <- VIIRS_values[is.finite(VIIRS_values)]
  
  data.frame(
    date = VIIRS_extract_date_from_filename(h5_file),
    file = basename(h5_file),
    viirs_snow = ifelse(length(VIIRS_values) > 0, mean(VIIRS_values / 100, na.rm = TRUE), NA_real_),
    n_valid = length(VIIRS_values))}


#================================================================================================================================================================================================
###SCHRITT F.3: VERARBEITUNG DER VIIRS-DATEN###
#================================================================================================================================================================================================

#F.3.1 Ergebnisliste vorbereiten
VALLEY_VIIRS_results <- vector("list", length(VALLEY_VIIRS_file_list))

#F.3.2 Progress-Bar für die Verarbeitung aller VIIRS-Dateien
VALLEY_VIIRS_progress <- txtProgressBar(
  min = 0,
  max = length(VALLEY_VIIRS_file_list),
  style = 3)

#F.3.3 Alle VIIRS-Dateien verarbeiten
for (file_index in seq_along(VALLEY_VIIRS_file_list)) {
  setTxtProgressBar(VALLEY_VIIRS_progress, file_index)
  
  VALLEY_VIIRS_results[[file_index]] <- VIIRS_process_daily_file(
    h5_file = VALLEY_VIIRS_file_list[file_index],
    target_vect = VALLEY_VIIRS_target_vect)}

close(VALLEY_VIIRS_progress)

#F.3.4 Tagesreihe der VIIRS-Schneebedeckung zusammenführen und speichern
VALLEY_VIIRS_daily_df <- bind_rows(VALLEY_VIIRS_results) %>%
  arrange(date)

write.csv(VALLEY_VIIRS_daily_df, VALLEY_VIIRS_output_csv, row.names = FALSE) #OPTIONAL!

#F.3.5 Kurze Kontrolle der Ergebnisse
print(head(VALLEY_VIIRS_daily_df, 10))
summary(VALLEY_VIIRS_daily_df$viirs_snow)

#================================================================================================================================================================================================
###SCHRITT F.4: ZUSAMMENFÜHREN MIT DER SIMULIERTEN FSCA-ZEITREIHE###
#================================================================================================================================================================================================

#F.4.1 Vergleichszeitraum der externen Plausibilisierung festlegen
VALLEY_VIIRS_compare_start <- as.Date("2023-05-01")
VALLEY_VIIRS_compare_end   <- as.Date("2023-09-30")

#F.4.2 Simulierte fSCA-Zeitreihe des Best-Runs mit der VIIRS-Tagesreihe zusammenführen
VALLEY_VIIRS_comparison_df <- VALLEY_calibration_best_run_df %>%
  dplyr::select(date, fSCA_mod) %>%
  dplyr::left_join(
    VALLEY_VIIRS_daily_df %>%
      dplyr::select(date, viirs_snow, n_valid),
    by = "date") %>%
  dplyr::filter(date >= VALLEY_VIIRS_compare_start, date <= VALLEY_VIIRS_compare_end) %>%
  dplyr::arrange(date)

#================================================================================================================================================================================================
###SCHRITT F.5: EINFACHE VISUELLE PLAUSIBILISIERUNG###
#================================================================================================================================================================================================

#F.5.1 Plot-Dataframe für den Vergleich von VIIRS und modellierter fSCA erzeugen
VALLEY_VIIRS_plot_df <- VALLEY_VIIRS_comparison_df %>%
  tidyr::pivot_longer(
    cols = c(viirs_snow, fSCA_mod),
    names_to = "series",
    values_to = "value") %>%
  dplyr::mutate(
    series = factor(
      series,
      levels = c("viirs_snow", "fSCA_mod"),
      labels = c("VIIRS Schneebedeckung", "Simulierter Best-Run")),
    panel = "Vergleich: modellierte fSCA vs. VIIRS Schneebedeckung | 2023-05-01 bis 2023-09-30")

#F.5.2 Plot des direkten Vergleichs zwischen VIIRS-Schneebedeckung und modellierter fSCA
ggplot(VALLEY_VIIRS_plot_df, aes(x = date, y = value, color = series)) +
  geom_line(linewidth = 1.0, na.rm = TRUE) +
  facet_wrap(~ panel, ncol = 1, scales = "fixed") +
  scale_color_manual(
    values = c(
      "VIIRS Schneebedeckung" = "blue",
      "Simulierter Best-Run"  = "red"),
    name = NULL) +
  labs(
    x = "Datum",
    y = "fraktionale Schneebedeckung (-)") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12))


#================================================================================================================================================================================================
###SCHRITT F.6: GÜTEMASSE DER EXTERNEN PLAUSIBILISIERUNG###
#================================================================================================================================================================================================

#F.6.1 Nur Tage mit gültigen VIIRS- und Modellwerten für die Gütemaße verwenden
VALLEY_VIIRS_evaluation_df <- VALLEY_VIIRS_comparison_df %>%
  dplyr::filter(!is.na(fSCA_mod), !is.na(viirs_snow))

#F.6.2 Gütemaße zwischen modellierter fSCA und VIIRS-Schneebedeckung berechnen
VALLEY_VIIRS_score_summary <- data.frame(
  n_obs = nrow(VALLEY_VIIRS_evaluation_df),
  RMSE  = VALLEY_calibration_rmse(
    obs = VALLEY_VIIRS_evaluation_df$viirs_snow,
    sim = VALLEY_VIIRS_evaluation_df$fSCA_mod),
  MBE   = VALLEY_calibration_mbe(
    obs = VALLEY_VIIRS_evaluation_df$viirs_snow,
    sim = VALLEY_VIIRS_evaluation_df$fSCA_mod),
  r     = VALLEY_calibration_pearson_r(
    obs = VALLEY_VIIRS_evaluation_df$viirs_snow,
    sim = VALLEY_VIIRS_evaluation_df$fSCA_mod),
  KGE   = VALLEY_calibration_kge(
    obs = VALLEY_VIIRS_evaluation_df$viirs_snow,
    sim = VALLEY_VIIRS_evaluation_df$fSCA_mod))

print(VALLEY_VIIRS_score_summary)
