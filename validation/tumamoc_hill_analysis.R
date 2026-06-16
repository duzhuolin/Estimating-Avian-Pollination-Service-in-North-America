# ============================================================
# Tumamoc Hill — Hummingbird-Pollinated Plant Density Trends
# Purpose: Compare long-term plant demographic trends (1906-2012)
#          against Aridlands biome AFV results
# Data:   Rodriguez-Buritica et al. (2013) Ecological Archives E094-083
# ============================================================

# ---- 1. Download data (run once) ----
data_url <- "https://www.esapubs.org/archive/ecol/e094/083/CsvFiles.zip"
zip_file <- "validation/CsvFiles.zip"
if (!file.exists(zip_file)) {
  download.file(data_url, zip_file, mode = "wb")
  unzip(zip_file, exdir = "validation/tumamoc_data")
}

# ---- 2. Load and prepare ----
library(data.table)

# Load key tables
density  <- fread("validation/tumamoc_data/SMDensity.csv")
species  <- fread("validation/tumamoc_data/Species.csv")
plots    <- fread("validation/tumamoc_data/Plots.csv")

# ---- 3. Identify hummingbird/bird-pollinated species ----
# Based on known pollination syndromes of Sonoran Desert flora
bird_pollinated <- data.table(
  Code = c("FOSP", "CAGI", "FEWI", "ECFA", "CAER", "HYEM",
           "SACO", "PEPA", "SILO", "DESC", "LYEX", "LYBE", "LYFR", "HICO"),
  CommonName = c("Ocotillo", "Saguaro", "Barrel cactus", "Hedgehog cactus",
                 "Fairyduster", "Desert lavender", "Chia", "Parry's penstemon",
                 "Longflower tubetongue", "Larkspur", "Wolfberry", "Wolfberry",
                 "Wolfberry", "Desert rosemallow"),
  Pollinator = c("Hummingbird", "Bat/Bird", "Hummingbird", "Hummingbird",
                 "Hummingbird", "Hummingbird", "Hummingbird", "Hummingbird",
                 "Hummingbird", "Hummingbird", "Hummingbird/Bee",
                 "Hummingbird/Bee", "Hummingbird/Bee", "Hummingbird")
)

# Merge species info
bird_spp <- merge(bird_pollinated, species[, .(Code, AcceptedName, Habit, Duration)],
                  by = "Code", all.x = TRUE)
cat("\nBird-pollinated species in dataset:\n")
print(bird_spp[, .(Code, CommonName, AcceptedName, Habit)])

# ---- 4. Extract density data for bird-pollinated species ----
bird_density <- density[Code %in% bird_pollinated$Code]
bird_density <- merge(bird_density, bird_pollinated[, .(Code, CommonName, Pollinator)],
                      by = "Code")

# ---- 5. Key analysis: density trends per species (1970-2012) ----
# Focus on years overlapping with your study period
bird_density_study <- bird_density[Year >= 1970]

# Compute mean density per species per year (across all plots)
annual_density <- bird_density_study[, .(
  mean_density = mean(Density, na.rm = TRUE),
  sd_density   = sd(Density, na.rm = TRUE),
  n_plots      = .N,
  total_density = sum(Density, na.rm = TRUE)
), by = .(Code, CommonName, Year)]

# ---- 6. Directional trend per species (1970-2012) ----
# Simple linear trend
species_trends <- annual_density[, {
  if (.N >= 3) {
    fit <- lm(mean_density ~ Year)
    s <- summary(fit)
    list(
      slope         = coef(fit)[2],
      p_value       = s$coefficients[2, 4],
      r_squared     = s$r.squared,
      first_density = mean_density[which.min(Year)],
      last_density  = mean_density[which.max(Year)],
      n_years       = .N
    )
  }
}, by = .(Code, CommonName)]

species_trends[, direction := fifelse(slope > 0, "Increase", "Decrease")]
cat("\n--- Density trends of hummingbird-pollinated plants (1970-2012) ---\n")
print(species_trends[order(slope)])

# ---- 7. Aggregate: all bird-pollinated species combined ----
total_annual <- bird_density_study[, .(
  total_density = sum(Density, na.rm = TRUE),
  n_species     = uniqueN(Code)
), by = Year]

cat("\n--- Aggregate trend (all bird-pollinated species) ---\n")
if (nrow(total_annual) >= 3) {
  fit_total <- lm(total_density ~ Year, data = total_annual)
  cat(sprintf("Slope: %.3f, p = %.4f, R² = %.3f\n",
              coef(fit_total)[2],
              summary(fit_total)$coefficients[2, 4],
              summary(fit_total)$r.squared))
}

# ---- 8. Simple plot ----
pdf("validation/tumamoc_bird_pollinated_density_trends.pdf", width = 12, height = 8)

par(mfrow = c(2, 1))

# Panel A: Per-species trends
plot(NA, xlim = c(1970, 2012), ylim = range(annual_density$mean_density, na.rm = TRUE),
     xlab = "Year", ylab = "Mean density (individuals/plot)",
     main = "Hummingbird-pollinated perennials — Tumamoc Hill")

cols <- rainbow(length(unique(annual_density$Code)))
for (i in seq_along(unique(annual_density$Code))) {
  sp <- unique(annual_density$Code)[i]
  sp_data <- annual_density[Code == sp]
  lines(sp_data$Year, sp_data$mean_density, col = cols[i], lwd = 2, type = "b", pch = 16)
}
legend("topright", legend = unique(annual_density$CommonName),
       col = cols, lwd = 1, cex = 0.7, ncol = 2)

# Panel B: Total
plot(total_annual$Year, total_annual$total_density,
     type = "b", pch = 16, lwd = 2,
     xlab = "Year", ylab = "Total density (all plots combined)",
     main = "Aggregate density — all hummingbird-pollinated species")
if (nrow(total_annual) >= 3) {
  abline(fit_total, col = "red", lty = 2)
}

dev.off()

cat("\nPlot saved to validation/tumamoc_bird_pollinated_density_trends.pdf\n")
cat("\n=== DONE ===\n")
