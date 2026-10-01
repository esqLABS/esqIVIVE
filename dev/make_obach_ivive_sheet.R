# Build the import sheet used by vignettes/articles/clearance-ivive-check.qmd
#
# Source: Obach 1999, Drug Metab Dispos 27(11):1350-1359
# (Prediction of human clearance of twenty-nine drugs from hepatic microsomal
# intrinsic clearance data: an examination of in vitro half-life approach and
# nonspecific binding to microsomes). Data already curated in
# inst/extdata/Obach_1999_Clmicrosomes.csv.
#
# Output: inst/extdata/Obach1999_IVIVE_input.xlsx with sheets
#   - Compounds : physchem, in vivo clearance, in vitro half-life and fu_mic
#   - Scenarios : IVIVE / PBK options, one row per scenario (edit to add more)
#   - Columns   : description and unit of every column
#
# Run from the package root: source("dev/make_obach_ivive_sheet.R")

library(openxlsx)

obach <- read.csv(
  "inst/extdata/Obach_1999_Clmicrosomes.csv",
  check.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)

halogen <- function(x) ifelse(is.na(x), 0, x)

compounds <- data.frame(
  Drug = obach$Drug,
  Ionization = c(Basic = "base", Neutral = "neutral", Acidic = "acid")[obach$Ionization],
  pKa = obach$pKa,
  MW_gmol = obach$MW,
  Cl = halogen(obach$Cl),
  Br = halogen(obach$Br),
  I = halogen(obach$I),
  F = 0,
  LogP = obach$LogP37C,
  LogD = obach$LogD37C,
  fu_plasma = obach[["Fraction Unbound in Plasma (fu)"]],
  BP = obach$BP,
  CLp_obs_mLminkg = obach[["Plasma Clearance _mlminkg"]],
  CLb_obs_mLminkg = obach[["Blood clearance mLminkg"]],
  Cmic_mgml = obach$Microsomal_mgml,
  halflife_min = obach[["half-life_Invitro_min"]],
  halflife_sd_min = obach[["half-life_sd"]],
  CLint_invitro_mLminmg = obach$Cl_invitro_mlminmg,
  fu_mic_obs = obach$fu_mic,
  fu_mic_obs_sd = obach$fu_mic_SD,
  References_in_vivo = obach$References,
  row.names = NULL
)

scenarios <- data.frame(
  scenario = paste0("S", 1:6),
  description = c(
    "Standard IVIVE (fu_mic = 1), PK-Sim partitioning, QSAR permeability",
    "fu_mic All_literature, PK-Sim partitioning, QSAR permeability",
    "fu_mic All_literature, high cell permeability, PK-Sim partitioning",
    "fu_mic All_literature, high cell permeability, Wood 2017 scaling factors, PK-Sim partitioning",
    "fu_mic Rodgers & Rowland, high cell permeability, Rodgers & Rowland partitioning",
    "Reference: measured fu_mic (Obach 1999), high cell permeability, PK-Sim partitioning"
  ),
  fu_mic_method = c(
    "none", "All_literature", "All_literature", "All_literature",
    "Rodgers & Rowland + fu", "measured"
  ),
  fu_mic_fallback = c(NA, NA, NA, NA, "All_literature", NA),
  permeability = c("QSAR", "QSAR", "high", "high", "high", "high"),
  permeability_high_cmmin = c(NA, NA, 1000, 1000, 1000, 1000),
  empirical_scalar = c("No", "No", "No", "Yes", "No", "No"),
  blood_plasma_ratio = "observed",
  pkml = c(
    rep("single-iv-pksim.pkml", 4),
    "single-iv-rodgers-rowland.pkml",
    "single-iv-pksim.pkml"
  )
)

columns <- data.frame(
  sheet = c(rep("Compounds", ncol(compounds)), rep("Scenarios", ncol(scenarios))),
  column = c(names(compounds), names(scenarios)),
  description = c(
    "Drug name",
    "Ionization class used by esqIVIVE (acid, base, neutral)",
    "Most relevant pKa (0 for neutrals)",
    "Molecular weight",
    "Number of chlorine atoms (PK-Sim effective MW)",
    "Number of bromine atoms",
    "Number of iodine atoms",
    "Number of fluorine atoms (not reported by Obach, set to 0)",
    "log octanol/water partition coefficient at 37 C",
    "log distribution coefficient at pH 7.4, 37 C",
    "Fraction unbound in plasma",
    "Blood to plasma concentration ratio",
    "Observed in vivo plasma clearance",
    "Observed in vivo blood clearance",
    "Microsomal protein concentration in the incubation",
    "In vitro half-life of substrate depletion",
    "Standard deviation of the in vitro half-life",
    "In vitro intrinsic clearance = 0.693 / t1/2 / Cmic",
    "Measured fraction unbound in microsomes",
    "Standard deviation of measured fu_mic",
    "Original references for the in vivo clearance (as cited by Obach 1999)",
    "Scenario identifier",
    "Scenario description",
    "fu_mic used in IVIVE: none (=1), measured, or any calculate_fu_in_vitro() partition_qspr",
    "Method used when fu_mic_method fails for a compound (e.g. R&R for non strong bases)",
    "Cell permeability: QSAR (PK-Sim formula from MW_eff and LogP) or high",
    "Cell permeability value used when permeability = high",
    "Apply Wood 2017 empirical scaling factors in IVIVE_clearance() (Yes / No)",
    "Blood/plasma ratio: observed (set blood cell partitioning from BP) or model (PK-Sim formula from LogP)",
    "pkml file in inst/extdata/pkml4htpbk (defines the partition coefficient method)"
  ),
  unit = c(
    "", "", "", "g/mol", "", "", "", "", "Log Units", "Log Units", "", "",
    "mL/min/kg", "mL/min/kg", "mg/mL", "min", "min", "mL/min/mg protein", "", "", "",
    "", "", "", "", "", "cm/min", "", "", ""
  )
)

wb <- createWorkbook()
header <- createStyle(textDecoration = "bold", fgFill = "#DDEBF7", border = "bottom")
for (sheet in list(
  list("Compounds", compounds),
  list("Scenarios", scenarios),
  list("Columns", columns)
)) {
  addWorksheet(wb, sheet[[1]])
  writeData(wb, sheet[[1]], sheet[[2]], headerStyle = header)
  freezePane(wb, sheet[[1]], firstRow = TRUE)
  setColWidths(wb, sheet[[1]], cols = seq_len(ncol(sheet[[2]])), widths = "auto")
}
saveWorkbook(wb, "inst/extdata/Obach1999_IVIVE_input.xlsx", overwrite = TRUE)
