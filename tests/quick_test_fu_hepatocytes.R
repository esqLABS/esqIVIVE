devtools::load_all()

calculate_fu_in_vitro(
  partition_qspr = "All PK-Sim Standard",
  log_lipophilicity = 0.46,
  ionization = c("neutral", 0, 0),
  type_system = "hepatocytes",
  fetal_bovine_serum_fraction = 0,
  microplate_type = 96,
  volume_medium = 0.15,
  pka = c(0, 0, 0),
  henry_law_constant = 1E-6,
  concentration_cells = 1
)
