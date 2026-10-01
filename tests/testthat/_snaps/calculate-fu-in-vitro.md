# calculate_fu_in_vitro: rodgers_rowland_fu, strong bases only

    Code
      calculate_fu_in_vitro(method = "rodgers_rowland_fu", lipophilicity = 2,
        ionization = c("neutral", 0), system = "cells", fbs_fraction = 0,
        microplate_type = 96, volume_medium = 0.22, pka = c(0, 0),
        henry_law_constant = 1e-06, fu_plasma = 0.3, blood_plasma_ratio = 1,
        concentration_cells = 2)
    Condition
      Error in `calculate_fu_in_vitro()`:
      ! "rodgers_rowland_fu" is only available for strong bases: first ionization class "base" with a pKa above 7.
      i PK-Sim uses protein binding instead for acids, neutral compounds and weak bases, which is not available here yet.

# calculate_fu_in_vitro rejects invalid method and system

    Code
      calculate_fu_in_vitro(method = "bogus", lipophilicity = 3, ionization = c(
        "acid", 0), system = "cells", pka = c(6, 0), concentration_cells = 2)
    Condition
      Error in `calculate_fu_in_vitro()`:
      ! `method` must be one of "austin", "hallifax", "turner", "kilford", "poulin", "all_literature", "poulin_theil", "poulin_theil_fu", "berezhkovskiy", "berezhkovskiy_fu", "pksim_standard", "pksim_standard_fu", "rodgers_rowland_fu", "schmitt", or "schmitt_fu", not "bogus".

---

    Code
      calculate_fu_in_vitro(method = "pksim_standard", lipophilicity = 3, ionization = c(
        "acid", 0), system = "bogus", pka = c(6, 0), concentration_cells = 2)
    Condition
      Error in `calculate_fu_in_vitro()`:
      ! `system` must be one of "microsomes" or "cells", not "bogus".

# calculate_fu_in_vitro rejects a method not available for the system

    Code
      calculate_fu_in_vitro(method = "kilford", system = "microsomes", lipophilicity = 3,
        ionization = c("neutral", 0), concentration_microsomes = 1)
    Condition
      Error in `calculate_fu_in_vitro()`:
      ! "kilford" is only available for cells.

---

    Code
      calculate_fu_in_vitro(method = "hallifax", system = "cells", lipophilicity = 3,
        ionization = c("neutral", 0), concentration_cells = 1)
    Condition
      Error in `calculate_fu_in_vitro()`:
      ! "hallifax" is only available for microsomes.

# calculate_fu_in_vitro names the missing arguments

    Code
      calculate_fu_in_vitro(method = "austin", system = "cells", lipophilicity = 3,
        ionization = c("neutral", 0), concentration_microsomes = 1)
    Condition
      Error in `calculate_fu_in_vitro()`:
      ! "austin" needs `concentration_cells` for cells.

---

    Code
      calculate_fu_in_vitro(method = "schmitt_fu", system = "microsomes",
        lipophilicity = 3, ionization = c("neutral", 0), concentration_microsomes = 1,
        fbs_fraction = 0, microplate_type = 96, volume_medium = 0.22)
    Condition
      Error in `calculate_fu_in_vitro()`:
      ! "schmitt_fu" needs `fu_plasma` for microsomes.

---

    Code
      calculate_fu_in_vitro(method = "pksim_standard", system = "microsomes",
        lipophilicity = 3, ionization = c("neutral", 0), concentration_microsomes = 1)
    Condition
      Error in `calculate_fu_in_vitro()`:
      ! "pksim_standard" needs `fbs_fraction`, `microplate_type`, and `volume_medium` for microsomes.

---

    Code
      calculate_fu_in_vitro(method = "all_literature", system = "microsomes",
        lipophilicity = 3, ionization = c("base", 0), pka = c(8, 0),
        concentration_microsomes = 1)
    Condition
      Error in `calculate_fu_in_vitro()`:
      ! "all_literature" needs `fu_plasma` and `blood_plasma_ratio` for strong bases (first ionization class "base" with a pKa above 7).

# calculate_fu_in_vitro warns when the compound probably evaporates

    Code
      fu <- calculate_fu_in_vitro(method = "austin", system = "microsomes",
        lipophilicity = 3, ionization = c("neutral", 0), concentration_microsomes = 1,
        microplate_type = 12, volume_medium = 0.5, henry_law_constant = 1)
    Condition
      Warning:
      The compound probably evaporates from the well.

