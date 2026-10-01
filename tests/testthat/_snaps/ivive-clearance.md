# ivive_clearance rejects a unit that does not fit the value type

    Code
      ivive_clearance(value_type = "rate_constant", unit = "mL/minutes/millioncells",
        value = 18.27, system = "cells", fu_in_vitro = 0.5, concentration_cells = 0.5)
    Condition
      Error in `ivive_clearance()`:
      ! `unit` must be one of "/minutes", "/hours", or "/seconds", not "mL/minutes/millioncells".

---

    Code
      ivive_clearance(value_type = "intrinsic_clearance", unit = "mL/minutes/millioncell",
        value = 18.27, system = "cells", concentration_cells = 0.5)
    Condition
      Error in `ivive_clearance()`:
      ! `unit` must be one of "mL/minutes/millioncells", "uL/minutes/millioncells", "L/minutes/millioncells", "mL/hours/millioncells", "uL/hours/millioncells", "L/hours/millioncells", "mL/seconds/millioncells", "uL/seconds/millioncells", "L/seconds/millioncells", "mL/minutes/cell", "uL/minutes/cell", "mL/hours/cell", "uL/hours/cell", "mL/seconds/cell", "uL/seconds/cell", "mL/minutes/mg protein", "uL/minutes/mg protein", "L/minutes/mg protein", "mL/hours/mg protein", "uL/hours/mg protein", "L/hours/mg protein", "mL/seconds/mg protein", "uL/seconds/mg protein", "L/seconds/mg protein", "mL/minutes", "uL/minutes", "mL/seconds", "uL/seconds", "mL/hours", "uL/hours", "mL/minutes/kg", "uL/minutes/kg", "mL/hours/kg", or "uL/hours/kg", not "mL/minutes/millioncell".
      i Did you mean "mL/minutes/millioncells"?

# ivive_clearance rejects a missing concentration

    Code
      ivive_clearance(value_type = "half_life", unit = "minutes", value = 3.9,
        system = "microsomes")
    Condition
      Error in `ivive_clearance()`:
      ! `concentration_microsomes` is needed for microsomes when `value_type` is "half_life" and `unit` is "minutes".

---

    Code
      ivive_clearance(value_type = "intrinsic_clearance", unit = "uL/minutes", value = 3.9,
        system = "cells")
    Condition
      Error in `ivive_clearance()`:
      ! `concentration_cells` is needed for cells when `value_type` is "intrinsic_clearance" and `unit` is "uL/minutes".

# ivive_clearance rejects invalid options

    Code
      ivive_clearance(value_type = "half_life", unit = "minutes", value = 3.9,
        system = "microsomes", concentration_microsomes = 1, empirical_correction = "Yes")
    Condition
      Error in `ivive_clearance()`:
      ! `empirical_correction` must be `TRUE` or `FALSE`, not "Yes".

---

    Code
      ivive_clearance(value_type = "half_life", unit = "minutes", value = 3.9,
        system = "microsomes", concentration_microsomes = 1, fu_in_vitro = 0)
    Condition
      Error in `ivive_clearance()`:
      ! `fu_in_vitro` must be above 0 and at most 1, not 0.

---

    Code
      ivive_clearance(value_type = "half_life", unit = "minutes", value = 3.9,
        system = "microsomes", concentration_microsomes = 1, species = "dog",
        empirical_correction = TRUE)
    Condition
      Error in `ivive_clearance()`:
      ! The empirical correction is only available for human and rat, not "dog".

# ivive_clearance rejects an invalid species or tissue

    Code
      ivive_clearance(value_type = "half_life", unit = "hours", value = 3, system = "microsomes",
        concentration_microsomes = 0.5, species = "bogus")
    Condition
      Error in `ivive_clearance()`:
      ! `species` must be one of "human", "rat", or "dog", not "bogus".

---

    Code
      ivive_clearance(value_type = "half_life", unit = "hours", value = 3, system = "microsomes",
        concentration_microsomes = 0.5, tissue = "bogus")
    Condition
      Error in `ivive_clearance()`:
      ! `tissue` must be one of "brain", "gonads", "heart", "kidney", "gut", "liver", or "lung", not "bogus".

# ivive_clearance: cytosol needs its own concentration

    Code
      ivive_clearance(value_type = "half_life", unit = "minutes", value = 3.9,
        system = "cytosol", concentration_microsomes = 1)
    Condition
      Error in `ivive_clearance()`:
      ! `concentration_cytosol` is needed for cytosol when `value_type` is "half_life" and `unit` is "minutes".

# ivive_clearance: the empirical correction is not available for the cytosol

    Code
      ivive_clearance(value_type = "intrinsic_clearance", system = "cytosol", unit = "mL/minutes/mg protein",
        value = 0.05, concentration_cytosol = 1, empirical_correction = TRUE)
    Condition
      Error in `ivive_clearance()`:
      ! The empirical correction is only available for microsomes and cells, not "cytosol".

# ivive_clearance rejects the old hepatocytes system

    Code
      ivive_clearance(value_type = "half_life", unit = "hours", value = 3, system = "hepatocytes",
        concentration_cells = 0.5)
    Condition
      Error in `ivive_clearance()`:
      ! `system` must be one of "microsomes", "cells", or "cytosol", not "hepatocytes".

# ivive_clearance prints only the inputs when verbose

    Code
      cl <- ivive_clearance(value_type = "intrinsic_clearance", value = 18.27, unit = "mL/minutes/millioncells",
        system = "cells", concentration_cells = 0.5, verbose = TRUE)
    Output
      --- ivive_clearance ---
      Inputs:
        value_type = intrinsic_clearance
        value = 18.27
        unit = mL/minutes/millioncells
        system = cells
        fu_in_vitro = 1
        empirical_correction = FALSE
        tissue = liver
        species = human

