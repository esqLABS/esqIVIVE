# calculate_in_vitro_compartments rejects invalid inputs

    Code
      calculate_in_vitro_compartments("bogus", fbs_fraction = 0, microplate_type = 96,
        volume_medium = 0.2)
    Condition
      Error in `calculate_in_vitro_compartments()`:
      ! `system` must be one of "microsomes" or "hepatocytes", not "bogus".

---

    Code
      calculate_in_vitro_compartments("microsomes", fbs_fraction = 0,
        microplate_type = 96, volume_medium = 0.2)
    Condition
      Error in `calculate_in_vitro_compartments()`:
      ! `concentration_microsomes` is needed for microsomes.

---

    Code
      calculate_in_vitro_compartments("microsomes", fbs_fraction = 0,
        microplate_type = 384, volume_medium = 0.2, concentration_microsomes = 1)
    Condition
      Error in `calculate_in_vitro_compartments()`:
      ! `microplate_type` must be one of 96, 48, 24 or 12, not 384.

