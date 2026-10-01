# ivive_michaelis_menten rejects invalid inputs

    Code
      ivive_michaelis_menten(system = "hepatocytes", vmax = 2, km = 1, species = "bogus")
    Condition
      Error in `ivive_michaelis_menten()`:
      ! `species` must be one of "human", "rat", or "dog", not "bogus".

---

    Code
      ivive_michaelis_menten(system = "hepatocytes", vmax = 2, km = 1, tissue = "bogus")
    Condition
      Error in `ivive_michaelis_menten()`:
      ! `tissue` must be one of "brain", "gonads", "heart", "kidney", "gut", "liver", or "lung", not "bogus".

---

    Code
      ivive_michaelis_menten(system = "hepatocytes", vmax = 2, km = 1, fu_in_vitro = 1.2)
    Condition
      Error in `ivive_michaelis_menten()`:
      ! `fu_in_vitro` must be above 0 and at most 1, not 1.2.

