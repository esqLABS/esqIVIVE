# ivive_michaelis_menten: a tissue without a scaling factor gives NA vmax and a warning

    Code
      mm <- ivive_michaelis_menten(system = "cells", vmax = 2, km = 1, tissue = "brain")
    Condition
      Warning in `ivive_michaelis_menten()`:
      The scaling factor CellsGO (the cells per gram tissue) is not supported for species "human" and tissue "brain": it is `NA` in the scaling factor table.
      i The result is `NA`.

# ivive_michaelis_menten rejects invalid inputs

    Code
      ivive_michaelis_menten(system = "cells", vmax = 2, km = 1, species = "bogus")
    Condition
      Error in `ivive_michaelis_menten()`:
      ! `species` must be one of "human", "rat", "dog", or "beagle", not "bogus".

---

    Code
      ivive_michaelis_menten(system = "cells", vmax = 2, km = 1, tissue = "bogus")
    Condition
      Error in `ivive_michaelis_menten()`:
      ! `tissue` must be one of "brain", "gonads", "heart", "kidney", "gut", "liver", or "lung", not "bogus".

---

    Code
      ivive_michaelis_menten(system = "cells", vmax = 2, km = 1, fu_in_vitro = 1.2)
    Condition
      Error in `ivive_michaelis_menten()`:
      ! `fu_in_vitro` must be above 0 and at most 1, not 1.2.

