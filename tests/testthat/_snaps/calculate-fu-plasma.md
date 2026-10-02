# calculate_fu_plasma rejects an invalid species

    Code
      calculate_fu_plasma(partition_albumin = 1, partition_globulin = 1,
        partition_membrane_lipids = 1, partition_neutral_lipids = 1, species = "bogus")
    Condition
      Error in `calculate_fu_plasma()`:
      ! `species` must be one of "human", "rat", "dog", "monkey", "rabbit", or "mouse", not "bogus".

