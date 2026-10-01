# calculate_plasma_partitions rejects invalid inputs

    Code
      calculate_plasma_partitions(method = "bogus", lipophilicity = 2, pka = c(3, 0),
      ionization = c("acid", 0))
    Condition
      Error in `calculate_plasma_partitions()`:
      ! `method` must be one of "logp" or "pplfer", not "bogus".

---

    Code
      calculate_plasma_partitions(method = "pplfer", lipophilicity = 2, pka = c(3, 0),
      ionization = c("acid", 0), lfer_e = 1)
    Condition
      Error in `calculate_plasma_partitions()`:
      ! "pplfer" needs `lfer_b`, `lfer_a`, `lfer_s`, and `lfer_v`.

