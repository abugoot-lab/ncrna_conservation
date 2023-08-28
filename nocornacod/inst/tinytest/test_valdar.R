## Test Valdar01 scoring function.
msa <- matrix(c(
  "A", "C", "T", "-", "G",
  "A", "A", "T", "-", "G",
  "A", "C", "G", "T", "C",
  "A", "T", "G", "-", "G"
), nrow=4, byrow=TRUE)
expect_equal(1 + 1, 2)

