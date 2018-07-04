Sys.setenv("R_TESTS" = "")

library(testthat)
library(mimosa)

test_check("mimosa")
