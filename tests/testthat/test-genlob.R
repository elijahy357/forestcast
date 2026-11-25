test_that("base function works", {
  expect_equal(genlob(10, 10, 1000, 10, 10, 18), data.frame(runID = 1,
                                                            age = 10,
                                                            n = 1000,
                                                            ht = 10,
                                                            ba = 10
                                                            ))
})
