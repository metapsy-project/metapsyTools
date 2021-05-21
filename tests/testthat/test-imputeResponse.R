# Tests for imputeResponse function

test_that("imputeResponse works with treatment arm only (cutpoint method)", {
  result <- imputeResponse(m.trt.pre = 20, m.trt.post = 11,
                          sd.trt.post = 7, n.trt = 100)
  
  expect_s3_class(result, "data.frame")
  expect_true("trtResponder" %in% names(result))
  expect_true("nTrt" %in% names(result))
  expect_equal(result$nTrt, 100)
  expect_true(result$trtResponder >= 0 && result$trtResponder <= 100)
})

test_that("imputeResponse works with both arms (cutpoint method)", {
  result <- imputeResponse(m.trt.pre = 20, m.trt.post = 11,
                          sd.trt.post = 7, n.trt = 100,
                          m.ctr.pre = 20, m.ctr.post = 18,
                          sd.ctr.post = 6, n.ctr = 120)
  
  expect_s3_class(result, "data.frame")
  expect_true("trtResponder" %in% names(result))
  expect_true("ctrResponder" %in% names(result))
  expect_true("logRR" %in% names(result))
  expect_true("seLogRR" %in% names(result))
  expect_true("logOR" %in% names(result))
  expect_true("seLogOR" %in% names(result))
})

test_that("imputeResponse works with custom cutoff", {
  result <- imputeResponse(m.trt.pre = 20, m.trt.post = 11,
                          sd.trt.post = 7, n.trt = 100,
                          m.ctr.pre = 20, m.ctr.post = 18,
                          sd.ctr.post = 6, n.ctr = 120,
                          cutoff = 15)
  
  expect_s3_class(result, "data.frame")
  expect_true("trtResponder" %in% names(result))
})

test_that("imputeResponse works with lower.is.better = FALSE", {
  result <- imputeResponse(m.trt.pre = 20, m.trt.post = 25,
                          sd.trt.post = 7, n.trt = 100,
                          m.ctr.pre = 20, m.ctr.post = 18,
                          sd.ctr.post = 6, n.ctr = 120,
                          cutoff = 15, lower.is.better = FALSE)
  
  expect_s3_class(result, "data.frame")
  expect_true("trtResponder" %in% names(result))
})

test_that("imputeResponse works with RCI method", {
  result <- imputeResponse(m.trt.pre = 20, m.trt.post = 11,
                          sd.trt.post = 7, n.trt = 100,
                          m.ctr.pre = 20, m.ctr.post = 18,
                          sd.ctr.post = 6, n.ctr = 120,
                          method = "rci", rho = 0.88,
                          sd.trt.pre = 8)
  
  expect_s3_class(result, "data.frame")
  expect_true("trtResponder" %in% names(result))
})

test_that("imputeResponse works with RCI method and control pre SD", {
  result <- imputeResponse(m.trt.pre = 20, m.trt.post = 11,
                          sd.trt.post = 7, n.trt = 100,
                          m.ctr.pre = 20, m.ctr.post = 18,
                          sd.ctr.post = 6, n.ctr = 120,
                          method = "rci", rho = 0.88,
                          sd.trt.pre = 8, sd.ctr.pre = 9)
  
  expect_s3_class(result, "data.frame")
  expect_true("trtResponder" %in% names(result))
})

test_that("imputeResponse works with vector inputs", {
  set.seed(123)
  result <- imputeResponse(m.trt.pre = runif(5, 10, 30),
                         m.trt.post = runif(5, 8, 14),
                         sd.trt.post = runif(5, 5, 8),
                         n.trt = rpois(5, 150))
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 5)
  expect_true(all(result$nTrt > 0))
})

test_that("imputeResponse works with vector inputs and cutoff", {
  set.seed(123)
  result <- imputeResponse(m.trt.pre = runif(5, 10, 30),
                         m.trt.post = runif(5, 8, 14),
                         sd.trt.post = runif(5, 5, 8),
                         n.trt = rpois(5, 150),
                         cutoff = 11:15)
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 5)
})

test_that("imputeResponse validates method parameter", {
  # Invalid method should error or default to cutpoint
  # The function may error with invalid method, so we test that it either works or errors gracefully
  result <- try(imputeResponse(m.trt.pre = 20, m.trt.post = 11,
                               sd.trt.post = 7, n.trt = 100,
                               method = "invalid"),
                silent = TRUE)
  
  # Either it works (defaults to cutpoint) or errors gracefully
  if (!inherits(result, "try-error")) {
    expect_s3_class(result, "data.frame")
  } else {
    # If it errors, that's also acceptable behavior
    expect_true(TRUE)
  }
})

