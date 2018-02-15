library(soGGi)
context("Test rbind and c give ChIPprofile object")
data(chipExampleBig)
expect_that(class(rbind(chipExampleBig[[1]],chipExampleBig[[1]])) == "ChIPprofile",
            is_true()
)

expect_that(class(c(chipExampleBig[[1]],chipExampleBig[[1]])) == "ChIPprofile",
            is_true()
)