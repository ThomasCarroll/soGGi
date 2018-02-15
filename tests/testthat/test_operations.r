library(soGGi)
context("Arithmetic Operation")
data(chipExampleBig)
expect_that(all(assays((chipExampleBig[[1]]+chipExampleBig[[2]])/2)[[1]]==
               assays(mean(chipExampleBig[[1]],chipExampleBig[[2]]))[[1]]),is_true()
)
