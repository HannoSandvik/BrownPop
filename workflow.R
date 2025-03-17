
Sys.setlocale("LC_CTYPE", "nb_NO.utf8")

source("BrownPop.R")
source("BrownPVA.R")

krykkje.tall <- c(29498, 28207, 27086, 33589, NA, 25045, NA, 27355, 27070, 28827, 25314, 25789, 24634, 24318, 22466, 23843, 23859, 25662, 26169, 28257, 24792, 19128, 21628, 18574, 16976, 15505, 11724, 10996, 10395, 10081, 10717, 10525, 9661)

popmod <- BrownPop(krykkje.tall, 1980:2012, nboot = 100000, name = "Kittiwake")

pva <- BrownPVA(popmod, crash.prop = 0.25)


# ---------------------
# Eller f.eks.:

pva <- BrownPVA(popmod, crash.prop = c(0, 0.1, 0.25, 0.5))
# eller
pva <- BrownPVA(popmod, crash.prop = seq(0, 0.25, 0.01))

print(pva)

