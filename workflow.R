
Sys.setlocale("LC_CTYPE", "nb_NO.utf8")

source("aux.R")
source("brownianplot.R")
source("popmod.R")
source("pva.R")
source("pvaplot.R")


krykkje.tall <- c(29498, 28207, 27086, 33589, NA, 25045, NA, 27355, 27070, 28827, 25314, 25789, 24634, 24318, 22466, 23843, 23859, 25662, 26169, 28257, 24792, 19128, 21628, 18574, 16976, 15505, 11724, 10996, 10395, 10081, 10717, 10525, 9661)
krykkje.modell <- popmod(krykkje.tall, 1980:2012, nboot=1000)
pva(krykkje.modell, crash.prop=0.25, crash.year=1, n=1000)

