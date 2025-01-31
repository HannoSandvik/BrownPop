######################
# Auxiliary  functions
#  (Hjelpefunksjoner)
######################


# rydda, men ikke "behovsprøvd"


# Naturlig logaritme skal forkortes som "ln"!
ln <- function(x) log(x)


# Dekadisk logaritme skal forkortes som "lg"!
lg <- function(x) log10(x)


# Tester om argumentene er like - ikke følsom for avrundingsfeil!
"%=%" <- function(arg1, arg2) { 
  attributes(arg1) <- NULL
  attributes(arg2) <- NULL
  return(identical(all.equal(arg1, arg2), TRUE))
}


# Tester om argumentene er ulike - ikke følsom for avrundingsfeil!
"%!=%" <- function(arg1, arg2) !(arg1 %=% arg2)


# Limer sammen tekstvariabler
"%+%" <- function(string1, string2) paste0(string1, string2)


# Tester om en vektor er en delmengde av en annen
"%E%" <- function(arg1, arg2) { #¤ sjekk om jeg klarer meg uten!
  e <- rep(F, length(arg1))
  if ((length(arg1) > 0) & (length(arg2) > 0)) {
    i <- 1
    mer <- T
    while (mer) {
      j <- 1
      endamer <- T
      while (endamer) {
        if (arg1[i] %=% arg2[j]) {
          e[i] <- TRUE
        }
        endamer <- (j < length(arg2) & !e[i])
        j <- j + 1
      }
      mer <- (i < length(arg1) & e[i])
      i <- i + 1
    }
  } else {
    e <- arg1 %=% arg2
  }
  return(identical(all(e), TRUE))
}


# Tester om en vektor _ikke_ er en delmengde av en annen
"%!E%" <- function(arg1, arg2) !(arg1 %E% arg2) #¤ sjekk om jeg klarer meg uten!


# Tester om elementene i argumentene er like - ikke følsom for avrundingsfeil 
"%==%" <- function(arg1, arg2) {
  if(is.numeric(arg1) & is.numeric(arg2)) {
    (abs(arg1 - arg2) < 1e-12)
  } else {
    (arg1 == arg2)
  }
}


# Beregner snittmengden av to vektorer
"%A%" <- function(set1, set2)
  if (is.null(set1)) logical(0) else as.vector(na.omit(set1[set1 %in% set2]))


# Beregner unionen av to vektorer
"%V%" <- function(set1, set2) sort(unique(c(set1, set2)))


# Avrunding til nærmeste heltall (men alltid oppover ved halve)
int <- function(x) floor(x + 0.5) #¤ sjekk om jeg klarer meg uten!


# "Pen" avrunding
rund <- function(tall, siffer=0, tegn=0, komma=FALSE, plus=FALSE) {#¤ sjekk om jeg klarer meg uten!
  ### kan henge seg opp når tallet har attributter!!!
  rund <- character(0)
  if (missing(tall)) tall <- NULL
  if (!is.numeric(tall)) {
    if (all(is.number(tall))) {
      tall <- as.numeric(tall)
    } else {
      tall <- NULL
    }
  }
  if (length(siffer) %!=% 1 || !is.numeric(siffer) || is.na(siffer)) {
    siffer <- 0
  } else {
    siffer <- int(siffer)
  }
  if (siffer >  12) siffer <-  12
  if (siffer < -12) siffer <- -12
  if (length(tegn) %!=% 1 || !is.numeric(tegn) || is.na(tegn)) {
    tegn <- 0
  } else {
    tegn <- int(tegn)
  }
  if (komma %!=% TRUE) komma <- FALSE
  if (plus  %!=% TRUE) plus  <- FALSE
  for (j in min(1, length(tall)):length(tall)) {
    t <- tall[j]
    if (is.null(t) || is.na(t) || is.infinite(t)) {
      if (is.null(t)) run <- ""
      if (is.nan(t))  run <- "NaN"
      if (is.na(t))   run <- "NA"
      if (t %=% Inf)  run <- "Inf"
      if (t %=% -Inf) run <- "-Inf"
    } else {
      run <- ""
      r <- int(t * 10^siffer)
      minus <- (r %<% 0)
      r <- abs(r)
      if (r %=% 0) {
        if (siffer > 0) {
          if (komma) {
            run <- "0,"
          } else {
            run <- "0."
          }
          for (i in 1:siffer) run <- run %+% "0"
        } else {
          run <- "0"
        }
      } else {
        s <- floor(lg(r))
        if (s >= siffer) {
          i <- 1
          while (s >= siffer) {
            ru <- substr(as.character(r), i, i)
            i <- i + 1
            if (ru %=% ".") {
              ru <- substr(as.character(r), i, i)
              i <- i + 1
            }
            if (ru %E% c("", "e")) {
              ru <- "0"
              r <- 0
            }
            run <- run %+% ru
            s <- s - 1
          }
          if (siffer > 0) {
            run <- run %+% ifelse(komma, ",", ".")
          }
        } else {
          i <- s + 1
          run <- ifelse(komma, "0,", "0.")
          while (i < siffer) {
            run <- run %+% "0"
            i <- i + 1
          }
          i <- 1
        }
        while (s >= 0) {
          ru <- substr(as.character(r), i, i)
          i <- i + 1
          if (ru %=% ".") {
            ru <- substr(as.character(r), i, i)
            i <- i + 1
          }
          if (ru %E% c("", "e")) {
            ru <- "0"
            r <- 0
          }
          run <- run %+% ru
          s <- s - 1
        }
      }
      if (minus) {
        run <- "-" %+% run
      } else {
        if (plus) {
          run <- "+" %+% run
        }
      }
    }
    #is (is.null(siffer)) #¤?
    while (nchar(run) < tegn) run <- " " %+% run
    rund[j] <- run
  }
  return(rund)
}
