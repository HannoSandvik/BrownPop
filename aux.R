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



#¤¤¤¤¤

heltall <- function(x, m = 6, n = 4) {
  if (is.numeric(x)) x <- round(x)
  x <-  as.character(x)
  while(nchar(x) < m)           x <- " " %+%  x
  while(nchar(x) < m + 1 + n)   x <-  x  %+% " "
  return(x)
}

destall <- function(x, m = 6, n = 4, z = 3, des = ",") {
  x <- round(as.numeric(x), z)
  z <- min(z, n)
  x <- unlist(strsplit(as.character(x), ".", fixed = T))
  en <- x[1]
  while(nchar(en) < m) en <- " " %+% en
  to <- x[2]
  if (is.na(to)) {
    if (z %=% 0) {
      to <- ""
    } else {
      to <- des %+% paste(rep("0", z), collapse = "")
    }
  } else {
    while(nchar(to) < z)   to <- to %+% "0"
    to <- des %+% to
  }
  while(nchar(to) < n + 1) to <- to %+% " "
  return(en %+% to)
}

