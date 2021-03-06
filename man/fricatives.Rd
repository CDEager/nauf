% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{fricatives}
\alias{fricatives}
\title{Catalan and Spanish intervocalic alveolar fricatives.}
\format{A data frame with 1622 rows and 6 variables:
\describe{
  \item{dur}{The duration of the fricative, as determiend by the onset and
    offset of aperiodic energy in the acoustic signal (milliseconds).}
  \item{pvoi}{The proportion of the fricative which is voiced (measured
    including 30 milliseconds of each of the surrounding vowels).}
  \item{lang}{The language of the map task interview (Catalan or Spanish).}
  \item{wordpos}{Position of the fricative in the word (Final, Initial, or
    Medial).}
  \item{uvoi}{The underlying voicing of the fricative (Neutralized, Voiced,
    or Voiceless).}
  \item{speaker}{Speaker identifier (s01 through s42).}
}}
\usage{
fricatives
}
\description{
A dataset containing duration and voicing measures for /s/ from 16 speakers
of Spanish and for /s/, /z/, and /S/ (the variant that is neutralized in
underlying voicing) for 26 speakers of Catalan.  The data are from a corpus of
map task interviews in Spanish and Catalan and all fricatives are intervocalic.
If you analyze the \code{fricatives} dataset in a publication, please cite
Hualde and Prieto (2014) from the references section below.
}
\section{References}{

Hualde, J. I., & Prieto, P. (2014). Lenition of intervocalic alveolar
fricatives in Catalan and Spanish. Phonetica, 71(2), 109-127.
}
\keyword{datasets}

