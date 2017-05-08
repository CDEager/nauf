

#' Spanish intervocalic plosives.
#'
#' A dataset containing measures of plosive strength
#' for instances of intervocalic Spanish /p/, /t/, /k/, /b/, /d/ and /g/.
#' The data are taken from read speech and informal interviews of 30 speakers in
#' Cuzco, Peru and 8 speakers in Lima, Peru; and from 18 speakers from
#' Valladolid, Spain in the task dialogues in the Spanish portion of the
#' Glissando Corpus (Garrido et al. 2013). If you analyze the \code{plosives}
#' dataset in a publication, please cite Eager (2017) from the references section below.
#'
#' @format A data frame with 5281 rows and 21 variables:
#' \describe{
#'   \item{cdur}{Total plosive duration, measured from preceding vowel intensity
#'     maximum to following vowel intensity maximum, in milliseconds. Set to
#'     \code{0} for elided plosives.}
#'   \item{vdur}{Duration of the period of voicelessness in the
#'     vowel-consonant-vowel sequence in milliseconds. Set to \code{0} for
#'     fully voiced plosives and elided plosives.}
#'   \item{vpct}{Percentage of the consonant duration which was voiceless.
#'     For non-elided plosives, \code{vpct = vdur / cdur}, and for elided
#'     plosives, \code{vpct = 0}.}
#'   \item{intdiff}{The maximum intensity in the following vowel minus the
#'     minimum intensity in the plosive, in decibels.  Set to \code{0} for
#'     elided plosives.}
#'   \item{intvel}{The maximum rising velocity of the intensity contour between
#'     the consonant minimum intensity and following vowel maximum intensity,
#'     in decibels per millisecond.  Set to \code{0} for elided plosives.}
#'   \item{voicing}{The underlying voicing of the plosive (Voiced or
#'     Voiceless).}
#'   \item{place}{Place of articulation (Bilabial, Dental, or Velar).}
#'   \item{stress}{Syllabic stress context (Tonic, Post-Tonic, or Unstressed).}
#'   \item{prevowel}{Preceding vowel phoneme identity (a, e, i, o, or u).}
#'   \item{posvowel}{Following vowel phoneme identity (a, e, i, o, or u).}
#'   \item{wordpos}{Position of the plosive in the word (Initial or Medial).}
#'   \item{wordfreq}{Number of times the word containing the plosive occurs in
#'     the CREA corpus (Real Academia Espanola).}
#'   \item{speechrate}{Local speech rate around the consonant in nuclei per
#'     second (nuclei located using De Jong and Wempe's (2009) script).}
#'   \item{spont}{Whether the speech was spontaneous (TRUE) or read (FALSE).}
#'   \item{dialect}{The city where the speaker is from (Cuzco, Lima, or
#'     Valladolid).}
#'   \item{sex}{The speaker's sex (Female or Male).}
#'   \item{age}{The speaker's age group (Older or Younger) based on whether
#'     they were over 40 years old, or 40 years old or younger at the time of
#'     recording.}
#'   \item{ed}{The speaker's highest level of education (Secondary or
#'     University).}
#'   \item{ling}{The speaker's language background (Bilingual or Monolingual)
#'     based on whether they spoke only Spanish, or both Spanish and Quechua.}
#'   \item{speaker}{Speaker identifier (s01 through s56).}
#'   \item{item}{Read speech item identifier (i01 through i54). Set to \code{NA}
#'     for spontaneous speech.}
#' }
#'
#' @section Note: The \code{\link[standardize]{ptk}} dataset in the
#' \code{standardize} package is a subset of the \code{plosives} dataset, but
#' with the speakers renumbered:
#'
#' \preformatted{
#' d <- droplevels(subset(plosives,
#'   dialect == "Valladolid" & voicing == "Voiceless"))
#'
#' levels(d$speaker)  # s39 to s56
#' levels(ptk$speaker)  # s01 to s18
#'
#' levels(d$speaker) <- levels(ptk$speaker)
#' d <- d[, colnames(ptk)]
#' rownames(d) <- NULL
#'
#' all.equal(d, ptk)  # TRUE
#' }
#'
#' @section References:
#' Eager, Christopher D. (2017). Contrast preservation and constraints on
#' individual phonetic variation. Doctoral thesis. University of Illinois at
#' Urbana-Champaign.
#'
#' Garrido, J. M., Escudero, D., Aguilar, L., Cardenoso, V., Rodero, E.,
#' de-la-Mota, C., ... Bonafonte, A. (2013). Glissando: a corpus for
#' multidisciplinary prosodic studies in Spanish and Catalan. Language Resources
#' and Evaluation, 47(4), 945-971.
#'
#' Real Academia Espanola. Corpus de referencia del espanol actual (CREA). Banco
#' de Datos. http://www.rae.es
#'
#' De Jong, N. H., & Wempe, T. (2009). Praat script to detect syllable nuclei
#' and measure speech rate automatically. Behavior Research Methods, 41(2),
#' 385-390.
"plosives"


#' Catalan and Spanish intervocalic alveolar fricatives.
#'
#' A dataset containing duration and voicing measures for /s/ from 16 speakers
#' of Spanish and for /s/, /z/, and /S/ (the variant that is neutralized in
#' underlying voicing) for 26 speakers of Catalan.  The data are from a corpus of
#' map task interviews in Spanish and Catalan and all fricatives are intervocalic.
#' If you analyze the \code{fricatives} dataset in a publication, please cite
#' Hualde and Prieto (2014) from the references section below.
#'
#' @format A data frame with 1622 rows and 6 variables:
#' \describe{
#'   \item{dur}{The duration of the fricative, as determiend by the onset and
#'     offset of aperiodic energy in the acoustic signal (milliseconds).}
#'   \item{pvoi}{The proportion of the fricative which is voiced (measured
#'     including 30 milliseconds of each of the surrounding vowels).}
#'   \item{lang}{The language of the map task interview (Catalan or Spanish).}
#'   \item{wordpos}{Position of the fricative in the word (Final, Initial, or
#'     Medial).}
#'   \item{uvoi}{The underlying voicing of the fricative (Neutralized, Voiced,
#'     or Voiceless).}
#'   \item{speaker}{Speaker identifier (s01 through s42).}
#' }
#'
#' @section References:
#' Hualde, J. I., & Prieto, P. (2014). Lenition of intervocalic alveolar
#' fricatives in Catalan and Spanish. Phonetica, 71(2), 109-127.
"fricatives"

