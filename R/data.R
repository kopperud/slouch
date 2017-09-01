#' Artiodactyl brain data 
#'
#' Morphological data on mean neocortex area (mm^2), mean brain size (g) and mean body size (g) for 43 species, including estimates of within-species observational error. These standard errors are not based directly on the unbiased sample variance estimator. The sample variances were first calculated, then a "global sample variance" was estimated using a weighted average. Lastly, the global sample variance was divided by the respective within-species sample size for each species to obtain the squared standard error as reported in the dataset.
#'
#' @docType data
#'
#' @usage data(neocortex)
#'
#' @format An object of class \code{"data.frame"}.
#'
#' @keywords datasets
#'
#' @references 
#' 
#' The following literature includes details on the original data collection.
#' 
#' \itemize{
#' \item{Haarmann, K., & Oboussier, H. (1972). Morphologishce und quantitative Neocortexuntersuchungen bei Boviden, ein Beitrag zur Phylogenie dieser Familie. II. Formen geringen Körpergewichts (3kg - 25kg) aus den Subfamilien Cephalophinae und Antilopinae. Mitteilungen Aus Dem Hamburgischen Zoologischen Museum Und Institut, 68, 231–269.}
#' 
#' \item{Oboussier, H. (1972). Morphologische und quantitative Neocortexuntersuchungen bei Boviden, ein Beitrag zur Phylogenie dieser Familie III. Formen über 75 kg Körpergewicht. Mitteilungen Aus Dem Hamburgischen Zoologischen Museum Und Institut, 68, 271–292.}
#' 
#' \item{Oboussier, H. (1978). Zur Kenntnis des Bergnyalas (Tragelaphus buxtoni) und des Bongos (Taurotragus euryceros). Untersuchungen über den Körperbau und das Gehirn. Zeitschrift Für Säugetierkunde, 43, 114–125.}
#' 
#' \item{Oboussier, H., & Möller, G. (1971). Zur Kenntnis des Gehirns der Giraffidae (Pecora, Artiodactyla, Mammalia) - ein Vergleich der Neocortex-Oberflåachengrösse). Zeitschrift Für Säugetierkunde, 36, 291–296.}
#' 
#' \item{Ronnefeld, U. (1970). Morphologische und quantitative Neocortexuntersuchungen bei Boviden, ein Beitrag zur Phylogenie dieser Familie. I. Formen mittlerer Körpergrösse (25 kg bis 75 kg). Gegenbaurs Morphologische Jahrbuch, 161–230.}
#'
#' }
#' 
#' See also:
#' \itemize{
#' 
#' \item{Haarmann, K. (1975). Morphological and histological study of neocortex of bovides (Antilopinae, Cephalophinae) and Tragulidae with comments on evolutionary development. Journal Fur Hirnforschung, 16, 93–116.}
#' 
#' \item{Oboussier, H. (1979). Evolution of the brain and phylogenetic development of African Bovidae. South African Journal of Zoology, 14(3), 119–124. https://doi.org/10.1080/02541858.1979.11447660}
#' 
#' }
#'
#' @examples
#' data(neocortex)
#' plot(neocortex$brain_mass_g_log_mean, neocortex$neocortex_area_mm2_log_mean)
#' 
"neocortex"