#' Protein sequences of the algae species Galdieria sulphuraria
#'
#' Data obtained from PLAZA Diatoms. Only genes containing domains used for
#' TF family classification (with PlantTFDB's or TAPscan's schemes) were kept
#' for package size issues.
#'
#' @name gsu
#' @format An AAStringSet object as returned 
#' by \code{Biostrings::readAAStringSet()}.
#' @references
#' Osuna-Cruz, C. M., Bilcke, G., Vancaester, E., De Decker, S., Bones, A. M., 
#' Winge, P., ... & Vandepoele, K. (2020). The Seminavis robusta genome 
#' provides insights into the evolutionary adaptations of benthic diatoms. 
#' Nature communications, 11(1), 1-13.
#' @examples
#' data(gsu)
#' @usage data(gsu)
"gsu"


#' Data frame of PlantTFDB's TF family classification scheme
#'
#' The classification scheme is the same as the one used by PlantTFDB,
#' except for three family/subfamily names that were changed to match the
#' output of \code{classify_tfs()} ('M_type' to 'M-type', 'MYB_related' to
#' 'MYB-related', and 'LBD (AS2/LOB)' to 'LBD'). As in PlantTFDB, families
#' without subfamilies have the family name as subfamily.
#' 
#' @name planttfdb_scheme
#' @format A data frame with the following variables:
#' \describe{
#'  \item{Family}{TF family name.}
#'  \item{Subfamily}{TF subfamily name.}
#'  \item{DBD}{DNA-binding domain}
#'  \item{Auxiliary}{Auxiliary domain}
#'  \item{Forbidden}{Forbidden domain}
#' }
#' @references 
#' Jin, J., Tian, F., Yang, D. C., Meng, Y. Q., Kong, L., Luo, J., & 
#' Gao, G. (2016). PlantTFDB 4.0: toward a central hub for transcription 
#' factors and regulatory interactions in plants. 
#' Nucleic acids research, gkw982.
#' @examples 
#' data(planttfdb_scheme)
#' @usage data(planttfdb_scheme)
"planttfdb_scheme"


#' Domain annotation for the algae species Galdieria sulphuraria
#'  
#' The data set was created using the function \code{annotate_domains()}.
#'
#' @name gsu_annotation
#' @format A list of 2 data frames (\strong{PlantTFDB} and
#' \strong{TAPscan}) with domain annotation, as returned
#' by \code{annotate_domains()}.
#' @examples
#' data(gsu_annotation)
#' @usage data(gsu_annotation)
"gsu_annotation"


#' Classes of transcription-associated proteins (TAPs) in TAPscan
#'
#' Families of transcription-associated proteins (TAPs) in TAPscan v4, their
#' classes, descriptions, and references. Data were obtained from TAPscan's
#' website (file \emph{_data/import-tapinfo/tapinfo_v4.csv} in the
#' GitHub repository Rensing-Lab/TAPscan-v4-website). HTML tags in
#' descriptions and references were removed. Family MADS_MIKC, which is not
#' included in the original file, was manually added as a TF family (as
#' family MADS), without description and references.
#'
#' @name tap_class
#' @format A data frame with the following variables:
#' \describe{
#'  \item{TAP_family}{TAP family or subfamily, as in the variable
#'  \strong{TAPscan_subfamily} of the output of \code{classify_tfs()}.}
#'  \item{TAP_class}{TAP class. One of 'TF' (transcription factor),
#'  'TR' (transcriptional regulator), or 'PT' (putative TAP).}
#'  \item{Description}{Description of the TAP family.}
#'  \item{References}{References for the TAP family, separated by ' | '.}
#'}
#' @references
#' Petroll, R., Varshney, D., Hiltemann, S., Finke, H., Schreiber, M.,
#' de Vries, J., & Rensing, S. A. (2025). Enhanced sensitivity of TAPscan v4
#' enables comprehensive analysis of streptophyte transcription factor
#' evolution. The Plant Journal, 121(1), e17184.
#' @examples
#' data(tap_class)
#' @usage data(tap_class)
"tap_class"


#' Data frame of TAPscan's classification scheme
#'
#' TAPscan v4 classification rules, with one row per rule. Genes are
#' assigned to a family (and subfamily, if any) if they have all domains
#' in variable \strong{Should} and none of the domains in variable
#' \strong{Should_not}. For families with more than one rule (e.g., bZIP),
#' genes are assigned to the family if they satisfy any of the rules.
#'
#' @name tapscan_scheme
#' @format A data frame with the following variables:
#' \describe{
#'  \item{Family}{TAPscan family, as in variable \strong{TAPscan_family} of
#'  the output of \code{classify_tfs()}.}
#'  \item{Subfamily}{TAPscan subfamily, as in variable
#'  \strong{TAPscan_subfamily} of the output of \code{classify_tfs()}.
#'  Families without subfamilies have the family name as subfamily.}
#'  \item{Rule}{Name of the rule in TAPscan's rules file.}
#'  \item{Should}{Domains that genes should have, separated by ' and '.
#'  Numbers in parentheses (e.g., 'Myb_DNA-binding (2)') indicate the
#'  required number of copies of a domain.}
#'  \item{Should_not}{Domains that genes should not have, separated
#'  by ' or '.}
#'  \item{TAP_class}{TAP class. One of 'TF' (transcription factor),
#'  'TR' (transcriptional regulator), or 'PT' (putative TAP). NA for
#'  families without a TAP class in TAPscan (see \code{tap_class}).}
#'}
#' @references
#' Petroll, R., Varshney, D., Hiltemann, S., Finke, H., Schreiber, M.,
#' de Vries, J., & Rensing, S. A. (2025). Enhanced sensitivity of TAPscan v4
#' enables comprehensive analysis of streptophyte transcription factor
#' evolution. The Plant Journal, 121(1), e17184.
#' @examples
#' data(tapscan_scheme)
#' @usage data(tapscan_scheme)
"tapscan_scheme"


#' Correspondence table between TAPscan and Plant-TFClass
#'
#' Table used by \code{classify_tfs()} to infer Plant-TFClass superclasses,
#' classes, and families from TAPscan families and subfamilies.
#' TAPscan families/subfamilies with multiple rows correspond to more than
#' one Plant-TFClass family, and the PlantTFDB subfamily of each gene
#' (variable \strong{PlantTFDB_subfamily}) is used to choose among them.
#'
#' @name planttfclass_scheme
#' @format A data frame with the following variables:
#' \describe{
#'  \item{TAPscan_family}{TAPscan family.}
#'  \item{TAPscan_subfamily}{TAPscan subfamily (or family, for families
#'  without subfamilies).}
#'  \item{PlantTFClass_superclass}{Plant-TFClass superclass.}
#'  \item{PlantTFClass_class}{Plant-TFClass class.}
#'  \item{PlantTFClass_family}{Plant-TFClass family (or class, for classes
#'  without families).}
#'  \item{PlantTFDB_subfamily}{PlantTFDB subfamily used to choose among
#'  multiple Plant-TFClass families for the same TAPscan family/subfamily.
#'  NA for TAPscan families/subfamilies with a single Plant-TFClass family.}
#'}
#' @references
#' Blanc-Mathieu, R., Dumas, R., Turchi, L., Lucas, J., & Parcy, F. (2024).
#' Plant-TFClass: a structural classification for plant transcription
#' factors. Trends in Plant Science, 29(1), 40-51.
#' @examples
#' data(planttfclass_scheme)
#' @usage data(planttfclass_scheme)
"planttfclass_scheme"
