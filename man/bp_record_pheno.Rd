% Generated manually to document bp_record_pheno()
\name{bp_record_pheno}
\alias{bp_record_pheno}
\title{Record Phenotypes in Long Format}
\usage{
bp_record_pheno(
  state,
  cohort_id,
  stage,
  individual_id,
  traits = NULL,
  pheno_matrix,
  available_tick,
  n_loc,
  reps,
  environment = NULL,
  p_value = NA_real_,
  drop_na = TRUE
)
}
\arguments{
\item{state}{Program state.}
\item{cohort_id}{Output cohort id.}
\item{stage}{Output stage name.}
\item{individual_id}{Integer vector of individual IDs; length must match rows in \code{pheno_matrix}.}
\item{traits}{Trait labels, either integer trait codes or character names. If \code{NULL}, labels are inferred from \code{colnames(pheno_matrix)} when present, otherwise default to \code{trait1}, \code{trait2}, ...}
\item{pheno_matrix}{Numeric vector or matrix of phenotypes. If a vector is provided, it is treated as a one-column matrix.}
\item{available_tick}{Tick when these observations become available.}
\item{n_loc}{Number of locations represented by this record block.}
\item{reps}{Number of reps represented by this record block.}
\item{environment}{Environment id/label(s) for this record block. Can be a scalar or a vector. If a vector is provided, it is aligned to phenotype columns with the same recycling rules as \code{traits}. If \code{NULL}, defaults to \code{NA}.}
\item{p_value}{Optional environment value(s) for this record block. Can be a
scalar or a vector. If a vector is provided, it is aligned to phenotype
columns with the same recycling rules as \code{traits}.}
\item{drop_na}{Logical; if \code{TRUE}, rows with \code{NA} phenotype values are omitted from the appended log.}
}
\value{Updated program state.}
\description{
Append phenotype observations to \code{state$phenotype_log}. This helper accepts
a dense or sparse phenotype matrix with one row per individual and one column
per trait label in \code{traits}. Output is stored in long format (one row per
observed phenotype cell).
}
\details{
For sparse MET workflows, pass trait columns that encode trial/environment
identity (for example \code{yield_env1}, \code{yield_env2}, ...) and set
\code{drop_na = TRUE} so missing cells are not logged.
}
\examples{
# Example sketch:
# state <- bp_record_pheno(
#   state = state,
#   cohort_id = "cohort_0001",
#   stage = "PYT",
#   individual_id = c(101, 102, 103),
#   traits = c("yield_env1", "yield_env2"),
#   pheno_matrix = matrix(c(4.1, NA, 3.9, 4.0, NA, 4.2), ncol = 2),
#   available_tick = state$time$tick + 4,
#   n_loc = 2,
#   reps = 1,
#   environment = NA_integer_,
#   p_value = NA_real_,
#   drop_na = TRUE
# )
}
