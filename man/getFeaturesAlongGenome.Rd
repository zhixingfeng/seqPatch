\name{getFeaturesAlongGenome}
\alias{getFeaturesAlongGenome}
\title{Get IPD for each position of reference genome}
\description{Get IPD for each position of reference genome}
\usage{
getFeaturesAlongGenome (alnsF, alnsIdx)

}
\arguments{
\item{alnsF}{object returned by \code{getAlignmentsWithFeatures} in \code{pbh5} package}
\item{alnsIdx}{object returned by \code{alnIndex} in \code{pbh5} package}

}
\value{
This function a list that contains IPD of each position. 

}
\references{

}

\author{Zhixing Feng \email{fzx08@mails.thu.edu.cn} }


\seealso{
\code{\link{detectModification}},
\code{\link{getContextEffectByPos}},
\code{\link{getGenomeSeq}}
}
\examples{

}
\keyword{SMRT sequencing}
\keyword{polymerase kinetic}
\keyword{DNA modification}
\keyword{hierarchical model}
\keyword{Empirical Bayes}

