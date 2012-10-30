\name{detectModification}
\alias{seqPatch}
\title{Detecting DNA modification from SMRT sequencing data}
\description{Detect DNA modification by hierarchical model}
\usage{
detectModification (genomeF.native, genomeF.ctrl, genomeSeq, context.effect=NULL, left.len=6, right.len=1,
                 is.filter.sample=FALSE, is.trim.cvg = TRUE, is.trim.pos = TRUE, min.cvg=5, min.pos=4,  method='hieModel')

}
\arguments{
\item{genomeF.native}{genomeF object of the native sample, returned by \code{getFeaturesAlongGenome}}
\item{genomeF.ctrl}{genomeF object of the control sample, returned by \code{getFeaturesAlongGenome}}
\item{genomeSeq}{genomeSeq object of reference genome, returned by \code{getGenomeSeq}}
\item{context.effect}{context.effect object that records context effect of historical data, returned by \code{getContextEffectByPos}}
\item{left.len}{Number of upstream bases in sequence context}
\item{right.len}{Number of downstream bases in sequence context}
\item{is.filter.sample}{Whether to filter outliers for native and control sample}
\item{is.filter.context.effect}{Whehter to filter outliers in context.effect}
\item{is.trim.cvg}{Whether to remove low coverage positions in context.effect}
\item{is.trim.pos}{Whether to remove context that contains too few positions}
\item{min.cvg}{Minimal coverage in context.effect.Positions whose coverage is smaller than min.cvg will be removed if is.trim.cvg is TRUE, and it will not be used if is.trim.cvg is FALSE}
\item{min.pos}{Minimal number of poisitions of each sequence context in context.effect. A sequence context will be removed if number of positions with that context is samller than min.pos. It will not be used if is.trim.pos is FALSE}
\item{method}{method could be 'hieModel' or 'CC'. Hierarchical model will be used if method='hieModel', and case-control method will be used if method='CC'}

}
\details{


}
\value{
The function returns a list. There are 4 elements in the list, which are \code{pos}, \code{neg}, \code{genome.start.pos} and \code{genome.start.neg}. 
\code{pos} and \code{neg} contain a list of objects for positive and negative strand respectively, each element of the list contains test statistic Z and estimated parameters of each position in a genome, and names of each object is name of reference genome(the string after > in the fasta file of the reference genome). The elements of each object are
\item{mu_native}{Sample mean of each position in native sample}
\item{sigma2_native}{Sample variance of each position in native sample}
\item{sampleSize_native}{Coverage of each position in native sample}

\item{mu_ctrl}{For case-control method, it is sample mean of each position in control sample.}
\item{sigma2_ctrl}{For case-control method, it is sample variance of each position in control sample.}
\item{sampleSize_ctrl}{Coverage of each position in control sample}

\item{theta}{Estimated hyperparameter of each position, \eqn{\theta}}
\item{kappa}{Estimated hyperparameter of each position, \eqn{\kappa}}
\item{upsilon}{Estimated hyperparameter of each position, \eqn{\upsilon}}
\item{tau2}{Estimated hyperparameter of each position, \eqn{\tau^2}}

\item{is_findContext}{A logical vector, indicates whether there are enough homologous positions are found for each position. It is TRUE if number of effective homologous position >=min.pos, otherwise it is false.(A homologous position is an effective homologous positions if and only if its coverage >= min.cvg)}

\item{t.stat}{test statistic Z. For case-controle method, Z = (mu_native - mu_ctrl)/(sqrt(sigma2_native/sampleSize_native + sigma2_ctrl/sampleSize_ctrl)), and for hierarchical model Z = (mu_native - mu_ctrl)/(sqrt(sigma2_native/sampleSize_native + sigma2_ctrl/(sampleSize_ctrl + kappa))).  }
\item{ipd.ratio}{IPD ratio of each position}

All of aboved objects (t.stat etc) are vectors, and the corresponding genome positions of their first element is in $genome.start.pos$ and $genome.start.neg$. (The leftmost position of genome is 1)

}
\references{

}

\author{Zhixing Feng \email{fzx08@mails.thu.edu.cn} }


\seealso{
\code{\link{getFeaturesAlongGenome}},
\code{\link{getContextEffectByPos}},
\code{\link{getGenomeSeq}}
}
\examples{
library(seqPatch)
### load data
load(paste(system.file(package = "seqPatch"), "/data/test.Rdata",sep = ""))

### load reference genome
genomeSeq <- getGenomeSeq(paste(system.file(package = "seqPatch"), "/data/refgenome.fasta",sep = ""))

### normalize 
alnsF.native <- normalizeByMovie(alnsF.native,alnsIdx.native)
alnsF.ctrl <- normalizeByMovie(alnsF.ctrl,alnsIdx.ctrl)

### get IPD for each position( logarithm transform IPD and filter outlirer automatically)
genomeF.native <- getFeaturesAlongGenome(alnsF.native, alnsIdx.native )
genomeF.ctrl <- getFeaturesAlongGenome(alnsF.ctrl, alnsIdx.ctrl)

### detect DNA modifications
detection <- detectModification(genomeF.native, genomeF.ctrl, genomeSeq, context.effect, method='hieModel')
}
\keyword{SMRT sequencing}
\keyword{polymerase kinetic}
\keyword{DNA modification}
\keyword{hierarchical model}
\keyword{Empirical Bayes}

