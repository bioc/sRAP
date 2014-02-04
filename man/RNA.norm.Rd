\name{RNA.norm}
\alias{RNA.norm}

\title{
Normalization for RNA-Seq Data
}
\description{
Takes a table of RPKM (Read Per Kilobase per Million reads [1]) gene expression values.
Rounds RPKM values based upon RPKM.cutoff(to avoid bias from low-coverage genes),
and then performs a log2 transformation of the data (so that the data more closely
follows a normal distribution).  The efficacy of this protocol is described in [2].

Output files will be created in the "Raw_Data" subfolder.
}
\usage{
RNA.norm(input.file, project.name, project.folder, RPKM.cutoff = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input.file}{
	Table of RPKM expression values.  Genes are represented in columns.
	Samples are represented in rows.
}
  \item{project.name}{
	Name for sRAP project.  This determines the names for output files.
}
  \item{project.folder}{
	Folder for sRAP output files
}
  \item{RPKM.cutoff}{
	Cut-off for rounding RKPM expression values.  If the default of 0.1 is used,
	genes with expression values consistently below 0.1 will essentially be ignored.
}
}
\value{
	Data frame of normalized expression values on a log2 scale.
	
	Just like the input table, genes are represented on columns,
	samples are represented in rows.
	
	This data frame is used for quality control and differential expression analysis.
}
\references{

[1] Mortazavi A, Williams BA, McCue K, Schaeffer L, and Wold B. (2008).
Mapping and quantifying mammalian transcriptomes by RNA-Seq.
\emph{Nat Meth},
5:621-628.

[2] Warden CD, Yuan Y-C, and Wu X. (2013).
Optimal Calculation of RNA-Seq Fold-Change Values.
\emph{in preparation},
}
\author{
Charles Warden <cwarden@coh.org>
}

\examples{
	
library("sRAP")

dir <- system.file("extdata", package="sRAP")
expression.table <- file.path(dir,"MiSeq_cufflinks_genes_truncate.txt")
sample.table <- file.path(dir,"MiSeq_Sample_Description.txt")
project.folder <- getwd()
project.name <- "MiSeq"

expression.mat <- RNA.norm(expression.table, project.name, project.folder)
}

\seealso{
sRAP goes through an entire analysis for an example dataset provided with the sRAP package.

Please post questions on the sRAP discussion group: http://sourceforge.net/p/bdfunc/discussion/srap/
}