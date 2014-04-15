\name{RNA.qc}
\alias{RNA.qc}

\title{
Gene Expression Quality Control Statistics
}
\description{
Provides descriptive statistics (median, top/bottom quartiles, mininum,maximum),
sample histograms and box-plot, sample dendrogram, principal component analysis plot.

Output files will be created in the "QC" subfolder.
}
\usage{
RNA.qc(sample.file, expression.table, project.name, project.folder, plot.legend = TRUE, color.palette = c("red", "blue", "green", "orange", "purple", "cyan", "pink", "maroon", "yellow", "grey", "black", colors()))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sample.file}{
	Tab-delimited text file providing group attributions for all samples considered for analysis.
}
  \item{expression.table}{
	Data frame with genes in columns and samples in rows.  Data should be log2 transformed.
	The RNA.norm function automatically creates this file.
}
  \item{project.name}{
	Name for sRAP project.  This determines the names for output files.
}
  \item{project.folder}{
	Folder for sRAP output files
}
  \item{plot.legend}{
	A logical value: Should legend be plotted within QC figures?
}
  \item{color.palette}{
	Colors for primary variable (specified in the second column of the sample file).
	If the primary variable is a continuous variable, this parameter is ignored.
}
}

\author{
Charles Warden <cwarden45@gmail.com>
}

\examples{
	
library("sRAP")
library("WriteXLS")

dir <- system.file("extdata", package="sRAP")
expression.table <- file.path(dir,"MiSeq_cufflinks_genes_truncate.txt")
sample.table <- file.path(dir,"MiSeq_Sample_Description.txt")
project.folder <- getwd()
project.name <- "MiSeq"

expression.mat <- RNA.norm(expression.table, project.name, project.folder)

RNA.qc(sample.table, expression.mat, project.name, project.folder, plot.legend=FALSE, color.palette=c("green","orange"))
}

\seealso{
sRAP goes through an entire analysis for an example dataset provided with the sRAP package.

Please post questions on the sRAP discussion group: http://sourceforge.net/p/bdfunc/discussion/srap/
}

