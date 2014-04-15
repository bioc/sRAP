\name{RNA.deg}
\alias{RNA.deg}

\title{
Differential Expression Statistics
}
\description{
Provides a table of differenitally expressed genes (in .xlsx format) as well
as differential expression statistics for all genes (in .xlsx format as well
as returned data frame).  Function automatically creates a heatmap for
differentially expressed genes and user can optionally also create box-plots
for each individual differentially expressed gene.  The efficacy of this
 protocol is described in [1].

Output files will be created in the "DEG" and "Raw_Data" subfolders.
}
\usage{
RNA.deg(sample.file, expression.table, project.name, project.folder, log2.fc.cutoff = 0.58, pvalue.cutoff = 0.05, fdr.cutoff = 0.05, box.plot = TRUE, ref.group = FALSE, ref = "none", method = "lm", color.palette = c("red", "blue", "green", "orange", "purple", "cyan", "pink", "maroon", "yellow", "grey", "black", colors()))
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
  \item{log2.fc.cutoff}{
	If the primary variable contains two groups with a specified reference, 
	this is the cut-off to define differentially expressed genes (default = 1.5,
	on a linear scale).  Otherwise, this variable is ignored
}
  \item{pvalue.cutoff}{
	Minimum p-value to define differentially expressed genes
}
  \item{fdr.cutoff}{
	Minimum false discovery rate (FDR) to define differentially expressed genes.
}
  \item{box.plot}{
	A logical value: Should box-plots be created for all differenitally expressed genes?
	If TRUE, then box-plots will be created in a separate subfolder.
}
  \item{ref.group}{
	A logical value: Is the primary variable 2 groups, with a reference group?
}
  \item{ref}{
	If the primary variable contains two groups (indicated by  ref.group = FALSE),
	this is the reference used to calculate fold-change values (so,
	the mean expression for the	reference group is substracted from
	the treatment group).
	
	Otherwise, this variable is ignored
}
  \item{method}{
	Method for calculating p-values:
	
	"lm" (Default) = linear regression
	"aov" = ANOVA
}
  \item{color.palette}{
	Colors for primary variable (specified in the second column of the sample file).
	If the primary variable is a continuous variable, this parameter is ignored.
}
}
\value{
	Data frame containing differential expression statistics.
	
	First column contains gene name.
	
	If the primary variable contains two groups (with a specified reference),
	then fold-change values are provided in the second column.
	
	P-values and FDR values are provided for each variable in subsequent columns,
	starting with the primary variable.
}
\references{
[1] Warden CD, Yuan Y-C, and Wu X. (2013).
Optimal Calculation of RNA-Seq Fold-Change Values.
\emph{Int J Comput Bioinfo In Silico Model},
2(6): 285-292
}
\author{
Charles Warden <cwarden45@gmail.com>
}

\examples{
	
library("sRAP")

dir <- system.file("extdata", package="sRAP")
expression.table <- file.path(dir,"MiSeq_cufflinks_genes_truncate.txt")
sample.table <- file.path(dir,"MiSeq_Sample_Description.txt")
project.folder <- getwd()
project.name <- "MiSeq"

expression.mat <- RNA.norm(expression.table, project.name, project.folder)

stat.table <- RNA.deg(sample.table, expression.mat, project.name, project.folder, box.plot=FALSE, ref.group=TRUE, ref="scramble",method="aov", color.palette=c("green","orange"))

}

\seealso{
sRAP goes through an entire analysis for an example dataset provided with the sRAP package.

Please post questions on the sRAP discussion group: http://sourceforge.net/p/bdfunc/discussion/srap/
}
