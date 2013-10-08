\name{RNA.bdfunc.signal}
\alias{RNA.bdfunc.signal}

\title{
Functional Enrichment for a Table of Normalized Gene Expression Values
}
\description{
Bi-Directional FUNCtional enrichment [1] compares expression values for up- and down-regulated genes are compared for at least
one gene set, using normalized expression values.  Gene sets are already 
defined for human and mouse gene symbols.  All other gene sets must be specified by the user.
The user can optionally output box-plots to visualze enrichment scores across samples
in different groups.
}
\usage{
RNA.bdfunc.signal(expression.table, sample.file, project.name, project.folder, species = NULL, enrichment.file = NULL, p.method = "t-test", p.adjust.method = "fdr", plot.flag = TRUE, color.palette = c("red", "blue", "green", "orange", "purple", "cyan", "pink", "maroon", "yellow", "grey", "black", colors()))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expression.table}{
	Data frame with genes in columns and samples in rows.  Data should be log2 transformed.
	The RNA.norm function automatically creates this file.
}
  \item{sample.file}{
	Tab-delimited text file providing group attributions for all samples considered for analysis.
}
  \item{project.name}{
	Name for sRAP project.  This determines the names for output files.
}
  \item{project.folder}{
	Folder for sRAP output files
}
  \item{species}{
	Name for species used for analysis.  If species is set to "human" or "mouse,"
	then pre-defined gene lists provided by the sRAP package are used.

	The default human gene list is created from gene ontology [2] and MSigDB [3]
	databases.  The default mouse gene list is created from the gene ontology [2] database.
}
  \item{enrichment.file}{
	Table of gene lists including up- and down-regulated genes.
	
	This is only necessary when defining a custom species.  This parameter is
	ignored when the species is set to "human" or "mouse".
}
  \item{p.method}{
	Method for calculating p-values
	
	"t-test" (Default) = t-test between up-regulated and down-regualted genes
	"mann-whitney" = Non-parametric Mann-Whitney U test between up-and down-regulated genes
	"ks" = Kolmogorov-Smirnov test between up- and down-regulated genes
}
  \item{p.adjust.method}{
	Method for calculating false discovery rate (FDR):
	
	"fdr" (Default)= B-H "Step-Up" FDR [4]
	"q-value" = Storey q-value [5]
	"none" = use unadjusted p-value without multiple hypothesis correction
}
  \item{plot.flag}{
	Logical value: Should box-plots be created for all gene sets?
	
	If primary variable is two groups called "positive" and "negative",
	this value also determines if ROC plot will be created.
}
  \item{color.palette}{
	Colors for primary variable (specified in the second column of the sample file).
	
	If method is set to "t-test," this variable is ignored.  In this special case,
	groups with an average t-test statistic above 2 are colored red, groups with
	an average t-test statistic below -2 are colored green, and all other groups
	are colored grey.
}
}
\references{

[1] Warden CD, Kanaya N, Chen S, and Yuan Y-C. (2013).
BD-Func: A Streamlined Algorithm for Predicting Activation and Inhibition of Pathways.
\emph{peerJ},
1:e159

[2] Ashburner M, Ball CA, Blake JA, Botstein D, Butler H, Cherry JM,
 Davis AP, Dolinski K, Dwight SS, Eppig JT, Harris MA, Hill DP, Issel-Tarver L,
 Kasarskis A, Lewis S, Matese JC, Richardson JE, Ringwald M, Rubin GM,
 and Sherlock G.(2000).
Gene Ontology: tool for the unification of biology
\emph{Nat Genet},
25:25-29

[3] Liberzon A, Subramanian A, Pinchback R, Thorvaldsdottir H, Tamayo P, and Mesirov JP.(2011).
Molecular signatures database (MSigDB) 3.0.
\emph{Bioinformatics},
27:1739-1740.

[4] Benjamini Y, and Hochberg Y.(1995).
Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing.
\emph{Journal of the Royal Statistical Society Series B},
57:289-300.

[5] Storey JD, and Tibshirani R. (2003).
Statistical significance for genomewide studies.
\emph{Proceedings of the National Academy of Sciences},
100:9440-9445.

}
\author{
Charles Warden <cwarden@coh.org>
}

\examples{
	
library("sRAP")
library("WriteXLS")
library("qvalue")
library("ROCR")

dir <- system.file("extdata", package="sRAP")
expression.table <- file.path(dir,"MiSeq_cufflinks_genes_truncate.txt")
sample.table <- file.path(dir,"MiSeq_Sample_Description.txt")
project.folder <- getwd()
project.name <- "MiSeq"

expression.mat <- RNA.norm(expression.table, project.name, project.folder)

stat.table <- RNA.deg(sample.table, expression.mat, project.name, project.folder, box.plot=FALSE, ref.group=TRUE, ref="scramble",method="aov", color.palette=c("green","orange"))

RNA.bdfunc.signal(expression.mat, sample.table, plot.flag=FALSE, project.name, project.folder, species="human")
}

\seealso{
sRAP goes through an entire analysis for an example dataset provided with the sRAP package.
}
