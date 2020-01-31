![searching under the
lamppost](https://source.unsplash.com/y41FEMqdJ3A/400x300)

Why this project?
-----------------

Like the proverbial man [searching for his lost keys under the lamp
post](https://www.matthewhirschey.com/articles/exploratory-mind) because
the light shines there, searching for biological truths often occurs
under ‘lamp posts’ because that’s where scientists can see. But what if
your keys are not under the light? Or your gene is totally unknown? What
do you do?

The scientific method has guided scientific minds for hundreds of years,
starting with a question, followed by a hypothesis, and then an
experimental path to test the prediction. While hypotheses are the
bedrock of science, the volume, complexity, and sophistication of modern
science necessitate new methods to generate hypotheses.

New tools in Data Science – a combination of computer programming, math
& statistics, and topical expertise – combined with the rapid adoption
of open science and data sharing allow scientists to access publicly
available datasets and interrogate these data *before* performing any
experiments.

Imagine having strong data to support your new hypothesis *before*
testing it. Welcome to data-driven hypothesis.

What is this project?
---------------------

The overall goal of the data-driven hypothesis (DDH) project is to use
new tools in Data Science to generate hypotheses supported by data that
can be tested in the lab.

Several high-quality, functional genomic datasets are published online
and made available with Creative Commons Attribution 4.0 International
[(CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/) licenses.
Functional genomics is a field of molecular biology that aims to
understand the function of all genes and proteins in a genome – a stated
goal of much basic science research. In functional genomics,
experimental strategies generally involve high-throughput, genome-wide
approaches rather than a more traditional “gene-by-gene” approach. The
advent and rapid adoption of data-sharing platforms, such as
[figshare.com](https://figshare.com) have provided high-quality data
sets for public interrogation. The DDH project aims to integrate
functional genomics data and holds tremendous promise to generate
hypotheses, data, and knowledge in order to provide a deep understanding
of the dynamic properties of an organism.

How does it work?
-----------------

This project began as an extension of a simple and common concept in
molecular biology called gene co-expression analysis. When a gene of
unknown function is identified, one strategy to learn something about
the new gene is to identify shared patterns of expression with other
genes. If unknown Gene X is expressed with known genes A, B, and C, then
you can infer that Gene X might be part of a functional module with A,
B, C. This approach is particularly powerful when genes A, B, and C are
part of a known biological pathway, which leads to the hypothesis that
Gene X might also be part of that same pathway.

![](methods_files/figure-markdown_strict/gene_coexpression-1.png)

Following on this idea, we set out to map genes to common functional
pathways based on dependence of a pathway for cell viability. Project
Achilles is a systematic effort by the [Broad
Institute](https://www.broadinstitute.org) as part of a larger [‘DepMap’
project](http://www.depmap.org) aimed at identifying and cataloging gene
essentiality across hundreds of well-characterized cancer cell lines
using highly standardized pooled genome-scale loss-of-function screens.
This project uses lentiviral-based pooled RNAi or CRISPR/Cas9 libraries
to systematically knock-out each gene in the genome, which allows for
the stable suppression/ablation of each gene individually in a subset of
cells within a pooled format allowing for genome wide interrogation of
gene essentiality. Using computational modeling, a normalized value of
gene essentiality is given for each gene in a single cell line. A lower
score means that a gene is more likely to be essential in a given cell
line. A score of -1 corresponds to the median of all common essential
genes, whereas a score of 0 is equivalent to a gene that is not
essential; a positive score indicates a gain in fitness and often
identifies tumor suppressor genes.

It is well-known that human cancer cell lines rely on different pathways
for their viability. Indeed this is the entire rationale for
personalized, precision medicine in cancer. The overall goal of the
‘DepMap’ project is to identify all essential genes in 2000 cell lines
over the 5-year project period to identify new therapeutic targets in
various cancers. Despite not knowing the mechanistic basis for why some
cell lines require specific genes while other cell lines do not, we
reasoned that intrinsic reliance of a cell on a pathway might allow
unbiased detection of novel genes participating in specific pathways.

![](methods_files/figure-markdown_strict/cell_dependency-1.png)

What did I do?
--------------

Essential gene data from Project Achilles were downloaded from the
DepMap portal at: [depmap.org](https://depmap.org/portal/download/). The
19Q4 release contains gene essentiality scores for 17440 genes across
689 cell lines, and was used for this project.

![](methods_files/figure-markdown_strict/dep_scores-1.png)

#### Patterns

To find patterns in gene dependencies across cell lines, we generated a
Pearson correlation matrix of all genes by all genes. This analysis
generated gene-gene correlation values that matched values published on
[depmap.org](https://depmap.org), validating the first step in our
analysis. High levels of gene expression are often thought to be
indicative of key genes for a given cell type. Thus, we next compared
dependency values to gene expression values. The [Cancer Cell Line
Encyclopedia](https://portals.broadinstitute.org/ccle/about) project is
a collaboration between the Broad Institute, and the Novartis Institutes
for Biomedical Research and its Genomics Institute of the Novartis
Research Foundation, which together conduct detailed genetic and
pharmacologic characterization of a large panel of human cancer models.
As of the CCLE 2019 release, 1249 cell lines have been characterized for
gene expression. In the 19Q4 DepMap release, 684 of the 689 cell lines
have gene expression data. Using these two datasets, we compared the
essentiality of a gene to its expression value.

![](methods_files/figure-markdown_strict/depVexp-1.png)

We predicted a V-shaped curve, with stronger dependencies as gene
expression increases. Surprisingly, we saw no relationship between gene
expression and gene essentiality, where genes with both low and high
expression displayed both gains and losses in fitness. The overall
observation from this dataset shows baseline gene expression levels are
poor indicators of the essentiality of a gene. This analysis also
highlighted that several genes were binned on the x-axis, i.e. could
have no measurable expression levels, but have assigned dependency
scores. Across 689 cell lines in the Achilles project, 16.6% of all gene
expression values are zero, confirming this notion.

#### Noise Reduction

Given cells do not express all genes, but might receive a dependency
sore in this experimental paradigm, we sought to remove dependency
scores for gene-cell line pairs that have an expression value of zero
under basal conditions. Of the 689 cell lines for which gene
essentiality data is collected, 684 have genome-wide gene expression
data. From these cell lines, we removed dependency scores for genes from
cell line that have a corresponding gene expression value of zero.

![](methods_files/figure-markdown_strict/expression_0-1.png)

For some genes expressed in highly specific and restricted cell types,
this operation removed many dependency values. After removing these
values, we found that highly specialized genes in discrete cell types
have too few cells with both gene expression values and gene
essentiality values to assign a meaningful correlation value. Thus, if a
gene was absent from too many cell lines, we omitted it to prevent
assigned values from relying on too few data points.

![](methods_files/figure-markdown_strict/dep_cleaning-1.png)

We set a threshold of no more than 589 zeros, meaning that if a gene had
fewer than 100 cell lines with dependency values, the correlation
pattern of a gene would be meaningless, and that gene was therefore
removed. This process removed 0 genes that had too few cells with
expression and dependency data. These ‘cleaned’ dependency data had
17440 remaining gene-dependency pairs, which were then used to generate
correlation matrix.

How does it work?
-----------------

To identify genes that shared similar patterns of essentiality with
other genes, thereby placing genes in functional pathways, we generated
a Pearson correlation matrix on these prioritized data in order to
quantify the similarity in dependency patterns and annotate genes in
functional pathways.

![](methods_files/figure-markdown_strict/r2-1.png)

This process generated approximately 304 million correlation values,
with a distribution centered around zero. This output produced a range
of maximum correlation values for each gene.

![](methods_files/figure-markdown_strict/achilles_max-1.png)

#### Statistics

Rather than setting an arbitrary threshold for the r^2 value that would
be considered a low, medium, or high correlation between two genes, we
performed a permutation test on the correlated data. A permutation test
involves permuting one or more variables in a data set before performing
the test, in order to break any existing relationships and simulate the
null hypothesis. In this case, we broke the relationship between
gene-gene pairs and the correlation values. We can then compare the true
statistic (mean correlation) to the generated distribution of null
statistics (fake means), along with standard deviations of these sampled
data. This strategy will give a better idea of where to draw a threshold
of a “significant correlation” for these analyses. We sampled 20,000 r^2
values from all gene-gene pairs without replacement simulating a virtual
Achilles dataset for a single cell. We then repeated this process 1000
times mimicking 1000 discrete cell lines.

This statistical analysis produced the following data:  
**Mean:** 0.0032181  
**Standard Deviation:** 0.0684781

Using a standard deviation threshold of 3, we calculated the boundaries
of r^2 values to be greater than 0.21 or lower than -0.2 for negative
correlations.

#### Pathway Analyses

To identify clusters of genes with shared relationships, we performed
gene set enrichment analysis. Enrichment analysis is a computational
method for inferring knowledge about a target gene set by comparing it
to annotated gene sets representing prior biological knowledge.
Enrichment analysis determines whether an input set of genes
significantly overlaps with annotated gene sets. For each gene in our
matrix, we determined the number of genes that were greater than or less
than 3 standard deviations away from the permuted mean. This target gene
list was then queried against 108 gene sets across a [broad range of
curated data](https://amp.pharm.mssm.edu/Enrichr/#stats). By leveraging
the [Enrichr](https://amp.pharm.mssm.edu/Enrichr/) resource from the
[Ma’ayan Laboratory](http://labs.icahn.mssm.edu/maayanlab/), we
determined the top ranked pathways, processes, drugs, cell lines,
tissues, or diseases, and ranked by p-value. In this setting, the
p-value is computed using a standard statistical method used by most
enrichment analysis tools: Fisher’s exact test or the hypergeometric
test. This is a binomial proportion test that assumes a binomial
distribution and independence for probability of any gene belonging to
any set. [See here for more information about how Enrichr computes its
associations](https://amp.pharm.mssm.edu/Enrichr/help#background).

How do I use this?
------------------

In order to identify the functional annotation of a single gene, begin
with a query. Entering a single gene in the search box produces a series
of tables and plots that identifies a functional map of the processes
that gene *might* be involved in. In some cases, a querying gene with
known functions will identify gene with well-established connections to
the query gene; in other cases, new genes and new biological process
might be identified, suggesting there is more to discover for well-known
pathways. Querying unknown genes is especially powerful, as the
associated genes and pathways provide a starting point for an otherwise
difficult problem to prioritize experimentally.

#### 1. Query YFG (your favorite gene)

As an example, we will query the protein P53. Typing in the given
protein name “P53” produces an error, because the official gene symbol
needs to be entered. Querying TP53 (official gene symbol) generates a
short summary of the gene, its name and list of aliases, and an Entrez
gene summary paragraph when available.

#### Summary

**Gene**: TP53  
**Name**: Tumor protein p53  
**aka**: p53, LFS1  
**Entrez ID**: 7157

This gene encodes a tumor suppressor protein containing transcriptional
activation, DNA binding, and oligomerization domains. The encoded
protein responds to diverse cellular stresses to regulate expression of
target genes, thereby inducing cell cycle arrest, apoptosis, senescence,
DNA repair, or changes in metabolism. Mutations in this gene are
associated with a variety of human cancers, including hereditary cancers
such as Li-Fraumeni syndrome. Alternative splicing of this gene and the
use of alternate promoters result in multiple transcript variants and
isoforms. Additional isoforms have also been shown to result from the
use of alternate translation initiation codons from identical transcript
variants (PMIDs: 12032546, 20937277). \[provided by RefSeq, Dec 2016\]

#### 2. Dependencies

The first plot shows the distribution of dependency scores across 689
cell lines ranked from lowest (strongest dependencies) to highest (no
dependency or inverse). Each of the 689 cell lines is represented by a
single point on the plot. Generally, values below -1 indicate the gene
of interest (TP53 in this example) is essential in that cell line;
values between -1 and 0, mean that cells lose fitness, but the gene is
not essential; values hovering around zero indicate that ablation of
TP53 has little effect on cell growth; values above 1, indicate that
knocking-out the gene leads to a fitness advantage. In the case of TP53,
several cells have a fitness advantage in its absence, consistent with
its role as a tumor suppressor.

The second plot is a histogram of dependency scores, showing the
distribution of scores for TP53. While the majority of cells have little
change in cellular growth when TP53 is absent (the histogram is centered
around zero), some cells require TP53 for growth (cells scoring below
-1), whereas in other cells TP53 functions as a tumor suppressor (cells
with a score above 1).

![](methods_files/figure-markdown_strict/dep_plots-1.png)

To identify the cells at the tails of the plots, a dependency table will
show the ranked cells by dependency score, with cell lineage information
appended. In some cases, specific cell types or lineages will show
consistent patterns of dependency on a gene.

#### Cells with strong TP53 genetic dependencies:

<table style="width:69%;">
<colgroup>
<col style="width: 16%" />
<col style="width: 34%" />
<col style="width: 18%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Cell Line</th>
<th style="text-align: left;">Lineage</th>
<th style="text-align: right;">Dependency Score</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">NMCG1</td>
<td style="text-align: left;">central_nervous_system</td>
<td style="text-align: right;">4.6</td>
</tr>
<tr class="even">
<td style="text-align: left;">DKMG</td>
<td style="text-align: left;">central_nervous_system</td>
<td style="text-align: right;">4.89</td>
</tr>
<tr class="odd">
<td style="text-align: left;">TUHR4TKB</td>
<td style="text-align: left;">kidney</td>
<td style="text-align: right;">3.65</td>
</tr>
<tr class="even">
<td style="text-align: left;">KMRC1</td>
<td style="text-align: left;">kidney</td>
<td style="text-align: right;">4.13</td>
</tr>
<tr class="odd">
<td style="text-align: left;">TTC642</td>
<td style="text-align: left;">soft_tissue</td>
<td style="text-align: right;">4.03</td>
</tr>
<tr class="even">
<td style="text-align: left;">BIN67</td>
<td style="text-align: left;">ovary</td>
<td style="text-align: right;">3.55</td>
</tr>
</tbody>
</table>

#### Cells with low or inverse TP53 genetic dependencies

<table style="width:56%;">
<colgroup>
<col style="width: 16%" />
<col style="width: 19%" />
<col style="width: 19%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Cell Line</th>
<th style="text-align: left;">Lineage</th>
<th style="text-align: right;">Dependency Score</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">BL70</td>
<td style="text-align: left;">lymphocyte</td>
<td style="text-align: right;">-1.36</td>
</tr>
<tr class="even">
<td style="text-align: left;">HCC15</td>
<td style="text-align: left;">lung</td>
<td style="text-align: right;">-1.29</td>
</tr>
<tr class="odd">
<td style="text-align: left;">KASUMI1</td>
<td style="text-align: left;">blood</td>
<td style="text-align: right;">-1.06</td>
</tr>
<tr class="even">
<td style="text-align: left;">HCC1143</td>
<td style="text-align: left;">breast</td>
<td style="text-align: right;">-0.93</td>
</tr>
<tr class="odd">
<td style="text-align: left;">CME1</td>
<td style="text-align: left;">soft_tissue</td>
<td style="text-align: right;">-0.89</td>
</tr>
<tr class="even">
<td style="text-align: left;">KMS34</td>
<td style="text-align: left;">plasma_cell</td>
<td style="text-align: right;">-0.89</td>
</tr>
</tbody>
</table>

Understanding the shape of the curve and distribution of the raw data
underlying the patterns is important for interpreting the results.

#### 3. Similar Patterns

Positive correlations of dependency scores are ranked for each gene.
Recall that these genes show similar patterns of dependencies in the
same cell lines. More simply, the cells that care about TP53 deletion
also care about deletion of these genes, implying a functional
relationship. In the Dependency Score Example heatmap schematic above,
TP53 is gene X, and genes with similar patterns would be genes A, B, and
C. The 117 genes that show a similar genetic dependencies as TP53 and
are above 3 standard deviations away from the resampled mean are
displayed.

<table style="width:78%;">
<colgroup>
<col style="width: 13%" />
<col style="width: 54%" />
<col style="width: 9%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Gene</th>
<th style="text-align: left;">Gene Name</th>
<th style="text-align: right;">R^2</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">TP53BP1</td>
<td style="text-align: left;">Tumor protein p53 binding protein 1</td>
<td style="text-align: right;">0.73</td>
</tr>
<tr class="even">
<td style="text-align: left;">CDKN1A</td>
<td style="text-align: left;">Cyclin dependent kinase inhibitor 1a</td>
<td style="text-align: right;">0.69</td>
</tr>
<tr class="odd">
<td style="text-align: left;">USP28</td>
<td style="text-align: left;">Ubiquitin specific peptidase 28</td>
<td style="text-align: right;">0.65</td>
</tr>
<tr class="even">
<td style="text-align: left;">CHEK2</td>
<td style="text-align: left;">Checkpoint kinase 2</td>
<td style="text-align: right;">0.64</td>
</tr>
<tr class="odd">
<td style="text-align: left;">ATM</td>
<td style="text-align: left;">Atm serine/threonine kinase</td>
<td style="text-align: right;">0.61</td>
</tr>
<tr class="even">
<td style="text-align: left;">XPO7</td>
<td style="text-align: left;">Exportin 7</td>
<td style="text-align: right;">0.44</td>
</tr>
</tbody>
</table>

These 117 genes were queried for gene set enrichment, and the gene sets
and pathways with the strongest statistical significance are shown.
Simply stated, these are the pathways that best represent the list of
genes that share similar genetic dependencies, and suggest that the
query gene is part of these pathways.

<table style="width:72%;">
<caption>Table continues below</caption>
<colgroup>
<col style="width: 36%" />
<col style="width: 36%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Gene Set</th>
<th style="text-align: left;">Gene List</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">KEA 2015</td>
<td style="text-align: left;">ATM</td>
</tr>
<tr class="even">
<td style="text-align: left;">Jensen DISEASES</td>
<td style="text-align: left;">Ataxia telangiectasia</td>
</tr>
<tr class="odd">
<td style="text-align: left;">WikiPathways 2019 Human</td>
<td style="text-align: left;">Integrated Cancer Pathway WP1971</td>
</tr>
<tr class="even">
<td style="text-align: left;">WikiPathways 2019 Human</td>
<td style="text-align: left;">ATM Signaling Pathway WP2516</td>
</tr>
<tr class="odd">
<td style="text-align: left;">SubCell BarCode</td>
<td style="text-align: left;">HCC827 Nuclear N4 30609389</td>
</tr>
<tr class="even">
<td style="text-align: left;">Jensen DISEASES</td>
<td style="text-align: left;">Nijmegen breakage syndrome</td>
</tr>
<tr class="odd">
<td style="text-align: left;">BioCarta 2016</td>
<td style="text-align: left;">RB Tumor Suppressor/Checkpoint Signaling in response to DNA damage Homo sapiens h rbPathway</td>
</tr>
<tr class="even">
<td style="text-align: left;">KEA 2015</td>
<td style="text-align: left;">MAPK8</td>
</tr>
<tr class="odd">
<td style="text-align: left;">SysMyo Muscle Gene Sets</td>
<td style="text-align: left;">MGS440 Human left ventricle idiopathic cardiomyopathy v healthy GSE1145 down</td>
</tr>
<tr class="even">
<td style="text-align: left;">SILAC Phosphoproteomics</td>
<td style="text-align: left;">down 5nM dasatinib vs ctrl K562 (Human) [19651622]</td>
</tr>
</tbody>
</table>

<table style="width:14%;">
<colgroup>
<col style="width: 13%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: right;">Overlap</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: right;">11/161</td>
</tr>
<tr class="even">
<td style="text-align: right;">6/30</td>
</tr>
<tr class="odd">
<td style="text-align: right;">6/44</td>
</tr>
<tr class="even">
<td style="text-align: right;">6/40</td>
</tr>
<tr class="odd">
<td style="text-align: right;">20/981</td>
</tr>
<tr class="even">
<td style="text-align: right;">5/22</td>
</tr>
<tr class="odd">
<td style="text-align: right;">4/13</td>
</tr>
<tr class="even">
<td style="text-align: right;">10/225</td>
</tr>
<tr class="odd">
<td style="text-align: right;">11/239</td>
</tr>
<tr class="even">
<td style="text-align: right;">11/352</td>
</tr>
</tbody>
</table>

#### 4. Dissimilar

Like the analysis for genes that share similar patterns, this analysis
can be used to find genes that share distinctly dissimilar patterns;
that is, genes that have an inverse correlation of dependences. Simply
stated, the cells that care about TP53 deletion *do not* care about
deletion of these genes, implying an inverse or opposing relationship.
In the Dependency Score Example heatmap schematic above, TP53 is gene X,
and genes with dissimilar patterns would be genes D, E, and F. The 181
genes that show inverse genetic dependencies to TP53 and are below 3
standard deviations away from the resampled mean are:

<table style="width:86%;">
<colgroup>
<col style="width: 11%" />
<col style="width: 63%" />
<col style="width: 11%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Gene</th>
<th style="text-align: left;">Gene Name</th>
<th style="text-align: right;">R^2</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">MDM2</td>
<td style="text-align: left;">Mdm2 proto-oncogene</td>
<td style="text-align: right;">-0.7</td>
</tr>
<tr class="even">
<td style="text-align: left;">PPM1D</td>
<td style="text-align: left;">Protein phosphatase, mg2+/mn2+ dependent 1d</td>
<td style="text-align: right;">-0.59</td>
</tr>
<tr class="odd">
<td style="text-align: left;">MDM4</td>
<td style="text-align: left;">Mdm4 regulator of p53</td>
<td style="text-align: right;">-0.47</td>
</tr>
<tr class="even">
<td style="text-align: left;">USP7</td>
<td style="text-align: left;">Ubiquitin specific peptidase 7</td>
<td style="text-align: right;">-0.45</td>
</tr>
<tr class="odd">
<td style="text-align: left;">PPM1G</td>
<td style="text-align: left;">Protein phosphatase, mg2+/mn2+ dependent 1g</td>
<td style="text-align: right;">-0.43</td>
</tr>
<tr class="even">
<td style="text-align: left;">TERF1</td>
<td style="text-align: left;">Telomeric repeat binding factor 1</td>
<td style="text-align: right;">-0.38</td>
</tr>
</tbody>
</table>

These 181 genes were also queried for gene set enrichment, and the gene
sets and pathways with the strongest statistical significance are shown.
Simply stated, these are the pathways that best represent the list of
genes that have inverse genetic dependencies.

<table>
<colgroup>
<col style="width: 23%" />
<col style="width: 63%" />
<col style="width: 12%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Gene Set</th>
<th style="text-align: left;">Gene List</th>
<th style="text-align: right;">Overlap</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">TF Perturbations Followed by Expression</td>
<td style="text-align: left;">ADAR OE HEK293 HUMAN GSE87198 RNASEQ UP</td>
<td style="text-align: right;">20/276</td>
</tr>
<tr class="even">
<td style="text-align: left;">Jensen TISSUES</td>
<td style="text-align: left;">Cervical carcinoma cell</td>
<td style="text-align: right;">88/4876</td>
</tr>
<tr class="odd">
<td style="text-align: left;">GTEx Tissue Sample Gene Expression Profiles up</td>
<td style="text-align: left;">GTEX-S95S-0002-SM-3NM8K blood male 60-69 years</td>
<td style="text-align: right;">77/3958</td>
</tr>
<tr class="even">
<td style="text-align: left;">GTEx Tissue Sample Gene Expression Profiles up</td>
<td style="text-align: left;">GTEX-XUYS-0002-SM-47JXL blood male 50-59 years</td>
<td style="text-align: right;">75/3854</td>
</tr>
<tr class="odd">
<td style="text-align: left;">GTEx Tissue Sample Gene Expression Profiles up</td>
<td style="text-align: left;">GTEX-UPK5-0003-SM-3NMDI blood male 40-49 years</td>
<td style="text-align: right;">77/4003</td>
</tr>
<tr class="even">
<td style="text-align: left;">GTEx Tissue Sample Gene Expression Profiles up</td>
<td style="text-align: left;">GTEX-T6MN-0002-SM-3NMAH blood male 50-59 years</td>
<td style="text-align: right;">78/4038</td>
</tr>
<tr class="odd">
<td style="text-align: left;">GTEx Tissue Sample Gene Expression Profiles up</td>
<td style="text-align: left;">GTEX-UPIC-0002-SM-3NMDC blood female 20-29 years</td>
<td style="text-align: right;">75/3954</td>
</tr>
<tr class="even">
<td style="text-align: left;">GTEx Tissue Sample Gene Expression Profiles up</td>
<td style="text-align: left;">GTEX-U3ZH-0002-SM-3NMDD blood male 30-39 years</td>
<td style="text-align: right;">75/4012</td>
</tr>
<tr class="odd">
<td style="text-align: left;">GTEx Tissue Sample Gene Expression Profiles up</td>
<td style="text-align: left;">GTEX-X638-0003-SM-47JZ1 blood female 70-79 years</td>
<td style="text-align: right;">73/3894</td>
</tr>
<tr class="even">
<td style="text-align: left;">GTEx Tissue Sample Gene Expression Profiles up</td>
<td style="text-align: left;">GTEX-SIU7-0001-SM-3NMAW blood male 50-59 years</td>
<td style="text-align: right;">71/3711</td>
</tr>
<tr class="odd">
<td style="text-align: left;">GTEx Tissue Sample Gene Expression Profiles up</td>
<td style="text-align: left;">GTEX-UPJH-0001-SM-3NMDE blood male 50-59 years</td>
<td style="text-align: right;">73/3891</td>
</tr>
<tr class="even">
<td style="text-align: left;">huMAP</td>
<td style="text-align: left;">RPL5</td>
<td style="text-align: right;">14/137</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Enrichr Submissions TF-Gene Coocurrence</td>
<td style="text-align: left;">RAD51</td>
<td style="text-align: right;">19/299</td>
</tr>
<tr class="even">
<td style="text-align: left;">RNA-Seq Disease Gene and Drug Signatures from GEO</td>
<td style="text-align: left;">PRMT9 HeLa knockdown GSE63953 down</td>
<td style="text-align: right;">24/498</td>
</tr>
<tr class="odd">
<td style="text-align: left;">GTEx Tissue Sample Gene Expression Profiles up</td>
<td style="text-align: left;">GTEX-V1D1-0003-SM-3NMDP blood male 50-59 years</td>
<td style="text-align: right;">69/3617</td>
</tr>
</tbody>
</table>

How to interpret these genes and pathways is more variable than the
positively correlated genes and pathways. In some cases, a negative
regulator of a gene has a negative correlation with that gene, such as
in this example with TP53 and its negative regulator MDM2. In other
cases, opposing *pathways* are shown, contrasting the TP53 enriched
pathway term “MAPK8” with the dissimilar enriched pathway term “RPL5”,
revealing two opposing biological pathways.

#### 5. Graph

Identifying genes that share similar patterns of dependency to a queried
unknown gene generates strong hypotheses about new functional
annotations and maps to new pathways. However, the strength of the
hypothesis cannot be fully inferred from single gene list. If a new gene
is associated with the queried gene, then you might infer a functional
relationship. However, if you inspect the top 10 genes with the queried
gene, then inspect the top 10 genes of each of those, building a
functional network graph of the top related genes might reveal a
stronger association of the new gene with your queried gene *and* its
top ranked genes.

![](methods_files/figure-markdown_strict/static_graph-1.png)

The data presented on datadrivenhypothesis.org is an interactive image
that can be zoomed, dragged, and manipulated for data exploration; the
static image above is a representative snapshot.

Where do I get more information?
--------------------------------

Code is available on the Hirschey Lab github account, including links to
download the raw data, and run the analyses in R from scratch.
Furthermore, the Broad Institute has a lot of information on their
[website](http://www.broadinstitute.org) and the website dedicated to
the [Dependency Map project](http://www.depmap.org) about how the raw
data were generated, and provide a list of references.

#### Code Availability

[Generate
Data](https://github.com/hirscheylab/ddh/blob/master/code/generate_depmap_data.R)  
[Statistical
Analyses](https://github.com/hirscheylab/ddh/blob/master/code/generate_depmap_stats.R)  
[Table
Generator](https://github.com/hirscheylab/ddh/blob/master/code/generate_depmap_tables.R)  
[Pathway
Generator](https://github.com/hirscheylab/ddh/blob/master/code/generate_depmap_pathways.R)

Methods updated January 31, 2020
