Introduction 
------------

Uveal melanoma (UM) is the most common intraocular cancer afflicting 7
in 1,000,000 individuals in the U.S and UK. 90% of uveal melanoma cases
have mutations in the GNAQ/11 genes. These genes encode for G-protein
coupled receptor (GPCR) subunits found within the inner leaflet of the
plasma membrane that act as an on/off switch for a wide range of
cellular processes. Mutations in the GNAQ/GNA11 genes lead to
constitutively active GPCRs, which causes increased cell growth and
proliferation. In the paper referenced below the authors demonstrate
that FR900359 (FR) prevents GDP/GTP exchange in GRPCs; a crucial step in
signal transduction, via small molecule allosteric inhibition.
Inactivation of the primary mutated pathway found in most UM cases leads
to arrest of cell proliferation, inhibition of secondary signaling, and
reinstated melanocyte differentiation leading to potential therapeutic
uses.

In this report I plan to reproduce the MA plot in Fig.6 B from Onken et
al. 2018 with RNA expression data published with the paper. This figure
depicts the fold change vs. average log expression of RNA between FR and
non FR treated cells. Identification of gene clusters suppressed by FR,
or not, helps us determine whether FR will be a useful treatment in UM.
Once I've recreated the MA plot depicted below I plan on using my own
RNAseq data collected from cells treated with a B-Raf inhibitor. B-Raf
is a proto-oncogenic serine/threonine kinase that is commonly expressed
in FR resistant cells. Expression profiles of FR treated cells, and
B-Raf inhibitor treated cells, allows us to compare their individual
efficacy while also illuminating a potential multidrug treatment.

Reference
---------

*Onken MD, Makepeace CM, Kaltenbronn KM, et al. Targeting nucleotide
exchange to inhibit constitutively active G protein α subunits in cancer
cells. Science Signaling. 2018;11(546):eaao6852.
doi:[[10.1126/scisignal.aao6852]{.underline}](https://doi.org/10.1126/scisignal.aao6852)*

![](media/image1.png){width="3.9739588801399823in"
height="3.7801071741032373in"}

Fig. 6. FR represses expression of differentiation genes by restoring
function of the PRC2. (A)Gaq-mutant 92.1 UM cells were treated with FR
or vehicle, and RNA was collected 1 and 3 days (1d and 3d, respectively)
after treatment for RNAseq analysis. Results of a multidimensional gene
expression analysis that compares the relative patterns of expression of
all genes across all samples and groups genes with similar patterns. The
graph shows samples positioned by their relative gene expression values
within each pattern. Dimension 1 (x axis; the most

represented pattern) shows separation based on vehicle treatment (red
balls) versus FR treatment (blue balls), whereas dimension 2 (y axis;
the second most represented pattern) shows separation based on time in
culture (indicated by 1d or 3d on balls). (B)MA plot (M, log ratio; A,
mean average) comparing gene expression between FR- and vehicle-treated
92.1 samples identifies a group of significantly reduced genes (circled;
fold change, \>2; FDR, q \<0.01) associated with FR treatment. (C) GO
analysis of the FR-repressed gene set \[circled in (B) with arrow\]. (D)
FR-repressed genes \[circled in (B) with arrow\] identified as targets
of the polycomb repressive complex2(PRC2)byGSEA.EGF,epidermalgrowth
factor; BMP2, bone morphogenetic protein; hESC, human embryonic stem
cells. (E) Effect of the EZH1/2 inhibitor GSK503 on morphological
differentiation elicited by FR. Representative fields are shown from one
of three experiments of 92.1 UM cells treated for 7 days with GSK503 and
for 3 days with FR and then imaged by phase-contrast microscopy. Scale
bar, 100 μm. (F) Effect of GSK503 on pigmentation of FR-treated cells,
visualized by macroscopic inspection. 92.1 cells were treated for 7 days
with GSK503 and for 3 days with FR and pelleted; representative images
from one of three experiments. (G) PRC2 inhibition by GSK503.
Immunoblots of 92.1 cells treated for 7 days with GSK503 show reduced
histone H3K27 trimethylation. Plot shows relative fraction of
trimethyl-histone H3K27 compared to dimethyl sulfoxide (DMSO) control
and normalized to total histone H3 from densitometry data from three
independent experiments. \*P \< 0.01 by t test; significance was
confirmed using q \< 0.01 by the FDR method of Benjamini and Hochberg.

Materials and Methods 
---------------------

To collect the data seen in Figure 6 B the authors grew cells in 100 nM
FR or DMSO in RPMI growth medium. After 1 and 3 days of treatment RNA
was collected from each treatment. HiSeq2500 was used to generate FastQ
raw data. This data was aligned to the whole genome using Bioconductor
in EdgeR. Once aligned, a 2.0 fold change in expression was considered
significant comparing FR treated to non-FR treated cells. This data was
then plotted as an MA plot, M = log ratio and A = mean average.

Since I have a CSV file containing the log ratio and mean averages, I'll
make the volcano plot using Rstudio. Once that is completed I will work
my way back through the pipeline described in their methods section.
I'll start with aligning my own data set using Bioconductor, then I'll
use SAM (significance analysis of microarrays) to rank gene expression.
