# lps_macrophages
Analyzing Affymetrix microarray data using R. The goal here is to determine the extent of macrophage polarization to the M1 subtype following treatment with a bacterial molecule called lipopolysaccharide (LPS) in cells from two genetically different types of mice. 

The "m1.R" file analyzes M1 macrophage marker expression from two different strains of mice (AKR/J and C3H/HeJ) following LPS treatment. The source of this microarray data is the Gene Expression Omnibus (GEO) repository (https://www.ncbi.nlm.nih.gov/geo/), from the GSE38705 series. The individual sample files from this series can be accessed using the accession numbers provided in the .R file.

Briefly, the following steps are done in the "m1.R" file:
1. Four samples of microarray data are acquired: LPS-treated AKR/J mice, control-treated AKR/J mice, LPS-treated C3H/HeJ mice, and control-treated C3H/HeJ mice. Each sample contains expression values for 22,416 probes.
2. NA values are removed from the data.
3. Differences in expression level between LPS vs. control is calculated for each strain.
4. The data is subsetted for analysis of just genes representing M1 macrophage markers.
5. The resulting M1 marker expression heat map from the "m1.R" file is shown below. These two strains show differences in gene expression following LPS treatment, thus providing an example of a gene-environment interaction.

<img src="images/m1_heat_map.png" width = 700>
