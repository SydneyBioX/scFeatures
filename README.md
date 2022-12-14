# scFeatures: Multi-view representations of single-cell and spatial data for disease outcome prediction

<img src="https://github.com/SydneyBioX/scFeatures/blob/master/inst/sticker.png?raw=true" align="right" width="200">

scFeatures is a tool that generates multi-view representations of single-cell and spatial data through the construction of a total of 17 feature types belonging to the following six categories. 

1. cell type proportions
2. cell type specific gene expressions
3. cell type specific pathway expressions
4. cell type specific cell-cell interaction (CCI) scores
5. overall aggregated gene expressions
6. spatial metrics

![Overview](https://github.com/SydneyBioX/scFeatures/blob/master/inst/overview.png?raw=true)

##  Installation 

The latest scFeatures can be installed using devtools: 

 ```
library(devtools)
devtools::install_github("SydneyBioX/scFeatures")
 ```
 
##  Vignettes

A number of vignettes are provided in the Vignettes folder of this repo. These can also be viewed here by clicking the link below.  
     
### Feaeture generation using scFeatures

* [Brief vignette](https://sydneybiox.github.io/scFeatures/articles/scFeatures_summary.html) - The brief vignette is a quick starting point and provides a checklist of the functions to apply on single-cell and spatial data.   
* [Detailed vignette with case studies](https://sydneybiox.github.io/scFeatures/articles/scFeatures_detail.html) - This provides an in-depth description for each feature types and explanation of the output. It contains a step-by-step case study to the application of scFeatures on i) scRNA-seq data, ii) spatial proteomics data and iii) spatial transcriptomics data.    


### Association study of the features with conditions

We provide a function `run_association_study_report(scfeatures_result, output_folder)` that automatically runs association study and generates an html output report. The report visualises the features most associated with conditions (eg, diseased vs non-diseased).   

*  [Association study vignette](https://sydneybiox.github.io/scFeatures/articles/scFeatures_associationstudy.html) - This explains the usage of the function.

## Reference

Cao, Y., Lin, Y., Patrick, E., Yang, P., & Yang, J. Y. H. (2022). scFeatures: multi-view representations of single-cell and spatial data for disease outcome prediction. In O. Vitek (Ed.), Bioinformatics (Vol. 38, Issue 20, pp. 4745–4753). Oxford University Press (OUP). https://doi.org/10.1093/bioinformatics/btac590 
