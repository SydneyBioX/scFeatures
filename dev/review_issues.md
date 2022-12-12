Hi @ycao6928,

Thank you for submitting **scFeatures** to Bioconductor.

I think the package needs a little more work, particularly to improve the documentation (man pages and vignettes), before it can be accepted.
It also feels a bit awkward for a Bioconductor package to not have first-class support for Bioconductor data structures (*SingleCellExperiment* and *SpatialExperiment*).
Yes, there's the `makeSeurat()` converter function but the main functionality requires a *Seurat* object as input.
It would good in future versions to allow a user to pass a _SingleCellExperiment_/_SpatialExperiment_ object as input (perhaps doing the conversion on the fly, if necessary, or, better yet, directly extracting the relevant data from the _SingleCellExperiment_/_SpatialExperiment_ object)

In my review below I have separated the issues into Required and Recommended points that I would ask you to address before the package can be accepted.
Would you please provide line-by-line comments to my initial review so that I know what changes I'm looking for in my re-review.

Cheers,
Pete

## Required

- [x] Most of the code in the vignettes is unevaluated and so is untested by `R CMD check`. Bioconductor expects that vignettes contain non-trivial evaluated code chunks and static vignettes are generally not acceptable.
  - When I tried running the vignette code line-by-line (like a new user would) there were several examples that returned an error. Such errors would be caught by `R CMD check` if the vignette code was evaluated. Some examples of where I encountered errors when running the vignette code:

```r
> feature_gene_prop_celltype  <- run_gene_prop_celltype(data_remove_mito, genes = genes_of_interest)
Error: BiocParallel errors
  1 remote errors, element index: 1
  1 unevaluated and other errors
  first remote error:
Error: Under current subsetting parameters, the default assay will be removed. Please adjust subsetting parameters or change default assay.
```
  - **Fixed**: The subsetted data included was causing issues.
    - TODO: Double check this is working properly.
```r
> feature_CCI <- run_CCI(data , species = \"Homo sapiens\" )
Error in run_CCI(data, species = \"Homo sapiens\") : 
  could not find function \"run_CCI\"
```

  - **Fixed**: use to include cell-cell interaction via CellChat, but it is not on CRAN/BioC, so it was removed.
  Now we implemented this function using the method SingleCellSignalR, which is also a top performing method
  as shown in `Evaluation of cell-cell interaction methods by integrating single-cell RNA sequencing data with spatial information`.
  and is available as a Bioconductor package. 

```r
genes_of_interest <- c(\"TIGIT\", \"PDCD1\")
feature_gene_prop_bulk <- run_gene_prop(data, genes = genes_of_interest )
Error: Under current subsetting parameters, the default assay will be removed. Please adjust subsetting parameters or change default assay.
```

```r
scfeatures_result <- scFeatures(data, 
             # input selective feature types to generate.      
             feature_types = c(\"proportion_raw\", \"pathway_gsva\", \"L_stats\", \"gene_mean_celltype\",    \"gene_prop_aggregated\") , 
             # by default assumes data is \"scrna\" (single-cell RNA-seq), this is now set to \"spatial_p\" (spatial proteomics).   
             type = \"spatial_p\" ,
             # by default assumes \"Homo sapiens\", this is now set to `Mus musculus`. 
             species  = \"Mus musculus\"  ,   
             # by default uses the 50 hallmark genes,  now set to user specified pathways.   
             geneset = list(\"pathway_a\" = c(\"SMA\" ,  \"Vimentin\"  ),  \"pathway_b\" = c(\"B7H3\" ,\"FoxP3\"   )) ,
             # by default uses top variable genes to generate the celltype specific gene expression feature category,now set to user defined genes. 
             celltype_genes = data.frame(celltype = c(\"Macrophages\"  , \"Mono/Neu\" ) , marker = c(\"CD3\" , \"p53\")) ,
             # by default uses top variable genes to generate the overall aggregated gene expression feature category, now set to user defined genes.    
             aggregated_genes =  c(\"CD3\" , \"p53\") , 
             # When passing as SingleCellExperiment or SpatialExperiment, by default we use the assay stored in \"logcount\" 
             assay = \"norm\",
            # By default we look for the sample info in \"sample\" column and the celltype info in \"celltype\" column
             sample = \"imageID\",
             celltype = \"cellType\",
      #  If users want to construct features from the spatial metrics category, by default we look for the \"x_cord\" and \"y_cord\" column 
             spatialCoords = c(\"x\", \"y\"),
             # by default uses single core
             ncores = 8)
Error: Cannot find 'cellType' in this Seurat object
```
 - **Fixed**: Code was not evaluated because of very long runtimes, but this has been changed and broken chunks have been fixed.

- [x] The man pages (documentation) is too terse and imprecise.
  - Each man page should include well-written:
    - Title, Description and Details (see https://r-pkgs.org/man.html#title-description-details)
    - Arguments (see https://r-pkgs.org/man.html#arguments)
    - Return value (see https://r-pkgs.org/man.html#return-value)
    - Examples (see https://r-pkgs.org/man.html#sec-man-examples)
  - To take one example, the documentation of `process_data()`:
    - The 'Description' just repeats the title, which doesn't explain what pre-processing actually is.
    - There's no real description of arguments (what sort of object is `data`, what sort of normalisation is applied when `normalise = TRUE`, etc.)
    - The 'Value' section is incorrect (the return value is actually a _Seurat_ object) and it should be explained what the returned object  contains.

  - You might take a look at popular Bioconductor packages such as **scater** or **edgeR** for inspiration for the documentation (although it needn't be that detailed).
  - Some of the required documentation detail is in the `scFeatures_detail` vignette, but it properly belongs in the man pages. The vignette should show users how to use the package as a whole and the man pages how to use individual functions.
  
  **Fixed**: All exported functions now have Title, Description, Arguments, Return value.  
  Non - exported functions now have description. 
    
- [x] The 'Value' section of many man pages is incorrect. E.g., `run_proportion_raw()` returns a _data.frame_ but is documented to return 'a _matrix_ of samples x features'. These details matter when it comes to documenting package code.
  
  **Fixed**: Corrected the return value type. 

- [x] The 'Introduction' of the each vignette, especially the first or main vignette (which I take to be `scFeatures_summary`) should read like an abstract; see https://contributions.bioconductor.org/docs.html#vignette-introduction
 **Fixed**: edited the introduction in `scFeatures_summary`

- [x] Please give the vignettes informative titles and filenames; see https://contributions.bioconductor.org/docs.html#vignettes
  - The titles and names seems sufficiently informative.
- [x] Please use [BiocStyle](https://www.bioconductor.org/packages/release/bioc/html/BiocStyle.html) for vignette formatting.
  - All vigettes now use `BiocStyle`.
- [x] Simply removing **ClassifyR** from `DESCRIPTION` is insufficient because the `scFeatures_detail` vignette still uses it to demonstrate the utility of **scFeatures**. Now that **ClassifyR** is passing builds again, please re-add **ClassifyR** to the `DESCRIPTION` if using it in the vignettes or elsewhere in **scFeatures**.
  - Re-added `ClassifyR`.
- [x] **RhpcBLASctl** should not be needed or used in a package. Why is it needed (and hidden) in the vignette?
  - It was added to fix a bug that is now forgotten, I removed all references.
- [x] Why is **inline** needed in the vignettes?
  - Same as above. It is no longer needed.
- [ ] The formatting in the `scFeatures_detail` vignette is a bit messy because it includes lots of verbose output (progress bars). Please revise.
  - Messages and warnings have been removed from the output.
- [x] In the `scFeatures_detail` vignette you are subsetting a vector outside its bounds; please fix:
  - This has been fixed by including different data with the package.

- [x] In the `scFeatures_associationstudy` vignette, please use a temporary directory via `tempdir()`, rather than the current working directory, to demonstrate report generation.
  - This vignette is not intended to be run as it simulates report generation and should be done with one's own results.
- [x] That said, the report generation fails when using `tempdir()` but should work.

```r
> output_folder <- tempdir()
> run_association_study_report(scfeatures_result, output_folder )
Error in abs_path(input) : The file 'output_report.Rmd' does not exist.
In addition: Warning message:
In normalizePath(path, winslash = winslash, mustWork = mustWork) :
  path[1]=\"output_report.Rmd\": No such file or directory
```
 **Fixed**: This is now fixed 


- [x] The **org.Hs.eg.db** package is an unspecified dependency used in the report generation demonstrated in `scFeatures_associationstudy`. Please include it as a dependency in the `DESCRIPTION`.
  - Added `org.Hs.eg.db` to suggests.
- [x] There are hardcoded paths in the report demonstrated in `scFeatures_associationstudy` (https://github.com/SydneyBioX/scFeatures/blob/80034c0ec37c206d63a08198e22dda351126b805/inst/extdata/output_report.Rmd#L228). This means it won't work on anyone else's computer.
  - Removed hard-coded path in favour of `system.file` call.
- [x] Please add a `BugReports` field to the `DESCRIPTION` (usually a link to the Issues page of the GitHub repo).
- [ ] Bioconductor requires documentation of `.rds`/`.Rdata` files in `inst/extdata` in an `inst/script/ directory`. See [data documentation](https://contributions.bioconductor.org/docs.html#doc-inst-script).
  - [ ] `inst/extdata` is usually used for 'raw' data, so these data might properly belong under `data/` rather than `inst/extdata/`; see https://contributions.bioconductor.org/data.html.
- [x] Please add a table of contents to each vignette.
- [x] All man pages should have runnable examples (see https://contributions.bioconductor.org/docs.html#examples)
Nick to check - is there are reason that we include /dontrun in some of the examples?
I saw this is stated as not recommended in the webpage. 
- [ ] What are the `dev` and `docs` folders? Please justify or remove them from the main branch of the git repo (Bioconductor requests that \"Any files or directories for other applications (Github Actions, devtools, etc.) should ideally be in a different branch and not submitted to the Bioconductor version of the package.\" (see https://contributions.bioconductor.org/general.html?q=unnec#undesirable-files)
  - [x] Remove `docs` folder.
  - [ ] Rremove `dev` folder.

## Recommended

- [ ] It is strongly recommended to add unit tests. I would suggest starting with the individual `scFeatures::run_*()` functions that underpin the wrapper `scFeatures::scFeatures()` function.
- [x] This may be personal preference, but the function outputs feel the wrong way around (with samples as rows and features as columns). This is the opposite of how 'rectangular' data are usually stored in Bioonductor, e.g., *SummarizedExperiment*, where rows are features (e.g., genes) and samples are columns. A side effect, is that this leads to very 'wide' objects that don't display very nicely when printed (at least that was my experience with the example data). I recommend at least documenting why you choose to return the results in this orientation.
  - This is one of the reasons why sfFeatures was created.' Rectangular' data are not suitable for statistical modelling, which requires standard sample by feature [design matrices](https://en.wikipedia.org/wiki/Design_matrix). The reason why the output is as such is that the purpose of the features is statisitcal modelling.
- [x] `makeSeurat()` also accepts Seurat objects (https://github.com/SydneyBioX/scFeatures/blob/709f075578bf01b1823dc39fe1d5617472c3f888/R/wrapper_run_scfeatures.R#L225-L237), but this isn't documented. Please document when a user would need this functionality.
  **Fixed**: Added in the function description. 

- [ ] Please try to cite relevant literature. E.g., in the `scFeatures_detail` vignette you write, \"the L values between the pairs of proteins are calculated using the L function defined in literature\" but no reference is given to the relevant literature.
- [x] In  the `scFeatures_detail` vignette is the advice, \"This can be obtained from performing cell type prediction using reference data, for example, using SCTransform from Seurat (see https://satijalab.org/seurat/articles/spatial_vignette.html)\". However, to my understanding, **SCTransform** is a normlization method, not a cell type prediction method, and the link points doesn't point to how to actually do the cell type prediction as best I can tell. Please clarify.
  **Fixed**: Update the link to https://satijalab.org/seurat/articles/spatial_vignette.html#integration-with-single-cell-data-1, which 
  corresponds to the section using SCTransform for "probabilistic transfer of annotations from a reference to a query set". 

- [ ] Please consider NOTES raised by `BiocCheck::BiocCheck()`. Code styling notes can be regarded as suggestions, but other points should be followed or reasons given for not following them. 
  
```
    * NOTE: Update R version dependency from 4.2.0 to 4.3.0.
    * NOTE: Consider adding the maintainer's ORCID iD in 'Authors@R' with
      'comment=c(ORCID=\"...\")'
    * NOTE: Avoid 1:...; use seq_len() or seq_along()
    * NOTE: Use accessors; don't access S4 class slots via '@' in
      examples/vignettes.
    * NOTE: Consider adding runnable examples to man pages that document exported
      objects.
    * NOTE: Consider adding unit tests. We strongly encourage them. See
      https://contributions.bioconductor.org/tests.html
```

- [ ] Consider adding a top-level README file.
- [ ] Recommend adding a NEWS file (see https://contributions.bioconductor.org/news.html)
- [x] Recommend adding a `inst/CITATION` file (see https://contributions.bioconductor.org/citation.html)
 **Fixed**: This is now added
 
- [ ] Consider adding a package-level man page (see https://contributions.bioconductor.org/docs.html#package-level-documentation).
- [ ] It's strongly recommended to avoid direct slot access with `@` or `slot()` of S4 objects. Instead, use accessor functions. That said, my understanding is that **Seurat** authors have not followed this advice, so it's unavoidable when interacting with **Seurat** objects, but please ensure you are doing this when interacting with Bioconductor objects.
- [ ] Why does `remove_mito()` remove more than just mitochondrial genes? I strongly recommend choosing a more precise name for the function or split the functionality up into `remove_mito()`, `remove_ribo()`, etc.


# Work Done

- Nick: Down-sample data to 550 samples in order to create a dataset <3MB
- Yue: Finished writing the documentation. All exported functions now have Title, Description, Arguments, Return value.  
Non - exported functions now have description. 
Fixed `scFeatures_associationstudy` vignette. 
Added back the cell cell communication features, now implemented using SingleCellSignalR instead of CellChat. 
Removed the real spatial transcriptomics and spatial proteomics data. 
For vignette, generate the spatial transcriptomics and spatial proteomics by adding fake x, y coordinates onto the single-cell RNA-seq data that's included in the package. 

