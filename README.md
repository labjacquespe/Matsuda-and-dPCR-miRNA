
```bash
Rscript  ~/def_link/extract_variables.R data/query.txt data/traits.tsv data/Gen3G_id_list.txt

Rscript explore_phenotype.R 
Rscript run_analysis.R quantification/alldata.tsv main_results_non_missing_scale T0
Rscript mirseq_analysis.R data/mirseq_plasma.cpm.tsv main_results_non_missing_scale/mirnaseq_cpm

Rscript forest_plot.R 

```