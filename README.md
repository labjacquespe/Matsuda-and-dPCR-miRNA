
To use the scripts, you will need a "data/" directory including dPCR data. 


### Process dPCR data 

```bash
Rscript process_dPCR_data.R
```

### Process sequencing data  

```bash
Rscript process_seq_data.R
```

### Run all associations and produce forestplots and boxplots.

```bash
Rscript associations.R <outdir>
```

### Replication
To reproduce results from the paper using the supplementary table:

```bash
Rscript test_associations_from_supp.R <outdir>
```
