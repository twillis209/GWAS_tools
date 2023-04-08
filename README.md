# GWAS_tools

This is a fork of [`GWAS_tools`](github.com/GRealesM/GWAS_tools), restructuring what was a bash program into a `snakemake` pipeline depending on a series of R scripts. 99% of the credit for this belongs to the original author; translating things into `snakemake` and `R` took only a few hours.

## Dependencies

The `environment.yaml` file can be used to create a `conda` environment called `gwas_tools` in which `snakemake` is available:

```
conda env create --file environment.yaml
```

The pipeline also expects the `liftOver` executable to be on the path ([download here](https://genome-store.ucsc.edu/)).

Lastly, the pipeline's R scripts also require the following packages:
- `data.table`
- `stringr`
- `magrittr`

## Testing the pipeline

The `workflow/rules/test_rules.smk` file contains a target which can be used to run the pipeline on a rheumatoid arthritis data set from [Okada et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3944098/). 

From the root directory of the repository, run:

```
snakemake --profile . test_pipeline
```
