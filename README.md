# `GWAS_tools`, a pipeline for processing GWAS summary statistics into a uniform format

This is a fork of [`GWAS_tools`](https://github.com/GRealesM/GWAS_tools), restructuring what was a `bash` program into a `snakemake` pipeline depending on a series of R scripts.

All but one dependency is managed using the accompanying `docker` image. `snakemake` will pull this from Dockerhub for you, but for the sake of transparency I've also included the `Dockerfile` at `docker/Dockerfile` if you would like to examine its contents. I've not yet found a stable download URL for the `liftOver` executable and I think copying my own instance may violate some sort of copyright, so for now you'll have to put it on your own path ([download here](https://genome-store.ucsc.edu/)).

This pipeline also depends on the [`1kGP_pipeline`, also hosted here on GitHub](https://github.com/twillis209/1kGP_pipeline). This workflow downloads and processes the 1000 Genomes Project Phase 3 data and its use should be entirely transparent as `snakemake` will download it and make use of it through the `module` statement.

## Testing the pipeline

The `workflow/rules/test_rules.smk` file contains a target which can be used to run the pipeline on a rheumatoid arthritis data set from [Okada et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3944098/).

From the root directory of the repository, run:

```
snakemake test_pipeline
```

You'll want to pass your own profile using the `--profile` CLI argument, too.
