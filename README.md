# `GWAS_tools`, a pipeline for processing GWAS summary statistics into a uniform format

This is a fork of [`GWAS_tools`](https://github.com/GRealesM/GWAS_tools), restructuring what was a `bash` program into a `snakemake` pipeline depending on a series of R scripts.

All but one dependency is managed using the accompanying `docker` image. `snakemake` will pull this from Dockerhub for you, but for the sake of transparency I've also included the `Dockerfile` at `docker/Dockerfile` if you would like to examine its contents. I've not yet found a stable download URL for the `liftOver` executable and I think copying my own instance may violate some sort of copyright, so for now you'll have to put it on your own path ([download here](https://genome-store.ucsc.edu/)).

This pipeline also depends on the [`1kGP_pipeline`, also hosted here on GitHub](https://github.com/twillis209/1kGP_pipeline). This workflow downloads and processes the 1000 Genomes Project Phase 3 data and its use should be entirely transparent as `snakemake` will download it and make use of it through the `module` statement.

## A version without containers

At present I'm not able to run this pipeline in my current HPC environment due to some issues with `apptainer` permissions, so have had to roll back the containerisation in the dedicated `conda` branch (`master` contains the `docker` version). This environment does not container `liftOver` and all dependencies are listed in `workflow/envs/global.yaml`.

## Importing this workflow in another `snakemake` workflow using `module` (and the role of the config)

I wrote this workflow so it could be plugged into others with use of the `module` statement. I do this like so in the [`igad_paper_pipeline` worfklow](https://github.com/twillis209/igad_paper_pipeline):

```
module gwas_pipeline:
    snakefile: github("twillis209/GWAS_tools", path = "workflow/Snakefile", branch = "conda")
    config: config["GWAS_tools"]
```

Overwriting the workflow's config with the `config` statement is necessary to customise the `GWAS_tools` workflow's handling of the data sets in your importing workflow, typically to tell it which column contains the `ref` allele and which the `alt` allele, e.g.
```
 scz:
    ref: "^A1$"
    alt: "^A2$"
```
You can also list `extraneous_columns` which you would like `GWAS_tools` to drop and `snps_to_retain` which do not appear in the `1kG` reference panel used to prune SNPs in the later part of the pipeline. See [the configfile for `igad_paper_pipeline`](https://github.com/twillis209/igad_paper_pipeline/blob/master/config/config.yaml) for an example of customising the workflow's behaviour in this way.

## Testing the pipeline

The `workflow/rules/test_rules.smk` file contains a target which can be used to run the pipeline on a rheumatoid arthritis data set from [Okada et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3944098/).

From the root directory of the repository, run:

```
snakemake test_pipeline
```

You'll want to pass your own profile using the `--profile` CLI argument, too.
