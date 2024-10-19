rule download_test_dataset:
    output:
        ensure("resources/gwas_pipeline/ra.tsv.gz", sha256 = "c4e32a29251af1ecd86436664ead74220516f3f6dd39d1d9ee31296eda6df7bc")
    params:
        url = "http://plaza.umin.ac.jp/~yokada/datasource/files/GWASMetaResults/RA_GWASmeta_European_v2.txt.gz"
    shell: "wget -O {output} {params.url}"

rule test_pipeline:
    input:
        "results/gwas_pipeline/reheaded/fixed_alleles/min_col_set/recalc_sumstats/liftover/realigned_alleles/ra.tsv.gz"
