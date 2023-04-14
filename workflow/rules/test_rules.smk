rule download_test_dataset:
    output:
        "resources/gwas_pipeline/ra.tsv.gz"
    run:
        shell("wget -O {output} http://plaza.umin.ac.jp/~yokada/datasource/files/GWASMetaResults/RA_GWASmeta_European_v2.txt.gz")

        md5sum_actual = shell("md5sum {output}", read = True).split(' ')[0]

        if md5sum_actual != "ca5fe1dd8400a6bd620d3a76d8d9bf18":
            raise Exception("md5sums do not match")

rule test_pipeline:
    input:
        "results/gwas_pipeline/reheaded/fixed_alleles/min_col_set/recalc_sumstats/liftover/ra-hg38.tsv.gz"
