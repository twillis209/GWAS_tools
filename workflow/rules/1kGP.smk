import pandas as pd

def get_hg19_vcf_sha256(w):
    daf = pd.read_csv(f"resources/1kG/hg19/sha256.tsv", sep = '\t', names = ['sha256', 'File'])

    return daf[daf.File.str.match(f'{w.chr}.vcf.gz')].sha256.values[0]

def get_hg38_vcf_sha256(w):
    daf = pd.read_csv(f"resources/1kG/hg38/sha256.tsv", sep = '\t', names = ['sha256', 'File'])

    return daf[daf.File.str.match(f'{w.chr}.vcf.gz')].sha256.values[0]

rule download_hg19_variant_reference:
    output:
        "resources/genome_reference/hg19.tsv.gz"
    group: "1kG"
    shell:
        """
        wget -O {output} ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz
        """

rule download_hg19_reference_sequence:
    output:
        ensure("resources/genome_reference/hg19.fa.zst", sha256 = "32f65df649ae46813bad00fee998542c7fd121aa9d01659e950ac307f2502693")
    params:
        compressed = "resources/genome_reference/hg19.fa.gz",
        uncompressed = "resources/genome_reference/hg19.fa"
    group: "1kG"
    shell:
        """
        wget -O {params.compressed} ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
        gunzip {params.compressed}
        zstd --rm {params.uncompressed}
        """

rule download_hg38_reference_sequence:
    output:
        ensure("resources/genome_reference/hg38.fa.zst", sha256 = "abf5670bf3de11a1e39176bc6596fb0597071e2be33ecae269855e502e1cdc53")
    resources:
        runtime = 30
    group: "1kG"
    shell:
        """
        wget -O {output} https://www.dropbox.com/s/xyggouv3tnamh0j/GRCh38_full_analysis_set_plus_decoy_hla.fa.zst?dl=1
        """

rule download_1kG_hg38_manifest:
    output:
        temp("resources/1kG/hg38/manifest.tsv")
    group: "1kG"
    shell:
        "wget -O {output} http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/phased-manifest_July2021.tsv"

rule process_1kG_hg38_manifest:
    input:
        "resources/1kG/hg38/manifest.tsv"
    output:
        "results/1kG/hg38/processed_manifest.tsv"
    group: "1kG"
    run:
        daf = pd.read_csv(input[0], sep = '\t', names = ['File', 'Byte', 'Checksum'])

        daf = daf[daf.File.str.match('.+\.vcf\.gz$')]

        daf = daf.assign(Chr=daf.File.str.extract('.+chr(\w+)\.filtered.+'))

        daf.to_csv(output[0], sep = '\t', index = False)

rule download_1kG_hg38_genotype_data:
    output:
        protected(ensure("resources/1kG/hg38/{chr}.vcf.gz", sha256 = get_hg38_vcf_sha256))
    resources:
        runtime = 60
    retries: 3
    group: "1kG"
    run:
        if wildcards.chr == 'chrX':
            shell("wget -O resources/1kG/hg38/chrX.vcf.gz http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz")
        else:
            shell("wget -O resources/1kG/hg38/{wildcards.chr}.vcf.gz http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_{wildcards.chr}.filtered.shapeit2-duohmm-phased.vcf.gz")

rule download_1kG_hg38_sample_metadata:
     output:
        "resources/1kG/hg38/ped.txt"
     group: "1kG"
     shell:
         "wget -O resources/1kG/hg38/ped.txt http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt"

rule download_1kG_hg19_genotype_data:
    output:
        protected(ensure("resources/1kG/hg19/{chr}.vcf.gz", sha256 = get_hg19_vcf_sha256))
    resources:
        runtime = 15
    group: "1kG"
    run:
        if wildcards.chr == 'chrX':
            shell("wget -O {output} http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz")
        elif wildcards.chr == 'chrY':
            shell("wget -O {output} http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz")
        else:
            shell("wget -O {output} http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.{wildcards.chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz")

rule vcf_to_pgen:
    input:
        "resources/1kG/{assembly}/{chr}.vcf.gz"
    output:
        temp(multiext("results/1kG/{assembly}/{chr}", ".pgen", ".pvar.zst", ".psam"))
    params:
        out = "results/1kG/{assembly}/{chr}",
        id_format = "@:#:\$r:\$a",
        max_allele_len = 20
    threads: 8
    resources:
        mem_mb=get_mem_mb
    group: "1kG"
    shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --vcf {input} --make-pgen vzs --out {params.out} --set-all-var-ids {params.id_format} --max-alleles 2 --new-id-max-allele-len {params.max_allele_len} truncate"

rule pgen_to_hap_and_legend:
    input:
        multiext("results/1kG/{assembly}/{chr}", ".pgen", ".pvar.zst", ".psam")
    output:
        temp(multiext("results/1kG/{assembly}/{chr}", ".haps", ".legend", ".sample"))
    params:
        stem = "results/1kG/{assembly}/{chr}",
    threads: 8
    resources:
        mem_mb=get_mem_mb
    group: "1kG"
    shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --pfile {params.stem} vzs --export hapslegend --out {params.stem}"

rule concatenate_legend_files:
    input:
        expand("results/1kG/{{assembly}}/{chr}.legend", chr = [f"chr{x}" for x in range(1,23)]+["chrX"])
    output:
        temp("results/1kG/{assembly}/combined.legend.gz")
    params:
        uncompressed_output = "results/1kG/{assembly}/combined.legend"
    resources:
        runtime = 20
    group: "1kG"
    shell:
        """
        for x in {input}; do
            tail -n +2 $x >>{params.uncompressed_output}
        done

        gzip {params.uncompressed_output}
        """
