# Collection of rules for preparing data for the pipeline


# get phenotype metadata
def get_pheno_metadata(x):
    import re
    with open(f"resources/metadata/phenotypes/{x.pheno}.yaml", "r") as f:
        metadata = yaml.safe_load(f)
    gcst = metadata["gcst_id"]
    match = re.search(r"GCST(\d{8})", gcst)
    if match:
        num = int(match.group(1))
        start = f"GCST{num - (num % 1000) + 1:08d}"
        end = f"GCST{num - (num % 1000) + 1000:08d}"
    else:
        raise ValueError("Invalid format: Expected a string starting with 'GCST' followed by 8 digits.")

    url = f"ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/{start}-{end}/{gcst}/harmonised"
    return url

rule dl_gwas_catalog:
    """
    Download harmonised GWAS summary statistics from GWAS catalog
    """
    input: "resources/metadata/phenotypes/{pheno}.yaml"
    output: directory("resources/gwas/gwas_catalog/{pheno}")
    params: url = get_pheno_metadata
    shell:
        """
        mkdir -p {output}
        wget -P {output} -r {params.url} --no-parent --no-directories
        """

rule prep_bioheart_geno:
    """
    Prepare BioHeart genotype data for the pipeline
    """
    input:
        bfile = expand("resources/genotypes/tenk10k_phase1/chr{{chr}}.{ext}", ext=["bed", "bim", "fam"]),
        extract = "resources/misc/genotype_extract/bioheart975.txt"
    output: expand("resources/genotypes/bioheart975/chr{{chr}}.{ext}", ext=["bed", "bim", "fam"])
    params:
        input_dir = "resources/genotypes/tenk10k_phase1",
        output_dir = "resources/genotypes/bioheart975",
        min_maf = 0.01
    shell:
        """
        plink \
            --bfile {params.input_dir}/chr{wildcards.chr} \
            --keep <(awk -vOFS='\\t' '{{print '0', $0}}' {input.extract}) \
            --make-bed \
            --maf {params.min_maf} \
            --out {params.output_dir}/chr{wildcards.chr}
        """

rule prep_freq_geno:
    input:
        expand("resources/genotypes/{{study}}/chr{{chr}}.{ext}", ext=["bed", "bim", "fam"])
    output:
        "resources/genotypes/{study}/chr{chr}.frq"
    params:
        prefix = "resources/genotypes/{study}/chr{chr}"
    shell:
        """
        plink \
            --bfile {params.prefix} \
            --freq \
            --out {params.prefix}
        """

