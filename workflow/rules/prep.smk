# Collection of rules for preparing data for the pipeline


# get phenotype metadata
def get_pheno_metadata(x):
    with open("resources/metadata/phenotypes/{x.pheno}.yaml", "r") as f:
        metadata = yaml.safe_load(f)
    return metadata["phenotypes"]

rule dl_gwas_catalog:
    """
    Download harmonised GWAS summary statistics from GWAS catalog
    """
    input: "resources/metadata/phenotypes/{pheno}.yaml"
    output: directory("resources/gwas/gwas_catalog/{pheno}")
    params: pheno=lambda wildcards: wildcards.pheno

