# run magma for scdrs based on harmonised GWAS catalog sumstats

# workflow to run magma
# URL: https://cncr.nl/research/magma/


rule magma_gene_annot:
    input:
        snp_loc = "resources/genotypes/{geno_set}/chr{chr}.bim",
        gene_loc = "resources/magma/geneanno.loc"
    output:
        "results/magma/{geno_set}/chr{chr}.genes.annot"
    params:
        gene_windows_bp=0,
        prefix_out = lambda x: f"results/magma/{x.geno_set}/chr{x.chr}"
    log: 
        "logs/magma/{geno_set}/chr{chr}.genes.annot.log"
    shell:
        """
        magma --annotate window={params.gene_windows_bp} \
            --snp-loc {input.snp_loc} \
            --gene-loc {input.gene_loc} \
            --out {params.prefix_out}
        """

# rule magma_mk_synonyms:
#     input:
#         snp_frq = "resources/genotypes/{geno_set}/chr{chr}.frq",
#         gene_loc = "resources/magma/geneanno.loc"
#     output:
#         "results/magma/{geno_set}/chr{chr}.genes.annot"
#     params:
#         gene_windows_bp=0,
#         prefix_out = lambda x: f"results/magma/{x.geno_set}/chr{x.chr}"
#     shell:
#         """
#         magma --annotate window={params.gene_windows_bp} \
#             --snp-loc {input.snp_loc} \
#             --gene-loc {input.gene_loc} \
#             --out {params.prefix_out}
#         """

def get_pheno_metadata(x):
    """
    Get the number of samples for a given pheno from metadata file
    """
    file = list(Path(f"resources/gwas/gwas_catalog/{x.pheno}/").glob("*.h.tsv.gz"))[0]
    file_meta = f'{str(file)}-meta.yaml'
    with open(file_meta, 'r') as f:
        meta = yaml.safe_load(f)
    n_sample = meta['samples'][0]['sample_size']
    return {'file': file, 'n_sample': int(n_sample)}

rule format_gwas_catalog:
    """
    Format GWAS catalog files for the pipeline
    """
    input:
        gwas_file = lambda x: get_pheno_metadata(x)['file'],
        frq_files = expand("resources/genotypes/{{geno_set}}/chr{chr}.frq", chr=range(1, 23))
    output:
        formatted_gwas = "resources/gwas/formatted/{geno_set}/{pheno}.formatted.tsv"
    log:
        "logs/gwas/{geno_set}/{pheno}.formatted.log"
    params:
        max_maf_diff = 0.2,
    conda: "renv"
    script: "snakescripts/format_gwas_catalog.R"

rule magma_calc_gene_assoc:
    input:
        gwas = rules.format_gwas_catalog.output.formatted_gwas,
        bfile = expand("resources/genotypes/{{geno_set}}/chr{{chr}}.{ext}", \
                        ext = ["bed", "bim", "fam"]),
        genes_annot = rules.magma_gene_annot.output
    output:
        raw = "results/magma/output/{pheno}/{geno_set}/chr{chr}.genes.raw",
        out = "results/magma/output/{pheno}/{geno_set}/chr{chr}.genes.out"
    log:
        "logs/magma/output/{pheno}/{geno_set}/chr{chr}.genes.log"
    params:
        prefix_out = lambda x, output: Path(str(output[0])).with_suffix('').with_suffix('').as_posix(),
        prefix_bfile = lambda x, input: Path(str(input.bfile[0])).with_suffix('').as_posix(),
        # n = lambda x: get_pheno_metadata(x)['n_sample']
    shell:
        """
        magma \
	        --bfile {params.prefix_bfile} \
	        --gene-annot {input.genes_annot} \
	        --pval {input.gwas} ncol=N \
	        --gene-model snp-wise=mean \
	        --out {params.prefix_out}
        """
    
# # aggregate results
rule magma_agg_gene_assoc:
    input:
        out = expand("results/magma/output/{{pheno}}/{{geno_set}}/chr{chr}.genes.out", \
               chr = range(1,23)),
        raw = expand("results/magma/output/{{pheno}}/{{geno_set}}/chr{chr}.genes.raw", \
               chr = range(1,23)),
    output:
        raw = "results/magma/aggregate/{geno_set}/{pheno}.genes.raw",
        out = "results/magma/aggregate/{geno_set}/{pheno}.genes.out"
    log:
        "logs/magma/aggregate/{geno_set}/{pheno}.genes.log"
    shell:
        """
        awk 'NR == FNR {{print; next}} FNR > 1' {input.out} > {output.out}
        awk 'NR == FNR {{print; next}} FNR > 2' {input.raw} > {output.raw}
		"""

# # update gene out with hgnc_symbols & fdr correction (with qvalue)
rule magma_format_output:
    input:
        out = rules.magma_agg_gene_assoc.output.out,
        gene_loc = "resources/magma/geneanno.loc"
    output:
        "results/magma/aggregate/{geno_set}/{pheno}.magma.tsv"
    conda: "renv"
    script:
        "snakescripts/magma_format_output.R"

# # run results for all trait
def get_trait_input(x):
    """
    Get the input for all traits
    """
    phenos = [f.with_suffix('').name for f in Path("resources/metadata/phenotypes").glob("*.yaml")]
    
    files = [f"results/magma/aggregate/{x.geno_set}/{p}.magma.tsv" for p in phenos]
    return files

rule magma_all_trait:
    input: get_trait_input
    output: "results/aggregate/{geno_set}.magma.parquet.gz"
    conda: "renv"
    log: "logs/aggregate/{geno_set}.magma.log"
    script: "snakescripts/aggregate/magma.R"