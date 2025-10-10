# Rules to run scDRS

wildcard_constraints:
    # exclude forward slash
    geno_set="[^/]+",
    phenotype="[^/]+"

rule magma_to_zscore:
    """
    Create zscore file for scDRS from collection of magma file
    format: https://martinjzhang.github.io/scDRS/file_format.html#pval-file-zscore-file
    """
    input: rules.magma_all_trait.output
    output: "resources/scdrs/zscore/{geno_set}.zscore.tsv"
    conda: "renv"
    log: "logs/scdrs/{geno_set}.magma_to_zscore.log"
    script: "snakescripts/scdrs/magma_to_zscore.R"

rule scdrs_munge_gs:
    """
    Create .gs file from .zscore file for scDRS
    """
    input: "resources/scdrs/zscore/{geno_set}.zscore.tsv"
    output: "resources/scdrs/gs/{geno_set}.gs"
    conda: "scverse"
    log: "logs/scdrs/{geno_set}.scdrs_munge_gs.log"
    params:
        fdr = 0.05,
        n_max = 1000
    shell:
        """
        scdrs munge-gs \
            --out-file {output} \
            --zscore-file {input} \
            --weight zscore \
            --n-max {params.n_max} \
            --fdr {params.fdr}
        """

checkpoint chunk_gs:
    """
    Chunk .gs file into smaller files (per phenotype) for scDRS
    """
    input:
        gs = "resources/scdrs/gs/{geno_set}.gs",
        config = "resources/scdrs/config/{geno_set}.yaml"
    output:
        directory("resources/scdrs/gs_chunked/{geno_set}")
    conda: "scverse"
    log: "logs/scdrs/{geno_set}.chunk_gs.log"
    shell:
        """
        mkdir -p {output}
        awk -v DIR={output} \
            'NR == 1 {{HEAD = $0; next}}
             {{TRAIT = $1;
              print HEAD RS $0 > DIR "/" TRAIT ".gs"}}' \
            {input.gs}
        """

rule scdrs_prep_h5ad_cov:
    """
    Prepare .h5ad file for scDRS
    """
    input:
        h5ad = "resources/scdrs/h5ad/{geno_set}.h5ad",
        sample_covar = "resources/scdrs/sample_covar/{geno_set}.sample_covar.csv",
        config = "resources/scdrs/config/{geno_set}.yaml",
        cell_meta = "resources/scdrs/cell_meta/{geno_set}.csv"
    output:
        prep_h5ad = "resources/scdrs/h5ad/{geno_set}.prep.h5ad",
        cov = "resources/scdrs/cov/{geno_set}.cov.tsv"
    conda: "scverse"
    resources:
        ncpus = 8,
        mem =  "200G"
    log: "logs/scdrs/{geno_set}.prep_h5ad_cov.log"
    script: "snakescripts/scdrs/prep/{wildcards.geno_set}.py"


rule scdrs_regress_h5ad_cov:
    """
    Prepare regressed .h5ad file for scDRS
    """
    input:
        prep_h5ad = "resources/scdrs/h5ad/{geno_set}.prep.h5ad",
        # gs = "resources/scdrs/gs/{geno_set}.gs",
        cov = "resources/scdrs/cov/{geno_set}.cov.tsv",
        config = "resources/scdrs/config/{geno_set}.yaml"
    output:
        # save as pickle as cannot save with standard .h5ad
        "resources/scdrs/h5ad/{geno_set}.reg.h5ad.pkl"
    conda: "scverse"
    resources:
        ncpus = 8,
        mem =  "200G"
    log: "logs/scdrs/{geno_set}.regress_h5ad_cov.log"
    script: "snakescripts/scdrs/regress_h5ad_cov.py"

rule scdrs_compute_score_api:
    """
    Compute scDRS score using the python API per phenotype
    """
    input:
        reg_pkl = "resources/scdrs/h5ad/{geno_set}.reg.h5ad.pkl",
        gs = "resources/scdrs/gs_chunked/{geno_set}/{phenotype}.gs",
        config = "resources/scdrs/config/{geno_set}.yaml"
    output:
        cell = "results/scdrs/cell_score/{geno_set}/{phenotype}.cell_score.tsv.parquet.gz",
        celltype = "results/scdrs/cell_type_stats/{geno_set}/{phenotype}.cell_type_stats.tsv"
    conda: "scverse"
    resources:
        time = "24:00:00",
        ncpus = 8,
        mem =  "200G"
    log: "logs/scdrs/{geno_set}-{phenotype}.compute_score.log"
    script: "snakescripts/scdrs/compute_score.py"

# rule scdrs_compute_group_stats_api:
#     """
#     Compute scDRS cell type stats using the python API per phenotype
#     """
#     input:
#         reg_pkl = "resources/scdrs/h5ad/{geno_set}.reg.h5ad.pkl",
#         gs = "resources/scdrs/gs_chunked/{geno_set}/{phenotype}.gs",
#         config = "resources/scdrs/config/{geno_set}.yaml",
#         cell = "results/scdrs/cell_score/{geno_set}/{phenotype}.cell_score.tsv.parquet.gz"
#     output:
#         celltype = "results/scdrs/cell_type_stats/{geno_set}/{phenotype}.cell_type_stats.tsv"
#     conda: "scverse"
#     resources:
#         ncpus = 8,
#         mem =  "8G"
#     log: "logs/scdrs/{geno_set}-{phenotype}.compute_group_stats.log"
#     script: "snakescripts/scdrs/compute_group_stats.py"

rule scdrs_get_top_score:
    """
    Calculate counts of cells annotate with specific cell types
    within the top nth quantile based on the normalised scDRS score
    """
    input:
        reg_pkl = "resources/scdrs/h5ad/{geno_set}.reg.h5ad.pkl",
        cell_score = "results/scdrs/cell_score/{geno_set}/{phenotype}.cell_score.tsv.parquet.gz",
        config = "resources/scdrs/config/{geno_set}.yaml"
    output: "results/scdrs/cell_type_top/{geno_set}/{phenotype}.cell_type_top.tsv"
    conda: "scverse"
    params:
        percentiles = [90, 95,99]
    resources:
        time = "24:00:00",
        ncpus = 8,
        mem =  "200G"
    log: "logs/scdrs/{geno_set}-{phenotype}.get_top_score.log"
    script: "snakescripts/scdrs/get_top_score.py"


# aggregate scdrs per phenotype
def scdrs_aggregate_stat(x):
    """
    Aggregate scDRS cell type stats per phenotype
    """
    dir_gs = Path(str(checkpoints.chunk_gs.get(**x).output[0]))
    gs_files = dir_gs.glob("*.gs")
    traits = [str(gs_file.stem) for gs_file in gs_files]

    celltype_stats_files = [f"results/scdrs/cell_type_stats/{x.geno_set}/{p}.cell_type_stats.tsv" for p in traits]
    celltype_top_files = [f"results/scdrs/cell_type_top/{x.geno_set}/{p}.cell_type_top.tsv" for p in traits]
    return {'cell_type_stats':celltype_stats_files,
            'cell_type_top': celltype_top_files}

rule scdrs_aggregate_stat:
    input: unpack(scdrs_aggregate_stat)
    output:
        cell_type_stats = "results/aggregate/{geno_set}.scdrs.cell_type_stats.tsv",
        cell_type_top = "results/aggregate/{geno_set}.scdrs.cell_type_top.tsv"
    resources:
        time = "24:00:00",
        ncpus = 8,
        mem =  "200G"
    log: "logs/aggregate/{geno_set}.scdrs_cell_type_stat.log"
    conda: "renv"
    script: "snakescripts/aggregate/scdrs_cell_type.R"
    

def scdrs_aggregate_score(x):
    """
    Aggregate scDRS cell-level scores per phenotype
    """
    dir_gs = Path(str(checkpoints.chunk_gs.get(**x).output[0]))
    gs_files = dir_gs.glob("*.gs")
    traits = [str(gs_file.stem) for gs_file in gs_files]

    cell_score_files = [f"results/scdrs/cell_score/{x.geno_set}/{p}.cell_score.tsv.parquet.gz" for p in traits]
    return cell_score_files

rule scdrs_aggregate_score:
    """
    Get umap coordinates (and annotations) of sampled cells
    and cell-level scores for plotting
    """
    input:
        cell_score = scdrs_aggregate_score,
        prep_h5ad = "resources/scdrs/h5ad/{geno_set}.prep.h5ad",
        config = "resources/scdrs/config/{geno_set}.yaml"
    output:
        scores  = "results/aggregate/{geno_set}.scdrs.cell_score.tsv.parquet.gz",
        mcp  = "results/aggregate/{geno_set}.scdrs.cell_mcp.tsv.parquet.gz"
    conda: "scverse"
    resources:
        time = "30:00:00",
        ncpus = 12,
        mem =  "240G"
    log: "logs/aggregate/{geno_set}.scdrs_cell_score.log"
    script: "snakescripts/aggregate/scdrs_cell_score.py"

