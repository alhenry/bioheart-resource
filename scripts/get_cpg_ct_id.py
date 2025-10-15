# get cpg id and cell index from h5ad file

import scanpy as sc
import pandas as pd
from pathlib import Path
import numpy as np
import h5py

from anndata.io import read_elem

STUDY = "bioheart975-full"
ADATA = f"resources/scdrs/h5ad/{STUDY}.prep.h5ad"

# get unique ct_id and cpg_id combinations
with h5py.File(ADATA) as f:
    df = pd.DataFrame(
        {'ct_id': read_elem(f["obs/ct_id"]),
         'cpg_id': read_elem(f["obs/cpg_id"])}
    ).drop_duplicates(subset=['ct_id', 'cpg_id'])


df.to_csv(f"resources/{STUDY}.ctid_cpgid.tsv", sep="\t", index=False)