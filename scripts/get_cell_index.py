# get cpg id and cell index from h5ad file

import scanpy as sc
import pandas as pd
from pathlib import Path
import numpy as np
import h5py

from anndata.io import read_elem


STUDY = "bioheart975-full"
ADATA = f"resources/scdrs/h5ad/{STUDY}.prep.h5ad"

with h5py.File(ADATA) as f:
    df = pd.DataFrame(
        {'index': read_elem(f["obs/_index"]),
         'cpg_id': read_elem(f["obs/cpg_id"])}
    )

df.to_csv(f"resources/{STUDY}.cell_cpgid.tsv", sep="\t", index=False)