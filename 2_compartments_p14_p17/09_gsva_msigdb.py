# ---------------------------------------------------------------------------
# GSVA via the python implementation, MSigDB v6.2.
#
# ESTABLISHES
#     Gene-set score matrices over the full (or subsampled) expression matrix
#     from step 07, using the python GSVA package rather than the Bioconductor
#     one. Both ran; both outputs are on the cluster.
#
# WHY A SECOND IMPLEMENTATION
#     The R GSVA steps (08, 08b) score a senescence-focused collection at P17.
#     These score the broad MSigDB collections across all cells, which was
#     tractable in the python implementation's parallel mode.
#
# KNOWN ODDITY - not a transcription error
#     The imports include pysnow, flask and flask_snow, which have nothing to do
#     with gene-set scoring. They are template leftovers and were present in the
#     scripts as run. They are left in place; note that the scripts will not
#     import without those packages installed.
#
# REFERENCE RELEASE
#     MSigDB v6.2 here, against v7.1 in the P5-P17 analysis. Do not substitute.
#
# INPUTS / OUTPUTS
#     see the paths in the body; outputs are the gsva.exprs.Full_eset.* csv
#     files under GSVA/NotImputed/ (150 MB to 2.9 GB each).
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 0. parameters - edit PROJ_DIR to point at your own copy
# ---------------------------------------------------------------------------
PROJ_DIR = "/project/def-jsjoyal/gaelcge"
GMT_DIR  = PROJ_DIR + "/Gene_list/Gmt.file"

import pandas as pd
from GSVA import gsva, gmt_to_dataframe
import pysnow
from flask import Flask, jsonify
from flask_snow import Snow
import os

if __name__ == '__main__':
    in_file  = 'exprs.Full_eset.NonImputed.scaled.csv'
    fname  = GMT_DIR + '/sup_h.c2.c5.c6.c7.all.v6.2.symbols.gmt'
    curr_dir=os.getcwd()
    #fname  = GMT_DIR + '/h.all.v6.2.symbols.gmt'
    out_file = 'gsva.exprs.Full_eset.NonImputed.sup_h.c2.c5.c6.c7.csv'
    # ex_matrix is a DataFrame with gene names as column names
    expression_df = pd.read_csv(in_file,index_col=0)
    expression_df.iloc[0:5,0:5]
    #expression_df = expression_df.iloc[:,0:2]
    # tf_names is read using a utility function included in Arboreto
    geneset_df = gmt_to_dataframe(fname)
    geneset_df.head()
    # compute the GSVA
    df = gsva(expression_df, geneset_df=geneset_df, method='gsva', kcdf='Gaussian', abs_ranking=False, min_sz=1, max_sz=None, parallel_sz=14, parallel_type='SOCK', mx_diff=True, tau=None, ssgsea_norm=True, verbose=False, tempdir=curr_dir)
    # write the GRN to file
    df.to_csv(out_file)
    exit()

