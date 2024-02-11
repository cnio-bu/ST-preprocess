import sys
sys.stderr = open(snakemake.log[0], "a+")
sys.stdout = open(snakemake.log[0], "a+")

import random
import numpy as np
import scanpy as sc
import paste as pst
import pandas as pd
import session_info

## SNAKEMAKE I/O ##
target_counts = snakemake.input["raw_counts"]
target_coords = snakemake.input["coordinates"]
ref_counts = snakemake.input["ref_raw_counts"]
ref_coords = snakemake.input["ref_coordinates"]

## FUNCTIONS ##
def spatial_scanpy(counts, coords):
    "Returns a scanpy spatial object. Requires raw and already filtered counts."
    scanpyobj = sc.read_csv(counts, delimiter = "\t", first_column_names = True)
    scanpyobj = scanpyobj.transpose()
    scanpyobj.obsm["spatial"] = np.genfromtxt(coords, delimiter = "\t", \
                                              skip_header = 1, usecols = (1, 2))
    return scanpyobj

def generalized_procrustes_analysis(X, Y, pi, output_params = False, matrix = False):
    """
    From https://github.com/raphael-group/paste visualization.py.
    Finds and applies optimal rotation between spatial coordinates of two layers (may also do a reflection).
    Args:
        X: np array of spatial coordinates (ex: sliceA.obs['spatial'])
        Y: np array of spatial coordinates (ex: sliceB.obs['spatial'])
        pi: mapping between the two layers output by PASTE
        output_params: Boolean of whether to return rotation angle and translations along with spatial coordinates.
        matrix: Boolean of whether to return the rotation as a matrix or an angle.
    Returns:
        Aligned spatial coordinates of X, Y, rotation angle, translation of X, translation of Y.
    """
    assert X.shape[1] == 2 and Y.shape[1] == 2

    tX = pi.sum(axis=1).dot(X)
    tY = pi.sum(axis=0).dot(Y)
    X = X - tX
    Y = Y - tY
    H = Y.T.dot(pi.T.dot(X))
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T.dot(U.T)
    Y = R.dot(Y.T).T
    if output_params and not matrix:
        M = np.array([[0,-1],[1,0]])
        theta = np.arctan(np.trace(M.dot(H))/np.trace(H))
        return X,Y,theta,tX,tY
    elif output_params and matrix:
        return X, Y, R, tX, tY
    else:
        return X,Y

## CODE ##
# Random seed
random.seed(1)

# Create Scanpy objects
target = spatial_scanpy(target_counts, target_coords)
ref = spatial_scanpy(ref_counts, ref_coords)

# Pair-wise alignment of slides
pis = pst.pairwise_align(target, ref)

# Compute aligned coordinates
Target, Ref, angle, tTarget, tRef = generalized_procrustes_analysis(target.obsm["spatial"], ref.obsm["spatial"], pis, \
                                                                    output_params = True, matrix = True)
Ref = Ref.dot(np.linalg.inv(angle.T))

# Create pandas dataframes
colnames = ["adj_x", "adj_y"]
target_df = pd.DataFrame(Target, index = target.obs_names, columns = colnames)
ref_df = pd.DataFrame(Ref, index = ref.obs_names, columns = colnames)

# Save dataframes
target_df.to_csv(snakemake.output["coordinates_aligned"], sep = "\t", header = True, index = True)
ref_df.to_csv(snakemake.output["ref_coordinates_aligned"], sep = "\t", header = True, index = True)

# Session info
session_info.show()