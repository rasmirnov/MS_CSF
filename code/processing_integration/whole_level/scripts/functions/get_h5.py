import h5py
import scipy
from anndata import AnnData
import numpy as np

def get_h5(adata: AnnData, outfile: str):
    counts = scipy.sparse.csc_matrix(adata.X)
    with h5py.File(outfile, 'w') as hf:
        x = hf.create_group('X')
        x.create_dataset('indptr', data=counts.indptr, compression="gzip", compression_opts=9,
                         chunks=True, dtype='int32')
        x.create_dataset('indices', data=counts.indices, compression="gzip", compression_opts=9,
                         chunks=True, dtype='int32')
        x.create_dataset('data', data=counts.data, compression="gzip",
                         compression_opts=9, chunks=True, dtype='float64')
        hf['X'].attrs['shape'] = np.array(counts.shape, dtype='int32')
