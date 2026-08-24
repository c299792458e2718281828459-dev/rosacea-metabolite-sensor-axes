#!/usr/bin/env python3
"""03_mebocost_receiver_side.py — M3 module
Metabolite-sensor communication scoring on scRNA-seq, receiver-side only
(microbial metabolites are exogenous senders).

Design (documented in Methods):
  - Database: MEBOCOST human met_sensor DB (793 pairs; kaifuchenlab/MEBOCOST,
    data/mebocost_db/human/human_met_sensor_update_Oct21_2025.tsv)
  - Candidates: M2 candidate metabolites (genus_metabolite_sensor_evidence.tsv
    metabolites with sensors)
  - Input: scRNA-seq expression matrix (genes x cells) + cell annotations
    (HRA006167, Chen et al. 2024, pending controlled-access approval)
  - Score for each (metabolite, sensor, cell_type):
      comm = mean_log_expr(sensor in cell type) * pct.exp(sensor in cell type)
             * specificity (log2 fold vs other cell types)
  - Significance: permutation test — shuffle sensor-metabolite links 1000x,
    FDR (BH) on empirical P values.
Output: output/m3_mebocost/comm_scores.tsv, comm_significant.tsv
"""
import sys
import argparse
import numpy as np
import pandas as pd
from scipy.sparse import issparse

def load_expression(path):
    """Load gene x cell matrix from h5ad, mtx dir, or tsv."""
    if path.endswith('.h5ad'):
        import anndata as ad
        a = ad.read_h5ad(path)
        X = a.X
        genes = list(a.var_names)
        cells = list(a.obs_names)
        return X, genes, cells, a.obs
    if path.endswith('.tsv') or path.endswith('.txt'):
        df = pd.read_csv(path, sep='\t', index_col=0)
        return df.values, list(df.index), list(df.columns), None
    raise ValueError('unsupported matrix format; use h5ad or tsv')

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--matrix', required=True, help='scRNA-seq expression matrix (h5ad or genes x cells tsv)')
    ap.add_argument('--cell-annotation', required=True, help='tsv with columns: cell, cell_type')
    ap.add_argument('--sensor-db', default='data/MEBOCOST-main/data/mebocost_db/human/'
                     'human_met_sensor_update_Oct21_2025.tsv')
    ap.add_argument('--candidates', default='output/m2_microbiome/'
                    'genus_metabolite_sensor_evidence.tsv')
    ap.add_argument('--n-perm', type=int, default=1000)
    ap.add_argument('--min-pct', type=float, default=10.0)
    ap.add_argument('--outdir', default='output/m3_mebocost')
    ap.add_argument('--seed', type=int, default=20260824)
    args = ap.parse_args()

    import os
    os.makedirs(args.outdir, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    # --- Load data
    X, genes, cells, obs = load_expression(args.matrix)
    if issparse(X):
        X = X.toarray()
    ann = pd.read_csv(args.cell_annotation, sep='\t')
    ann = ann.set_index('cell')
    ann = ann.loc[[c for c in cells if c in ann.index]]
    keep_cells = list(ann.index)
    X = X[:, [cells.index(c) for c in keep_cells]]
    cells = keep_cells
    cell_types = sorted(ann['cell_type'].unique())
    print(f'cells: {len(cells)}, cell types: {cell_types}')

    # --- Sensor DB + candidates
    db = pd.read_csv(args.sensor_db, sep='\t')
    cand = pd.read_csv(args.candidates, sep='\t')
    cand = cand[cand['metabolite_sensors_MEBOCOST'].notna()]
    # candidate metabolite names -> DB HMDB-matched names (name-level match)
    db['syn'] = db['metName'].str.lower()
    pairs = []
    for _, r in cand.iterrows():
        met = r['metabolite']
        sensor = r['metabolite_sensors_MEBOCOST']
        base = met.split(' (')[0].lower()
        dbr = db[db['syn'].str.contains(base, regex=False, na=False)]
        for _, dr in dbr.iterrows():
            pairs.append((met, sensor, r['sensor_type']))
    pairs = sorted(set(pairs))
    print(f'candidate (metabolite, sensor) pairs: {len(pairs)}')

    # --- Sensor expression per cell type
    gidx = {g.upper(): i for i, g in enumerate(genes)}
    out_rows = []
    for met, sensor, stype in pairs:
        if sensor.upper() not in gidx:
            continue
        i = gidx[sensor.upper()]
        gx = X[i, :]
        for ct in cell_types:
            m = np.array([c == ct for c in ann['cell_type']])
            if m.sum() < 3:
                continue
            x = gx[m]
            pct = 100 * (x > 0).mean()
            avg = np.log1p(x).mean()
            others = gx[~m]
            avg_other = np.log1p(others).mean() if len(others) else 0
            spec = avg - avg_other
            if pct >= args.min_pct:
                out_rows.append([met, sensor, stype, ct, pct, avg, spec,
                                 pct * avg * spec])
    comm = pd.DataFrame(out_rows, columns=['metabolite', 'sensor', 'sensor_type',
                                           'cell_type', 'pct_exp', 'avg_log1p',
                                           'specificity', 'comm_raw'])
    if len(comm) == 0:
        print('no communications passing min_pct'); sys.exit(0)

    # --- Permutation test: shuffle sensor-metabolite mapping
    # Observed comm per row; permuted comm distribution per row
    pvals = []
    for _, r in comm.iterrows():
        obs = r['comm_raw']
        # permute: random sensor gene from DB (not the true one) on same cell type
        rand_sensors = db['Gene_name'].str.upper().unique()
        perm_vals = []
        n_avail = 0
        for _ in range(args.n_perm):
            s = rng.choice(rand_sensors)
            if s not in gidx:
                continue
            m = np.array([c == r['cell_type'] for c in ann['cell_type']])
            x = X[gidx[s], m]
            pct = 100 * (x > 0).mean()
            avg = np.log1p(x).mean()
            others = X[gidx[s], ~m]
            avg_other = np.log1p(others).mean() if len(others) else 0
            perm_vals.append(pct * avg * (avg - avg_other))
            n_avail += 1
            if n_avail >= 200:
                break
        perm_vals = np.array(perm_vals)
        pvals.append((1 + (perm_vals >= obs).sum()) / (1 + len(perm_vals)))
    comm['p_value'] = pvals
    comm['fdr'] = comm['p_value'] * len(comm) / comm['p_value'].rank(method='first')
    comm['fdr'] = comm['fdr'].clip(upper=1)
    comm = comm.sort_values('fdr')

    comm.to_csv(f'{args.outdir}/comm_scores.tsv', sep='\t', index=False)
    sig = comm[comm['fdr'] < 0.05]
    sig.to_csv(f'{args.outdir}/comm_significant.tsv', sep='\t', index=False)
    print(f'total scored: {len(comm)}, significant (FDR<0.05): {len(sig)}')
    if len(sig):
        print(sig.head(15).to_string())

if __name__ == '__main__':
    main()
