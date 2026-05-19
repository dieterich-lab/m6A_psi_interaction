#!/usr/bin/env python3

import numpy as np
import pandas as pd
from argparse import ArgumentParser
from pathlib import Path
from joblib import Parallel, delayed


SEED = 1984


def mean_mutual_nn(a, b):
    """
    Mean mutual nearest-neighbor distance between two sets.
    More stable than min distance.
    """
    dist = np.abs(a[:, None] - b[None, :])
    a_to_b = dist.argmin(axis=1)
    b_to_a = dist.argmin(axis=0)

    mnnd = []
    for i, j in enumerate(a_to_b):
        if b_to_a[j] == i:
            mnnd.append(dist[i, j])
            
    if len(mnnd) == 0:
        return np.nan

    return np.mean(mnnd)


def mean_mutual_nn_self(a):
    """
    Mean mutual nearest-neighbor distance within one set.
    Used for homotypic clustering.
    """
    dist = np.abs(a[:, None] - a[None, :]).astype(float)
    np.fill_diagonal(dist, np.inf)

    nn = dist.argmin(axis=1)

    mnnd = []
    for i, j in enumerate(nn):
        if nn[j] == i:
            mnnd.append(dist[i, j])

    if len(mnnd) == 0:
        return np.nan

    return np.mean(mnnd)


def get_observed_and_null(read_id, g, mod_a, mod_b,
                 n_perm, min_sites, min_sites_each, seed):
    """Find per-read MEAN distance of mutual NN modification sites."""
    rng = np.random.default_rng(seed + hash(read_id) % 10**6)

    g_a = g[g["call_code"] == mod_a]
    g_b = g[g["call_code"] == mod_b]

    pos_a_all = g_a["forward_read_position"].to_numpy()
    p_a = g_a["call_prob"].to_numpy()

    pos_b_all = g_b["forward_read_position"].to_numpy()
    p_b = g_b["call_prob"].to_numpy()

    # observed (hard threshold)
    # use of modkit "fail" should be equivalent...
    keep_a = p_a > 0.5
    keep_b = p_b > 0.5

    pos_a = pos_a_all[keep_a]
    pos_b = pos_b_all[keep_b]

    n_a = pos_a.size
    n_b = pos_b.size

    is_self = (mod_a == mod_b)
    
    if is_self:
        # filter out "sparse" reads
        if n_a < min_sites_each:
            return None
        if n_a < min_sites:
            return None
        obs = mean_mutual_nn_self(pos_a)
    else:
        # filter out "sparse" reads
        if n_a < min_sites_each or n_b < min_sites_each:
            return None
        if (n_a + n_b) < min_sites:
            return None
        obs = mean_mutual_nn(pos_a, pos_b)
    if np.isnan(obs):
        return None

    # null (Bernoulli sampling)
    # P(U<p)=p ~ Bernouilli(p) since U is uniform [0,1)
    null = np.full(n_perm, np.nan)
    valid = 0
    max_tries = n_perm * 5
    tries = 0
    while valid < n_perm and tries < max_tries:
        tries += 1
        # Bernoulli draws
        samp_a = rng.random(p_a.size) < p_a
        samp_b = rng.random(p_b.size) < p_b
        # skip invalid draws
        if samp_a.sum() < 1 or samp_b.sum() < 1:
            continue
        
        samp_pos_a = pos_a_all[samp_a]
        samp_pos_b = pos_b_all[samp_b]
        
        samp_n_a = samp_pos_a.size
        samp_n_b = samp_pos_b.size
        
        if is_self:
            if samp_n_a < min_sites_each:
                continue
            if samp_n_a < min_sites:
                continue
            null_val = mean_mutual_nn_self(samp_pos_a)
        else:
            if samp_n_a < min_sites_each or samp_n_b < min_sites_each:
                continue
            if (samp_n_a + samp_n_b) < min_sites:
                continue
            null_val = mean_mutual_nn(samp_pos_a, samp_pos_b) 
        null[valid] = null_val
        valid += 1
        
    null_full = null.copy()
    null = null[~np.isnan(null)]
    # arbitrary... what value is meaningfull?
    if len(null) < int(n_perm/2):
        return None
    # remove degenerate nulls
    if np.std(null, ddof=1) == 0:
        return None

    # empirical p-values
    p_left = (1 + np.sum(null <= obs)) / (1 + len(null))
    p_right = (1 + np.sum(null >= obs)) / (1 + len(null))
    p_two = 2 * min(p_left, p_right)
    p_two = min(p_two, 1.0)

    return {
        "read_id": read_id,
        "n_a": n_a,
        "n_b": n_b,
        "read_len": int(g["read_length"].iloc[0]),
        "obs_distance": obs,
        "null_distance": null_full,
        "null_mean": null.mean(),
        "null_median": np.median(null),
        "null_sd": null.std(ddof=1),
        "p_left": p_left,
        "p_right": p_right,
        "p_two": p_two,
        "n_perm": valid
    }


def bh_fdr(pvals):
    n = pvals.size
    order = np.argsort(pvals)
    ranked = pvals[order]
    qvals = np.empty(n)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        q = ranked[i] * n / (i + 1)
        q = min(q, prev)
        prev = q
        qvals[i] = q
    out = np.empty(n)
    out[order] = qvals
    return out


def main():
    parser = ArgumentParser()
    parser.add_argument("--calls", required=True, 
                        help="Selected calls from modkit extract calls")
    parser.add_argument("--codes", nargs=2, type=str, default=["a", "17802"],
                        help='Two modification codes')
    parser.add_argument("--output", required=True,
                        help='Output directory')
    parser.add_argument('--name', type=str, required=True,
                        help='Name')
    parser.add_argument('--plot-name', type=str, default=None,
                        help='name of plot (w/o extension)')
    parser.add_argument('--min-sites', type=int, default=10,
                        help='Min. number of sites (confounded) per read')
    parser.add_argument('--min-sites-each', type=int, default=10,
                        help='Min. number of sites (for each modification separately) per read. Value is included.')
    parser.add_argument("--n-perm", type=int, default=500)
    parser.add_argument("--n-jobs", type=int, default=24)
    args = parser.parse_args()
   
    if args.plot_name is not None:
        args.name = f"{args.name}_{args.plot_name}"
        
    if args.min_sites < 2:
        print("Setting [--min-sites] to 2")
        args.min_sites = 2
    
    mod_a, mod_b = args.codes
    is_self = (mod_a == mod_b)

    calls = pd.read_csv(args.calls, sep="\t")
    groups = list(calls.groupby("read_id", sort=False))
    results = Parallel(n_jobs=args.n_jobs)(
        delayed(get_observed_and_null)(
            read_id, group, mod_a, mod_b,
            args.n_perm,
            args.min_sites,
            args.min_sites_each,
            SEED
        )
        for read_id, group in groups
    )
    results = [r for r in results if r is not None]
    null_mat = np.vstack([r["null_distance"] for r in results])
    results = [{k:v for k,v in r.items() if k != "null_distance"} for r in results]
    df = pd.DataFrame(results)
    
    # FDR 2-sided
    df["fdr_two"] = bh_fdr(df["p_two"].values)
    
    columns = {
        "n_a": f"n_{mod_a}",
        "n_b": f"n_{mod_b}"
    }
    if is_self:
        columns.pop("n_b")
        df.drop(columns=["n_b"], inplace=True)
        args.name = f"{args.name}_self{mod_a}"
    df.rename(columns=columns, inplace=True)

    df.to_csv(Path(args.output, f"{args.name}_cooccurence_motif.tsv.gz"), sep="\t", index=False)
    
    np.savez(
        Path(args.output, f"{args.name}_cooccurence_motif_nulls.npz"),
        nulls=null_mat,
    )


if __name__ == "__main__":
    main()
