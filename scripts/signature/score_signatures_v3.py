#!/usr/bin/env python3

from __future__ import annotations

import os
import pickle
from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd


_SHORT_COLS = [
    "read_id",
    "forward_read_position",
    "ref_position",
    "chrom",
    "ref_strand",
    "read_length",
    "call_prob",
    "call_code",
    "base_qual",
    "ref_kmer",
    "query_kmer",
    "canonical_base"
]

_MAP = {
    "m5C": "m",
    "Y": "17802",
}

# signature convention
_SIG = [
    ("none", ()),
    ("C2U", ("score_C2U",)),
    ("m5C", ("score_m5C",)),
    ("Y", ("score_Y",)),
    ("C2U+m5C", ("score_C2U", "score_m5C")),
    ("C2U+Y", ("score_C2U", "score_Y")),
    ("m5C+Y", ("score_m5C", "score_Y")),
    ("C2U+m5C+Y", ("score_C2U", "score_m5C", "score_Y")),
]

_SIG_NAMES = [s for s, _ in _SIG]

_RC = str.maketrans("ACGTUacgtu", "TGCAAtgcaa")

# bin sizes
_c_bins = 4
_u_bins = 4
_len_bins = 4
_min_size = 50

# chunk size for processing
_chunksize = 1_000_000


def phred_to_prob(q):
    """
    Convert Phred Q score to probability of correct base.
    Q = -10 log10(p_err)
    """
    q = np.asarray(q, dtype=float)
    return 1.0 - np.power(10.0, -q / 10.0)


def log_ratio_score(p, p0):
    """
    Positive log evidence relative to background threshold.
    Returns max(log(p/p0), 0)
    """
    p = np.asarray(p, dtype=float)
    p = np.clip(p, 1e-12, 1.0)
    return np.maximum(np.log(p / p0), 0.0)


def window_local_normalized(sum_arr, count_arr, window_size):
    """
    Max over windows of (sum_evidence / ref_base_count)
    """
    sum_arr = np.asarray(sum_arr, dtype=float)
    count_arr = np.asarray(count_arr, dtype=float)

    if sum_arr.size == 0:
        return 0.0

    if window_size <= 1 or sum_arr.size <= window_size:
        denom = count_arr.sum()
        return float(sum_arr.sum() / denom) if denom > 0 else 0.0

    csum = np.cumsum(np.insert(sum_arr, 0, 0.0))
    ccount = np.cumsum(np.insert(count_arr, 0, 0.0))

    win_sum = csum[window_size:] - csum[:-window_size]
    win_count = ccount[window_size:] - ccount[:-window_size]

    valid = win_count > 0
    if not np.any(valid):
        return 0.0

    scores = np.zeros_like(win_sum)
    scores[valid] = win_sum[valid] / win_count[valid]

    return float(scores.max())


def iter_chunks(path: str | Path, chunksize: int, no_header: bool):
    if no_header:
        yield from pd.read_csv(
            path,
            sep="\t",
            header=None,
            names=_SHORT_COLS,
            chunksize=chunksize,
        )
    else:
        yield from pd.read_csv(
            path,
            sep="\t",
            usecols=_SHORT_COLS,
            chunksize=chunksize,
        )
        
        
def motif_mask(kmer: pd.Series, mode: str) -> pd.Series:
    if mode == "DRACH":
        # actually the motif id RACH*
        drach = kmer.str.match(r"^[AG]AC[ACT][ACGT]$")
        return ~drach.fillna(False)
    if mode == "YC":
        yc = kmer.str.match(r"^[ACGT][CT]C[ACGT][ACGT]$")
        return yc.fillna(False)
    return pd.Series(True, index=kmer.index)


def summarize_chunk(
    chunk: pd.DataFrame,
    mod_prob_thresh: float,
    c2u_basequal_thresh: float,
    c2u_score_mode: str,
    window_size: int,
    m5c_motif_mode: str
) -> pd.DataFrame:
    chunk = chunk.copy()

    chunk["call_prob"] = pd.to_numeric(chunk["call_prob"], errors="coerce")
    chunk["base_qual"] = pd.to_numeric(chunk["base_qual"], errors="coerce")
    chunk["read_length"] = pd.to_numeric(chunk["read_length"], errors="coerce")
    chunk["forward_read_position"] = pd.to_numeric(chunk["forward_read_position"], errors="coerce")

    # keep rows ordered for cluster mode
    chunk = chunk.sort_values(["read_id", "forward_read_position"], kind="mergesort")

    # orient ref kmer to read orientation
    ref = chunk["ref_kmer"].astype(str)
    qry = chunk["query_kmer"].astype(str)
    neg = chunk["ref_strand"].eq("-")
    # replace values where the condition is False!
    ref_oriented = ref.where(~neg, ref.map(lambda s: s.translate(_RC)[::-1]))

    ref_center = ref_oriented.str[2]
    qry_center = qry.str[2]

    # ref base for possible signatures
    is_C_ref = ref_center.eq("C")
    is_U_ref = ref_center.eq("T")

    code = chunk["call_code"].astype(str)
    is_m5c_site = code.eq(_MAP["m5C"]) & motif_mask(ref_oriented, m5c_motif_mode)
    is_y_site = code.eq(_MAP["Y"])

    # C->U proxy: C in reference, T in read, and no direct mod call at same site.
    direct_here = is_m5c_site | is_y_site
    is_c2u_site = is_C_ref & qry_center.eq("T") & (~direct_here)

    # Site-level evidence scores.
    c2u_prob = phred_to_prob(chunk["base_qual"].to_numpy())
    c2u_evidence = np.where(
        is_c2u_site.to_numpy(),
        log_ratio_score(c2u_prob, phred_to_prob(c2u_basequal_thresh)),
        0.0,
    )
    m5c_evidence = np.where(
        is_m5c_site.to_numpy(),
        log_ratio_score(chunk["call_prob"].to_numpy(), mod_prob_thresh),
        0.0,
    )
    y_evidence = np.where(
        is_y_site.to_numpy(),
        log_ratio_score(chunk["call_prob"].to_numpy(), mod_prob_thresh),
        0.0,
    )

    tmp = pd.DataFrame(
        {
            "read_id": chunk["read_id"].astype(str),
            "forward_read_position": chunk["forward_read_position"],
            "read_length": chunk["read_length"],
            "n_C_ref": is_C_ref.astype(np.int32),
            "n_U_ref": is_U_ref.astype(np.int32),
            "n_C2U": is_c2u_site.astype(np.int32),
            "n_m5C": is_m5c_site.astype(np.int32),
            "n_Y": is_y_site.astype(np.int32),
            "sum_C2U": c2u_evidence,
            "sum_m5C": m5c_evidence,
            "sum_Y": y_evidence,
        }
    )

    agg = tmp.groupby("read_id", as_index=False).agg(
        read_length=("read_length", "max"),
        n_C_ref=("n_C_ref", "sum"),
        n_U_ref=("n_U_ref", "sum"),
        n_C2U=("n_C2U", "sum"),
        n_m5C=("n_m5C", "sum"),
        n_Y=("n_Y", "sum"),
        sum_C2U=("sum_C2U", "sum"),
        sum_m5C=("sum_m5C", "sum"),
        sum_Y=("sum_Y", "sum"),
    )

    if c2u_score_mode == "window":
        c2u_score = (
            tmp.groupby("read_id", sort=False)
            .apply(
                lambda df: window_local_normalized(
                    df["sum_C2U"].to_numpy(),
                    df["n_C_ref"].to_numpy(),
                    window_size,
                ),
                include_groups=False
            )
            .reset_index(name="score_C2U")
        )
        agg = agg.merge(c2u_score, on="read_id", how="left")

    return agg

        
def summarize_by_read(
    path: str | Path,
    mod_prob_thresh: float,
    c2u_basequal_thresh: float,
    c2u_score_mode: str,
    window_size: int,
    m5c_motif_mode: str,
    chunksize: int,
    n_jobs: int,
    no_header: bool,
) -> pd.DataFrame:
    chunks: list[pd.DataFrame] = []
    worker = partial(
        summarize_chunk,
        mod_prob_thresh=mod_prob_thresh,
        c2u_basequal_thresh=c2u_basequal_thresh,
        c2u_score_mode=c2u_score_mode,
        window_size=window_size,
        m5c_motif_mode=m5c_motif_mode,
    )

    chunk_iter = iter_chunks(path, chunksize=chunksize, no_header=no_header)

    if n_jobs and n_jobs > 1:
        with ProcessPoolExecutor(max_workers=n_jobs) as ex:
            for out in ex.map(worker, chunk_iter, chunksize=1):
                chunks.append(out)
    else:
        for chunk in chunk_iter:
            chunks.append(worker(chunk))

    concat_chunks = pd.concat(chunks, ignore_index=True)
    reads = concat_chunks.groupby("read_id", as_index=False).agg(
        read_length=("read_length", "max"),
        n_C_ref=("n_C_ref", "sum"),
        n_U_ref=("n_U_ref", "sum"),
        n_C2U=("n_C2U", "sum"),
        n_m5C=("n_m5C", "sum"),
        n_Y=("n_Y", "sum"),
        sum_C2U=("sum_C2U", "sum"),
        sum_m5C=("sum_m5C", "sum"),
        sum_Y=("sum_Y", "sum"),
    )
    # always global, normalized
    reads["score_Y"] = np.where(reads["n_U_ref"] > 0, reads["sum_Y"] / reads["n_U_ref"], 0.0)
    reads["score_m5C"] = np.where(reads["n_C_ref"] > 0, reads["sum_m5C"] / reads["n_C_ref"], 0.0)
        
    if c2u_score_mode == "global":
        # normalize by reference base count i.e. "signature opportunity"
        reads["score_C2U"] = np.where(reads["n_C_ref"] > 0, reads["sum_C2U"] / reads["n_C_ref"], 0.0)
    else:
        # max over read scores across chunks
        agg = concat_chunks.groupby("read_id", as_index=False).agg(
            score_C2U=("score_C2U", "max"),
        )
        reads = reads.merge(agg, on="read_id", how="left")

    # if any chunk missing for any reason... fill with zeros
    for col in ["score_C2U", "score_m5C", "score_Y"]:
        if col not in reads.columns:
            reads[col] = 0.0

    return reads


def qcut_edges(values: pd.Series, n_bins: int) -> np.ndarray:
    """
    Returns monotonic bin edges with -inf/+inf guards.
    Falls back to a single bin if the data are too sparse.
    """
    v = pd.to_numeric(values, errors="coerce").dropna().to_numpy()
    if v.size < 2 or np.unique(v).size < 2:
        return np.array([-np.inf, np.inf], dtype=float)

    qs = np.unique(np.quantile(v, np.linspace(0, 1, n_bins + 1)))
    if qs.size < 3:
        return np.array([-np.inf, np.inf], dtype=float)

    qs[0] = -np.inf
    qs[-1] = np.inf
    return qs


def cut_bins(values: pd.Series, edges: np.ndarray) -> pd.Series:
    if edges.size <= 2:
        return pd.Series(["all"] * len(values), index=values.index, dtype="string")
    b = pd.cut(pd.to_numeric(values, errors="coerce"), bins=edges, include_lowest=True)
    return b.astype("string").fillna("NA")


def make_bins(control_reads: pd.DataFrame, len_bins: int, c_bins: int, u_bins: int) -> dict:
    return {
        "len_edges": qcut_edges(control_reads["read_length"], len_bins),
        "c_edges": qcut_edges(control_reads["n_C_ref"], c_bins),
        "u_edges": qcut_edges(control_reads["n_U_ref"], u_bins),
    }


def assign_bins(reads: pd.DataFrame, model: dict) -> pd.DataFrame:
    x = reads.copy()
    x["len_bin"] = cut_bins(x["read_length"], model["len_edges"])
    x["c_bin"] = cut_bins(x["n_C_ref"], model["c_edges"])
    x["u_bin"] = cut_bins(x["n_U_ref"], model["u_edges"])
    x["bin_key"] = list(zip(x["len_bin"].astype(str), x["c_bin"].astype(str), x["u_bin"].astype(str)))
    return x


def add_signatures(reads: pd.DataFrame) -> pd.DataFrame:
    """
    Compute raw subset scores and assign signatures.
    Significance is handled entirely by the empirical control null.
    """
    x = reads.copy()

    # ensure these exist in both score modes
    if "score_C2U" not in x.columns:
        x["score_C2U"] = 0.0
    if "score_m5C" not in x.columns:
        x["score_m5C"] = 0.0
    if "score_Y" not in x.columns:
        x["score_Y"] = 0.0

    x["score_none"] = 0.0
    x["score_C2U+m5C"] = x["score_C2U"] + x["score_m5C"]
    x["score_C2U+Y"] = x["score_C2U"] + x["score_Y"]
    x["score_m5C+Y"] = x["score_m5C"] + x["score_Y"]
    x["score_C2U+m5C+Y"] = x["score_C2U"] + x["score_m5C"] + x["score_Y"]

    score_cols = [f"score_{s}" for s in _SIG_NAMES]
    scores = x[score_cols].to_numpy(dtype=float)

    best_idx = scores.argmax(axis=1)
    best_sorted = np.sort(scores, axis=1)

    x["signature"] = [_SIG_NAMES[i] for i in best_idx]
    x["score"] = scores[np.arange(len(x)), best_idx]
    x["score_margin"] = x["score"] - best_sorted[:, -2]

    return x


def fit_null_model(
    control_reads: pd.DataFrame,
    len_bins: int,
    c_bins: int,
    u_bins: int,
    min_bin_size: int,
    thr_quantile: float,
) -> dict:
    """
    Empirical null for the *best subset score* from the control reads.
    """
    ctrl = add_signatures(control_reads)
    model = make_bins(ctrl, len_bins=len_bins, c_bins=c_bins, u_bins=u_bins)
    ctrl = assign_bins(ctrl, model)

    null_bins = {}
    for (key, sig), sub in ctrl.groupby(["bin_key", "signature"], dropna=False):
        if len(sub) < min_bin_size:
            continue
        null_bins[(str(key), str(sig))] = np.sort(sub["score"].to_numpy(dtype=float))

    global_null = {
        sig: np.sort(ctrl.loc[ctrl["signature"].eq(sig), "score"].to_numpy(dtype=float))
        for sig in ctrl["signature"].unique()
    }

    pooled = np.concatenate([v for v in global_null.values() if len(v) > 0]) if any(len(v) > 0 for v in global_null.values()) else np.array([], dtype=float)
    if pooled.size == 0:
        pooled = np.array([0.0], dtype=float)

    fallback_thr = float(np.quantile(pooled, thr_quantile))
    candidate_thresholds = {
        sig: float(np.quantile(vals, thr_quantile)) if len(vals) > 0 else fallback_thr
        for sig, vals in global_null.items()
    }
    for sig in _SIG_NAMES:
        candidate_thresholds.setdefault(sig, fallback_thr)

    model["candidate_thresholds"] = candidate_thresholds
    model["candidate_threshold_fallback"] = fallback_thr
    model["null_bins"] = null_bins
    model["global_null"] = global_null
    model["global_null_all"] = np.sort(pooled)
    model["min_bin_size"] = min_bin_size
    model["thr_quantile"] = thr_quantile
    return ctrl, model


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """
    Benjamini-Hochberg FDR.
    """
    p = np.asarray(pvals, dtype=float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(1, n + 1))
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty_like(ranked)
    q[order] = np.clip(ranked, 0.0, 1.0)
    return q


def empirical_tail_p(sorted_null: np.ndarray, obs: np.ndarray) -> np.ndarray:
    """
    Tail probability with +1 smoothing:
      p = (1 + #{null >= obs}) / (n + 1)
    """
    if sorted_null.size == 0:
        return np.ones_like(obs, dtype=float)
    obs = np.asarray(obs, dtype=float)
    n = sorted_null.size
    idx = np.searchsorted(sorted_null, obs, side="left")
    return (n - idx + 1) / (n + 1)


def score_reads(query_reads: pd.DataFrame, model: dict) -> pd.DataFrame:
    x = add_signatures(query_reads)
    x = assign_bins(x, model)

    pvals = np.empty(len(x), dtype=float)

    for (key, sig), idx in x.groupby(["bin_key", "signature"]).groups.items():
        idx = np.asarray(list(idx))

        null = model["null_bins"].get(
            (str(key), str(sig)),
            model["global_null"].get(str(sig), model["global_null_all"]),
        )
        pvals[idx] = empirical_tail_p(
            null,
            x.loc[idx, "score"].to_numpy(dtype=float),
        )

    x["p_value"] = pvals

    # two-stage candidate filter
    x["candidate_threshold"] = x["signature"].map(model["candidate_thresholds"]).fillna(model["candidate_threshold_fallback"])
    x["candidate"] = x["score"] > x["candidate_threshold"]

    # FDR only on candidates
    x["fdr"] = 1.0
    cand_idx = x["candidate"].to_numpy()
    if cand_idx.sum() > 0:
        x.loc[cand_idx, "fdr"] = bh_fdr(x.loc[cand_idx, "p_value"].to_numpy())
    x["significant_fdr_0_05"] = x["candidate"] & (x["fdr"] < 0.05)

    return x


def signature_summary(scored: pd.DataFrame) -> pd.DataFrame:
    return (
        scored.groupby("signature", dropna=False)
        .agg(
            n_reads=("read_id", "size"),
            median_score=("score", "median"),
            median_p=("p_value", "median"),
            median_fdr=("fdr", "median"),
            significant_reads=("significant_fdr_0_05", "sum"),
            median_margin=("score_margin", "median"),
        )
        .reset_index()
        .sort_values(["significant_reads", "n_reads"], ascending=False)
    )


def main():
    parser = ArgumentParser()
    parser.add_argument("--control", required=True, help="Extract-calls table for control")
    parser.add_argument("--condition", required=True, help="Extract-calls table for condition")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--name", type=str, default="dflt", help="Name prefix")

    parser.add_argument("--min-prob", type=float, default=0.85, help="Min. call prob for m5C/Y scoring")
    parser.add_argument("--min-qual", type=int, default=30, help="Min. base qual for C->U scoring")
    parser.add_argument("--thr-quantile", type=float, default=0.995, help="Candidate threshold quantile (FDR)")

    parser.add_argument("--c2u-mode", choices=["global", "window"], default="global",
                        help="C2U scoring model")
    parser.add_argument("--window-size", type=int, default=51,
                        help="Sliding window size for window mode (only relevant for C2U model)")
    
    parser.add_argument("--motif-mode", choices=["all", "DRACH", "YC"], default="all",
                        help="""Optional motif filtering: all - no filtering, DRACH - ignore m5C sites on DRACH 
                        motifs (C position), or YC - keep m5C sites on YC motifs (Y=C/T).""")
    
    parser.add_argument("--n-jobs", type=int, default=1,
                        help="Parallel chunk processing; use >1 only if RAM allows")
    parser.add_argument("--chunksize", type=int, default=_chunksize,
                        help="CSV chunk size")
    parser.add_argument("--no-header", action="store_true", help="Input files have no header (hard-coded columns)")

    args = parser.parse_args()

    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Control reads: {args.control}.")
    print(f"Condition reads: {args.condition}.")
    print(f"Writing \"{args.name}\" results to: {args.output}.")
    print(f"Used for scoring: min_prob = {args.min_prob:.4f}, min_qual = {args.min_qual} ({phred_to_prob(args.min_qual):.4f}).")
    print(f"Threshold (adaptive FDR correction): quantile = {args.thr_quantile:.4f}.")
    print(f"C2U mode: {args.c2u_mode} (window size = {args.window_size} only used in window mode).")
    print(f"Motif mode: {args.motif_mode}.")

    control_reads = summarize_by_read(
        args.control,
        args.min_prob,
        args.min_qual,
        args.c2u_mode,
        args.window_size,
        args.motif_mode,
        args.chunksize,
        args.n_jobs,
        args.no_header
    )

    condition_reads = summarize_by_read(
        args.condition,
        args.min_prob,
        args.min_qual,
        args.c2u_mode,
        args.window_size,
        args.motif_mode,
        args.chunksize,
        args.n_jobs,
        args.no_header
    )

    control_reads_scored, model = fit_null_model(
        control_reads,
        _len_bins,
        _c_bins,
        _u_bins,
        _min_size,
        args.thr_quantile,
    )

    scored = score_reads(condition_reads, model)

    scored.to_csv(outdir / f"{args.name}_scores.tab.gz", sep="\t", index=False, compression="gzip")
    signature_summary(scored).to_csv(outdir / f"{args.name}_summary.tab", sep="\t", index=False)
    
    control_reads_scored.to_csv(outdir / f"{args.name}_control.tab.gz", sep="\t", index=False, compression="gzip")

    print("Significant reads (FDR < 0.05):", int(scored["significant_fdr_0_05"].sum()))
    print("Signatures:")
    print(scored["signature"].value_counts().reindex(_SIG_NAMES, fill_value=0).to_string())

    with open(outdir / f"{args.name}_null_model.pkl", "wb") as fh:
        pickle.dump(model, fh)


if __name__ == "__main__":
    main()
