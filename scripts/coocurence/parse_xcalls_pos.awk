BEGIN {
    FS = "\t"
    OFS = "\t"

    # column indices
    read_id_idx = 1
    forward_read_position_idx = 2
    read_length_idx = 12
    call_prob_idx = 13
    call_code_idx = 14
    query_kmer_idx = 17
    fail_idx = 20
    inferred_idx = 21
    within_alignment_idx = 22

    # ---- motif lists ----
    split("GGACT GGACA GAACT AGACT GGACC TGACT AAACT GAACA AGACA AGACC GAACC TGACA TAACT AAACA TGACC TAACA AAACC TAACC", m6A_motifs, " ")
    for (i in m6A_motifs) motif["a", m6A_motifs[i]] = 1

    split("GTTCA GTTCC GTTCG GTTCT TGTAG", Y_base, " ")
    for (i in Y_base) motif["17802", Y_base[i]] = 1

    # generated motifs for 17802
    split("A C G T", p1, " ")
    split("A G", p2, " ")
    split("A G", p4, " ")
    split("A C G T", p5, " ")

    for (i in p1)
        for (j in p2)
            for (k in p4)
                for (l in p5)
                    motif["17802", p1[i] p2[j] "T" p4[k] p5[l]] = 1

}

NR==1 {print "read_id", "forward_read_position", "read_length", "call_code", "call_prob", "query_kmer", "fail"; next}

{
    # keep only aligned, non-inferred
    if ($within_alignment_idx != "true")
        next
    if ($inferred_idx != "false")
        next
    # keep only a or 17802
    if ($call_code_idx != "a" && $call_code_idx != "17802")
        next
    # motif filter
    if (!(($call_code_idx SUBSEP $query_kmer_idx) in motif))
        next 
    
    print $read_id_idx, $forward_read_position_idx, $read_length_idx, $call_code_idx, $call_prob_idx, $query_kmer_idx, $fail_idx
    
}
