import argparse
from Bio import SeqIO
from Bio.Align import PairwiseAligner
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

def align_one_query(query_id, query_seq, reference_sequences, aligner_settings):

    aligner = PairwiseAligner()
    aligner.mode = aligner_settings["mode"]
    aligner.match_score = aligner_settings["match_score"]
    aligner.mismatch_score = aligner_settings["mismatch_score"]
    aligner.open_gap_score = aligner_settings["open_gap_score"]
    aligner.extend_gap_score = aligner_settings["extend_gap_score"]
    aligner.query_left_open_gap_score   = aligner_settings["query_left_open_gap_score"]
    aligner.query_left_extend_gap_score = aligner_settings["query_left_extend_gap_score"]
    aligner.query_right_open_gap_score  = aligner_settings["query_right_open_gap_score"]
    aligner.query_right_extend_gap_score= aligner_settings["query_right_extend_gap_score"]
    aligner.target_left_open_gap_score  = aligner_settings["target_left_open_gap_score"]
    aligner.target_left_extend_gap_score= aligner_settings["target_left_extend_gap_score"]
    aligner.target_right_open_gap_score = aligner_settings["target_right_open_gap_score"]
    aligner.target_right_extend_gap_score=aligner_settings["target_right_extend_gap_score"]

    best_score_for_query = float("-inf")
    best_alignment_tuple = None

    for ref_id, ref_seq in reference_sequences.items():
        all_alns = aligner.align(query_seq, ref_seq)

        local_best_score = float("-inf")
        local_best_data = None

        for aln in all_alns:
            q_blocks = aln.aligned[0]  # query blocks
            r_blocks = aln.aligned[1]  # reference blocks

            # If no matched blocks, skip
            if len(q_blocks) == 0 or len(r_blocks) == 0:
                continue

            q_start_index = q_blocks[0][0]
            q_end_index   = q_blocks[-1][1]
            s_start_index = r_blocks[0][0]
            s_end_index   = r_blocks[-1][1]

            q_start = q_start_index + 1
            s_start = s_start_index + 1

            virtual_qstart = q_start - (s_start - 1)
            virtual_qend = virtual_qstart + len(ref_seq) - 1

            # Build the corrected query substring (or add gaps if needed)
            if virtual_qstart < 1:
                num_gaps_left = 1 - virtual_qstart
                padded_query = "-" * num_gaps_left + query_seq
                qseq_corrected = padded_query[:len(ref_seq)]
            else:
                qseq_corrected = query_seq[virtual_qstart - 1 : virtual_qstart - 1 + len(ref_seq)]

            if len(qseq_corrected) < len(ref_seq):
                qseq_corrected += "-" * (len(ref_seq) - len(qseq_corrected))

            sseq_corrected = ref_seq

            # Skip left padding for identity calculation if virtual_qstart < 1
            if virtual_qstart < 1:
                padding = 1 - virtual_qstart
                effective_qseq = qseq_corrected[padding:]
                effective_sseq = sseq_corrected[padding:]
            else:
                effective_qseq = qseq_corrected
                effective_sseq = sseq_corrected

            matches = sum(
                1 for (qq, ss) in zip(effective_qseq, effective_sseq)
                if qq == ss and qq != "-" and ss != "-"
            )
            mismatches_raw = sum(
                1 for (qq, ss) in zip(effective_qseq, effective_sseq)
                if qq != ss and qq != "-" and ss != "-"
            )

            identity_pct = 0
            if (matches + mismatches_raw) > 0:
                identity_pct = matches / (matches + mismatches_raw) * 100
                identity_pct = round(identity_pct, 2)

            # Gather mutation strings
            mutations = []
            for i in range(len(sseq_corrected)):
                q_res = qseq_corrected[i]
                s_res = sseq_corrected[i]
                if q_res != s_res and q_res != "-" and s_res != "-":
                    mut_position = virtual_qstart + i
                    mutations.append(f"{s_res}{mut_position}{q_res}")

            mutation_count = len(mutations)
            mutations_str = "-".join(mutations) if mutation_count else "No mutations"

            alignment_score = aln.score

            # Keep top scoring alignment for (query, ref)
            if alignment_score > local_best_score:
                local_best_score = alignment_score
                local_best_data = (
                    query_id,
                    ref_id,
                    identity_pct,
                    len(aln.query),   
                    mutation_count,
                    virtual_qstart,
                    virtual_qend,
                    1,               
                    len(ref_seq),    
                    qseq_corrected,
                    sseq_corrected,
                    mutations_str,
                )

        # Compare best alignment for this ref to best overall for this query
        if local_best_data and local_best_score > best_score_for_query:
            best_score_for_query = local_best_score
            best_alignment_tuple = local_best_data

    return best_alignment_tuple

def main():
    parser = argparse.ArgumentParser(description="Align sequences and compute mutations (semi-global) in parallel.")
    parser.add_argument("--query", type=str, required=True, help="Path to the query FASTA file")
    parser.add_argument("--ref", type=str, required=True, help="Path to the reference AA FASTA file")
    parser.add_argument("--output", type=str, required=True, help="Path for output")
    parser.add_argument("--threads", type=int, default=4, help="Number of parallel worker processes")
    args = parser.parse_args()

    # Read FASTA
    query_sequences = {rec.id: str(rec.seq) for rec in SeqIO.parse(args.query, "fasta")}
    reference_sequences = {rec.id: str(rec.seq) for rec in SeqIO.parse(args.ref, "fasta")}

    # semi-global alignment settings
    aligner_settings = {
        "mode": "global",
        "match_score": 2,
        "mismatch_score": -1,
        "open_gap_score": -10,
        "extend_gap_score": -0.5,

        # Zero out all end-gap open/extend scores for both query and target
        "query_left_open_gap_score":   0,
        "query_left_extend_gap_score": 0,
        "query_right_open_gap_score":  0,
        "query_right_extend_gap_score":0,
        "target_left_open_gap_score":  0,
        "target_left_extend_gap_score":0,
        "target_right_open_gap_score": 0,
        "target_right_extend_gap_score":0,
    }

    results = []

    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = []
        for q_id, q_seq in query_sequences.items():
            fut = executor.submit(
                align_one_query,
                q_id, q_seq,
                reference_sequences,
                aligner_settings
            )
            futures.append(fut)

        for fut in futures:
            alignment_tuple = fut.result()  
            if alignment_tuple:
                results.append(alignment_tuple)

    columns = [
        "qseqid", "sseqid", "pident", "length", "mismatch",
        "qstart", "qend", "sstart", "send",
        "qseq", "sseq", "mutations"
    ]
    df_results = pd.DataFrame(results, columns=columns)
    df_results.to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
    main()
