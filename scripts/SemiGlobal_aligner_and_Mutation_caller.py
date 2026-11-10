#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio.Align import PairwiseAligner
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

WORKER_REFERENCES = None
WORKER_ALIGNER = None

def build_aligner():
    
    aligner = PairwiseAligner()

    aligner.mode = "global"

    aligner.match_score = 2
    aligner.mismatch_score = -1

    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    aligner.query_left_open_gap_score = 0
    aligner.query_left_extend_gap_score = 0

    aligner.query_right_open_gap_score = 0
    aligner.query_right_extend_gap_score = 0

    aligner.target_left_open_gap_score = 0
    aligner.target_left_extend_gap_score = 0

    aligner.target_right_open_gap_score = 0
    aligner.target_right_extend_gap_score = 0

    return aligner


def initialize_worker(reference_sequences):
    """
    Called once when each worker starts.

    This prevents us from passing the full reference dictionary
    into every individual query job.
    """

    global WORKER_REFERENCES
    global WORKER_ALIGNER

    WORKER_REFERENCES = reference_sequences
    WORKER_ALIGNER = build_aligner()


def align_one_query(query_item):
    """
    For one query:

        1. SCORE query against all germlines.
        2. Identify the best-scoring germline.
        3. Perform ONE full alignment against that germline.
        4. Identity mutations.
    """

    query_id, query_seq = query_item

    global WORKER_REFERENCES
    global WORKER_ALIGNER


    # ========================================================
    # STEP 1:
    # Find best germline using score() only
    # ========================================================

    best_score = float("-inf")

    best_ref_id = None
    best_ref_seq = None


    for ref_id, ref_seq in WORKER_REFERENCES:

        score = WORKER_ALIGNER.score(
            query_seq,
            ref_seq
        )

        # IMPORTANT:
        # Use >, NOT >=.
        # if two references tie, keep the first one encountered.

        if score > best_score:

            best_score = score
            best_ref_id = ref_id
            best_ref_seq = ref_seq


    # No usable reference
    if best_ref_id is None:
        return None


    # ========================================================
    # STEP 2:
    # Full alignment ONLY against the best germline
    # ========================================================

    alignments = WORKER_ALIGNER.align(
        query_seq,
        best_ref_seq
    )

    if len(alignments) == 0:
        return None

    aln = alignments[0]

    q_blocks = aln.aligned[0]
    r_blocks = aln.aligned[1]


    if len(q_blocks) == 0 or len(r_blocks) == 0:
        return None

    q_start_index = q_blocks[0][0]
    q_end_index = q_blocks[-1][1]

    s_start_index = r_blocks[0][0]
    s_end_index = r_blocks[-1][1]

    q_start = q_start_index + 1
    s_start = s_start_index + 1


    # ========================================================
    # Reconstruct query position relative to full
    # germline sequence
    # ========================================================

    virtual_qstart = q_start - (s_start - 1)

    virtual_qend = (
        virtual_qstart +
        len(best_ref_seq) -
        1
    )


    # ========================================================
    # Build corrected query sequence
    # ========================================================

    if virtual_qstart < 1:

        num_gaps_left = 1 - virtual_qstart

        padded_query = (
            "-" * num_gaps_left +
            query_seq
        )

        qseq_corrected = padded_query[
            :len(best_ref_seq)
        ]

    else:

        qseq_corrected = query_seq[
            virtual_qstart - 1:
            virtual_qstart - 1 + len(best_ref_seq)
        ]


    # Pad right side if query is shorter than reference
    if len(qseq_corrected) < len(best_ref_seq):

        qseq_corrected += (
            "-" *
            (len(best_ref_seq) - len(qseq_corrected))
        )


    sseq_corrected = best_ref_seq

    if virtual_qstart < 1:

        padding = 1 - virtual_qstart

        effective_qseq = qseq_corrected[
            padding:
        ]

        effective_sseq = sseq_corrected[
            padding:
        ]

    else:

        effective_qseq = qseq_corrected
        effective_sseq = sseq_corrected


    matches = sum(

        1
        for qq, ss in zip(
            effective_qseq,
            effective_sseq
        )

        if (
            qq == ss
            and qq != "-"
            and ss != "-"
        )
    )


    mismatches_raw = sum(

        1
        for qq, ss in zip(
            effective_qseq,
            effective_sseq
        )

        if (
            qq != ss
            and qq != "-"
            and ss != "-"
        )
    )


    if (matches + mismatches_raw) > 0:

        identity_pct = (
            matches /
            (matches + mismatches_raw) *
            100
        )

        identity_pct = round(
            identity_pct,
            2
        )

    else:

        identity_pct = 0

    # ========================================================
    # Mutation calls
    # ========================================================

    mutations = []


    for i in range(len(sseq_corrected)):

        q_res = qseq_corrected[i]
        s_res = sseq_corrected[i]


        if (
            q_res != s_res
            and q_res != "-"
            and s_res != "-"
        ):

            mut_position = (
                virtual_qstart + i
            )

            mutations.append(
                f"{s_res}{mut_position}{q_res}"
            )


    mutation_count = len(mutations)


    if mutation_count > 0:

        mutations_str = "-".join(
            mutations
        )

    else:

        mutations_str = "No mutations"


    return (

        query_id,

        best_ref_id,

        identity_pct,

        len(aln.query),

        mutation_count,

        virtual_qstart,

        virtual_qend,

        1,

        len(best_ref_seq),

        qseq_corrected,

        sseq_corrected,

        mutations_str
    )

# ============================================================   
# Main
# ============================================================   

def main():

    parser = argparse.ArgumentParser(

        description=(
            "Optimized semi-global germline alignment "
            "and mutation caller for SHMART."
        )
    )


    parser.add_argument(
        "--query",
        required=True,
        help="Query amino-acid FASTA"
    )


    parser.add_argument(
        "--ref",
        required=True,
        help="Reference germline amino-acid FASTA"
    )


    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV file"
    )


    parser.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Number of worker processes (default: 8)"
    )


    parser.add_argument(
        "--chunksize",
        type=int,
        default=50,
        help=(
            "Number of queries sent to each worker batch "
            "(default: 50)"
        )
    )


    args = parser.parse_args()


    # ========================================================
    # Read references
    # ========================================================

    print("Reading reference germlines...")


    reference_sequences = [

        (
            rec.id,
            str(rec.seq)
        )

        for rec in SeqIO.parse(
            args.ref,
            "fasta"
        )
    ]


    if len(reference_sequences) == 0:

        raise ValueError(
            f"No reference sequences found in {args.ref}"
        )


    print(
        f"Reference germlines: "
        f"{len(reference_sequences)}"
    )


    # ========================================================
    # Read query sequences
    # ========================================================

    print("Reading query sequences...")


    query_sequences = [

        (
            rec.id,
            str(rec.seq)
        )

        for rec in SeqIO.parse(
            args.query,
            "fasta"
        )
    ]


    if len(query_sequences) == 0:

        raise ValueError(
            f"No query sequences found in {args.query}"
        )


    print(
        f"Query sequences: "
        f"{len(query_sequences)}"
    )

    print(
        f"Workers: "
        f"{args.threads}"
    )

    print(
        f"Chunksize: "
        f"{args.chunksize}"
    )


    # ========================================================
    # Run alignments
    # ========================================================

    columns = [

        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "qstart",
        "qend",
        "sstart",
        "send",
        "qseq",
        "sseq",
        "mutations"
    ]


    results = []


    print(
        "\nStarting germline alignment..."
    )


    with ProcessPoolExecutor(

        max_workers=args.threads,

        initializer=initialize_worker,

        initargs=(
            reference_sequences,
        )

    ) as executor:


        mapped_results = executor.map(

            align_one_query,

            query_sequences,

            chunksize=args.chunksize
        )


        for i, alignment_tuple in enumerate(
            mapped_results,
            start=1
        ):

            if alignment_tuple is not None:

                results.append(
                    alignment_tuple
                )


            # Progress update
            if i % 500 == 0:

                print(
                    f"Processed "
                    f"{i}/{len(query_sequences)} "
                    f"queries",
                    flush=True
                )


    # ========================================================
    # Output
    # ========================================================

    print(
        "\nCreating output table..."
    )


    df_results = pd.DataFrame(

        results,

        columns=columns
    )


    df_results.to_csv(

        args.output,

        sep="\t",

        index=False
    )


    print(
        f"Completed queries: "
        f"{len(results)}"
    )

    print(
        f"Output: "
        f"{args.output}"
    )


# ============================================================
# Run
# ============================================================

if __name__ == "__main__":
    main()
