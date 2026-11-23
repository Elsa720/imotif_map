"""Determine sequence strands and normalize FASTA sequences to C-rich orientation."""


import argparse

import pandas as pd
from Bio import SeqIO

def get_strand_and_unify(fa_input, bed_output, fa_output):
    records = list(SeqIO.parse(fa_input, "fasta"))
    unified_records = []
    bed_data = []

    for rec in records:
        seq_str = str(rec.seq).upper()
        c_count = seq_str.count('C')
        g_count = seq_str.count('G')
        
        # Parse coordinates from an ID such as chr1:100-200.
        header = rec.id
        chrom = header.split(':')[0]
        start, end = header.split(':')[1].split('-')
        
        # If G > C, the C-rich sequence is on the negative strand.
        if g_count > c_count:
            strand = "-"
            # Reverse-complement the sequence to make it C-rich.
            new_seq = rec.seq.reverse_complement()
        else:
            strand = "+"
            new_seq = rec.seq
            
        # Store the BED fields.
        bed_data.append([chrom, start, end, header, ".", strand])
        
        # Store the normalized FASTA record.
        rec.seq = new_seq
        rec.id = header  # Keep the ID unchanged for downstream matching.
        rec.description = f"strand={strand}"
        unified_records.append(rec)

    # Write the output files.
    pd.DataFrame(bed_data).to_csv(bed_output, sep='\t', index=False, header=False)
    SeqIO.write(unified_records, fa_output, "fasta")
    print(
        "Processing complete!\n"
        f"BED file saved to: {bed_output}\n"
        f"Normalized FASTA file saved to: {fa_output}"
    )


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Determine the strand of each FASTA sequence and normalize it "
            "to the C-rich orientation."
        )
    )
    parser.add_argument(
        "fa_input",
        help="Input FASTA file (sequence IDs must use the chr1:100-200 format)",
    )
    parser.add_argument("bed_output", help="Output BED file")
    parser.add_argument("fa_output", help="Output normalized FASTA file")
    return parser.parse_args()


def main():
    args = parse_args()
    get_strand_and_unify(args.fa_input, args.bed_output, args.fa_output)


if __name__ == "__main__":
    main()



