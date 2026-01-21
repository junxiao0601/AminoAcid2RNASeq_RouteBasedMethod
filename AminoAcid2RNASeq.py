import pandas as pd

# Mapping from amino acids to their corresponding RNA codons
amino_acid_to_codons = {
    'A': ['GCU', 'GCC', 'GCA', 'GCG'],
    'C': ['UGU', 'UGC'],
    'D': ['GAU', 'GAC'],
    'E': ['GAA', 'GAG'],
    'F': ['UUU', 'UUC'],
    'G': ['GGU', 'GGC', 'GGA', 'GGG'],
    'H': ['CAU', 'CAC'],
    'I': ['AUU', 'AUC', 'AUA'],
    'K': ['AAA', 'AAG'],
    'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
    'M': ['AUG'],
    'N': ['AAU', 'AAC'],
    'P': ['CCU', 'CCC', 'CCA', 'CCG'],
    'Q': ['CAA', 'CAG'],
    'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
    'T': ['ACU', 'ACC', 'ACA', 'ACG'],
    'V': ['GUU', 'GUC', 'GUA', 'GUG'],
    'W': ['UGG'],
    'Y': ['UAU', 'UAC'],
    '*': ['UAA', 'UAG', 'UGA']  # Stop codons
}

condons_to_amino_acid = {
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'UGU': 'C', 'UGC': 'C',
    'GAU': 'D', 'GAC': 'D',
    'GAA': 'E', 'GAG': 'E',
    'UUU': 'F', 'UUC': 'F',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    'CAU': 'H', 'CAC': 'H',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I',
    'AAA': 'K', 'AAG': 'K',
    'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'AUG': 'M',
    'AAU': 'N', 'AAC': 'N',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAA': 'Q', 'CAG': 'Q',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 
    'AGA': 'R', 'AGG': 'R',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'AGU': 'S', 'AGC': 'S',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'UGG': 'W',
    'UAU': 'Y', 'UAC': 'Y',
    'UAA': '*', 'UAG': '*', 'UGA': '*'
}


def amino_acid_list_to_rna_seq_list(amino_acid_seq):
    """
    Convert an amino acid sequence to a list of possible RNA sequences.

    Args:
        amino_acid_seq (str): A string representing the amino acid sequence.
    Returns:
        list: A list of strings representing possible RNA sequences.
    """
    rna_seq_list = []
    for i in range(len(amino_acid_seq)):
        aa = amino_acid_seq[i]
        condons = amino_acid_to_codons[aa]

        if i == 0:
            rna_seq_list = condons
        else:
            candidate_rna_seq_list = []
            check_not_match = 0
            for rna_seq in rna_seq_list:
                for condon in condons:
                    if rna_seq[-2:] == condon[:2]:
                        candidate_rna_seq_list.append(rna_seq + condon[2])
                    else:
                        check_not_match += 1
            if check_not_match == len(rna_seq_list) * len(condons):
                print(f"Warning: There is no possible overlap between amino acid {amino_acid_seq[i-1]} at position {i-1} and {aa} at position {i}.")
            rna_seq_list = candidate_rna_seq_list

    return rna_seq_list


def from_rna_to_aa_seq(rna_seq):
    """
    Convert an RNA sequence to its corresponding amino acid sequence.

    Args:
        rna_seq (str): A string representing the RNA sequence.
    Returns:
        str: A string representing the corresponding amino acid sequence.
    """
    aa_seq = ""
    for i in range(len(rna_seq)-2):
        codon = rna_seq[i:i+3]
        aa = condons_to_amino_acid[codon]
        aa_seq += aa

    return aa_seq


def evaluate_conversion(rna_seq_list, original_aa_seq, index):
    """
    Evaluate the conversion from RNA sequences to amino acid sequences.

    Args:
        rna_seq_list (list): A list of RNA sequences.
        original_aa_seq (str): The original amino acid sequence.
    Returns:
        dict: A dictionary with RNA sequences as keys and boolean values indicating correctness.
    """
    evaluation_result = []
    for rna_seq in rna_seq_list:
        converted_aa_seq = from_rna_to_aa_seq(rna_seq)
        if (converted_aa_seq == original_aa_seq):
            evaluation_result.append([rna_seq, True])
    if evaluation_result == []:
        print("No valid RNA sequences can be generated. Sorry!")
    else:
        evaluation_result = pd.DataFrame(evaluation_result, columns=['RNA_Sequence', 'Correctness'])
        csv_name = f'AminoAcid2RNASeq_Evaluation_{index}.csv'
        evaluation_result.to_csv(csv_name, index=False)
        print(f"The RNA sequences with their entropy and score have been saved to 'AminoAcid2RNASeq_Output_{index}.csv'.")

    return evaluation_result


def main():
    i = 1
    while True:
        test_seq = input("Enter an amino acid sequence (type 'q' to exit): ").strip().upper()
        if test_seq.lower() == 'q':
            break
        rna_seqs = amino_acid_list_to_rna_seq_list(test_seq)
        rna_seqs_table = pd.DataFrame(rna_seqs, columns=['RNA_Sequence'])
        print("------------------------------------------------------------------------------")
        print("Input Amino Acid Sequence:", test_seq)
        print(f"Generated {len(rna_seqs)} RNA sequences: Here are the complete list:")
        print(rna_seqs_table)
        evaluate_conversion(rna_seqs, test_seq, i)
        if rna_seqs != []:
            i += 1
        print("------------------------------------------------------------------------------")


if __name__ == "__main__":
    main()