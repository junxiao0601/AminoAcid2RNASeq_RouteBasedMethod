# AminoAcid2RNASeq_RouteBasedMethod
Hi,

**Problem**: 
    An amino acid can come from one or more RNA codons, meaning one amino acid sequence can be coded by one or more RNA sequence(s).

**Mission**: 
   These codes are to convertan  amino acid sequence (alphabet sequence) to a list of possible RNA sequence(s) coding that amino acid sequence.


**How to Use**: Download the whole folder and run the code in either Main.py or AminoAcid2RNASeq.py in the downloaded folder.
   *Stop = "*" in the sequence

   *These codes output a list of possible RNA sequence(s) in a csv.
   *Evaluation: "correctness" = whether the output RNA sequence can return the input amino acid sequence.


**Method**:
   Inductive overlapping of codons between two adjacent amino acids.


**Contributor**: Jun Xiao
