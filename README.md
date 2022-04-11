# Assembly fasta reads
Iterative algorithm for assembling fasta reads.

# Algorithm

Algorithm in each iteration is going through all loaded sequences, finding the longest prefix or suffix, then merging them in place.

First run is performed finding exact -fix of specified minimum length without any miss-matches.

In every next iteration, similarity function is used to finding the longest -fix, loosening similarity threshold at the same time in each run.

# Usage 

Script requires only Python >= 3.8 and does not use any external libraries.
