# ParsePedMap4PairwiseSNPMatrix
The script will parse the ped and map file with group ids to generate pairwise SNP matrix.

# Note
If no groups are defined, provide list of ids followed by 1 (tab-delimited)

USAGE: perl ParsePedMap4PairwiseSNPMatrix.pl SNP.plink.ped SNP.plink.map group.txt

# Input
1) Plink ped and map files
2) Group file (group.txt)

group.txt

IS1004  1

IS10302 1

IS1212  1

IS1219  2

IS1233  1

IS12697 1

IS12937 2

IS12945 2

IS13294 2

IS13444 2

IS13549 1

IS13809 3

# Output
1) Genotype.fasta
2) inputfile.Matrix.txt

# Matrix file

  IS1004	IS10302	IS1212	IS1219

IS1004	0	6247	7668	7075

IS10302	6247	0	6514	6109

IS1212	7668	6514	0	7206

IS1219	7075	6109	7206	0
