PhyloHAP

Computer programmes used for manuscript: Can parallel ecological speciation be detected with phylogenetic analyses? by Noelia Pérez-Pereira, Humberto Quesada and Armando Caballero

FILES:

PhyloHAP.c : General model simulation programme
PhyloHAPsym.c : Littorina model for multiple origin in sympatry
PhyloHAPallo.c : Littorina model for single origin in allopatry
PhyloDis.c : Programme to handle with distances and trees
genlib.c : Genetic routines
libhdr : Declaration of external variables for genlib.c
OTHER NECESSARY PROGRAMMES: dnadist, gendist and fitch from PHYLIP package
OTHER NECESSARY FILE: seedfile (file with initial seed for random number generator)

COMPILATION:

time cc -o PhyloHAP -w -O PhyloHAP.c genlib.c -lm
time cc -o PhyloHAPsym -w -O PhyloHAPsym.c genlib.c -lm
time cc -o PhyloHAPallo -w -O PhyloHAPallo.c genlib.c -lm
time cc -o PhyloDis -w -O PhyloDis.c genlib.c -lm

RUNNING

General model
qsub script_PhyloHAP1.sh 0.0005 0.0005 50000 10000 1000 16 1 100 500 &

Littorina model
qsub script_PhyloHAPalloA.sh 0.000037 0.000036 50000 10000 100 16 0 22 20000 &
qsub script_PhyloHAPsymA.sh 0.000037 0.000036 50000 10000 100 16 0 22 20000 &
