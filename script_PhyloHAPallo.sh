#!/bin/bash
#$ -cwd

rm script_PhyloHAPallo.sh.*

#script_full_reps_new.sh <MECO> <MLOC> <GENS> <TGEN> <REPS> <NLOCI> <TYPE> <POP> <Ne>
#To be used only in a machine with /scratch directory

#Check number of arguments
if [ $# -ne 9 ]  
then
	echo "Usage: $0 <MECO> <MLOC> <GENS> <TGEN> <REPS> <NLOCI> <TYPE> <POP> <Ne>" 
	exit 1
fi

#Set arguments
MECO=$1
MLOC=$2
GENS=$3
TGEN=$4
REPS=$5
NLOCI=$6
TYPE=$7
POP=$8
Ne=$9

#Working directory
WDIR=$PWD 

#Output directory
mkdir -p $WDIR/A11_$MECO.$MLOC.$TYPE.$POP.$Ne 

#Scratch directory
mkdir -p /state/partition1/scratch11/$SLURM_JOBID/

#Copy all files in scratch directory
cp PhyloHAPallo /state/partition1/scratch11/$SLURM_JOBID/
cp PhyloDis /state/partition1/scratch11/$SLURM_JOBID/
cp seedfile /state/partition1/scratch11/$SLURM_JOBID/
cp gendist /state/partition1/scratch11/$SLURM_JOBID
cp fitch /state/partition1/scratch11/$SLURM_JOBID
cp dnadist /state/partition1/scratch11/$SLURM_JOBID

#File with information of node and directory
touch $WDIR/$SLURM_JOBID.`hostname`.`date +%HH%MM`

#Move to scratch directory
cd /state/partition1/scratch11/$SLURM_JOBID

#Run the programme

START=$(date +%s)

#################### START OF REPLICATES ####################

for ((a=1; a<=$REPS; a++))
do

time ./PhyloHAPallo>>out<<@
0
-99
$POP	pop
$Ne	nind
$MECO	mrECO
$MLOC	mrLOC
10	NCRO
$NLOCI	 NLOCI
0.00000002	u
0.0	L (morgans)
0.99	s
100000	PREGEN
$GENS	GEN
$TGEN	TREEGEN
100	NSAM
$TYPE	INITYPE (0=random, 1=habitat, 2=locality)
@

#################### START TREES FOR GENERATIONS AND REPLICATES ####################

mv phylfreqfile1.txt phylfreqfile0.txt

g=0

for x in ./phylfreqfile*.txt
do

cp $x phylfreqfile.txt

#################### FREQ TREE ####################

./gendist <<EOF
phylfreqfile.txt
y
EOF

mv outfile distances
cat distances >> DISTANCES_FREQ_$g
rm distances

g=$[g + TGEN]

done


#################### HAPLO TREE ####################

mv phylfile1.txt phylfile0.txt

g=0

for x in ./phylfile*.txt
do

c=1
d=5

#### b<=NSAM or less
for ((b=1; b<=10; b++))
do 

echo $c $d >> kk

sed -n "$c,${d}p" < $x > phylfile.txt

./dnadist <<EOF
phylfile.txt
D
D
I
y
EOF

mv outfile distances
cat distances >> DISTANCES_HAPLO_$g
rm distances

c=$[c + 5]
d=$[d + 5]

done

g=$[g + TGEN]

done

echo "Rep  $a" >> repfile

#################### RUN PhyloDis ####################

time ./PhyloDis>>out<<@
0
4	pop
$GENS	GEN
$TGEN	TREEGEN
10	NSAM (haplotypes)
$a	REPLICATES
@

cp -r /state/partition1/scratch11/$SLURM_JOBID/TREEFILE $WDIR/A11_$MECO.$MLOC.$TYPE.$POP.$Ne 
cp -r /state/partition1/scratch11/$SLURM_JOBID/repfile $WDIR/A11_$MECO.$MLOC.$TYPE.$POP.$Ne 
cp -r /state/partition1/scratch11/$SLURM_JOBID/datafile* $WDIR/A11_$MECO.$MLOC.$TYPE.$POP.$Ne 

#################### END OF REPLICATES ####################
done

############### COPY FILES INTO HOME DIRECTORY ################

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "Run took  $DIFF seconds" >> timefile


cp -r /state/partition1/scratch11/$SLURM_JOBID/datafile* $WDIR/A11_$MECO.$MLOC.$TYPE.$POP.$Ne 
# cp -r /state/partition1/scratch11/$SLURM_JOBID/phylfile* $WDIR/A11_$MECO.$MLOC.$TYPE.$POP.$Ne 
# cp -r /state/partition1/scratch11/$SLURM_JOBID/phylfreqfile* $WDIR/A11_$MECO.$MLOC.$TYPE.$POP.$Ne 

rm DISTANCES_FREQ_$[GENS + TGEN]
rm DISTANCES_HAPLO_$[GENS + TGEN]

# cp -r /state/partition1/scratch11/$SLURM_JOBID/tree* $WDIR/A11_$MECO.$MLOC.$TYPE.$POP.$Ne 
# cp -r /state/partition1/scratch11/$SLURM_JOBID/DISTANCES* $WDIR/A11_$MECO.$MLOC.$TYPE.$POP.$Ne 
cp -r /state/partition1/scratch11/$SLURM_JOBID/TREEFILE $WDIR/A11_$MECO.$MLOC.$TYPE.$POP.$Ne 
cp -r /state/partition1/scratch11/$SLURM_JOBID/timefile $WDIR/A11_$MECO.$MLOC.$TYPE.$POP.$Ne 
cp -r /state/partition1/scratch11/$SLURM_JOBID/out $WDIR/A11_$MECO.$MLOC.$TYPE.$POP.$Ne 
cp -r /state/partition1/scratch11/$SLURM_JOBID/salfile $WDIR/A11_$MECO.$MLOC.$TYPE.$POP.$Ne 
cp -r /state/partition1/scratch11/$SLURM_JOBID/seedfile $WDIR/

################### CLEANING OF SCRATCH ####################

rm -r /state/partition1/scratch11/$SLURM_JOBID/
rm $WDIR/$SLURM_JOBID.*
