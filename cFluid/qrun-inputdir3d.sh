#!/bin/bash

#This script submits a job for each of the input files
#that are prefixed by "in-" in the directory given as
#the first argument with the number of cores specified
#by the second argument.

#Must be run in the directory of the executable

#example:
#           ./qrun-inputdir3d.sh input-rt3d 2 2 2
#
#       Will run the executable for every input file
#       prefixed with "in-" in the directory input-rt3d
#       with 2 x partitions, 2 y partitions, and 2 z partitions
#       with 2*2*2 = 8 total cores.
        

#executable
EXENAME="cFluid"

#input directory 
INDIR=${1%/}

#domain partition
PX=$2
PY=$3
PZ=$4

#total cores
NP=$((${PX}*${PY}*${PZ}))

if [ "$#" -ne 4 ]; then
    echo "Error: incorrect number of input arguments"
    exit 1;
elif [ ! -d "$INDIR" ]; then
    echo "Error: input directory does not exist"
    exit 2
elif [ $PX -le 0 -o $PY -le 0 -o $PZ -le 0 ]; then
    echo "Error: each partition must be specified by a positive integer"
    exit 3
fi


DIR=`pwd`
EXE="$DIR/$EXENAME"
QRUNDIR="$INDIR/qrun"
PARENTOUTDIR="$DIR/out-$INDIR-p${PX}${PY}${PZ}"

echo -e "path to executable:\n\t $EXE\n"
echo -e "directory containing the input files:\n\t $INDIR\n"
echo -e "directory where qrun files will be generated:\n\t $QRUNDIR\n"
echo -e "directory to write output of each run:\n\t $PARENTOUTDIR\n"
mkdir -p $QRUNDIR
mkdir -p $PARENTOUTDIR

QRUNHEADER="$QRUNDIR/qrunheader.txt"
echo '#!/bin/bash' > $QRUNHEADER
echo '#$ -q all.q' >> $QRUNHEADER
echo '#$ -R y' >> $QRUNHEADER
echo '#$ -pe mpi' "$NP" >> $QRUNHEADER
echo -e '#$ -cwd\n' >> $QRUNHEADER
echo 'source /etc/profile.d/modules.sh' >> $QRUNHEADER
echo -e 'source ~/.bashrc\n' >> $QRUNHEADER

for infile in "$INDIR"/in-* ; do
    infile=`basename "$infile"`
    outdir="$PARENTOUTDIR/out-$infile"; mkdir -p $outdir
    qrunfile="$QRUNDIR/qrun-$infile"
    jobname="$EXENAME-$infile-np$NP"
    eval cat $QRUNHEADER > $qrunfile
    echo -e "\nmpirun -np $NP $EXE -d 3 -p $PX $PY $PZ -i $INDIR/$infile -o $outdir" >> $qrunfile
    chmod 700 $qrunfile
    eval qsub -N $jobname -o $outdir/qsubout -e $outdir/qsuberror $qrunfile
done



