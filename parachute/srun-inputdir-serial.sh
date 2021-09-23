#!/bin/bash

DIR=$(pwd)
EXENAME="parachute"
EXE="$DIR/$EXENAME"

if [ "$#" -ne 1 ]; then
    echo "Error: must provide input directory name as argument"
    exit -1;
fi

#input directory
INDIR=${1%/}

if [ ! -d "$INDIR" ]; then
    echo "Error: input directory does not exist"
    exit -1
fi

SRUNDIR="$INDIR/srun-serial"
PARENTOUTDIR="$DIR/out-$INDIR-serial"

echo -e "path to executable:\n\t $EXE\n"
echo -e "directory containing the input files:\n\t $INDIR\n"
echo -e "directory where qrun files will be generated:\n\t $SRUNDIR\n"
echo -e "directory to write output of each run:\n\t $PARENTOUTDIR\n"
mkdir -p $SRUNDIR
mkdir -p $PARENTOUTDIR

SRUNHEADER="$SRUNDIR/srunheader.txt"
echo -e '#!/bin/bash\n' > $SRUNHEADER

echo '#SBATCH --ntasks=1                     # Run on single cpu' >> $SRUNHEADER
echo -e '#SBATCH --mem-per-cpu=6G   # memory per cpu-core\n' >> $SRUNHEADER

echo 'module load mvapich2/2.3.3-gcc-8.3.1' >> $SRUNHEADER
echo '#export LD_LIBRARY_PATH="/opt/sundials/sundials-2.7.0-opt/lib:/opt/petsc/petsc-3.13.4-opt/lib:/usr/lib64:$LD_LIBRARY_PATH"' >> $SRUNHEADER


for infile in "$INDIR"/in-* ; do
    infile=`basename "$infile"`
    outdir="$PARENTOUTDIR/out-$infile"; mkdir -p $outdir
    srunfile="$SRUNDIR/srun-$infile"
    jobname="$EXENAME-$infile"
    cat $SRUNHEADER > $srunfile
    echo -e "\n$EXE -d 3 -i $INDIR/$infile -o $outdir" >> $srunfile
    chmod 700 $srunfile
    sbatch -J $jobname -o $outdir/slurm.log $srunfile
done

