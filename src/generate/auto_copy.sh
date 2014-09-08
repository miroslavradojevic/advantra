# script makes several folders with prefix (happens with synth data)
# list dirs using prefix and store run.txt in each (run.txt has folder name)
# list dirs using prefix and copy det.ijm in each (det.ijm is alway the same)

# create folders
rm -rf PAT.*	
mkdir PAT.1
mkdir PAT.2
mkdir PAT.3
mkdir PAT.4

# create det.ijm
rm -rf det.ijm # del just in case
echo "bullcrap" >> det.ijm




DIR=$(cd $(dirname "$0"); pwd)

echo ""
echo "current dir: "
echo $DIR
echo ""

# list dirs
find $DIR -type d -iname '*PAT*' -exec sh -c 'echo $0' {} \;

#DIR=$0; echo $DIR >> $DIR"/run.txt";

# list dirs using prefix and store run.txt in each (run.txt has folder name)
find $DIR -type d -iname '*PAT*' -exec sh -c 'DIR_FOUND=$0; echo "example command; run_command -arg "$DIR_FOUND"/det.ijm" $DIR_FOUND >> $DIR_FOUND/run.txt; ' {} \;

# list dirs using prefix and copy det.ijm in each (det.ijm is alway the same)
find $DIR -type d -iname '*PAT*' -exec sh -c 'DIR_FOUND=$0; cp det.ijm $DIR_FOUND;' {} \;
# find $DIR -type d -iname '*PAT*' -exec sh -c 'cp det.ijm $0;' {} \;

# example on how to generate string that would start all the scripts
#   find . -iname '*.txt' -exec echo "start script command "{}";" \;
# more examples
# delete file type
# find /scratch/mradojevic/critpoint_tests/synth_ext/ -type f -iname '*run.sh*' -exec rm -rf {} \;
# print file type
# find /scratch/mradojevic/critpoint_tests/synth_ext/ -type f -iname '*run.sh*' -exec echo {} \;
# generate file type
# find /scratch/mradojevic/critpoint_tests/synth_ext/ -type d -iname '*SYNGEN*' -exec sh -c 'DIR_FOUND=$0; echo "module load java;export _JAVA_OPTIONS=\"-Xms4G -Xmx6G\";/scratch/mradojevic/critpoint_tests/fiji/ImageJ-linux64 --headless -batch "$DIR_FOUND"/detect.ijm "$DIR_FOUND"/" >> $DIR_FOUND/run.sh; ' {} \;
# generate start job commands on cluster
# find /scratch/mradojevic/critpoint_tests/synth_ext/ -iname 'run.sh' -exec echo "qsub -q day -l h_vmem=20G  "{}";" \; >> start_jobs.txt

