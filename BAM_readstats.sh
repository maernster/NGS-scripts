#!bin/bash

#this is a bash script that calculates the read number and average read length of individual BAM files
#Individuals in question need to be supplemented as arguments (i.e. bash bamstats.sh indivA indivB indivC ...). Without the ".bam" at the end.
#BAM files have to be stored in the same directory as this script and should be named after the individuals (i.e. indivA.bam indivB.bam indivC.bam ...)
#stats are stored in a stats.csv file

#FOR WORKING ON UPPMAX, THE FOLLOWING MODULES NEED TO BE LOADED:
#module load bioinfo-tools bwa/0.7.17 samtools/1.10

stats_file='stats.csv'
rm -f $stats_file
touch $stats_file
echo "library,number_reads,mean_read_length,mean_cov,std_cov" >> $stats_file
for bamfile in "$@"
do
  TOTAL_NUMBER_READS=$(samtools view -c ${bamfile}.bam)
  TOTAL_READ_LENGTH=$(samtools view ${bamfile}.bam | awk '{print length($10)}' | awk '{ sum += $1 } END { print sum }')
  AVERAGE_READ_LENGTH=$(echo $((TOTAL_READ_LENGTH/TOTAL_NUMBER_READS)))
  TOTAL_COVERAGE=$(samtools depth -aa -Q 20 ${bamfile}.bam | awk '{print $3}' | awk '{ sum += $1 } END { print sum }')
  TOTAL_NUMBER_SITES=$(samtools depth ${bamfile}.bam | wc -l)
  AVERAGE_READ_DEPTH=$(echo $((TOTAL_COVERAGE/TOTAL_NUMBER_SITES)))
  STD_READ_DEPTH=$(samtools depth -aa ${bamfile}.bam | awk '{print $3}' | awk '{sum+=$1; sumsq+=$1*$1}END{print sqrt(sumsq/NR - (sum/NR)**2)}')
  str=$bamfile","$TOTAL_NUMBER_READS","$AVERAGE_READ_LENGTH","$AVERAGE_READ_DEPTH","$STD_READ_DEPTH
  echo $str >> $stats_file
done
