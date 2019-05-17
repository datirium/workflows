#!/bin/bash

echo "StepI: clean input"

cut -f 6,7,8 $1 | grep -v "start\|chrM" > peak1.bed
cut -f 6,7,8 $2 | grep -v "start\|chrM" > peak2.bed

bamToBed -i $3 |sed 's/\s$//g' | awk -v var=$5 'BEGIN {OFS="\t"}
     {if ($1 !="chrM" && $6=="+" && $1 !~/random/  && $2>0 && $3>0)
          print $1,$2+var,$3+var>"read1.bed";
      else if ($1 !="chrM" && $6=="-" && $1 !~/random/   && $2>var && $3>var)
          print $1,$2-var,$3-var>"read1.bed";
      else 
          print $0 > "/dev/null"}' &
bamToBed -i $4 |sed 's/\s$//g' | awk -v var=$6 'BEGIN {OFS="\t"}
     {if ($1 !="chrM" && $6=="+" && $1 !~/random/   && $2>0 && $3>0)
          print $1,$2+var,$3+var>"read2.bed";
      else if ($1 !="chrM" && $6=="-" && $1 !~/random/   && $2>var && $3>var)
          print $1,$2-var,$3-var>"read2.bed";
      else 
          print $0 > "/dev/null"}' &

wait

echo "StepII: classify common or unique peaks"

intersectBed -a peak1.bed -b peak2.bed -u | sort -k1,1 -k2,2n -k3,3n > common_peak1.bed &
intersectBed -a peak2.bed -b peak1.bed -u | sort -k1,1 -k2,2n -k3,3n > common_peak2.bed &
intersectBed -a peak1.bed -b peak2.bed -v | sort -k1,1 -k2,2n -k3,3n > unique_peak1.bed &
intersectBed -a peak2.bed -b peak1.bed -v | sort -k1,1 -k2,2n -k3,3n > unique_peak2.bed &

wait 

cat common_peak1.bed common_peak2.bed |sort -k1,1 -k2,2n | mergeBed -i - > common_peak.bed

echo "StepIII: count peak read"

if [ -f MAnorm.bed ];
then
rm MAnorm.bed
fi
coverageBed -a unique_peak1.bed -b read1.bed | sort -k1,1 -k2,2n -k3,3n  | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"unique_peak1" >> "MAnorm.bed"; print $4 > "unique_peak1_count_read1"}' &
coverageBed -a unique_peak1.bed -b read2.bed | sort -k1,1 -k2,2n -k3,3n  | awk '{print $4 > "unique_peak1_count_read2"}' &
wait 
coverageBed -a common_peak1.bed -b read1.bed | sort -k1,1 -k2,2n -k3,3n  | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"common_peak1" >> "MAnorm.bed";print $4 > "common_peak1_count_read1"}' &
coverageBed -a common_peak1.bed -b read2.bed | sort -k1,1 -k2,2n -k3,3n  | awk '{print $4 > "common_peak1_count_read2"}' &
wait
coverageBed -a common_peak2.bed -b read1.bed | sort -k1,1 -k2,2n -k3,3n  | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"common_peak2"  >> "MAnorm.bed";print $4 > "common_peak2_count_read1"}' &
coverageBed -a common_peak2.bed -b read2.bed |sort -k1,1 -k2,2n -k3,3n  |  awk '{print $4 > "common_peak2_count_read2"}' &
wait
coverageBed -a unique_peak2.bed -b read1.bed | sort -k1,1 -k2,2n -k3,3n  | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"unique_peak2">> "MAnorm.bed";print $4 > "unique_peak2_count_read1"}' &
coverageBed -a unique_peak2.bed -b read2.bed | sort -k1,1 -k2,2n -k3,3n  | awk '{print $4 > "unique_peak2_count_read2"}' &

wait 

cat common_peak1_count_read1 common_peak2_count_read1 > common_peak_count_read1
cat common_peak1_count_read2 common_peak2_count_read2 > common_peak_count_read2

if [ -f MAnorm_merge.bed ];
then
rm MAnorm_merge.bed
fi

cat  unique_peak1.bed | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"unique_peak1" >> "MAnorm_merge.bed"}'
coverageBed -a common_peak.bed -b read1.bed | sort -k1,1 -k2,2n -k3,3n  | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"merged_common_peak" >> "MAnorm_merge.bed"; print $4 > "merge_common_read1"}' &
coverageBed -a common_peak.bed -b read2.bed | sort -k1,1 -k2,2n -k3,3n  | awk '{print $4 > "merge_common_read2"}' &
wait
cat  unique_peak2.bed | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"unique_peak2" >> "MAnorm_merge.bed"}'

cat unique_peak1_count_read1 merge_common_read1  unique_peak2_count_read1 > merge_common_peak_count_read1
cat unique_peak1_count_read2 merge_common_read2  unique_peak2_count_read2 > merge_common_peak_count_read2

echo "SetpIV: normalize using common peaks"

manorm.R