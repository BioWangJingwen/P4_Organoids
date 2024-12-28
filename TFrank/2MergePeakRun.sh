nohup_number=0

dir="/chenfeilab/Avocado/P4_Organoids/P4_8TFNet/P4_8_1GetNet/P4_8_1_2_MergePeakRun/"
mkdir $dir
cd $dir

Refdir="/chenfeilab/Avocado/P4_Organoids/P4_8TFNet"

# index=`seq 1 10 131` 

while read Sam
do
  echo $Sam

  sed "s/RefSam/$Sam/g" $Refdir/P4_8_1_2MergePeakRef.R > $dir/P4_8_1_2MergePeak$Sam.R

  R CMD BATCH P4_8_1_2MergePeak$Sam.R &
  sleep 12
   
  nohup_number=`echo $nohup_number+1 | bc`
  if [[ $nohup_number -gt 30 ]]
  then
    echo "waiting..."
    wait
    nohup_number=0
  fi
  
done < /chenfeilab/Avocado/P4_Organoids/P4_8TFNet/P4_8_1GetNet/P4_8_1_2_MergePeakEach/P4_8_1_2_samples.txt









