#!/bin/bash

ulimit -c unlimited
#rm core* *json f*
#rm -rf ABC_case_unittest_log

declare -A func_map=()
func_map["0"]="1_6"

loop=5

for food in 50
do
   for sd in 30 # 50 #200
   do
     let lit="$sd*5"
     for limit in $sd
     do
       let iter="$food*$sd"
       for mcn in 100 #2000 #$iter 20000 200
       do
         for f in 0
	     do 
           str=${func_map[$f]}
    	   lb=`echo $str | cut -d '_' -f 1`
    	   up=`echo $str | cut -d '_' -f 2`
    	   #let step="$up/1"
    	   step=$up
    	   
    	   file='test_case_'$f'_'$food'_'$sd'_'$limit'_'$mcn'.json'
    	   if [ -f $file ]; then
    	     continue
    	     #rm $file
           fi
    	   touch $file
    	   echo '{}' > $file
    	   echo "AAA $f $lb $up $file"
           for alg in 3 4 5 6 7
           do
                ./ABC_case_unittest \
                        -config_path=/mnt/WXRG0274_ssd/yyang4/ficus2/experient/test_case/ \
                        -FoodNumber=$food \
                        -SD=$sd \
                        -limit=$limit \
                        -MCN=$mcn \
                        -algorithm_type=$alg \
                        -lb=$lb \
                        -up=$up \
                        -func=$f \
                        -output=$file \
                        -loop=$loop \
                        -step=$step \
                        --alsologtostderr \
                        -test_case_time=case_time.txt \
#                       -v 4
           done
           continue
	       let m="$food/10"
	       echo "BBB $m"
           ./SFLA_func_unittest \
                   -config_path=/mnt/WXRG0274_ssd/yyang4/ficus2/experient/func_my2/ \
                   -P=$food \
                   -SD=$sd \
                   -miter=10 \
                   -MCN=$mcn \
                   -algorithm_type=6 \
                   -lb=$lb \
                   -up=$up \
                   -func=$f \
                   -output=$file \
                   -loop=$loop \
                   -step=$step \
                   -M=$m \
                   -I=10 \
                   --alsologtostderr \
#                   -v 1
         done 
       done
     done
  done
done

