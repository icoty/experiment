#!/bin/bash

ulimit -c unlimited
#rm core* *json f*
#rm -rf ABC_func_unittest_log SFLA_func_unittest_log

declare -A func_map=()
func_map["1"]="-100_100"
func_map["2"]="-10_10"
func_map["3"]="-100_100"
func_map["4"]="-10_10"
func_map["5"]="-5.12_5.12"
func_map["6"]="-600_600"
func_map["7"]="-32_32"
func_map["8"]="-10_10"

loop=10

for food in 100
do
   for sd in 200
   do
     for limit in 200
     do
       for mcn in 2500
       do
         for f in 8 #5 7 # 1 2 4 6 # 5 7 8 #3
	 do 
           str=${func_map[$f]}
	   lb=`echo $str | cut -d '_' -f 1`
	   up=`echo $str | cut -d '_' -f 2`
	   #let step="$up/1"
	   step=$up
	   
	   file='f'$f'_'$food'_'$sd'_'$limit'_'$mcn'.json'
	   if [ -f $file ]; then
	     continue
	     #rm $file
           fi
	   touch $file
	   echo '{}' > $file
	   echo "AAA $f $lb $up $file"
           for alg in 0 1 2 3 4 #5
           do
                ./ABC_func_unittest \
                        -config_path=/mnt/WXRG0274_ssd/yyang4/ficus2/func_my1/ \
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
#                        -v 4

           done

	   let m="$food/10"
	   echo "BBB $m"
           ./SFLA_func_unittest \
                   -config_path=/mnt/WXRG0274_ssd/yyang4/ficus2/func_my1/ \
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

