#!/bin/bash

ulimit -c unlimited
rm core*
#rm *json *jpg
#rm -rf ABC_unittest_log SFLA_unittest_log
loop=10
for food in 80 #50 #80 100 50
do
   for d in 20 #30 40 50
   do
     let total="$food*$d"
     if [ "$total" -gt "2500" ]; then
        continue 
     fi 

     let lit="$d*5"
     for limit in $lit # 100 200
     do
       for mcn in 2000
       do
	  for data in qws2.txt
	  do 
              suffix=`echo $data | cut -d '.' -f 1`
	      file=$suffix'_'$food'_'$d'_'$limit'_'$mcn'.json'
	      if [ -f $file ]; then
	        continue
	        #rm $file
              fi 
	      touch $file
	      echo '{}' > $file
	      let up="$d-1"
              for alg in 0 1 2 3 4 # 5
              do
                    ./ABC_unittest \
                            -config_path=/mnt/WXRG0274_ssd/yyang4/ficus2/my1/ \
                            -test_data=$data \
                            -FoodNumber=$food \
                            -D=$d \
                            -K=4 \
                            -SD=$food \
                            -limit=$limit \
                            -MCN=$mcn \
                            -algorithm_type=$alg \
                            -lb=0 \
                            -up=$up \
                            -output=$file \
                            -loop=$loop \
                            -step=$up \
                            --alsologtostderr \
#                            -v 4

              done
	      let m="$food/10"
              ./SFLA_unittest \
                      -config_path=/mnt/WXRG0274_ssd/yyang4/ficus2/my1/ \
                      -test_data=qws2.txt \
                      -P=$food \
                      -M=$m \
                      -I=10 \
                      -D=$d \
                      -K=4 \
                      -SD=$food \
                      -miter=$limit \
                      -MCN=$mcn \
                      -step=$up \
                      -algorithm_type=6 \
                      -lb=0 \
                      -up=$up \
                      -output=$file \
                      -loop=$loop \
                      --alsologtostderr \
                     #-v 1
          done
       done
     done
  done
done

