!# /bin/bash

nvals=(200 300 400 600 800 1000 1500 2000 2500 3000 4000 5000)
nvals=(200 300)
reps=3

for ((ii=0; ii<${#nvals[@]} ;ii++))
do
    n=${nvals[$ii]}
    fnm="logging/progress-$n.out"
    echo $fnm
    nohup nice R CMD BATCH --no-save --no-restore "--args $n $reps " labor_script.R $fnm &
done
