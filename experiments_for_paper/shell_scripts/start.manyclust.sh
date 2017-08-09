#! /bin/bash

betavals=(1 2 3 4)
propvals=(3)
nvals=(800)
pvals=(4000)
epsvals=(0.1 0.25)
cvals=(3)
extravals=(2)
exp=(1)

reps=1000

for ((i1=0; i1<${#betavals[@]} ;i1++))
do
for ((i2=0; i2<${#propvals[@]} ;i2++))
do
for ((i3=0; i3<${#nvals[@]} ;i3++))
do
for ((i4=0; i4<${#pvals[@]} ;i4++))
do
for ((i5=0; i5<${#epsvals[@]} ;i5++))
do
for ((i6=0; i6<${#cvals[@]} ;i6++))
do
for ((i7=0; i7<${#extravals[@]} ;i7++))
do
    beta=${betavals[$i1]}
    prop=${propvals[$i2]}
    n=${nvals[$i3]}
    p=${pvals[$i4]}
    eps=${epsvals[$i5]}
    c=${cvals[$i6]}
    extra=${extravals[$i7]}

    fnm="logging/progress-$beta-$prop-$n-$p-$eps-$c-$extra-$reps-$exp.out"
    echo $fnm

    nohup nice R CMD BATCH --no-save --no-restore "--args $beta $prop $n $p $eps $c $extra $reps $exp" main.R $fnm &
done
done
done
done
done
done
done
