cat $0

exe=../../slmm/slmmir

ctr=0
function run {
    ctr=$(expr $ctr + 1)
    cmd="OMP_NUM_THREADS=48 KMP_AFFINITY=balanced $exe -method $method -ic gaussianhills -ic cosinebells -ic correlatedcosinebells -ic slottedcylinders -we 0 -rit -dmc eh -T $(expr 12 \* $ncycle) -d2c -ode $ode -ne $ne -np $np -nsteps $nstep -prefine $prefine -mono $cdrglb -lim $cdrlcl -timeint $timeint"
    grepcmd='grep "^C \|^L \|^M "'
    echo "cmd> $ctr $cmd"
    eval "$cmd | $grepcmd"
}

cdrglb=caas-node
cdrlcl=caas
for ne in 5 10 20 40 80; do
    for nstepfac in 1; do
        for ode in divergent; do
            nstep=$(expr $ne \* 6)
            nstep=$(expr $nstep \* $nstepfac)
            timeint=interp
            prefine=5
            ncycle=10
            for method in pcslu pcsl; do
                for np in 4 6 8 9 12; do
                    run
                done
                ncycle=100
            done
        done
    done
done
