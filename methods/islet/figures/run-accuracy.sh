cat $0

exe=../../slmm/slmmir

mkdir -p tmp

ctr=0
function run {
    ctr=$(expr $ctr + 1)
    cmd="OMP_NUM_THREADS=48 KMP_AFFINITY=balanced $exe -method pcsl -ic gaussianhills -ic cosinebells -ic correlatedcosinebells -ic slottedcylinders -we 0 -rit -dmc eh -T 12 -d2c -lauritzen -midpoint-check -ode $ode -ne $ne -np $np -nsteps $nstep -prefine $prefine -mono $cdrglb -lim $cdrlcl -timeint $timeint"
    grepcmd='grep "^C \|^L \|^M "'
    echo "cmd> $ctr $cmd"
    eval "$cmd | $grepcmd"
}

cdrglbs=(none caas-node);
cdrlcls=(none caas     );
for ne in 5 10 20 40 80; do
    for nstepfac in 1 5; do
        for ode in rotate nondivergent divergent; do
            nstep=$(expr $ne \* 6)
            nstep=$(expr $nstep \* $nstepfac)
            timeint=exact
            for icdr in 1; do
                cdrglb=${cdrglbs[$icdr]}
                cdrlcl=${cdrlcls[$icdr]}
                prefine=0
                np=4
                run
            done
            for icdr in 0 1; do
                cdrglb=${cdrglbs[$icdr]}
                cdrlcl=${cdrlcls[$icdr]}
                prefine=0
                for np in 4; do
                    cdrglb=${cdrglb:0:4}
                    run
                done
            done
            timeint=interp
            for icdr in 0 1; do
                cdrglb=${cdrglbs[$icdr]}
                cdrlcl=${cdrlcls[$icdr]}
                prefine=5
                for np in $(seq 5 13) 16; do
                    run
                done
            done
            timeint=exact
            for icdr in 0; do
                cdrglb=${cdrglbs[$icdr]}
                cdrlcl=${cdrlcls[$icdr]}
                prefine=0
                for np in $(seq 5 13) 16; do
                    run
                done                
            done
        done
    done
done
