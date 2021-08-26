cat $0

exe=../../slmm/slmmir

ctr=0
function run {
    ctr=$(expr $ctr + 1)
    cmd="OMP_NUM_THREADS=48 KMP_AFFINITY=balanced $exe -method pcsl -ic gaussianhills -ic cosinebells -ic correlatedcosinebells -ic slottedcylinders -we 0 -rit -dmc eh -T 12 -d2c -lauritzen -lauritzen-io -ode $ode -ne $ne -np $np -nsteps $nstep -prefine $prefine -mono $cdrglb -lim $cdrlcl -timeint $timeint -o mixing-0/${ode}-$timeint-nsteps${nstep}-prefine${prefine}-$cdrglb-$cdrlcl-ne$ne-np$np"
    grepcmd='grep "^C \|^L \|^M "'
    echo "cmd> $ctr $cmd"
    eval "$cmd | $grepcmd"
}

cdrglbs=(none caas-node);
cdrlcls=(none caas     );
for ne in 20 40; do
    for nstepfac in 1 5; do
        for ode in nondivergent; do
            nstep=$(expr $ne \* 6)
            nstep=$(expr $nstep \* $nstepfac)
            icdr=1
            timeint=exact
            prefine=0
            cdrglb=${cdrglbs[$icdr]}
            cdrlcl=${cdrlcls[$icdr]}
            np=4
            run
            timeint=interp
            prefine=5
            for np in 6 8 9 12; do
                run
            done
            cdrglb=${cdrglb:0:4}
            timeint=exact
            prefine=0
            np=4
            run
        done
    done
done
