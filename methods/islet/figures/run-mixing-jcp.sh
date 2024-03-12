cat $0

exe=slmmir

ctr=0
function run {
    ctr=$(expr $ctr + 1)
    cmd="OMP_NUM_THREADS=1 KMP_AFFINITY=balanced ../$exe -method pcsl -ic gaussianhills -ic cosinebells -ic correlatedcosinebells -ic slottedcylinders -we 0 -rit -T 12 -d2c -lauritzen -lauritzen-io -ode $ode -ne $ne -np $np -nsteps $nstep -prefine $prefine -mono $cdrglb -lim $cdrlcl -timeint $timeint -rotate-grid -o mixing-jcp/${ode}-$timeint-nsteps${nstep}-prefine${prefine}-$cdrglb-$cdrlcl-ne$ne-np$np"
    grepcmd='grep "^C \|^L \|^M "'
    echo "cmd> $ctr $cmd"
    eval "$cmd | $grepcmd"
}

cdrglbs=(none caas-node);
cdrlcls=(none caas);
for tne in 20 40; do
    tgtres=$(($tne * 3))
    for nstepfac in 1 5; do
        for ode in nondivergent; do
            icdr=0
            cdrglb=${cdrglbs[$icdr]}
            cdrlcl=${cdrlcls[$icdr]}
            timeint=exact
            prefine=0
            for np in 4 6 8 11 13; do
                ne=$(($tgtres / ($np - 1)))
                nstep=$(($nstepfac * 2 * ($np - 1) * $ne))
                echo ">>> $ne $np $nstep"
                run
            done
        done
    done
done
