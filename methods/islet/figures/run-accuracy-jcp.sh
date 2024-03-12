cat $0

exe=slmmir

ctr=0
function run {
    ctr=$(expr $ctr + 1)
    cmd="OMP_NUM_THREADS=64 KMP_AFFINITY=balanced ../$exe -method pcsl -ic gaussianhills -ic cosinebells -ic correlatedcosinebells -ic slottedcylinders -we 0 -rit -T 12 -d2c -ode $ode -ne $ne -np $np -nsteps $nstep -prefine 0 -mono none -lim none -timeint exact -rotate-grid"
    grepcmd='grep "^C \|^L \|^M "'
    echo "cmd> $ctr $cmd"
    eval "$cmd | $grepcmd"
}

for tne in 5 10 20 40 80 160; do
    for nstepfac in 1 5; do
        for ode in nondivergent; do
            tgtres=$(($tne * 3))
            for np in $(seq 4 13) 16; do
                ne=$(($tgtres / ($np - 1)))
                nstep=$(($nstepfac * 2 * ($np - 1) * $ne))
                if [[ $ne == 1 ]]; then
                    continue
                fi
                run
            done
        done
    done
done
