cat $0

exe=slmmir

ctr=0
function run {
    ctr=$(expr $ctr + 1)
    cmd="OMP_NUM_THREADS=64 KMP_AFFINITY=balanced ../$exe -method pcsl -ic vortextracer -we 0 -rit -T 24 -d2c -ode movingvortices -ne $ne -np $np -nsteps $nstep -prefine $prefine -mono $cdrglb -lim $cdrlcl -timeint $timeint -rotate-grid"
    grepcmd='grep "^C \|^L \|^M "'
    echo "cmd> $ctr $cmd"
    eval "$cmd | $grepcmd"
}

cdrglb=none
cdrlcl=none
timeint=exact
prefine=0
for nstepfac in 1 5; do
    for tne in 5 10 20 40 80 160; do
        tgtres=$(($tne * 3))
        for np in $(seq 4 13) 16; do
            ne=$(($tgtres / ($np - 1)))
            nstep=$(($nstepfac * 2 * ($np - 1) * $ne))
            if [[ $ne == 1 ]]; then
                continue
            fi
            echo ">>> $ne $np $nstep"
            run
        done
    done
done
