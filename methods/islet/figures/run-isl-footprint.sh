cat $0

exe=../../slmm/slmmir

ctr=0
function run {
    ctr=$(expr $ctr + 1)
    cmd="OMP_NUM_THREADS=48 KMP_AFFINITY=balanced $exe -method pcsl $ics -ic gaussianhills -we 0 -rit -dmc eh -T 12 -d2c -ode $ode -ne $ne -np $np -nsteps $nstep -prefine $prefine -mono $cdrglb -lim $cdrlcl -timeint $timeint -footprint"
    grepcmd='grep "^C \|^L \|^M \|footprint>"'
    echo "cmd> $ctr $cmd"
    eval "$cmd | $grepcmd"
}

cdrlcl=caas
ncycle=1
ode=nondivergent
for ne in 30; do
    for nstepfac in 1 5; do
        nstep=$(expr $ne \* 6)
        nstep=$(expr $nstep \* $nstepfac)
        for ode in nondivergent; do
            timeint=exact
            prefine=0
            cdrglb=caas
            np=4
            run
            timeint=interp
            prefine=5
            cdrglb=caas-node
            for np in 6 8 12; do
                run
            done
        done
    done
done
