cat $0

exe=../../slmm/slmmir

ctr=0
function run {
    ctr=$(expr $ctr + 1)
    ics1="-ic gaussianhills -ic cosinebells -ic slottedcylinders"
    ics2="-ic zero -ic zero -ic zero"
    for pg in 2 $(expr $np - 2) $np; do
        cmd="OMP_NUM_THREADS=48 KMP_AFFINITY=balanced $exe -method pcsl $ics1 -ic toychem1 -ic toychem2 $ics2 -we 0 -rit -dmc eh -T 12 -d2c -ode $ode -ne $ne -np $np -pg $pg -nsteps $nstep -prefine $prefine -mono $cdrglb -lim $cdrlcl -timeint $timeint -midpoint-check"
        grepcmd='grep "^C \|^L \|^M "'
        echo "cmd> $ctr $cmd"
        eval "$cmd | $grepcmd"
    done
}

cdrglb=caas-node
cdrlcl=caas
for ne in 5 10 20 40 80; do
    for nstepfac in 1 5; do
        for ode in nondivergent; do #divergent; do
            nstep=$(expr $ne \* 6)
            nstep=$(expr $nstep \* $nstepfac)
            timeint=exact
            prefine=0
            np=4
            run
            timeint=interp
            prefine=5
            for np in 6 8 9 12; do
                run
            done
        done
    done
done
