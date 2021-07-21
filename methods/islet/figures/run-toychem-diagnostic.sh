cat $0

exe=../../slmm/slmmir

ctr=0
function run {
    ctr=$(expr $ctr + 1)
    for pg in $(expr $np - 2) $np; do
        cmd="OMP_NUM_THREADS=48 KMP_AFFINITY=balanced $exe -method pcsl $ics -ic gaussianhills -ic toychem1 -ic toychem2 -we 0 -rit -dmc eh -T $(expr 12 \* $ncycle) -d2c -ode $ode -ne $ne -np $np -pg $pg -nsteps $nstep -prefine $prefine -mono $cdrglb -lim $cdrlcl -timeint $timeint"
        grepcmd='grep "^C \|^L \|^M \|toy "'
        echo "cmd> $ctr $cmd"
        eval "$cmd | $grepcmd"
    done
}

cdrlcl=caas
ncycle=10
nstep=576
for ne in 30; do
    for ode in nondivergent; do
        timeint=exact
        prefine=0
        cdrglb=caas
        np=4
        run
        cdrglb=caas-node
        run
        timeint=interp
        prefine=5
        for np in 6 8 9 12; do
            run
        done
    done
done
