cat $0
exe=../../slmm/slmmir

mkdir -p filament-imgs

function run {
    cmd="$exe -method pcsl -ode nondivergent -ic slottedcylinders -T 12 -ne $ne -nsteps $nsteps -timeint $timeint -nonunimesh 0 -np $np -dmc eh -mono $cdrglb -lim caas -lauritzen -we $we -io internal -o filament-imgs/ne$ne-np$np-nstep$nsteps-$timeint-$cdrglb-pr$prefine -res 256 -rit -prefine $prefine $d2c"
    echo "cmd> $cmd"
    eval "OMP_NUM_THREADS=48 KMP_AFFINITY=balanced $cmd"
}

d2c="-d2c -io-nodss"
for ne in 20 40; do
    for nstepfac in 1 5; do
        nsteps=$(expr $ne \* 6)
        nsteps=$(expr $nsteps \* $nstepfac)
        we=$(expr $nsteps / 2)
        timeint=exact
        prefine=0
        cdrglb=caas
        np=4
        run
        cdrglb=caas-node
        run
        timeint=interp
        prefine=5
        cdrglb=caas-node
        for np in 6 8; do
            run
        done
        if [[ $ne == 20 ]]; then
            np=12
            run
        fi
    done
done
