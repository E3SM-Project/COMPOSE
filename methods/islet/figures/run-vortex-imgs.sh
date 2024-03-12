cat $0

exe=slmmir

ctr=0
function run {
    ctr=$(expr $ctr + 1)
    cmd="OMP_NUM_THREADS=1 KMP_AFFINITY=balanced ../$exe -method pcsl -ic vortextracer -rit -T 24 -d2c -ode movingvortices -ne $ne -np $np -nsteps $nstep -prefine $prefine -mono $cdrglb -lim $cdrlcl -timeint $timeint -rotate-grid -we $nstep -io-type internal -io-recon $iorecon -res 1024 -o vortex-imgs/ne${ne}np${np}nstep${nstep}${iorecon}"
    grepcmd='grep "^C \|^L \|^M "'
    echo "cmd> $ctr $cmd"
    eval "$cmd | $grepcmd"
}

cdrglb=none
cdrlcl=none
timeint=exact
prefine=0
iorecon=bilin

ne=3
np=13
nstep=$((12 * 6))
run
iorecon=constant
run
iorecon=bilin

ne=6
np=13
nstep=$((24 * 6))
run

ne=12
np=13
nstep=$((48 * 6))
run
