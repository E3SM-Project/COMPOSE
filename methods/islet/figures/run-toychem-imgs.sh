exe=../../slmm/slmmir
datadir=toychem-imgs

ne=30
nstep=$(expr 48 \* 12)
dt=dt30min
we=$(expr 48 \* 6)

np=8
glbcdr=caas-node
pg=$np
name=toychem-nondiv-ne${ne}pr5np${np}pg${pg}-${glbcdr}-caas-$dt
cmd="$exe -method pcsl -ode nondivergent -ic gaussianhills -ic toychem1 -ic toychem2 -T 12 -nsteps $nstep -timeint interp -ne $ne -np ${np} -dmc eh -d2c -mono $glbcdr -lim caas -we $we -io internal -res 256 -o $datadir/$name -rit -prefine 5 -pg $pg"
echo "cmd> $cmd"
eval "$cmd"

glbcdr=caas
np=4
pg=0
name=toychem-nondiv-ne${ne}pr0np${np}pg${pg}-${glbcdr}-caas-$dt
cmd="$exe -method pcsl -ode nondivergent -ic gaussianhills -ic toychem1 -ic toychem2 -T 12 -nsteps $nstep -timeint exact -ne $ne -np ${np} -dmc eh -d2c -mono $glbcdr -lim caas -we $we -io internal -res 256 -o $datadir/$name -rit -prefine 0 -pg $pg"
echo "cmd> $cmd"
eval "$cmd"
