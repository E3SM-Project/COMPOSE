exe=../../slmm/slmmir

ne=20
np=6
nstep=$(($ne * 6))
we=110

glbcdrs=(none caas-node none caas-node)
bases=(Gll Gll GllNodal GllNodal)
ncycles=(1 1 1 1)

for i in $(seq 0 3); do
    ncycle=${ncycles[$i]}
    glbcdr=${glbcdrs[$i]}
    basis=${bases[$i]}
    echo $i $ncycle $glbcdr $basis
    name=instab-nondiv-ne${ne}np${np}-${glbcdr}-cyc${ncycle}-$basis
    time=$(($ncycle * 12))
    cmd="OMP_PROC_BIND=false OMP_NUM_THREADS=8 $exe -method pcsl -ode nondivergent -ic gaussianhills -ic cosinebells -ic slottedcylinders -T $time -nsteps $nstep -timeint exact -ne $ne -np $np -dmc eh -d2c -mono $glbcdr -lim caas -we $we -io-type internal -io-start-cycle $ncycle -res 256 -o $name -rit -prefine 0 -basis $basis -rhot0"
    echo "cmd> $cmd"
    eval "$cmd"
done
