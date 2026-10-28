Build:
  Modify make.inc
  make

Run:
  mpas_standalone is the main driver.

To remake dependencies,
  bash make-depends.sh

Example run script:

$ cat runseq.sh
prob=...  # Fill this in
nstep=20
chalo=15
np=64
for qu in 480 240 120 60 30 15; do
    f=/path/to/mesh.QU.${qu}km.151026
    time srun -n $np ./mpas_standalone -m $f.nc -g $f.graph.info -p $f.graph.info.part.$np -nstep $nstep -chalo $chalo -o $prob
    nstep=$(( $nstep * 2 ))
done
$ bash runseq.sh

The $prob variable is a string that encodes the problem and algorithms to run.
See optstr2options in mpas_advect_ir_test.cpp.

mpas_standalone requires a graph file and a parmetis partition. mpas_python.cpp
some utilities for interacting with this code in Python. For graph files in
particular, call
    set_mesh
and then
    make_graph_file.
Then make a partition file:
    gpmetis graph-info-file-name num-procs
mpas.hy demonstrates how to use mpas_python.cpp, but you will likely want to
translate mpas.hy to python.
