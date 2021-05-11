nthr=1
for nremap in 1; do
    for limiter in none caas; do
        for np in 4 5 8; do
            for nphys in 2 3 4 6; do
                for fv2gll in idem elrecon; do
                    if [[ $np -lt $nphys ]]; then
                        continue
                    fi
                    if [[ $nphys == 6 && $np -lt 8 ]]; then
                        continue;
                    fi
                    if [[ $fv2gll == elrecon && $nphys -lt 3 ]]; then
                        continue;
                    fi
                    for ne in 5 10 20 40 80 160 320; do
                        for ic in gaussianhills cosinebells slottedcylinders; do
                            OMP_NUM_THREADS=$nthr OMP_PROC_BIND=false OMP_PLACES=cores ../physgrid -lim $limiter -ne $ne -np $np -nphys $nphys -nremap $nremap -rho constant -q $ic -fv2gll $fv2gll -omitrho
                        done
                    done
                done
            done
        done
    done
done
