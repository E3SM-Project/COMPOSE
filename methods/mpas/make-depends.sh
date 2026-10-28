for i in *.cpp; do
    g++ -MM $i -DMPAS_HAVE_NETCDF
done > make.depends
