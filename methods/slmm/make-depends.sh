for i in *.cpp; do
    g++ -I../../siqk -MM $i
done > make.depends
