for i in *.cpp; do
    g++ -MM $i
done > make.depends
