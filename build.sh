cd src
rm -rf build/
rm *.c *.html *.so
easycython GWSWEX.pyx
mv *.so ../ && clear && echo "Built and Moved .so file" || echo "Build failed!"