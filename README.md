Instrucción de compilación aplicación:  
g++ tests/EKF_GEOS3.cpp src/*cpp -lm -std=c++23 -o bin/main.exe
cd bin
main.exe
pause  

Instrucción de compilación  para los test unitarios:  
g++ tests/tests.cpp src/*.cpp -lm -std=c++23 -o bin/tests.exe
cd bin
tests.exe
pause
