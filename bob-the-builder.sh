g++ -std=c++17 -O2 -c main.cpp      
g++ -std=c++17 -O2 -c multigrid.cpp 
g++ -std=c++17 -O2 -c grid.cpp 
g++ -std=c++17 -O2 -c field.cpp 

g++ -std=c++17 -O2 main.o multigrid.o grid.o field.o -o polly


