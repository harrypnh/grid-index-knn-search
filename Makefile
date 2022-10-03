all : getResults makeIndex clean

getResults : getResults.o GridIndex.o
	g++ -std=c++17 -o getResults getResults.o GridIndex.o
makeIndex : makeIndex.o GridIndex.o
	g++ -std=c++17 -o makeIndex makeIndex.o GridIndex.o

getResults.o : getResults.cpp
	g++ -std=c++17 -c getResults.cpp
makeIndex.o : makeIndex.cpp
	g++ -std=c++17 -c makeIndex.cpp
GridIndex.o : GridIndex.cpp
	g++ -std=c++17 -c GridIndex.cpp
clean : 
	rm *.o
