CC = gcc
CXX = g++
CFLAGS = -fopenmp -O3 -DNDEBUG -static-libgcc -static-libstdc++ -static -lpthread -lm -lgomp -finline-functions -fomit-frame-pointer -fno-operator-names -march=x86-64 -msse4.2 -mfpmath=sse -ffast-math -fexcess-precision=fast -ffinite-math-only #-DPARALLEL#-DFLACH #-DCOMPARE_TIMING_ENGINES #-I. #-fopenmp 
CXXFLAGS = -fopenmp -O3 -DNDEBUG -static-libgcc -static-libstdc++ -static -lpthread -lm -lgomp -finline-functions -fomit-frame-pointer -fno-operator-names -march=x86-64 -msse4.2 -mfpmath=sse -ffast-math -fexcess-precision=fast -ffinite-math-only # -DPARALLEL #-DFLACH #-DCOMPARE_TIMING_ENGINES #-I. #-fopenmp 

SRC     = parser_helper.cpp timer_interface.cpp Circuit.cpp global.cpp Vcell.cpp
OBJ     = $(SRC:.cpp=.o) $(SRCflute:.c=.o)

all: sizerbin

sizerbin: main.cpp ${OBJ}  
	$(CXX) $(CXXFLAGS) -o sizerbin main.cpp ${OBJ} -lm

parser_helper.o: parser_helper.h parser_helper.cpp
	$(CC) $(CFLAGS) -c -o parser_helper.o parser_helper.cpp

timer_interface.o: timer_interface.h timer_interface.cpp
	$(CC) $(CFLAGS) -c -o timer_interface.o timer_interface.cpp

Circuit.o: Circuit.h Circuit.cpp
	$(CC) $(CFLAGS) -c -o Circuit.o Circuit.cpp

Vcell.o: Vcell.h Vcell.cpp
	$(CC) $(CFLAGS) -c -o Vcell.o Vcell.cpp

global.o: global.h global.cpp
	$(CC) $(CFLAGS) -c -o global.o global.cpp

clean:
	rm -f *.o core* sizerbin
