# CPP := g++
# CPPFLAGS := -O3 -std=c++11 -DNDEBUG
# LFLAGS :=

CPP = g++

CFLAGS = -g -c -O3

# LFLAGS = -lgsl -lgslcblas -lm

LFLAGS =

DEPENDPATH = -I.

objects = src/main.o \
		  src/nicole.o \
		  src/dust.o \
		  src/utils.o \
		  src/lsoda_cpp/lsoda.o \
		  src/lsoda_cpp/linear.o 

run: $(objects)
	$(CPP) -o run $(objects) $(LFLAGS)

src/main.o: src/main.cpp
	(cd src; $(CPP) $(CFLAGS) $(DEPENDPATH) main.cpp)
	
src/nicole.o: src/nicole.hpp src/nicole.cpp
	(cd src; $(CPP) $(CFLAGS) $(DEPENDPATH) nicole.cpp)

src/dust.o: src/dust.hpp src/dust.cpp 
	(cd src; $(CPP) $(CFLAGS) $(DEPENDPATH) dust.cpp)

src/utils.o: src/utils.hpp src/utils.cpp 
	(cd src; $(CPP) $(CFLAGS) $(DEPENDPATH) utils.cpp)

src/lsoda_cpp/linear.o: src/lsoda_cpp/linear.hpp src/lsoda_cpp/linear.cpp 
	(cd src/lsoda_cpp; $(CPP) $(CFLAGS) $(DEPENDPATH) linear.cpp)

src/lsoda_cpp/lsoda.o: src/lsoda_cpp/lsoda.hpp src/lsoda_cpp/lsoda.cpp 
	(cd src/lsoda_cpp; $(CPP) $(CFLAGS) $(DEPENDPATH) lsoda.cpp)

clean:
	rm run $(objects)