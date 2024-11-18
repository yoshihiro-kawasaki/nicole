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
		  src/lsoda.o \
		  src/linear.o 

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

src/linear.o: src/linear.hpp src/linear.cpp 
	(cd src; $(CPP) $(CFLAGS) $(DEPENDPATH) linear.cpp)

src/lsoda.o: src/lsoda.hpp src/lsoda.cpp 
	(cd src; $(CPP) $(CFLAGS) $(DEPENDPATH) lsoda.cpp)

clean:
	rm run $(objects)