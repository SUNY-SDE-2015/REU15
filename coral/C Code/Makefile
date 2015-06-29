
CFLAGS =  -O3 
#CFLAGS =  -g 
CC = g++
cc = gcc
AR = ar
ARFLAGS = rv

MYLIB = 
LOCALLIB = 
LOCALINCLUDE = $(LOCALLIB:.cpp=.hpp)
LOCALOBJECTS = $(LOCALLIB:.cpp=.o)
LINK = -lm 
UTILITY = 
FILES =  CoralSimulation
.SUFFIXES: .c .cpp




.cpp.o:
	@echo 'Compiling $<'
	$(CC) $(CFLAGS) -c $<

.c.o:
	$(cc) $(CFLAGS) -c $<

all:	$(FILES)

CoralSimulation: CoralSimulation.o $(UTILITY)
	echo $@
	$(CC) $(CFLAGS) -o $@ $@.o $(LINK) $(UTILITY) 
	#$(CC) $(CFLAGS) -o $@ $@.o $(LINK) $(UTILITY) $(MYLIB).a 


#$(MYLIB).a:	$(LOCALLIB) $(LOCALOBJECTS) $(LOCALINCLUDE) 
#		rm -f $(MYLIB).a
#		$(AR) $(ARFLAGS) $(MYLIB).a $(LOCALOBJECTS)
#		@echo  Done


clean:
	rm -fr *.o *~ $(FILES) $(MYLIB).a 


