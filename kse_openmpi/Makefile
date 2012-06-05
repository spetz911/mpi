CFLAGS		= -g
HEADERS		= fvm.h kse.h
OBJECTS		= fvm.o kse.o
TARGETS		= kse

all : $(TARGETS)

kse : $(OBJECTS) main.o
	mpicc -o $@ $^ $(LDFLAGS)
	mkdir -p out/

#test : $(OBJECTS) test.o
#	mpicc -o $@ $^ $(LDFLAGS)


%.o : %.c $(HEADERS)
	gcc $(CXXFLAGS) $(CFLAGS) $(DEBUG) $(INC) -c $< -o $@

clean :
	rm -rf *.o* *~ $(TARGETS)


