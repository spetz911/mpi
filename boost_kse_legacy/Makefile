CFLAGS		= -g
HEADERS		= fvm.h kse.h
OBJECTS		= fvm.o kse.o
TARGETS		= kse
LDFLAGS		= -lboost_mpi -lboost_serialization

all : $(TARGETS)

kse : $(OBJECTS) main.o
	mpic++ -o $@ $^ $(LDFLAGS)
	mkdir -p out/

%.o : %.cc $(HEADERS)
	mpic++ $(CXXFLAGS) $(CFLAGS) $(DEBUG) $(INC) -c $< -o $@

clean :
	rm -rf *.o* *~ $(TARGETS) out/

run:
	mpirun -np 4 ./kse
