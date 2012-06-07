CXX         = g++
CXXLINKER   = g++
CXXFLAGS    = -c -I/opt/mpich/include
LINKERFLAGS = -lmpi
OBJ         = work_withCell.o functions.o main.o

all: kse

kse: $(OBJ)
	mpiCC -o kse $(OBJ) $(LINKERFLAGS)

main.o: main.cpp structures.h work_withCell.h functions.h const.h
	$(CXX) $(CXXFLAGS) $(INCLUDEPATH) main.cpp


run:
	mpirun -np 4 ./kse

work_withCell.o: work_withCell.cpp work_withCell.h
	$(CXX) $(CXXFLAGS) $(INCLUDEPATH) work_withCell.cpp

functions.o: functions.cpp functions.h
	$(CXX) $(CXXFLAGS) $(INCLUDEPATH) functions.cpp
	
clean:
	rm kse $(OBJ)
