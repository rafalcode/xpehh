
nooptf = -g -Wall -fopenmp
noopt = -g -Wall
optimization = -O3 
objects = ehhpop.o CmdLine.o 
BINARIES = *.o xpehh_omp xpehh hzy test
libs = -L/opt/local/lib -I/usr/include -lgsl -lgslcblas -lm

xpehh_omp: xpehh_omp.o $(objects)
	g++ $(libs) -o $@ $^ $(nooptf)

xpehh: xpehh.o $(objects)
	g++ $(libs) -o $@ $^ $(optimization)

hzy: hzy.o $(objects)
	g++ $(libs) -o $@ $^ $(optimization)
 
test: test.o $(objects)
	g++ $(libs) -o $@ $^ $(optimization)
 
# %.o: %.cpp
# 	g++ -c $< -o $@ $(noopt)
#
ehhpop.o: ehhpop.cpp
	g++ -c $< -o $@ $(noopt)
Cmdline.o: Cmdline.cpp
	g++ -c $< -o $@ $(noopt)
hzy.o: hzy.cpp
	g++ -c $< -o $@ $(noopt)
test.o: test.cpp
	g++ -c $< -o $@ $(noopt)
xpehh.o: xpehh.cpp
	g++ -c $< -o $@ $(noopt)
xpehh_omp.o: xpehh_omp.cpp
	g++ -c $< -o $@ $(nooptf)

clean:
	rm ${BINARIES}
	# rm *.o mpop # wtf is mpop?
