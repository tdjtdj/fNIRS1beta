SRC = main.cpp mcmc.cpp cholesky.cpp randgen.cpp mybspline.cpp hrf.cpp cpoly.cpp knots.cpp
OBJ = main.o mcmc.o cholesky.o randgen.o mybspline.o hrf.o cpoly.o knots.o
CFLAGS = -O3
LFLAGS = 
LIBS = -lm -lfftw3_threads -lfftw3 -lpthread
gp : ${OBJ} 
	g++ -o foo2 ${OBJ} ${LIBS}

${OBJ} : ${SRC}
	g++ -c ${CFLAGS}  ${SRC}

{OBJ} : cholesky.h LGCP.h randgen.h

clean:
	rm ${OBJ}
	rm *~
