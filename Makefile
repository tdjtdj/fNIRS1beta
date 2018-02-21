SRC = main.cpp mcmc.cpp cholesky.cpp randgen.cpp mybspline.cpp hrf.cpp knots.cpp statistics.cpp config_info.cpp dlm.cpp
OBJ = main.o mcmc.o cholesky.o randgen.o mybspline.o hrf.o knots.o statistics.o config_info.o dlm.o
CFLAGS = -O3
LFLAGS = 
LIBS = -lm -lfftw3
gp : ${OBJ} 
	g++ -o foo ${OBJ} ${LIBS}

${OBJ} : ${SRC}
	g++ -c ${CFLAGS}  ${SRC}

{OBJ} : cholesky.h LGCP.h randgen.h

clean:
	rm ${OBJ}
	rm *~
