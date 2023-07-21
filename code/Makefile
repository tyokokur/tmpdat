PROGRAM = tjygo
FILES.c = main.c setup.c report.c FFT_solve.c MDE.c
FILES.o = ${FILES.c:.c=.o}
UTILS.h = common.h 
L_FFTW3 = -L/share/apps/software/fftw-3.3.8/lib 

CC = gcc
WFLAGS = -Wall

all: ${PROGRAM}

${PROGRAM}: ${FILES.o}
	${CC} ${L_FFTW3} -o $@ ${WFLAGS} ${FILES.o} -lfftw3 -lm
main.o: ${UTILS.h} consts.h input.h setup.h MDE.h
setup.o:  ${UTILS.h}
report.o: ${UTILS.h}
FFT_solve.o: ${UTILS.h}
MDE.o: ${UTILS.h} FFT_solve.h
