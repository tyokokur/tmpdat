PROGRAM = tjygo
SRCS    = main.c setup.c FFT_solve.c DGC.c PB.c gets.c report.c

OBJDIR = ./bin/obj
SRCDIR = ./bin/src

CC      = gcc -g
OBJS    = ${foreach src, ${SRCS}, ${addprefix ${OBJDIR}/, ${src:.c=.o}}}
LIBS = -fopenmp -lm -lfftw3 -lfftw3_omp

# LINK 
${PROGRAM}: ${OBJS}
	${CC} -Wall ${OBJS} -o ${PROGRAM} ${LIBS}

# INCLUDE DEPENDENCY INFO FOR EXISTING .o FILES
-include ${foreach src, ${SRCS}, ${addprefix ${OBJDIR}/, ${src:.c=.d}}}

# COMPILE AND AUTOGENERATE DEPENDENCIES
${OBJDIR}/%.o: ${SRCDIR}/%.c Makefile
	@# COMPILE
	${CC} -c ${SRCDIR}/$*.c -o ${OBJDIR}/$*.o ${LIBS}

	@# Generate dependency
	@${CC} -MMD -MP -c $< -o $@ ${LIBS}

	@# Edit $.d to deal with SRCS in another folder than make
	@#mv -f ${DEPDIR}/$*.d ${DEPDIR}/$*.d.tmp
	@#sed -e 's|.*:|$*.o:|' < ${DEPDIR}/$*.d.tmp > ${DEPDIR}/$*.d

	@# Edit *.d to deal with renaming files (scottmcpeak.com/autodepend/autodepend.html)
	@#cp -f ${DEPDIR}/$*.d ${DEPDIR}/$*.d.tmp
	@#sed -e 's/.*://' -e 's/\\$$//' < ${DEPDIR}/$*.d.tmp | fmt -1 | \
	 #sed -e 's/^ *//' -e 's/$$/:/' >> ${DEPDIR}/$*.d
	@#rm -f ${DEPDIR}/$*.d.tmp


# CLEAN
.PHONY = clean clean-comp win
clean: 
	rm -f "("* ph*.dat it*.dat el*.dat printout.dat 
	rm -f ./logs/*
win:
	@mv W_a* W_in.dat 2>/dev/null; true
clean-comp:
	rm -f ${PROGRAM} ${OBJDIR}/*.o ${DEPDIR}/*.d ${DEPDIR}/*.d.tmp

