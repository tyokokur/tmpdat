PROGRAM = tjygo
SRCS    = main.c setup.c FFT_solve.c MDE.c gets.c report.c

OBJDIR = ./bin/obj
SRCDIR = ./bin/src
DEPDIR := ./bin/dep

CC      = gcc
OBJS    = ${foreach src, ${SRCS}, ${addprefix ${OBJDIR}/, ${src:.c=.o}}}
L_FFTW3 = -L/share/apps/software/fftw-3.3.8/lib 
LIBS    = -lm -lfftw3

# LINK 
${PROGRAM}: ${OBJS}
	${CC} ${L_FFTW3} ${OBJS} -o ${PROGRAM} ${LIBS}

# INCLUDE DEPENDENCY INFO FOR EXISTING .o FILES
-include ${DEPDIR}/${OBJS:.o=.d}

# COMPILE AND AUTOGENERATE DEPENDENCIES
${OBJDIR}/%.o: ${SRCDIR}/%.c
	@# COMPILE
	${CC} -c ${L_FFTW3} ${SRCDIR}/$*.c -o ${OBJDIR}/$*.o ${LIBS}

	@# Generate dependency
	${CC} -MM ${L_FFTW3} ${SRCDIR}/$*.c > ${DEPDIR}/$*.d ${LIBS}

	@# Edit $.d to deal with SRCS in another folder than make
	@mv -f ${DEPDIR}/$*.d ${DEPDIR}/$*.d.tmp
	@sed -e 's|.*:|$*.o:|' < ${DEPDIR}/$*.d.tmp > ${DEPDIR}/$*.d

	@# Edit *.d to deal with renaming files (scottmcpeak.com/autodepend/autodepend.html)
	@cp -f ${DEPDIR}/$*.d ${DEPDIR}/$*.d.tmp
	@sed -e 's/.*://' -e 's/\\$$//' < ${DEPDIR}/$*.d.tmp | fmt -1 | \
	  sed -e 's/^ *//' -e 's/$$/:/' >> ${DEPDIR}/$*.d
	@rm -f ${DEPDIR}/$*.d.tmp

# CLEAN
.PHONY = clean clean-comp win
clean: 
	rm -f "("* ph*.dat it*.dat el*.dat printout.dat
win:
	@mv W_a* W_in.dat 2>/dev/null; true
clean-comp:
	rm -f ${PROGRAM} ${OBJDIR}/*.o ${DEPDIR}/*.d

