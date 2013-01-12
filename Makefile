CC=			gcc
CXX=		g++
CFLAGS=		-g  -O2 
CXXFLAGS=	$(CFLAGS)
DFLAGS=		-DMAIN_INDEX#-DHAVE_PTHREAD #-DDEBUG#-D_FILE_OFFSET_BITS=64 
OBJS=		utils.o bwt.o bwtio.o  is.o \
			bntseq.o bwtmisc.o\
		bwt_gen.o  QSufSort.o\
		4bit_bntseq.o 4bit_bwt_gen.o hapmap.o   rbwt.o LookUpTable.o\
	        localPattern.o mixRef.o index1.o	
PROG=		index
INCLUDES=	
LIBS=		-lm -lz -lpthread 
SUBDIRS=	. 

.SUFFIXES:.c .o .cc

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@
.cc.o:
		$(CXX) -c $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)




index:$(OBJS) 
		$(CC) $(CFLAGS) $(DFLAGS) $(OBJS) -o $@ $(LIBS)


clean:
		rm -f gmon.out *.o a.out $(PROG) *~ *.a


