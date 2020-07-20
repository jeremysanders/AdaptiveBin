# objects to be made
objf = FITSGeneral.o FITSImage.o FITSFile.o FITSPosn.o
libf = FITSmm.a

# dependancy headers
fitsh   = FITSGeneral.h FITSImage.h FITSFile.h FITSPosn.h

CXX=g++

# c options to pass
CXXFLAGS = -Wall -g -O2 -I/home/jss/code/libs -I/home/jss/code/classes -I/usr/include/cfitsio -I.

FITSmm.a : $(objf)
	ar ruv FITSmm.a $?
	ranlib FITSmm.a

clean :
	rm -f $(objf) $libf
