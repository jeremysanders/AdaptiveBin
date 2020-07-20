# objects
CXX=g++
LD=g++
objFITS = FITSmm/FITSmm.a
objParammm = parammm/libparammm.a

programs = AdaptiveBin ABPixelCopy MakeMask AnnuliMap AdaptiveAnnuli AdaptiveContour AdaptiveBinT RayMap

all:	$(programs)

objMergeBinMap = MergeBinMap.o $(objFITS) $(objParammm)
objAdaptiveContour = AdaptiveContour.o $(objFITS) $(objParammm)
objAdaptiveBin = AdaptiveBin.o binmodule.o $(objFITS) $(objParammm)
objAdaptiveBlock = AdaptiveBlock.o SigCalc.o $(objFITS) $(objParammm)
objABPostSmooth = ABPostSmooth.o $(objFITS) $(objParammm)
objABPixelCopy = ABPixelCopy.o $(objFITS) $(objParammm)
objBinOnGrid = BinOnGrid.o $(objFITS) $(objParammm)
objRayMap = RayMap.o $(objFITS) $(objParammm)
objAnnuliMap = AnnuliMap.o $(objFITS) $(objParammm)
objMakeMask = MakeMask.o $(objFITS) $(objParammm)
objAdaptiveAnnuli = AdaptiveAnnuli.o $(objFITS) $(objParammm)
objAdaptiveBinT = AdaptiveBinT.o binmodule.o $(objFITS) $(objParammm)

# header files
headAdaptiveBin = Coord.hh
headAdaptiveBlock = Coord.hh SigCalc.hh

# c++ options
CXXFLAGS = -Wall -g -O2 -IFITSmm -I.


# object files
AdaptiveContour.o : version.hh
AdaptiveBin.o : binmodule.hh version.hh
SigCalc.o : $(headAdaptiveBlock)
AdaptiveBlock.o : $(headAdaptiveBlock)
ABPostSmooth.o :
ABPixelCopy.o : version.hh
binmodule.o: binmodule.hh
BinOnGrid.o :
MergeBinMap.o :
RayMap.o : version.hh
AnnuliMap.o : version.hh
MakeMask.o :
AdaptiveAnnuli:
AdaptiveBinT.o : binmodule.hh version.hh

# programs
AdaptiveAnnuli: $(objAdaptiveAnnuli) $(objFITS)
	g++ -o AdaptiveAnnuli $(objAdaptiveAnnuli) $(objFITS) -lm \
	-lcfitsio $(objParammm)
AnnuliMap : $(objAnnuliMap) $(objFITS)
	g++ -o AnnuliMap $(objAnnuliMap) $(objFITS) -lm -lcfitsio \
	$(objParammm)
RayMap : $(objRayMap) $(objFITS)
	g++ -o RayMap $(objRayMap) $(objFITS) -lm -lcfitsio \
	$(objParammm)
MergeBinMap : $(objMergeBinMap) $(objFITS)
	g++ -o MergeBinMap $(objMergeBinMap) $(objFITS) -lm \
	-lcfitsio
AdaptiveContour: $(objAdaptiveContour) $(objFITS)
	g++ -o AdaptiveContour $(objAdaptiveContour) $(objFITS) -lm \
	-lcfitsio
AdaptiveBin : $(objAdaptiveBin) $(objFITS)
	g++ -o AdaptiveBin $(objAdaptiveBin) $(objFITS) -lm -lcfitsio \
	$(objParammm)
AdaptiveBlock : $(objAdaptiveBlock) $(objFITS)
	g++ -o AdaptiveBlock $(objAdaptiveBlock) $(objFITS) -lm -lcfitsio
ABPostSmooth : $(objABPostSmooth) $(objFITS)
	g++ -o ABPostSmooth $(objABPostSmooth) $(objFITS) -lm -lcfitsio \
	-lngmath $(objParammm)
ABPixelCopy : $(objABPixelCopy) $(objFITS)
	g++ -o ABPixelCopy $(objABPixelCopy) $(objFITS) -lm -lcfitsio \
	$(objParammm)
BinOnGrid : $(objBinOnGrid) $(objFITS)
	g++ -o BinOnGrid $(objBinOnGrid) $(objFITS) -lm -lcfitsio \
	$(objParammm)
MakeMask : $(objMakeMask) $(objFITS)
	g++ -o MakeMask $(objMakeMask) $(objFITS) -lm -lcfitsio \
	$(objParammm)
AdaptiveBinT : $(objAdaptiveBinT) $(objFITS)
	g++ -o AdaptiveBinT $(objAdaptiveBinT) $(objFITS) -lm -lcfitsio \
	$(objParammm)

FITSmm/FITSmm.a:
	$(MAKE) -C FITSmm FITSmm.a

parammm/libparammm.a:
	$(MAKE) -C parammm libparammm.a

clean:
	-rm -f *.o parammm/*.o parammm/*.a FITSmm/*.o FITSmm/*.a $(programs)
