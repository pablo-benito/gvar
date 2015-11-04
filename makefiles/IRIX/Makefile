ARCH = `uname`
D = ../..

TARGET = gvar

CC = CC -O2 -fullwarn -OPT:Olimit=2026 -D_MSBF_ -D_msbf_
#CC = CC -g -fullwarn  -D_MSBF_ -D_msbf_
# add -D_MODIFY_LINE_NUMBERS_ to process goes10 from before 3 Feb 1998
INCLUDEDIR = /usr/include/CC 
#SYSLIBS = -lC -lc -L/usr/local/lib -lmmalloc -lm 
SYSLIBS = -lC -lc -L/usr/local/lib -lm 
cc = cc -D_MSBF_ -D_msbf_

# for gnu g++ uncomment the next three lines
#CC = gcc -O 
#INCLUDEDIR =  /usr/local/lib/g++-include
#SYSLIBS = -L/usr/local/lib -lstdc++ -lg++ -lm 
#cc = gcc

CCFLAGS = -I. -I$D -I $D/f2c -I $(INCLUDEDIR) -I/usr/include

F2CLIBS = -L $D/f2c/libF77 -lF77$(ARCH) \
	  -L $D/f2c/libI77 -lI77$(ARCH)

IMAGELIBS = \
	-L $D/external_libs/jpeg -ljpeg$(ARCH) \
	-L $D/external_libs/tiff -ltiff$(ARCH) 


LIBS = $(F2CLIBS) $(IMAGELIBS) $(SYSLIBS)	

SRC = \
$D/block0/attitudeAngle.C	\
$D/block0/clampMode.C \
$D/block0/grid.C \
$D/block0/imagerDoc.C \
$D/block0/iscan.C  \
$D/block11/sscan.C \
$D/block11/doc/sounderDoc.C \
$D/block11/imc/imc.C \
$D/block11/sadId.C \
$D/block11/spacelook/imager/ISblock1.C	\
$D/block11/spacelook/imager/imagerSpacelook.C \
$D/gvar.C \
$D/FrameList.C \
$D/gvarBlock/gvarBlock.C \
$D/gvarHeader/gvarHeader.C \
$D/block1-10/lineDoc.C \
$D/block1-10/writeImagerLine.C \
$D/utils/bitfield/bitfield.C \
$D/utils/cdaTime/cdaTime.C \
$D/utils/crcCheck/crcCheck.C \
$D/utils/fileBuf/fileBin.C \
$D/utils/fileBuf/fileBuf.C \
$D/utils/fileBuf/fileJpeg.C \
$D/utils/fileBuf/fileTiff.C \
$D/utils/fileBuf/fileMcIDAS.C \
$D/utils/fileBuf/genericFile.C \
$D/utils/options/options.C \
$D/utils/radiometric/radiometric.C \
$D/utils/radiometric/setRadiometric.C \
$D/utils/reSample/reSample.C \
$D/utils/selFloat/selFloat.C \
$D/utils/string/stringUtils.C	\
$D/utils/unpack/unpack10.C \
$D/utils/options/channelOptions.C \
$D/utils/options/sectorOptions.C \
$D/block1-10/skipIt.C \
$D/block0/biasStatistics.C \
$D/block0/calibration.C \
$D/Navigation.9411/gimloc3.C \
$D/Navigation.9411/Gimloc3.Om.f \
version.C

OBJ = \
ISblock1.o \
attitudeAngle.o \
bitfield.o \
cdaTime.o \
clampMode.o \
crcCheck.o \
fileBin.o \
fileBuf.o \
fileJpeg.o \
fileTiff.o \
fileMcIDAS.o \
genericFile.o \
grid.o \
gvar.o \
FrameList.o \
gvarBlock.o \
gvarHeader.o \
imagerDoc.o \
imagerSpacelook.o \
imc.o \
iscan.o  \
lineDoc.o \
options.o \
radiometric.o \
reSample.o \
sadId.o \
selFloat.o \
setRadiometric.o \
sscan.o \
sounderDoc.o \
stringUtils.o \
unpack10.o \
version.o \
channelOptions.o \
sectorOptions.o \
writeImagerLine.o \
biasStatistics.o \
calibration.o \
gimloc3.o \
Gimloc3.Om.o \
skipIt.o 


$(TARGET): $(OBJ)
	$(CC) $(CCFLAGS) -o $(TARGET) $(OBJ) $(LIBS)
	rm version.C

install: 
	make depend
	make $(TARGET)
	touch $D/bin/gvar.$(ARCH)
	mv $D/bin/gvar.$(ARCH) $D/bin/gvar.$(ARCH).old
	cp $(TARGET) $(HOME)/bin/gvar.$(ARCH)
clean: 
	rm -f $(TARGET) *.o *~* #* *.bak

depend:
	make version.C
	make readme.h
	make usage.h
	makedepend  -fDepend -- $(CCFLAGS) -- $(SRC) 

include Depend

Gimloc3.Om.o: $D/Navigation.9411/Gimloc3.Om.f
	$D/f2c/src/f2c.$(ARCH) $D/Navigation.9411/Gimloc3.Om.f
	$(cc) -I$D/f2c -c Gimloc3.Om.c 

gimloc3.o : $D/Navigation.9411/gimloc3.C 
	$(CC) $(CCFLAGS) -c  $D/Navigation.9411/gimloc3.C 

biasStatistics.o: $D/block0/biasStatistics.C
	$(CC) $(CCFLAGS) -c  $D/block0/biasStatistics.C	

calibration.o:$D/block0/calibration.C 
	$(CC) $(CCFLAGS) -c  $D/block0/calibration.C 

attitudeAngle.o:	$D/block0/attitudeAngle.C	
	$(CC) $(CCFLAGS) -c $D/block0/attitudeAngle.C	

channelOptions.o:	$D/utils/options/channelOptions.C
	$(CC) $(CCFLAGS) -c $D/utils/options/channelOptions.C

sectorOptions.o:	$D/utils/options/sectorOptions.C
	$(CC) $(CCFLAGS) -c $D/utils/options/sectorOptions.C

clampMode.o : $D/block0/clampMode.C
	$(CC) $(CCFLAGS) -c $D/block0/clampMode.C

imagerDoc.o : $D/block0/imagerDoc.C
	$(CC) $(CCFLAGS) -c $D/block0/imagerDoc.C

iscan.o : $D/block0/iscan.C 
	$(CC) $(CCFLAGS) -c $D/block0/iscan.C 

sscan.o : $D/block11/sscan.C 
	$(CC) $(CCFLAGS) -c $D/block11/sscan.C 

sounderDoc.o : $D/block11/doc/sounderDoc.C 
	$(CC) $(CCFLAGS) -c $D/block11/doc/sounderDoc.C 

imc.o : $D/block11/imc/imc.C
	$(CC) $(CCFLAGS) -c $D/block11/imc/imc.C

sadId.o: $D/block11/sadId.C
	$(CC) $(CCFLAGS) -c $D/block11/sadId.C

ISblock1.o: $D/block11/spacelook/imager/ISblock1.C	
	$(CC) $(CCFLAGS) -c $D/block11/spacelook/imager/ISblock1.C

imagerSpacelook.o: $D/block11/spacelook/imager/imagerSpacelook.C 
	$(CC) $(CCFLAGS) -c $D/block11/spacelook/imager/imagerSpacelook.C 

gvar.o:	$D/gvar.C 
	$(CC) $(CCFLAGS) -c $D/gvar.C 

FrameList.o: $D/FrameList.C
	$(CC) $(CCFLAGS) -c $D/FrameList.C 

writeImagerLine.o:	$D/block1-10/writeImagerLine.C 
	$(CC) $(CCFLAGS) -c $D/block1-10/writeImagerLine.C 

skipIt.o:	$D/block1-10/skipIt.C
	$(CC) $(CCFLAGS) -c $D/block1-10/skipIt.C

gvarBlock.o :	$D/gvarBlock/gvarBlock.C 
	$(CC) $(CCFLAGS) -c $D/gvarBlock/gvarBlock.C 

gvarHeader.o : $D/gvarHeader/gvarHeader.C 
	$(CC) $(CCFLAGS) -c $D/gvarHeader/gvarHeader.C 

lineDoc.o : 	$D/block1-10/lineDoc.C 
	$(CC) $(CCFLAGS) -c 	$D/block1-10/lineDoc.C 

bitfield.o : $D/utils/bitfield/bitfield.C 
	$(CC) $(CCFLAGS) -c $D/utils/bitfield/bitfield.C 

cdaTime.o : $D/utils/cdaTime/cdaTime.C
	$(CC) $(CCFLAGS) -c $D/utils/cdaTime/cdaTime.C

fileBuf.o : $D/utils/fileBuf/fileBuf.C 
	$(CC) $(CCFLAGS) -c $D/utils/fileBuf/fileBuf.C 

fileBin.o : $D/utils/fileBuf/fileBin.C
	$(CC) $(CCFLAGS) -c $D/utils/fileBuf/fileBin.C

fileJpeg.o : $D/utils/fileBuf/fileJpeg.C 
	$(CC) $(CCFLAGS) -c $D/utils/fileBuf/fileJpeg.C 

fileTiff.o : $D/utils/fileBuf/fileTiff.C 
	$(CC) $(CCFLAGS) -c $D/utils/fileBuf/fileTiff.C 

fileMcIDAS.o : $D/utils/fileBuf/fileMcIDAS.C 
	$(CC) $(CCFLAGS) -c $D/utils/fileBuf/fileMcIDAS.C 

genericFile.o : $D/utils/fileBuf/genericFile.C 
	$(CC) $(CCFLAGS) -c $D/utils/fileBuf/genericFile.C 

options.o : $D/utils/options/options.C $D/README $D/USAGE readme.h usage.h 
	$(CC) $(CCFLAGS) -c $D/utils/options/options.C 

radiometric.o : $D/utils/radiometric/radiometric.C 
	$(CC) $(CCFLAGS) -c $D/utils/radiometric/radiometric.C 

setRadiometric.o : $D/utils/radiometric/setRadiometric.C
	$(CC) $(CCFLAGS) -c $D/utils/radiometric/setRadiometric.C

selFloat.o : $D/utils/selFloat/selFloat.C 
	$(CC) $(CCFLAGS) -c $D/utils/selFloat/selFloat.C 

stringUtils.o : $D/utils/string/stringUtils.C	
	$(CC) $(CCFLAGS) -c $D/utils/string/stringUtils.C	

reSample.o : $D/utils/reSample/reSample.C 
	$(CC) $(CCFLAGS) -c $D/utils/reSample/reSample.C 

unpack10.o: $D/utils/unpack/unpack10.C 
	$(CC) $(CCFLAGS) -c $D/utils/unpack/unpack10.C 

grid.o: $D/block0/grid.C
	$(CC) $(CCFLAGS) -c $D/block0/grid.C

crcCheck.o: $D/utils/crcCheck/crcCheck.C 
	$(CC) $(CCFLAGS) -c $D/utils/crcCheck/crcCheck.C 

mkdoc: $D/utils/mkdoc/mkdoc.C
	$(CC) $D/utils/mkdoc/mkdoc.C -o mkdoc

readme.h: $D/README mkdoc
	mkdoc $D/README readme.h

version.C: mkdoc
	date>VERSION
	mkdoc VERSION version.C
	echo "const char * version() { return VERSION; }" | cat >> version.C

version.o: version.C 
	$(CC) $(CCFLAGS) -c version.C

usage.h: $D/USAGE mkdoc 
	mkdoc $D/USAGE usage.h

$D/external_libs/jpeg/libjpeg$(ARCH).a:
	cd $D/external_libs/jpeg; make ; mv libjpeg.a libjpeg$(ARCH).a

$D/external_libs/tiff/libtiff$(ARCH).a:
	cd $D/external_libs/tiff; make ; mv libtiff.a libtiff$(ARCH).a







