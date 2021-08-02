#Use the command line:
#qmake -o Makefile TissueFoldingUI.pro
#To generate the Makefile. Then you can build the project either with eclipse, or with
#the commands on the header of the Makefile


#the line necessary for the compilation of non-visual interface model is below. This should be in the makefile inside TissueFolding/Debug/:
#g++ -L/usr/include/  -L/opt/Qt5.2.1/5.2.1/gcc_64/lib  -o "TissueFolding" $(OBJS) $(USER_OBJS) $(LIBS) $(INCLUDEPATH) -lgsl  -lgslcblas -lpardiso500-GNU461-X86-64 -fopenmp -llapack -lpthread
	
#FOR LEGION, when you add a new source/header file to the sourcecode, you do need to update the *.mk files. 
#These include the files under "/home/melda/Documents/TissueFolding/TissueFolding/Debug/" and check legion page
# (https://wiki.rc.ucl.ac.uk/wiki/Development_Tools)
#                              "/home/melda/Documents/TissueFolding/TissueFolding/Debug/SourceCode"
#Install the qt module to use correct qmake first (below s for QT 4, QT 5 is different if you ever need it, :
#module load qt/4.8.6/gnu-4.9.2
#then qmake: 
#qmake -o Makefile TissueFoldingUI.pro
#The procedure is, you need to build the latest version (updates the *.mk files). Then go to the makefile of the non-visual version
# under  "/home/melda/Documents/TissueFolding/TissueFolding/Debug/". Change the compiler line with the line given below. Make the latest non-visual version.
# These steps will ensure the the subdirectory repositories are updated. Then copy all (as listed above) to Legion. The make file of Legion 
# is not the same as the makfile of Tethys. The compiler line is:
#   	g++  -o "TissueFolding" $(OBJS) $(USER_OBJS) $(LIBS) -L${OPENBLASROOT}/lib -lopenblas -lpardiso500-GNU481-X86-64  -fopenmp  -I/shared/ucl/apps/gsl/1.16/gcc/include  -lgsl -lgslcblas
# Then you can make on legion.

#curr path for ubuntu
CurrPath = /home/melda/Documents/TissueFolding/TissueFolding/

#curr path for MacOS:
#CurrPath = ./

#CurrPath for Legion:
#CurrPath = /home/ucgamto/Scratch/Projects/TissueFolding/TissueFolding/

TARGET = $$CurrPath/Debug/TissueFolding

QMAKE_CFLAGS_RELEASE += -fopenmp
QMAKE_CFLAGS_DEBUG += -fopenmp
QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS +=  -fopenmp
 
# Input
OBJECTS_DIR += $$CurrPath/Debug/

HEADERS       +=  	$$CurrPath/SourceCode/*.h

SOURCES += $$CurrPath/SourceCode/main.cpp \
	$$CurrPath/SourceCode/Prism.cpp \
	$$CurrPath/SourceCode/ReferenceShapeBase.cpp \
 	$$CurrPath/SourceCode/ShapeBase.cpp \
        $$CurrPath/SourceCode/Simulation.cpp \
        $$CurrPath/SourceCode/Node.cpp \
        $$CurrPath/SourceCode/ModelInputObject.cpp \
	$$CurrPath/SourceCode/RandomGenerator.cpp \
	$$CurrPath/SourceCode/NewtonRaphsonSolver.cpp \
	$$CurrPath/SourceCode/Analysis.cpp \
	$$CurrPath/SourceCode/CellMigration.cpp

#libs and includes for linux for independent license pardiso:
#LIBS += -L/usr/include -lgsl -lgslcblas -lpardiso500-GNU461-X86-64  -fopenmp -llapack  ---- Old pardiso
LIBS += -L/usr/include -lgsl -lgslcblas -lgomp -fopenmp -lpardiso600-GNU720-X86-64   -llapack 

#libs and includes for linux for MKL: --- clashes with gsl cblas!!!
#LIBS += -L/usr/include -lgsl  -lpardiso500-GNU461-X86-64  -fopenmp -llapack -DMKL_ILP64 -m64 -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed  -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
#INCLUDEPATH += ${MKLROOT}/include


# libs and includes for MacOS
# LIBS += -L/usr/include -L/usr/local/lib/ -lgsl -lgslcblas -L/usr/local/Cellar/boost/1.58.0/include -lpardiso500-MACOS-X86-64
# INCLUDEPATH += /usr/local/Cellar/boost/1.58.0/include /usr/local/include/


# install
target.path   =  $$CurrPath/
sources.files =  $$SOURCES $$HEADERS $$RESOURCES TissueFolding.pro
sources.path  =  $$CurrPath
INSTALLS += target sources
