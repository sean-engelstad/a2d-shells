include ../../Makefile.in
include ../../A2D_SHELLS.mk

CXX_OBJS = AMDInterface.o \
           TACSBVecDistribute.o \
		   TACSNodeMap.o \
		   TACSParallelMat.o \
		   KSM.o \
		   TACSMatDistribute.o \
		   TACSBlockCyclicMat.o \
		   TACSBVec.o \
		   TACSSchurMat.o \
		   TACSSerialPivotMat.o \
		   BCSCMatPivot.o \
		   BCSRMat.o \
		   BCSRMatFact.o \
			BCSRMatMult.o \
			BCSRMatFact1.o \
			BCSRMatMult1.o \
			BCSRMatFact2.o \
			BCSRMatMult2.o \
			BCSRMatFact3.o \
			BCSRMatMult3.o \
			BCSRMatFact4.o \
			BCSRMatMult4.o \
			BCSRMatFact5.o \
			BCSRMatMult5.o \
			BCSRMatFact6.o \
			BCSRMatMult6.o \
			BCSRMatFact8.o \
			BCSRMatMult8.o \
			BCSCMatPivot.o

default: ${CXX_OBJS}

clean:
	rm -f *.o
	rm -f **/*.o
	rm -f **/**/*.o