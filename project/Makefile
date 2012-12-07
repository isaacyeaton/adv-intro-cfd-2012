SRC = ./src
BIN = ./bin
CFILE = $(SRC)/cavity.template-to-students.c
FFILE = $(SRC)/cavity.template-to-students.f95
CPPFILE = $(SRC)/DrivenCavity.template-to-students.cpp

clean: 
	rm $(BIN)/*
	
c3: $(CFILE)
	gcc -O3 $(CFILE) -o $(BIN)/cavity_c3 -lm

f3: $(FFILE)
	#rm $(BIN)/*.mod
	gfortran -O3 $(FFILE) -o $(BIN)/cavity_f3
	mv *.mod $(BIN)

f: $(FFILE)
	#rm $(BIN)/*.mod
	gfortran $(FFILE) -o $(BIN)/cavity_f
	#clang $(FFILE) -o $(BIN)/cavity_f
	mv *.mod $(BIN)
