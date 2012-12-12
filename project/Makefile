SRC = ./src
BIN = ./bin
CBIN = $(BIN)/cavity_c
CFILE = $(SRC)/cavity-solver.c
FFILE = $(SRC)/cavity.template-to-students.f95
CPPFILE = $(SRC)/DrivenCavity.template-to-students.cpp

clean: 
	rm $(BIN)/*

cleantxt:
	rm *.dat
	rm *.out
	rm *.tec
	
c3: $(CFILE)
	clang -O3 $(CFILE) -o $(BIN)/cavity_c3 -lm

f3: $(FFILE)
	gfortran -O3 $(FFILE) -o $(BIN)/cavity_f3
	mv *.mod $(BIN)

c: $(CFILE)
	clang $(CFILE) -o $(CBIN) -lm

f: $(FFILE)
	gfortran $(FFILE) -o $(BIN)/cavity_f
	mv *.mod $(BIN)

run: c
	$(CBIN)

run3 :c3
	$(BIN)/cavity_c3
