CFILE = cavity.template-to-students.c
FFILE = cavity.template-to-students.f95
CPPFILE = DrivenCavity.template-to-students.cpp

compc: $(CFILE)
	gcc -O3 $(CFILE) -o cavity_c -lm

compf: $(FFILE)
	gfortran -O3 $(FFILE) -o cavity_f
