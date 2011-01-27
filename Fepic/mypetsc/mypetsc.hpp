#ifndef MYPETSC_HPP
#define MYPETSC_HPP

#include "petscsnes.h"
#include "petscksp.h"
#include <iostream>
#include <iomanip>



void inline MatrixInfo(Mat& K) {
	MatInfo info;
	double  mal, nz_a, nz_u;
	MatGetInfo(K,MAT_LOCAL,&info);
	mal  = info.nz_allocated;
	nz_a = info.nz_used;
	nz_u = info.nz_unneeded;
	std::cout << mal << " allocated\t\t" << nz_a << " used\t" << nz_u << " unneeded" << std::endl;
}

void inline View(Mat &K, const char* name = "matrix.m") {

	PetscViewer viewer;
	PetscViewerASCIIOpen(PETSC_COMM_WORLD, name, &viewer);
    PetscObjectSetName((PetscObject)viewer,"a");
	PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
	//PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_DENSE);
	MatView(K, viewer);
	PetscViewerDestroy(viewer);
}

void inline View(Vec &v, const char* name = "VETOR") {

    PetscViewer viewer;
	PetscViewerASCIIOpen(PETSC_COMM_WORLD, name, &viewer);
    PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB); 
	VecView(v, viewer);
	PetscViewerDestroy(viewer);
}

// OBSOLETO D+++
void inline ViewXPM(Mat &K, const char* name = "MATRIX") {
    std::ofstream file(name);
    int width, heigth;
    
    MatGetSize(K, &heigth, &width);
    
    file << "/* XPM */" << std::endl;
    file << "static char * blarg_xpm[] = {" << std::endl;
    file << "\"" <<  width << " " << heigth << " 3 1\",";
    file << "\"a c #000000\"," << std::endl;
    file << "\"b c #FFFFFF\"," << std::endl;
    file << "\"c c #FF0000\"," << std::endl;
    
    for (int i=0;i<heigth-1;i++) {
        file << "\"";
        double val;
        for (int j=0; j<width; j++) {
            MatGetValue(K,i,j,&val);
            if (val==1) file << "c";
            else if (val != 0) file << "b";
            else file << "a";
        }
        file << "\"," << std::endl;
    }
    
    file << "\"";
    for (int j=0; j<width; j++) {
        double val;
        MatGetValue(K,heigth-1,j,&val);
        if (val==1) file << "c";
        else if (val != 0) file << "b";
        else file << "a";
    }
    file << "\"" << std::endl;
    file << "}";
    
    file.close();
    
}

void inline Assembly(Mat &K) {
	MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);
};

void inline Assembly(Vec &v) {
	VecAssemblyBegin(v);
	VecAssemblyEnd(v);
}

void inline Destroy(Mat &K) {
	MatDestroy(K);
}

void inline Destroy(Vec &v) {
	VecDestroy(v);
}

void inline Destroy(PC &pc) {
	PCDestroy(pc);
}

void inline Destroy(KSP &ksp) {
	KSPDestroy(ksp);
}

void inline Destroy(SNES &snes) {
	SNESDestroy(snes);
}

#endif

