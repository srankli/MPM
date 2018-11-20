#include <cstdlib>
#include <vector>
#include <set>
#include <string>
#include <iostream>

#include "BoundaryCondition.h"

#include "Solver_1D_Mechanics_1D2_Explicit_FixedMem.h"
#include "Solver_1D_Hydromechanics_1D2_Explicit_FixedMem.h"
#include "Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem.h"

#include "TimeCurve.h"

int main(void)
{
	//extern void test_ReadTmpDataFile();
	//extern void test_WriteTmpDataFile();
	//test_WriteTmpDataFile();
	//test_ReadTmpDataFile();

	//extern void test_Solver_1D_Mechanics_1D2_Explicit_FixedMem(void);
	//test_Solver_1D_Mechanics_1D2_Explicit_FixedMem();

	//extern void test_Solver_1D_Hydromechanics_1D2_Explicit_FixedMem_2(void);
	//test_Solver_1D_Hydromechanics_1D2_Explicit_FixedMem_2();

	//extern void test_Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem_4(void);
	//test_Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem_4();

	//extern void test_forceBC(void);
	//test_forceBC();

	//extern void test_FileBuffer(void);
	//test_FileBuffer();

	extern void test_OutputRequest(void);
	test_OutputRequest();

	//extern void test_TmpDataToHdf5(void);
	//test_TmpDataToHdf5();

	system("pause");
	return 0;
}