#ifdef _DEBUG
#include <iostream>
#endif

#include "LinearElasticityModel.h"

// --------------------------- Linear Elasticity -----------------------------
void LinearElasticityModel::calStiffnessMatrix()
{
	// Elastic stiffness matrix
	double elastic_tmp;
	elastic_tmp = E / (1.0 + nu) / (1.0 - 2.0 * nu);
	De[0][0] = De[1][1] = De[2][2] = elastic_tmp * (1.0 - nu);
	De[0][1] = De[0][2] = De[1][0] = De[1][2] = De[2][0] = De[2][1] = elastic_tmp * nu;
	De[3][3] = De[4][4] = De[5][5] = E / (1 + nu); // 2G
	De[0][3] = De[0][4] = De[0][5] = 0.0;
	De[1][3] = De[1][4] = De[1][5] = 0.0;
	De[2][3] = De[2][4] = De[2][5] = 0.0;
	De[3][0] = De[3][1] = De[3][2] = De[3][4] = De[3][5] = 0.0;
	De[4][0] = De[4][1] = De[4][2] = De[4][3] = De[4][5] = 0.0;
	De[5][0] = De[5][1] = De[5][2] = De[5][3] = De[5][4] = 0.0;


}

void LinearElasticityModel::integration_explicit(
	const double dstrain[6], const double stress[6],
	double dstress[6], double destrain[6], double dpstrain[6])
{
	dstress[0] = De[0][0] * dstrain[0] + De[0][1] * dstrain[1] + De[0][2] * dstrain[2];
	dstress[1] = De[1][0] * dstrain[0] + De[1][1] * dstrain[1] + De[1][2] * dstrain[2];
	dstress[2] = De[2][0] * dstrain[0] + De[2][1] * dstrain[1] + De[2][2] * dstrain[2];
	dstress[3] = De[3][3] * dstrain[3];
	dstress[4] = De[4][4] * dstrain[4];
	dstress[5] = De[5][5] * dstrain[5];

	destrain[0] = dstrain[0];
	destrain[1] = dstrain[1];
	destrain[2] = dstrain[2];
	destrain[3] = dstrain[3];
	destrain[4] = dstrain[4];
	destrain[5] = dstrain[5];

	dpstrain[0] = 0.0;
	dpstrain[1] = 0.0;
	dpstrain[2] = 0.0;
	dpstrain[3] = 0.0;
	dpstrain[4] = 0.0;
	dpstrain[5] = 0.0;

};

void LinearElasticityModel::integration_implicit(
	const double dstrain[6], const double stress[6],
	double dstress[6], double destrain[6], double dpstrain[6],
	double tangentialMatrix[6][6])
{
	dstress[0] = De[0][0] * dstrain[0] + De[0][1] * dstrain[1] + De[0][2] * dstrain[2];
	dstress[1] = De[1][0] * dstrain[0] + De[1][1] * dstrain[1] + De[1][2] * dstrain[2];
	dstress[2] = De[2][0] * dstrain[0] + De[2][1] * dstrain[1] + De[2][2] * dstrain[2];
	dstress[3] = De[3][3] * dstrain[3];
	dstress[4] = De[4][4] * dstrain[4];
	dstress[5] = De[5][5] * dstrain[5];

	destrain[0] = dstrain[0];
	destrain[1] = dstrain[1];
	destrain[2] = dstrain[2];
	destrain[3] = dstrain[3];
	destrain[4] = dstrain[4];
	destrain[5] = dstrain[5];

	dpstrain[0] = 0.0;
	dpstrain[1] = 0.0;
	dpstrain[2] = 0.0;
	dpstrain[3] = 0.0;
	dpstrain[4] = 0.0;
	dpstrain[5] = 0.0;

	tangentialMatrix[0][0] = De[0][0];
	tangentialMatrix[0][1] = De[0][1];
	tangentialMatrix[0][2] = De[0][2];
	tangentialMatrix[0][3] = De[0][3];
	tangentialMatrix[0][4] = De[0][4];
	tangentialMatrix[0][5] = De[0][5];

	tangentialMatrix[1][0] = De[0][0];
	tangentialMatrix[1][1] = De[1][1];
	tangentialMatrix[1][2] = De[1][2];
	tangentialMatrix[1][3] = De[1][3];
	tangentialMatrix[1][4] = De[1][4];
	tangentialMatrix[1][5] = De[1][5];

	tangentialMatrix[2][0] = De[2][0];
	tangentialMatrix[2][1] = De[2][1];
	tangentialMatrix[2][2] = De[2][2];
	tangentialMatrix[2][3] = De[2][3];
	tangentialMatrix[2][4] = De[2][4];
	tangentialMatrix[2][5] = De[2][5];

	tangentialMatrix[3][0] = De[3][0];
	tangentialMatrix[3][1] = De[3][1];
	tangentialMatrix[3][2] = De[3][2];
	tangentialMatrix[3][3] = De[3][3];
	tangentialMatrix[3][4] = De[3][4];
	tangentialMatrix[3][5] = De[3][5];

	tangentialMatrix[4][0] = De[4][0];
	tangentialMatrix[4][1] = De[4][1];
	tangentialMatrix[4][2] = De[4][2];
	tangentialMatrix[4][3] = De[4][3];
	tangentialMatrix[4][4] = De[4][4];
	tangentialMatrix[4][5] = De[4][5];

	tangentialMatrix[5][0] = De[5][0];
	tangentialMatrix[5][1] = De[5][1];
	tangentialMatrix[5][2] = De[5][2];
	tangentialMatrix[5][3] = De[5][3];
	tangentialMatrix[5][4] = De[5][4];
	tangentialMatrix[5][5] = De[5][5];
};

#ifdef _DEBUG

void LinearElasticityModel::showParam()
{
	size_t i, j;

	std::cout << "E: " << E << "; nu: " << nu << std::endl;

	for (i = 0; i < 6; i++)
	{
		for (j = 0; j < 6; j ++)
		{
			std::cout << De[i][j] << " ";
		}
		std::cout << std::endl;
	}
}
#endif // _DEBUG