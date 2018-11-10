#ifndef __LINEARELASTICITY_H__
#define __LINEARELASTICITY_H__

#include "ConstitutiveModel.h"

// ------------------------------ User Specific Model ------------------------------
struct UserSpecificModelParam : public ConstitutiveModelParam
{
	UserSpecificModelParam() :
		ConstitutiveModelParam(ConstitutiveModelType::UserSpecific) {}
};

class UserSpecificModel : public ConstitutiveModel
{
public:
	UserSpecificModel(ConstitutiveModelParam *param) :
		ConstitutiveModel(param) {}
	void integration_explicit(
		const double dstrain[6], const double stress[6],
		double dstress[6], double destrain[6], double dpstrain[6]) {}
	void integration_implicit(
		const double dstrain[6], const double stress[6],
		double dstress[6], double destrain[6], double dpstrain[6],
		double tangentialMatrix[6][6]) {}
	void setParam(ConstitutiveModelParam *param)
	{
		UserSpecificModelParam *pa
			= static_cast<UserSpecificModelParam *>(param);
	}
};


// ------------------------------ LinearElasticity Model ------------------------------
struct LinearElasticityModelParam : public ConstitutiveModelParam
{
	LinearElasticityModelParam() :
		ConstitutiveModelParam(ConstitutiveModelType::LinearElasticity) {}
	double E;
	double nu;
};

class LinearElasticityModel : public ConstitutiveModel
{
protected:
	double E; // Young's modulus
	double nu; // Possion's ratio
	double De[6][6]; // Stiffness matrix

public:
	LinearElasticityModel(ConstitutiveModelParam *param) :
		ConstitutiveModel(param) { setParam(param); }
	~LinearElasticityModel() {}
	void integration_explicit(
		const double dstrain[6], const double stress[6],
		double dstress[6], double destrain[6], double dpstrain[6]);
	void integration_implicit(
		const double dstrain[6], const double stress[6],
		double dstress[6], double destrain[6], double dpstrain[6],
		double tangentialMatrix[6][6]);
	void setParam(ConstitutiveModelParam *param)
	{
		LinearElasticityModelParam *pa
			= static_cast<LinearElasticityModelParam *>(param);
		E = pa->E, nu = pa->nu;
		calStiffnessMatrix();
	}
#ifdef _DEBUG
	void showParam();
#endif // _DEBUG

protected:
	void calStiffnessMatrix();
};

// ------------------------------ VonMises Model ------------------------------
struct VonMisesModelParam : public ConstitutiveModelParam
{
	VonMisesModelParam() :
		ConstitutiveModelParam(ConstitutiveModelType::VonMises) {}
	double E;
	double nu;
};

class VonMisesModel : public ConstitutiveModel
{
public:
	VonMisesModel(ConstitutiveModelParam *param) :
		ConstitutiveModel(param) {}
	void integration_explicit(
		const double dstrain[6], const double stress[6],
		double dstress[6], double destrain[6], double dpstrain[6]) {}
	void integration_implicit(
		const double dstrain[6], const double stress[6],
		double dstress[6], double destrain[6], double dpstrain[6],
		double tangentialMatrix[6][6]) {}
	void setParam(ConstitutiveModelParam *param)
	{
		VonMisesModelParam *pa = static_cast<VonMisesModelParam *>(param);
	}
};


// ------------------------------ Tresca Model ------------------------------
struct TrescaModelParam : public ConstitutiveModelParam
{
	TrescaModelParam() :
		ConstitutiveModelParam(ConstitutiveModelType::Tresca) {}
	double E;
	double nu;
};

class TrescaModel : public ConstitutiveModel
{
public:
	TrescaModel(ConstitutiveModelParam *param) :
		ConstitutiveModel(param) {}
	void integration_explicit(
		const double dstrain[6], const double stress[6],
		double dstress[6], double destrain[6], double dpstrain[6]) {}
	void integration_implicit(
		const double dstrain[6], const double stress[6],
		double dstress[6], double destrain[6], double dpstrain[6],
		double tangentialMatrix[6][6]) {}
	void setParam(ConstitutiveModelParam *param)
	{
		TrescaModelParam *pa = static_cast<TrescaModelParam *>(param);
	}
};


// ------------------------------ MohrCoulomb Model ------------------------------
struct MohrCoulombModelParam : public ConstitutiveModelParam
{
	MohrCoulombModelParam() :
		ConstitutiveModelParam(ConstitutiveModelType::MohrCoulomb) {}
	double E;
	double nu;
};

class MohrCoulombModel : public ConstitutiveModel
{
public:
	MohrCoulombModel(ConstitutiveModelParam *param) :
		ConstitutiveModel(param) {}
	void integration_explicit(
		const double dstrain[6], const double stress[6],
		double dstress[6], double destrain[6], double dpstrain[6]) {}
	void integration_implicit(
		const double dstrain[6], const double stress[6],
		double dstress[6], double destrain[6], double dpstrain[6],
		double tangentialMatrix[6][6]) {}
	void setParam(ConstitutiveModelParam *param)
	{
		MohrCoulombModelParam *pa = static_cast<MohrCoulombModelParam *>(param);
	}
};

#endif