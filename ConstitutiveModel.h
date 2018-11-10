#ifndef __CONSTITUTIVEMODEL_H__
#define __CONSTITUTIVEMODEL_H__

// Add the name of all consitutive model
enum class ConstitutiveModelType : unsigned int
{
	Invalid = 0,
	UserSpecific = 1,
	LinearElasticity = 2,
	VonMises = 3,
	Tresca = 4,
	MohrCoulomb = 5
};

class ConstitutiveModel;

// Used to initialized cm data
struct ConstitutiveModelParam
{
	friend ConstitutiveModel;
protected:
	ConstitutiveModelType type;
public:
	ConstitutiveModelParam(
		ConstitutiveModelType tp = ConstitutiveModelType::Invalid) :
		type(tp) {}
	inline ConstitutiveModelType getType() noexcept { return type; };
};

// base class of constitutiveModel
class ConstitutiveModel
{
protected:
	// type of this model
	ConstitutiveModelType type;
	/*
	1. Stress - strain integration function for explicit algorithm
	does not generate tangential matrix.
	2. dstrain[] is the input parameter:
	dstrain[0]: strain11;
	dstrain[1]: strain22;
	dstrain[2]: strain33;
	dstrain[3]: strain12;
	dstrain[4]: strain23;
	dstrain[5]: strain31.
	
	typedef void (ConstitutiveModel::*integrationExplicit) (
		const double dstrain[6],
		const double stress[6],
		double dstress[6],
		double destrain[6],
		double dpstrain[6]);
	// Pointer of integration function, this pointer is for:
	//     1. Dynamic polymorphism of consitutive model class;
	//     2. Easier change of integration algorithm.
	integrationExplicit intExpFunc;
	
	1. Stress - strain integration function for explicit algorithm
	does not generate tangential matrix.
	2. dstrain[] is the input parameter:
	dstrain[0]: strain11;
	dstrain[1]: strain22;
	dstrain[2]: strain33;
	dstrain[3]: strain12;
	dstrain[4]: strain23;
	dstrain[5]: strain31.
	
	typedef void (ConstitutiveModel::*integrationImplicit) (
		const double dstrain[6],
		const double stress[6],
		double dstress[6],
		double destrain[6],
		double dpstrain[6],
		double tangentialMatrix[6][6]);
	// Pointer of integration function, this pointer is for:
	//     1. Dynamic polymorphism of consitutive model class;
	//     2. Easier change of integration algorithm.
	integrationImplicit intImpFunc;*/

public:
	ConstitutiveModel(ConstitutiveModelParam *param) :
		type(param->type) {}
	~ConstitutiveModel() {}
	inline ConstitutiveModelType getType() { return type; }
	// Set parameters of the constitutive model
	virtual void setParam(ConstitutiveModelParam *param) = 0;

	/*
		1. Stress - strain integration function for explicit
	algorithm (does not generate tangential matrix).
		2. Input parameter: dstrain[], stress[],
			1) dstrain[]: strain increment
				dstrain[0]: strain11, dstrain[1]: strain22,
				dstrain[2]: strain33, dstrain[3]: strain12,
				dstrain[4]: strain23, dstrain[5]: strain31.
			2) stress[]: stress state
				stress[0]: stress11, stress[1]: stress22,
				stress[2]: stress33, stress[3]: stress12,
				stress[4]: stress23, stress[5]: stress31.
	*/
	virtual void integration_explicit(
		const double dstrain[6], const double stress[6],
		double dstress[6], double destrain[6], double dpstrain[6]) = 0;
	/*
		1. Stress - strain integration function for implicit algorithm.
		2. Input parameter: dstrain[], stress[],
			1) dstrain[]: strain increment
				dstrain[0]: strain11, dstrain[1]: strain22,
				dstrain[2]: strain33, dstrain[3]: strain12,
				dstrain[4]: strain23, dstrain[5]: strain31.
			2) stress[]: stress state
				stress[0]: stress11, stress[1]: stress22,
				stress[2]: stress33, stress[3]: stress12,
				stress[4]: stress23, stress[5]: stress31.
	*/
	virtual void integration_implicit(
		const double dstrain[6], const double stress[6],
		double dstress[6], double destrain[6],
		double dpstrain[6], double tangentialMatrix[6][6]) = 0;
};

#endif