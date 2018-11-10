#include <new>

#include "CMManager.h"

/* -------------------------- CMManager --------------------------- */
ConstitutiveModel *CMManager::add_model(ConstitutiveModelParam *param)
{
	if (!param) return nullptr;
	
	switch (param->getType())
	{
	case ConstitutiveModelType::UserSpecific:
		return new(alloc(sizeof(UserSpecificModel))) UserSpecificModel(param);
	case ConstitutiveModelType::LinearElasticity:
		return new(alloc(sizeof(LinearElasticityModel))) LinearElasticityModel(param);
	case ConstitutiveModelType::VonMises:
		return new(alloc(sizeof(VonMisesModel))) VonMisesModel(param);
	case ConstitutiveModelType::Tresca:
		return new(alloc(sizeof(TrescaModel))) TrescaModel(param);
	case ConstitutiveModelType::MohrCoulomb:
		return new(alloc(sizeof(MohrCoulombModel))) MohrCoulombModel(param);
	default:
		return nullptr;
	}
	
	return nullptr;
}
