#ifndef __CMMANAGER_H__
#define __CMMANAGER_H__

#include "MemoryManager.h"

// Include all constitutive model here
#include "ConstitutiveModel.h"
#include "LinearElasticityModel.h"
//#include "VonMisesModel.h"
//#include "TrescaModel.h"
//#include "MohrCoulombModel.h"

class CMManager : public FlexibleSizeMemory
{
public:
	CMManager() {}
	~CMManager() {}
	
	ConstitutiveModel *add_model(ConstitutiveModelParam *param);
	inline ConstitutiveModel *get_first() noexcept
	{
		return static_cast<ConstitutiveModel *>(FlexibleSizeMemory::get_first());
	}
	inline ConstitutiveModel *get_next(ConstitutiveModel *prev_model) noexcept
	{
		return static_cast<ConstitutiveModel *>(FlexibleSizeMemory::get_next(prev_model));
	}
};

#endif