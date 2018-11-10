#include <new>

#include "BoundaryCondition.h"

int AccelerationBCManager::add_bc(AccelerationBCParam *param)
{
	void *ptr;
	AccelerationBCParam_Fixed *abc_param1;
	AccelerationBC_Fixed *abc1;
	AccelerationBCParam_Constant *abc_param2;
	AccelerationBC_Constant *abc2;
	AccelerationBCParam_Linear *abc_param3;
	AccelerationBC_Linear *abc3;
	AccelerationBCParam_CubicSmooth *abc_param4;
	AccelerationBC_CubicSmooth *abc4;

	if (!param) return -1;

	switch (param->getType())
	{
	case TimeCurveType::Zero:
		abc_param1 = static_cast<AccelerationBCParam_Fixed *>(param);
		ptr = alloc(sizeof(AccelerationBC_Fixed));
		if (!ptr) return -3;
		abc1 = new(ptr) AccelerationBC_Fixed;
		abc1->index = abc_param1->index;
		abc1->dof = abc_param1->dof;
		abc1->init();
		return 0;
	case TimeCurveType::Constant:
		abc_param2 = static_cast<AccelerationBCParam_Constant *>(param);
		ptr = alloc(sizeof(AccelerationBC_Constant));
		if (!ptr) return -3;
		abc2 = new(ptr) AccelerationBC_Constant;
		abc2->index = abc_param2->index;
		abc2->dof = abc_param2->dof;
		abc2->curve.x_const = abc_param2->a;
		abc2->init();
		return 0;
	case TimeCurveType::Linear:
		abc_param3 = static_cast<AccelerationBCParam_Linear *>(param);
		ptr = alloc(sizeof(AccelerationBC_Linear));
		if (!ptr) return -3;
		abc3 = new(ptr) AccelerationBC_Linear;
		abc3->index = abc_param3->index;
		abc3->dof = abc_param3->dof;
		abc3->curve.x_begin = abc_param3->a_begin;
		abc3->curve.x_end = abc_param3->a_end;
		abc3->curve.t_len = abc_param3->t_length;
		abc3->init();
		return 0;
	case TimeCurveType::CubicSmooth:
		abc_param4 = static_cast<AccelerationBCParam_CubicSmooth *>(param);
		ptr = alloc(sizeof(AccelerationBC_CubicSmooth));
		if (!ptr) return -3;
		abc4 = new(ptr) AccelerationBC_CubicSmooth;
		abc4->index = abc_param4->index;
		abc4->dof = abc_param4->dof;
		abc4->curve.x_begin = abc_param4->a_begin;
		abc4->curve.x_end = abc_param4->a_end;
		abc4->curve.t_len = abc_param4->t_length;
		abc4->init();
		return 0;
	default:
		break;
	}

	return -2;
}

int VelocityBCManager::add_bc(VelocityBCParam *param)
{
	void *ptr;
	VelocityBCParam_Fixed *vbc_param1;
	VelocityBC_Fixed *vbc1;
	VelocityBCParam_Constant *vbc_param2;
	VelocityBC_Constant *vbc2;
	VelocityBCParam_Linear *vbc_param3;
	VelocityBC_Linear *vbc3;
	VelocityBCParam_CubicSmooth *vbc_param4;
	VelocityBC_CubicSmooth *vbc4;
	
	if (!param) return -1;

	switch (param->getType())
	{
	case TimeCurveType::Zero:
		vbc_param1 = static_cast<VelocityBCParam_Fixed *>(param);
		ptr = alloc(sizeof(VelocityBC_Fixed));
		if (!ptr) return -3;
		vbc1 = new(ptr) VelocityBC_Fixed;
		vbc1->index = vbc_param1->index;
		vbc1->dof = vbc_param1->dof;
		vbc1->init();
		return 0;
	case TimeCurveType::Constant:
		vbc_param2 = static_cast<VelocityBCParam_Constant *>(param);
		ptr = alloc(sizeof(VelocityBC_Constant));
		if (!ptr) return -3;
		vbc2 = new(ptr) VelocityBC_Constant;
		vbc2->index = vbc_param2->index;
		vbc2->dof = vbc_param2->dof;
		vbc2->curve.x_const = vbc_param2->v;
		vbc2->init();
		return 0;
	case TimeCurveType::Linear:
		vbc_param3 = static_cast<VelocityBCParam_Linear *>(param);
		ptr = alloc(sizeof(VelocityBC_Linear));
		if (!ptr) return -3;
		vbc3 = new(ptr) VelocityBC_Linear;
		vbc3->index = vbc_param3->index;
		vbc3->dof = vbc_param3->dof;
		vbc3->curve.x_begin = vbc_param3->v_begin;
		vbc3->curve.x_end = vbc_param3->v_end;
		vbc3->curve.t_len = vbc_param3->t_length;
		vbc3->init();
		return 0;
	case TimeCurveType::CubicSmooth:
		vbc_param4 = static_cast<VelocityBCParam_CubicSmooth *>(param);
		ptr = alloc(sizeof(VelocityBC_CubicSmooth));
		if (!ptr) return -3;
		vbc4 = new(ptr) VelocityBC_CubicSmooth;
		vbc4->index = vbc_param4->index;
		vbc4->dof = vbc_param4->dof;
		vbc4->curve.x_begin = vbc_param4->v_begin;
		vbc4->curve.x_end = vbc_param4->v_end;
		vbc4->curve.t_len = vbc_param4->t_length;
		vbc4->init();
		return 0;
	default:
		break;
	}

	return -2;
}


int MassForceBCManager::add_bc(MassForceBCParam *param)
{
	void *ptr;
	MassForceBCParam_Constant *mfbc_param2;
	MassForceBC_Constant *mfbc2;
	MassForceBCParam_Linear *mfbc_param3;
	MassForceBC_Linear *mfbc3;
	MassForceBCParam_CubicSmooth *mfbc_param4;
	MassForceBC_CubicSmooth *mfbc4;

	if (!param) return -1;

	switch (param->getType())
	{
	case TimeCurveType::Constant:
		mfbc_param2 = static_cast<MassForceBCParam_Constant *>(param);
		ptr = alloc(sizeof(MassForceBC_Constant));
		if (!ptr) return -3;
		mfbc2 = new(ptr) MassForceBC_Constant;
		mfbc2->index = mfbc_param2->index;
		mfbc2->dof = mfbc_param2->dof;
		mfbc2->curve.x_const = mfbc_param2->massForce;
		mfbc2->init();
		return 0;
	case TimeCurveType::Linear:
		mfbc_param3 = static_cast<MassForceBCParam_Linear *>(param);
		ptr = alloc(sizeof(MassForceBC_Linear));
		if (!ptr) return -3;
		mfbc3 = new(ptr) MassForceBC_Linear;
		mfbc3->index = mfbc_param3->index;
		mfbc3->dof = mfbc_param3->dof;
		mfbc3->curve.x_begin = mfbc_param3->massForce_begin;
		mfbc3->curve.x_end = mfbc_param3->massForce_end;
		mfbc3->curve.t_len = mfbc_param3->t_length;
		mfbc3->init();
		return 0;
	case TimeCurveType::CubicSmooth:
		mfbc_param4 = static_cast<MassForceBCParam_CubicSmooth *>(param);
		ptr = alloc(sizeof(MassForceBC_CubicSmooth));
		if (!ptr) return -3;
		mfbc4 = new(ptr) MassForceBC_CubicSmooth;
		mfbc4->index = mfbc_param4->index;
		mfbc4->dof = mfbc_param4->dof;
		mfbc4->curve.x_begin = mfbc_param4->massForce_begin;
		mfbc4->curve.x_end = mfbc_param4->massForce_end;
		mfbc4->curve.t_len = mfbc_param4->t_length;
		mfbc4->init();
		return 0;
	default:
		break;
	}

	return -2;
}


int SurfaceForceBCManager::add_bc(SurfaceForceBCParam *param)
{
	void *ptr;
	SurfaceForceBCParam_Constant *sfbc_param2;
	SurfaceForceBC_Constant *sfbc2;
	SurfaceForceBCParam_Linear *sfbc_param3;
	SurfaceForceBC_Linear *sfbc3;
	SurfaceForceBCParam_CubicSmooth *sfbc_param4;
	SurfaceForceBC_CubicSmooth *sfbc4;

	if (!param) return -1;

	switch (param->getType())
	{
	case TimeCurveType::Constant:
		sfbc_param2 = static_cast<SurfaceForceBCParam_Constant *>(param);
		ptr = alloc(sizeof(SurfaceForceBC_Constant));
		if (!ptr) return -3;
		sfbc2 = new(ptr) SurfaceForceBC_Constant;
		sfbc2->index = sfbc_param2->index;
		sfbc2->dof = sfbc_param2->dof;
		sfbc2->curve.x_const = sfbc_param2->surfaceForce;
		sfbc2->init();
		return 0;
	case TimeCurveType::Linear:
		sfbc_param3 = static_cast<SurfaceForceBCParam_Linear *>(param);
		ptr = alloc(sizeof(SurfaceForceBC_Linear));
		if (!ptr) return -3;
		sfbc3 = new(ptr) SurfaceForceBC_Linear;
		sfbc3->index = sfbc_param3->index;
		sfbc3->dof = sfbc_param3->dof;
		sfbc3->curve.x_begin = sfbc_param3->surfaceForce_begin;
		sfbc3->curve.x_end = sfbc_param3->surfaceForce_end;
		sfbc3->curve.t_len = sfbc_param3->t_length;
		sfbc3->init();
		return 0;
	case TimeCurveType::CubicSmooth:
		sfbc_param4 = static_cast<SurfaceForceBCParam_CubicSmooth *>(param);
		ptr = alloc(sizeof(SurfaceForceBC_CubicSmooth));
		if (!ptr) return -3;
		sfbc4 = new(ptr) SurfaceForceBC_CubicSmooth;
		sfbc4->index = sfbc_param4->index;
		sfbc4->dof = sfbc_param4->dof;
		sfbc4->curve.x_begin = sfbc_param4->surfaceForce_begin;
		sfbc4->curve.x_end = sfbc_param4->surfaceForce_end;
		sfbc4->curve.t_len = sfbc_param4->t_length;
		sfbc4->init();
		return 0;
	default:
		break;
	}

	return -2;
}
