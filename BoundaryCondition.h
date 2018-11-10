#ifndef __BOUNDARYCONDITION_H__
#define __BOUNDARYCONDITION_H__

#include "TimeCurve.h"
#include "MemoryManager.h"

enum class DegreeOfFreedom : unsigned long long
{
	invalid = 0,
	x = 1,
	y = 2,
	z = 3,
	x_f = 4, // fluid phase
	y_f = 5, // fluid phase
	z_f = 6  // fluid phase
};

// ----------------- Acceleration Boundary Condition ------------------
// Parameters for initializing acceleration boundary conditions
struct AccelerationBCParam
{
protected:
	TimeCurveType type;
public:
	size_t index;
	DegreeOfFreedom dof;
	AccelerationBCParam(TimeCurveType tp = TimeCurveType::Invalid) : type(tp) {}
	inline TimeCurveType getType() noexcept { return type; }
};

struct AccelerationBCParam_Fixed : public AccelerationBCParam
{
	AccelerationBCParam_Fixed() : AccelerationBCParam(TimeCurveType::Zero) {}
};

struct AccelerationBCParam_Constant : public AccelerationBCParam
{
	double a;
	AccelerationBCParam_Constant() : AccelerationBCParam(TimeCurveType::Constant) {}
};

struct AccelerationBCParam_Linear : public AccelerationBCParam
{
	double a_begin;
	double a_end;
	double t_length;
	AccelerationBCParam_Linear() : AccelerationBCParam(TimeCurveType::Linear) {}
};

struct AccelerationBCParam_CubicSmooth : public AccelerationBCParam
{
	double a_begin;
	double a_end;
	double t_length;
	AccelerationBCParam_CubicSmooth() : AccelerationBCParam(TimeCurveType::CubicSmooth) {}
};

/*
Acceleration boundary conditions

Note: Acceleration boundary condition can be specified on nodes and particles.
*/
struct AccelerationBC
{
	size_t index;
	DegreeOfFreedom dof;
	virtual TimeCurveType getType(void) = 0;
	virtual double a(double t) noexcept = 0;
	virtual void init(void) noexcept {}
};

struct AccelerationBC_Fixed : public AccelerationBC
{
	TimeCurveZero_inline curve;
	TimeCurveType getType(void) noexcept { return curve.getType(); }
	double a(double t) noexcept { return curve.x(t); }
};

struct AccelerationBC_Constant : public AccelerationBC
{
	TimeCurveConstant_inline curve;
	TimeCurveType getType(void) noexcept { return curve.getType(); }
	double a(double t) noexcept { return curve.x(t); }
};

struct AccelerationBC_Linear : public AccelerationBC
{
	TimeCurveLinear_inline curve;
	TimeCurveType getType(void) noexcept { return curve.getType(); }
	double a(double t) noexcept { return curve.x(t); }
	void init(void) noexcept { return curve.init(); }
};

struct AccelerationBC_CubicSmooth : public AccelerationBC
{
	TimeCurveCubicSmooth_inline curve;
	TimeCurveType getType(void) noexcept { return curve.getType(); }
	double a(double t) noexcept { return curve.x(t); }
	void init(void) noexcept { return curve.init(); }
};

class AccelerationBCManager : public FlexibleSizeMemory
{
public:
	AccelerationBCManager() {}
	~AccelerationBCManager() {}

	int add_bc(AccelerationBCParam *param);
	inline AccelerationBC *get_first() noexcept
	{
		return static_cast<AccelerationBC *>(FlexibleSizeMemory::get_first());
	}
	inline AccelerationBC *get_next(AccelerationBC *prev_model) noexcept
	{
		return static_cast<AccelerationBC *>(FlexibleSizeMemory::get_next(prev_model));
	}
};

// ------------------ Velocity Boundary Condition ------------------
// Parameters for initializing velocity boundary conditions
struct VelocityBCParam
{
protected:
	TimeCurveType type;
public:
	size_t index;
	DegreeOfFreedom dof;
	VelocityBCParam(TimeCurveType tp = TimeCurveType::Invalid) : type(tp) {}
	inline TimeCurveType getType() noexcept { return type; }
};

// Zero
struct VelocityBCParam_Fixed : public VelocityBCParam
{
	VelocityBCParam_Fixed() : VelocityBCParam(TimeCurveType::Zero) {}
};

// Constant
struct VelocityBCParam_Constant : public VelocityBCParam
{
	double v;
	VelocityBCParam_Constant() : VelocityBCParam(TimeCurveType::Constant) {}
};

// Linear
struct VelocityBCParam_Linear : public VelocityBCParam
{
	double v_begin;
	double v_end;
	double t_length;
	VelocityBCParam_Linear() : VelocityBCParam(TimeCurveType::Linear) {}
};

// Cubic Spline
struct VelocityBCParam_CubicSmooth : public VelocityBCParam
{
	double v_begin;
	double v_end;
	double t_length;
	VelocityBCParam_CubicSmooth() : VelocityBCParam(TimeCurveType::CubicSmooth) {}
};

/*
Velocity Boundary Conditions

Note: Velocity boundary condition can be specified on nodes and particles.
*/
struct VelocityBC
{
	size_t index;
	DegreeOfFreedom dof;
	virtual TimeCurveType getType(void) = 0;
	virtual double v(double t) noexcept = 0;
	virtual double a(double t) noexcept = 0;
	virtual void init(void) noexcept {}
};

struct VelocityBC_Fixed : public VelocityBC
{
	TimeCurveZero_inline curve;
	TimeCurveType getType(void) noexcept { return curve.getType(); }
	double v(double t) noexcept { return curve.x(t); }
	double a(double t) noexcept { return curve.dx_dt(t); }
};

struct VelocityBC_Constant : public VelocityBC
{
	TimeCurveConstant_inline curve;
	TimeCurveType getType(void) noexcept { return curve.getType(); }
	double v(double t) noexcept { return curve.x(t); }
	double a(double t) noexcept { return curve.dx_dt(t); }
};

struct VelocityBC_Linear : public VelocityBC
{
	TimeCurveLinear_inline curve;
	TimeCurveType getType(void) noexcept { return curve.getType(); }
	double v(double t) noexcept { return curve.x(t); }
	double a(double t) noexcept { return curve.dx_dt(t); }
	void init(void) noexcept { return curve.init(); }
};

struct VelocityBC_CubicSmooth : public VelocityBC
{
	TimeCurveCubicSmooth_inline curve;
	TimeCurveType getType(void) noexcept { return curve.getType(); }
	double v(double t) noexcept { return curve.x(t); }
	double a(double t) noexcept { return curve.dx_dt(t); }
	void init(void) noexcept { return curve.init(); }
};


class VelocityBCManager : public FlexibleSizeMemory
{
public:
	VelocityBCManager() {}
	~VelocityBCManager() {}

	int add_bc(VelocityBCParam *param);
	inline VelocityBC *get_first() noexcept
	{
		return static_cast<VelocityBC *>(FlexibleSizeMemory::get_first());
	}
	inline VelocityBC *get_next(VelocityBC *prev_model) noexcept
	{
		return static_cast<VelocityBC *>(FlexibleSizeMemory::get_next(prev_model));
	}
};


/* --------------------- Mass Force Boundary Conditions --------------------- */
// body force per unit mass on particles
// Parameters for initializing mass force BC
struct MassForceBCParam
{
protected:
	const TimeCurveType type;
public:
	size_t index;
	DegreeOfFreedom dof;
	MassForceBCParam(TimeCurveType tp = TimeCurveType::Invalid) : type(tp) {}
	inline TimeCurveType getType() noexcept { return type; }
};

struct MassForceBCParam_Constant : public MassForceBCParam
{
	double massForce;
	MassForceBCParam_Constant() :
		MassForceBCParam(TimeCurveType::Constant) {}
};

struct MassForceBCParam_Linear : public MassForceBCParam
{
	double massForce_begin;
	double massForce_end;
	double t_length;
	MassForceBCParam_Linear() :
		MassForceBCParam(TimeCurveType::Linear) {}
};

struct MassForceBCParam_CubicSmooth : public MassForceBCParam
{
	double massForce_begin;
	double massForce_end;
	double t_length;
	MassForceBCParam_CubicSmooth() :
		MassForceBCParam(TimeCurveType::CubicSmooth) {}
};

// Mass force
// body force per unit mass on particles
struct MassForceBC
{
	size_t index;
	DegreeOfFreedom dof;
	virtual TimeCurveType getType(void) = 0;
	virtual double massForce(double t) noexcept = 0;
	virtual void init(void) noexcept {}
};

struct MassForceBC_Constant : public MassForceBC
{
	TimeCurveConstant_inline curve;
	TimeCurveType getType(void) noexcept { return curve.getType(); }
	double massForce(double t) noexcept { return curve.x(t); }
};

struct MassForceBC_Linear : public MassForceBC
{
	TimeCurveLinear_inline curve;
	TimeCurveType getType(void) noexcept { return curve.getType(); }
	double massForce(double t) noexcept { return curve.x(t); }
	void init(void) noexcept { return curve.init(); }
};

struct MassForceBC_CubicSmooth : public MassForceBC
{
	TimeCurveCubicSmooth_inline curve;
	TimeCurveType getType(void) noexcept { return curve.getType(); }
	double massForce(double t) noexcept { return curve.x(t); }
	void init(void) noexcept { return curve.init(); }
};

class MassForceBCManager : public FlexibleSizeMemory
{
public:
	MassForceBCManager() {}
	~MassForceBCManager() {}

	int add_bc(MassForceBCParam *param);
	inline MassForceBC *get_first() noexcept
	{
		return static_cast<MassForceBC *>(FlexibleSizeMemory::get_first());
	}
	inline MassForceBC *get_next(MassForceBC *prev_model) noexcept
	{
		return static_cast<MassForceBC *>(FlexibleSizeMemory::get_next(prev_model));
	}
};

/* --------------------- Surface Force Boundary Conditions --------------------- */
// Note that the surface force here is actually surface stress times area
// Parameters for initializing mass force BC
struct SurfaceForceBCParam
{
protected:
	const TimeCurveType type;
public:
	size_t index;
	DegreeOfFreedom dof;
	SurfaceForceBCParam(TimeCurveType tp = TimeCurveType::Invalid) : type(tp) {}
	inline TimeCurveType getType() noexcept { return type; }
};

struct SurfaceForceBCParam_Constant : public SurfaceForceBCParam
{
	double surfaceForce;
	SurfaceForceBCParam_Constant() :
		SurfaceForceBCParam(TimeCurveType::Constant) {}
};

struct SurfaceForceBCParam_Linear : public SurfaceForceBCParam
{
	double surfaceForce_begin;
	double surfaceForce_end;
	double t_length;
	SurfaceForceBCParam_Linear() :
		SurfaceForceBCParam(TimeCurveType::Linear) {}
};

struct SurfaceForceBCParam_CubicSmooth : public SurfaceForceBCParam
{
	double surfaceForce_begin;
	double surfaceForce_end;
	double t_length;
	SurfaceForceBCParam_CubicSmooth() :
		SurfaceForceBCParam(TimeCurveType::CubicSmooth) {}
};

// surface force
struct SurfaceForceBC
{
	size_t index;
	DegreeOfFreedom dof;
	virtual TimeCurveType getType(void) = 0;
	virtual double surfaceForce(double t) noexcept = 0;
	virtual void init(void) noexcept {}
};

struct SurfaceForceBC_Constant : public SurfaceForceBC
{
	TimeCurveConstant_inline curve;
	TimeCurveType getType(void) noexcept { return curve.getType(); }
	double surfaceForce(double t) noexcept { return curve.x(t); }
};

struct SurfaceForceBC_Linear : public SurfaceForceBC
{
	TimeCurveLinear_inline curve;
	TimeCurveType getType(void) noexcept { return curve.getType(); }
	double surfaceForce(double t) noexcept { return curve.x(t); }
	void init(void) noexcept { return curve.init(); }
};

struct SurfaceForceBC_CubicSmooth : public SurfaceForceBC
{
	TimeCurveCubicSmooth_inline curve;
	TimeCurveType getType(void) noexcept { return curve.getType(); }
	double surfaceForce(double t) noexcept { return curve.x(t); }
	void init(void) noexcept { return curve.init(); }
};

class SurfaceForceBCManager : public FlexibleSizeMemory
{
public:
	SurfaceForceBCManager() {}
	~SurfaceForceBCManager() {}

	int add_bc(SurfaceForceBCParam *param);
	inline SurfaceForceBC *get_first() noexcept
	{
		return static_cast<SurfaceForceBC *>(FlexibleSizeMemory::get_first());
	}
	inline SurfaceForceBC *get_next(SurfaceForceBC *prev_model) noexcept
	{
		return static_cast<SurfaceForceBC *>(FlexibleSizeMemory::get_next(prev_model));
	}
};

#endif