#ifndef __PARTICLE_1D_MECHANICS_H__
#define __PARTICLE_1D_MECHANICS_H__

#include "Particle.h"
#include "MemoryManager.h"

struct NodeVar_1D_Mechanics;
class ObjectByParticle_1D_Mechanics;
class Solver_1D_Mechanics_1D2_Explicit;
class Solver_1D_Mechanics_1D2_Explicit_FixedMem;
struct MassForceOnParticleParam_1D_Mechanics;
struct SurfaceForceOnParticleParam_1D_Mechanics;

// Used to initialize ParticleParam_1D_Mechanics
struct ParticleParam_1D_Mechanics
{
	// terms with *** MUST BE initialized
	double x; // location ***

	double mass; // ***
	double density; // ***
	double momentum1; // default = 0.0

	double stress11; // stress, take tension as positive, default = 0.0

	double strain11; // total strain of soil skeleton, default = 0.0
	double estrain11; // elastic strain of soil skeleton, default = 0.0
	double pstrain11; // plastic strain of soil skeleton, default = 0.0

	ParticleParam_1D_Mechanics() :
		momentum1(0.0),
		stress11(0.0), strain11(0.0),
		estrain11(0.0), pstrain11(0.0) {}
};


/* --------------------------------------------------------- */
struct Particle_1D_Mechanics : public Particle
{
	// Geometric properties
	double x;
	// Physical properties
	double mass;
	double density;
	double momentum1;
	double stress11;
	double strain11;
	double estrain11;
	double pstrain11;
};


// variables for calculation
struct ParticleVar_1D_Mechanics : public ParticleVar
{
	friend ObjectByParticle_1D_Mechanics;
public:
	// point to two node variables of nodes
	NodeVar_1D_Mechanics *nodeVar1, *nodeVar2;
	// natural coordinates
	double xi;
	// velocity of solid phase
	double v1;
	// increment of total strain
	double dstrain11;

	union // value of shape function
	{
		double N[2];
		struct
		{
			double N1;
			double N2;
		};
	};
	union
	{
		double dN_dx[2][1];
		struct
		{
			double dN1_dx;
			double dN2_dx;
		};
	};
};


class ObjectByParticle_1D_Mechanics : public ObjectByParticle
{
	friend Solver_1D_Mechanics_1D2_Explicit;
	friend Solver_1D_Mechanics_1D2_Explicit_FixedMem;
protected:
	Particle_1D_Mechanics *particles;
public:
	ObjectByParticle_1D_Mechanics() : 
		ObjectByParticle(sizeof(Particle_1D_Mechanics), SimulationType::Mechanics_1D),
		particles(nullptr) {}
	~ObjectByParticle_1D_Mechanics() {}

	int addParticle(ParticleParam_1D_Mechanics *pcl_param,
		            ConstitutiveModelParam *pcl_cm);
	void finish_init();

	inline Particle_1D_Mechanics *getParticleById(size_t id) noexcept
	{
		return id && id <= particleNum ? particles + id - 1 : nullptr;
	}
	inline Particle_1D_Mechanics *getParticles(void) noexcept
	{
		return particles;
	}
};

#endif