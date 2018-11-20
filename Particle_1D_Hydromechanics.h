#ifndef __PARTICLE_1D_HYDROMECHANICS_H__
#define __PARTICLE_1D_HYDROMECHANICS_H__

#include "Particle.h"
#include "MemoryManager.h"

struct NodeVar_1D_Hydromechanics;
class ObjectByParticle_1D_Hydromechanics;
class Solver_1D_Hydromechanics_1D2_Explicit;
class Solver_1D_Hydromechanics_1D2_Explicit_FixedMem;

// Used to initialize ParticleParam_1D_Hydromechanics
struct ParticleParam_1D_Hydromechanics
{
public:
	// terms with *** MUST BE initialized
	double x; // location ***
	
	// Momentum of solid - fluid mixture
	double momentum1_m; // ***

	// ----------- Solid phase ------------
	double mass_s; // ***
	double density_s; // ***

	double estress11; // effective stress, take tension as positive, default = 0.0
	double strain11; // total strain of soil skeleton, default = 0.0
	double estrain11; // elastic strain of soil skeleton, default = 0.0
	double pstrain11; // plastic strain of soil skeleton, default = 0.0

	// --------- fluid phase ---------
	double mass_f; // ***
	double density_f; // ***

	double p; // pore pressure, take compressive as positive, default = 0.0
	double Kf; // bulk modulus of fluid***

	/*
	Darcy's law:
	nv = k * (-grad pressure + fluid_density * bodyforce)
	note that k = kw / unit_weight of water
	*/
	double k; // permeability *** (unit: N*s/m4)

public:
	ParticleParam_1D_Hydromechanics() :
		momentum1_m(0.0),
		estress11(0.0), strain11(0.0),
		estrain11(0.0), pstrain11(0.0),
		p(0.0){}
};


// Used to initialize SurfaceForceOnParticle_1D_Hydromechanics
struct Particle_1D_Hydromechanics : public Particle
{
public:
	double x; // location
	
	// Physical properties
	double n; // porosity
	double stress11; // total stress = effctive stress - pore pressure

	// ----------- Solid phase ------------
	double mass_s;
	double density_s;

	double estress11; // effective stress, take tension as positive
	double strain11; // total strain of soil skeleton
	double estrain11; // elastic strain of soil skeleton
	double pstrain11; // plastic strain of soil skeleton
	
	// --------- fluid phase ---------
	double mass_f;
	double density_f;
	//double momentum1_f; // is not momentum strictly speaking

	// pore pressure, take compression as positive
	double p;
	// bulk modulus of fluid
	double Kf;

	double momentum1_m; // linear momentum of solid - fluid mixture 
	// permeability: refer to class ParticleParam_1D_Hydromechanics for its defination.
	double k;

public:
	Particle_1D_Hydromechanics() :
		x(0.0), n(0.0), stress11(0.0),
		mass_s(0.0), density_s(0.0), momentum1_m(0.0),
		estress11(0.0), strain11(0.0), estrain11(0.0), pstrain11(0.0),
		mass_f(0.0), density_f(0.0), //momentum1_f(0.0),
		p(0.0), Kf(0.0), k(0.0) {}
};


// variables for calculation
struct ParticleVar_1D_Hydromechanics : public ParticleVar
{
	friend ObjectByParticle_1D_Hydromechanics;
public:
	// point to two node variables of nodes
	NodeVar_1D_Hydromechanics *nodeVar1, *nodeVar2;

	double xi;// natural coordinate in element
	
	double volume; // volume represented by particle
	double avgdensity_s; // (1 - n) * density_s
	double avgdensity_f; // n * density_f

	double a1_s; // acceleration of soild phase
	double v1_s; // velocity of solid phase
	double dstrain11; // increment of total strain
	
	// relative apparent velocity of fluid phase relative to solid phase
	// w = n * (v1_f - v1_s)
	double w;
	double dw_dx; // "volumetric strain" of fluid

	union // value of shape function
	{
		double N[2];
		struct
		{
			double N1;
			double N2;
		};
	};
	union // one order derivative of shape function
	{
		double dN_dx[2][1];
		struct
		{
			double dN1_dx;
			double dN2_dx;
		};
	};
};


class ObjectByParticle_1D_Hydromechanics : public ObjectByParticle
{
	friend Solver_1D_Hydromechanics_1D2_Explicit;
	friend Solver_1D_Hydromechanics_1D2_Explicit_FixedMem;
protected:
	// Particles
	Particle_1D_Hydromechanics *particles;

public:
	ObjectByParticle_1D_Hydromechanics() :
		ObjectByParticle(sizeof(Particle_1D_Hydromechanics), SimulationType::Hydromechanics_1D),
		particles(nullptr) {}
	~ObjectByParticle_1D_Hydromechanics() {}
	
	int addParticle(ParticleParam_1D_Hydromechanics *pcl_param,
					ConstitutiveModelParam *pcl_cm);
	void finish_init(void);

	// These functions can only be called after finish_init()
	inline Particle_1D_Hydromechanics *getParticleById(size_t id) noexcept
	{
		return id <= particleNum && id ? particles + id - 1 : nullptr;
	}
	inline Particle_1D_Hydromechanics *getParticles(void) noexcept
	{
		return particles;
	}
};

#endif