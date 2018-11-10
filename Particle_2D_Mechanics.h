#ifndef __PARTICLE_2D_MECHANICS_H__
#define __PARTICLE_2D_MECHANICS_H__

#include "Particle.h"
#include "MemoryManager.h"

struct NodeVar_2D_Mechanics;
class ObjectByParticle_2D_Mechanics;
class Solver_2DPlaneStrain_Mechanics_R2D4_Explicit;
class Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem;
struct MassForceOnParticleParam_2D_Mechanics;
struct SurfaceForceOnParticleParam_2D_Mechanics;

// Used to initialize ParticleParam_1D_Mechanics
struct ParticleParam_2D_Mechanics
{
	// terms with *** MUST BE initialized
	// location ***
	union
	{
		double coords[2];
		struct { double x, y; };
	};
	double mass; // ***
	double density; // ***
	// momentum, default = 0.0
	union
	{
		double momentum[2];
		struct
		{
			double momentum1;
			double momentum2;
		};
	};
	// stress, take tension as positive, default = 0.0
	union
	{
		double stress[6];
		struct
		{
			double stress11;
			double stress22;
			double stress33;
			double stress12;
			double stress23;
			double stress31;
		};
	};
	// total strain of soil skeleton, default = 0.0
	union
	{
		double strain[3];
		struct
		{
			double strain11;
			double strain22;
			double strain12;
		};
	};
	// elastic strain of soil skeleton, default = 0.0
	union
	{
		double estrain[3];
		struct
		{
			double estrain11;
			double estrain22;
			double estrain12;
		};
	};
	// plastic strain of soil skeleton, default = 0.0
	union
	{
		double pstrain[3];
		struct
		{
			double pstrain11;
			double pstrain22;
			double pstrain12;
		};
	};
	ParticleParam_2D_Mechanics() :
		momentum1(0.0), momentum2(0.0),
		stress11(0.0), stress22(0.0), stress33(0.0),
		stress12(0.0), stress23(0.0), stress31(0.0),
		strain11(0.0), strain22(0.0), strain12(0.0),
		estrain11(0.0), estrain22(0.0), estrain12(0.0),
		pstrain11(0.0), pstrain22(0.0), pstrain12(0.0){}
};


// Dimension_PhysicalPhenomenon_BackgroundMeshType
struct Particle_2D_Mechanics : public Particle
{
public:
	// Geometric properties
	union
	{
		double coords[2];
		struct { double x, y; };
	};
	// Physical properties
	double mass;
	double density;
	// momentum
	union
	{
		double momentum[2];
		struct
		{
			double momentum1;
			double momentum2;
		};
	};
	// stress
	union
	{
		double stress[6];
		struct
		{
			double stress11;
			double stress22;
			double stress33;
			double stress12;
			double stress23;
			double stress31;
		};
	};
	// total strain
	union
	{
		double strain[3];
		struct
		{
			double strain11;
			double strain22;
			double strain12;
		};
	};
	// elastic strain
	union
	{
		double estrain[3];
		struct
		{
			double estrain11;
			double estrain22;
			double estrain12;
		};
	};
	// plastic strain
	union
	{
		double pstrain[3];
		struct
		{
			double pstrain11;
			double pstrain22;
			double pstrain12;
		};
	};
};


// variables for calculation
struct ParticleVar_2D_Mechanics : public ParticleVar
{
	friend ObjectByParticle_2D_Mechanics;
public:
	// point to two node variables of nodes
	union
	{
		NodeVar_2D_Mechanics *nodeVar[4];
		struct
		{
			NodeVar_2D_Mechanics *nodeVar1, *nodeVar2;
			NodeVar_2D_Mechanics *nodeVar3, *nodeVar4;
		};
	};
	// natural coordinate in element
	union
	{
		double naturalCoords[2];
		struct { double xi, eta; };
	};
	// velocity
	union
	{
		double velocity[2];
		struct { double v1, v2; };
	};
	// increment of total strain
	union
	{
		double dstrain[3];
		struct
		{
			double dstrain11;
			double dstrain22;
			double dstrain12;
		};
	};
	// W for Jaumann
	double dW12;
	union // value of shape function
	{
		double N[4];
		struct
		{
			double N1;
			double N2;
			double N3;
			double N4;
		};
	};
	union
	{
		double dN_dx[4][2];
		struct
		{
			double dN1_dx, dN1_dy;
			double dN2_dx, dN2_dy;
			double dN3_dx, dN3_dy;
			double dN4_dx, dN4_dy;
		};
	};
};


class ObjectByParticle_2D_Mechanics : public ObjectByParticle
{
	friend Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem;
protected:
	Particle_2D_Mechanics *particles;
public:
	ObjectByParticle_2D_Mechanics() :
		ObjectByParticle(sizeof(Particle_2D_Mechanics), SimulationType::Mechanics_2D),
		particles(nullptr) {}
	~ObjectByParticle_2D_Mechanics() {}

	int addParticle(ParticleParam_2D_Mechanics *pcl_param,
		ConstitutiveModelParam *pcl_cm);
	void finish_init();

	/*
	Create particles from triangular element
		1. x1, y1, x2, y2, x3, y3 are coordinates of nodes
		2. pcl_param->mass is unused
		3. characteristic_size = sqrt(area of elem).
		4. min_characteristic_size < 0 means no requirement
		   on maximum chararcteristic size of element.
	*/
	void addParticleFromTriElement(
		double x1, double y1,
		double x2, double y2,
		double x3, double y3,
		ParticleParam_2D_Mechanics *pcl_param,
		ConstitutiveModelParam *pcl_cm,
		double max_characteristic_size = -1.0);
	
	inline Particle_2D_Mechanics *getParticleById(size_t id) noexcept
	{
		return id && id <= particleNum ? particles + id - 1 : nullptr;
	}
	inline Particle_2D_Mechanics *getParticles(void) noexcept
	{
		return particles;
	}
};

#endif