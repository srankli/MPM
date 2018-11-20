#ifndef __PARTICLE_H__
#define __PARTICLE_H__

#include <string>
#include "CMManager.h"
#include "BoundaryCondition.h"

struct Particle;
class ObjectByParticle;

struct Element;
class ConstitutiveModel;
class Solver;


/*
Type of simulation currently supported by this program.
Seems that the only difference between plane
strain and plane stress is constitutive model.
Therefore they are not discerned here.
*/
enum class SimulationType : unsigned long long
{
	InvalidSimulation = 0,
	Mechanics_1D = 1,
	Hydromechanics_1D = 2,
	Mechanics_2D = 3,
	Hydromechanics_2D = 4
};

class SimulationTypeName
{
protected:
	static const char *type_name[];
	static const size_t type_num;
public:
	static const char *InvalidSimulation; // 0
	static const char *Mechanics_1D;      // 1
	static const char *Hydromechanics_1D; // 2
	static const char *Mechanics_2D;      // 3
	static const char *Hydromechanics_2D; // 4

	static const char *getName(size_t id);
};

// ------------------------ material points -----------------------
struct ParticleVar
{
	friend ObjectByParticle;
public:
	Particle *particle;
	/*
	Points to element in which the particle lies
	Should be null if particle does not lie in any elements
	*/
	Element *element;

protected: 
	// Used by ParticleVar List
	ParticleVar *next;
};

struct Particle
{
	size_t index;
	// The object that this particle belong
	ObjectByParticle *object;
	// Constitutive model
	ConstitutiveModel *cm;
	
	// Whether the particle is in mesh
	bool isInMesh;
	// Variables for calculation
	ParticleVar *particleVar;
};

// objects that do not merge with other during deformation
class ObjectByParticle
{
	friend Solver;
protected:
	// Object info
	static size_t curIndex;
	size_t index;
	std::string name;
	SimulationType stype;

protected:
	// Particles
	size_t curParticleIndex;
	size_t particleNum;
	FixedSizeMemeory particles_mem;
	// Constitutive models
	CMManager constitutiveModels_mem;
	// Boundary conditions on particles
	size_t massForceBCNum;
	MassForceBCManager massForceBCs_mem;
	size_t surfaceForceBCNum;
	SurfaceForceBCManager surfaceForceBCs_mem;
	
public:
	ObjectByParticle(size_t pcl_size,
		SimulationType stp = SimulationType::InvalidSimulation) :
		index(++curIndex), name(""), stype(stp),
		curParticleIndex(0), particleNum(0), particles_mem(pcl_size),
		constitutiveModels_mem(),
		massForceBCNum(0), surfaceForceBCNum(0)
	{
		reset_list();
	}
	~ObjectByParticle()
	{
		particles_mem.clear();
		constitutiveModels_mem.clear();
		massForceBCs_mem.clear();
		surfaceForceBCs_mem.clear();
	}

	inline size_t getIndex() noexcept { return index; }
	inline void setName(const char *na) noexcept { if (na) name = na; }
	inline const std::string &getName() noexcept { return name; }
	inline SimulationType getSimulationType() { return stype; }
	inline size_t getParticleNum() noexcept { return particleNum; }
	inline Particle *getFirstParticle(void) noexcept
	{
		return static_cast<Particle *>(particles_mem.get_first());
	}
	inline Particle *getNextParticle(Particle *iter) noexcept
	{
		return static_cast<Particle *>(particles_mem.get_next(iter));
	}

	// Maybe needed to be redefined
	virtual bool validateParticleId(size_t id) noexcept
	{
		return id && id <= particleNum ? true : false;
	}
	virtual int addMassForceBC(MassForceBCParam *mfp_param)
	{
		int res = -1;
		if (validateParticleId(mfp_param->index))
		{
			res = massForceBCs_mem.add_bc(mfp_param);
			if (!res) ++massForceBCNum;
		}
		return res;
	}
	virtual int addSurfaceForceBC(SurfaceForceBCParam *sfp_param)
	{
		int res = -1;
		if (validateParticleId(sfp_param->index))
		{
			res = surfaceForceBCs_mem.add_bc(sfp_param);
			if (!res) ++surfaceForceBCNum;
		}
		return res;
	}
	
	void finish_init()
	{
		Particle *pcl_iter;
		ConstitutiveModel *cm_iter;

		particles_mem.compress();
		constitutiveModels_mem.compress();
		massForceBCs_mem.compress();
		surfaceForceBCs_mem.compress();
		
		// Associate particles with constitutive model
		for (pcl_iter = static_cast<Particle *>(particles_mem.get_first()),
			 cm_iter = constitutiveModels_mem.get_first();
			 pcl_iter && cm_iter;
			 pcl_iter = static_cast<Particle *>(particles_mem.get_next(pcl_iter)),
			 cm_iter = constitutiveModels_mem.get_next(cm_iter))
			pcl_iter->cm = cm_iter;
	}

// ---------------- List of ParticleVars ---------------- 
protected: 	
	ParticleVar list_head;
	ParticleVar *list_tail;
	size_t list_num;
	ParticleVar *list_cur, *list_cur_prev;
public:
	/* List of ParticleVar:
	Note that this list does not
	allocate memory for ParticleVar.
	*/
	inline size_t get_num_list(void) noexcept { return list_num; }
	void add_list(ParticleVar *item)
	{
		item->next = nullptr;
		list_tail->next = item;
		list_tail = item;
		++list_num;
	}
	void reset_list(void)
	{
		list_tail = &list_head;
		list_head.next = nullptr;
		list_num = 0;
		list_cur = nullptr;
		list_cur_prev = nullptr;
	}
	inline void start_list(void) noexcept
	{
		list_cur_prev = &list_head;
		list_cur = list_head.next;
	}
	inline void next_list(void) noexcept
	{
		list_cur_prev = list_cur;
		list_cur = list_cur ? list_cur->next : nullptr;
	}
	// delete the item pointed by cur and move to the next item
	ParticleVar *del_cur_list(void)
	{
		ParticleVar *tmp;
		if (list_cur)
		{
			tmp = list_cur;
			list_cur_prev->next = list_cur->next;
			list_cur = list_cur->next;
			--list_num;
			return tmp;
		}
		return nullptr;
	}
	// This version of iteration cannot called del_cur()
	inline void start2_list(void) noexcept { list_cur = list_head.next; }
	inline void next2_list(void) noexcept { list_cur = list_cur ? list_cur->next : nullptr; }
	inline ParticleVar *get_cur_list(void) noexcept { return list_cur; }
};

#endif