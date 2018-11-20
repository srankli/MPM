#ifndef __SOLVER_TEST_H__
#define __SOLVER_TEST_H__

#include "Solver.h"
#include "Mesh_1D2.h"
#include "Particle_1D_Mechanics.h"
#include "OutputRequest.h"

class Solver_Test : public Solver
{
	Mesh_1D2 &mesh;
	std::vector<ObjectByParticle_1D_Mechanics> &pcl_objects;

public:
	Solver_Test(
		double time_step,
		Mesh_1D2 &mh,
		std::vector<ObjectByParticle_1D_Mechanics> &pcl_objs,
		OutputRequest &out,
		const char *na = "") :
		Solver(time_step, out, na),
		mesh(mh), pcl_objects(pcl_objs) {}
	// The copy function is used to get data from the previous step (solver)
	Solver_Test(
		double time_step,
		Solver_Test &prev_solver,
		const char *na = "") :
		Solver(time_step, prev_solver, na),
		mesh(prev_solver.mesh), pcl_objects(prev_solver.pcl_objects) {}
	~Solver_Test() {}

	int init(void)
	{
		std::cout << "step_id: " << index << std::endl;
		return 0;
	}

	int iteration(void)
	{
		std::cout << "  iter_id: " << iteration_index << std::endl;
		return 0;
	}
};

#endif