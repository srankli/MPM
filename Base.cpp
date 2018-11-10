#ifdef _DEBUG
#include <iostream>
#endif

#include "Mesh.h"
#include "Particle.h"
#include "Solver.h"
#include "OutputRequest.h"

/* --------------------------- MeshType ---------------------------- */
const char *MeshTypeName::type_name[] = {
	"InvalidMesh", // 0
	"Mesh_1D2",    // 1
	"Mesh_R2D4"   // 2
};

const char *MeshTypeName::InvalidMesh = MeshTypeName::type_name[0];
const char *MeshTypeName::Mesh_1D2 = MeshTypeName::type_name[1];
const char *MeshTypeName::Mesh_R2D4 = MeshTypeName::type_name[2];

const size_t MeshTypeName::type_num = 
	sizeof(MeshTypeName::type_name) / sizeof(MeshTypeName::type_name[0]);

const char *MeshTypeName::getName(size_t id)
{
	if (id < type_num) return type_name[id];
	else return type_name[0]; // Invalid Simulation
}


/* ------------------------ SimulationType ------------------------- */
const char *SimulationTypeName::type_name[] = {
	"InvalidSimulation", // 0
	"Mechanics_1D", // 1
	"Hydromechanics", // 2
	"Mechanics_2D", // 3 
	"Hydromechanics" // 4
};

const char *SimulationTypeName::InvalidSimulation = SimulationTypeName::type_name[0];
const char *SimulationTypeName::Mechanics_1D = SimulationTypeName::type_name[1];
const char *SimulationTypeName::Hydromechanics_1D = SimulationTypeName::type_name[2];
const char *SimulationTypeName::Mechanics_2D = SimulationTypeName::type_name[3];
const char *SimulationTypeName::Hydromechanics_2D = SimulationTypeName::type_name[4];

const size_t SimulationTypeName::type_num
	= sizeof(SimulationTypeName::type_name) / sizeof(SimulationTypeName::type_name[0]);

const char *SimulationTypeName::getName(size_t id)
{
	if (id < type_num) return type_name[id];
	else return type_name[0]; // Invalid Simulation
}

/* ------------------------ ObjectByParticle ------------------------ */
size_t ObjectByParticle::curIndex = 0;


/* --------------------------- Solver Error --------------------------- */
const std::string SolverError::errorCodeToExplanation[] = {
	std::string("No default error explanation.")
};

const size_t SolverError::errorExplanationNum
	= sizeof(SolverError::errorCodeToExplanation) /
	  sizeof(SolverError::errorCodeToExplanation[0]);

SolverError::SolverError(size_t cd, const char *msg) :
	code(cd), message(msg) {}

const std::string &SolverError::getExplaination(void)
{
	if (code >= errorExplanationNum) code = 0;
		
	return errorCodeToExplanation[code];
}

/* -------------------------------- Solver --------------------------------- */
int Solver::solve(double time_increment)
{
	double t_tol;

	init();

	t_cal = 0.0;
	output->output(iteration_index, t_cal, true);
	//std::cout << output->getStepTime() << std::endl;

	// The first iteration
	t_increment = time_increment;
	// At the first step, time increment of acceleration integration
	// should be half of that of velocity integration 
	t_increment_a = t_increment / 2.0;
	t_cal += t_increment;
	t_tol = 0.01 * t_increment;
	if (t_cal > (t_step - t_tol))
	{
		if ((t_cal - t_step) > 2.0 * t_tol)
		{
			t_increment -= (t_cal - t_step);
			t_increment_a = t_increment / 2.0;
			t_cal = t_step;
		}
		goto last_iteration;
	}
	++iteration_index;
	iteration();
	output->output(iteration_index, t_cal);

	// The other iteration
	t_increment_a = t_increment;
	while (true)
	{
		// t_increment may change with time and thus recalculated every iteration
		t_cal += t_increment;
		t_tol = 0.01 * t_increment;
		if (t_cal > (t_step - t_tol))
		{
			if ((t_cal - t_step) > 2.0 * t_tol)
			{
				t_increment -= (t_cal - t_step);
				t_increment_a = t_increment;
				t_cal = t_step;
			}
			goto last_iteration;
		}
		++iteration_index;
		iteration();
		// output result
		output->output(iteration_index, t_cal);
	}

	// The last iteration
last_iteration:
	++iteration_index;
	iteration();
	output->output(iteration_index, t_cal, true);
	output->complete();

	return 0;
}
