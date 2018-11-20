#ifdef _DEBUG
#include <iostream>

#include "Solver_1D_Hydromechanics_1D2_Explicit_FixedMem.h"
#include "TmpDataToHdf5.h"

// 1 point
void test_Solver_1D_Hydromechanics_1D2_Explicit_FixedMem_1(void)
{
	double time_step = 1.0;

	Mesh_1D2 mesh;
	double coords[] = { 0.0, 1.0 };
	mesh.initMesh(coords, sizeof(coords) / sizeof(coords[0]));
	// Velocity boundary conditions
	VelocityBCParam_Fixed vbc;
	vbc.index = 1;
	vbc.dof = DegreeOfFreedom::x;
	mesh.addVelocityBC(&vbc);
	vbc.dof = DegreeOfFreedom::x_f;
	mesh.addVelocityBC(&vbc);
	VelocityBCParam_Constant vbc2;
	vbc2.v = 0.0;
	vbc2.index = 2;
	vbc2.dof = DegreeOfFreedom::x;
	mesh.addVelocityBC(&vbc2);
	vbc2.v = -0.1;
	vbc2.dof = DegreeOfFreedom::x_f;
	mesh.addVelocityBC(&vbc2);
	mesh.finish_init();

	std::vector<ObjectByParticle_1D_Hydromechanics> objects;
	ParticleParam_1D_Hydromechanics pcl_param;
	LinearElasticityModelParam le_param;

	objects.resize(1);

	objects[0].setName("testObject1");
	// particle 1
	pcl_param.x = 0.5;
	pcl_param.mass_s = 1.0;
	pcl_param.density_s = 2.0;
	pcl_param.mass_f = 1.0;
	pcl_param.density_f = 2.0;
	pcl_param.Kf = 100.0;
	pcl_param.k = 0.01; // ??
	le_param.E = 50.0;
	le_param.nu = 0.0;
	objects[0].addParticle(&pcl_param, &le_param);
	/*
	// Boundary conditions
	MassForceBCParam_Constant mfp_param;
	// particle 1
	mfp_param.dof = DegreeOfFreedom::x;
	mfp_param.massForce = 0.01;
	mfp_param.index = 1;
	objects[0].addMassForceBC(&mfp_param);
	// Surface forces
	SurfaceForceBCParam_Constant sfp_param1;
	sfp_param1.dof = DegreeOfFreedom::x;
	sfp_param1.surfaceForce = 1.0;
	sfp_param1.index = 5;
	SurfaceForceBCParam_Linear sfp_param2;
	sfp_param2.dof = DegreeOfFreedom::x;
	sfp_param2.surfaceForce_begin = 0.0;
	sfp_param2.surfaceForce_end = 1.0;
	sfp_param2.t_length = time_step / 3.0 * 2.0;
	sfp_param2.index = 5;
	objects[0].addSurfaceForceBC(&sfp_param1);
	*/
	// complete initialization
	objects[0].finish_init();

	OutputRequest output_request("test_result");
	Output *out1;
	std::vector<size_t> pcl_id;
	std::vector<unsigned long long> field_id;

	pcl_id.push_back(1);
	field_id.push_back(static_cast<unsigned long long>(OutputFieldType_1D_Hydromechanics::x));
	field_id.push_back(static_cast<unsigned long long>(OutputFieldType_1D_Hydromechanics::n));
	field_id.push_back(static_cast<unsigned long long>(OutputFieldType_1D_Hydromechanics::p));
	field_id.push_back(static_cast<unsigned long long>(OutputFieldType_1D_Hydromechanics::estress11));

	out1 = output_request.addOutput("test_output1", 5);
	out1->addObject(&objects[0], &field_id, &pcl_id);

	output_request.outputMeshGeometry(&mesh);

	Solver_1D_Hydromechanics_1D2_Explicit_FixedMem *solver;
	solver = new Solver_1D_Hydromechanics_1D2_Explicit_FixedMem(time_step, mesh, objects, output_request);
	solver->solve(0.1);
	delete solver;

	TmpDataToHdf5 trans;
	trans.init(output_request);
	trans.generate();
}


// 1D consolidation, one side impermeable
void test_Solver_1D_Hydromechanics_1D2_Explicit_FixedMem_2(void)
{
	size_t i;

	double time_step = 1.0;

	Mesh_1D2 mesh;
	double coords[] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };
	mesh.initMesh(coords, sizeof(coords) / sizeof(coords[0]));
	
	// Velocity boundary conditions
	VelocityBCParam_Fixed vbc;
	vbc.index = 6;
	vbc.dof = DegreeOfFreedom::x;
	mesh.addVelocityBC(&vbc);
	vbc.dof = DegreeOfFreedom::x_f;
	mesh.addVelocityBC(&vbc);
	VelocityBCParam_Constant vbc2;
	vbc2.v = 0.0;
	vbc2.index = 2;
	vbc2.dof = DegreeOfFreedom::x;
	//mesh.addVelocityBC(&vbc2);
	vbc2.v = -0.1;
	vbc2.dof = DegreeOfFreedom::x_f;
	//mesh.addVelocityBC(&vbc2);
	mesh.finish_init();

	std::vector<ObjectByParticle_1D_Hydromechanics> objects;
	ParticleParam_1D_Hydromechanics pcl_param;
	LinearElasticityModelParam le_param;
	objects.resize(1);
	objects[0].setName("testObject1");
	// particle 1
	pcl_param.mass_s = 1.0;
	pcl_param.density_s = 4.0;
	pcl_param.mass_f = 1.0;
	pcl_param.density_f = 4.0;
	//pcl_param.p = 0.0;
	pcl_param.Kf = 100.0;
	pcl_param.k = 0.01;
	le_param.E = 10.0;
	le_param.nu = 0.0;
	for (i = 0; i < 5; i++)
	{
		pcl_param.x = 0.5 + (double)i;
		objects[0].addParticle(&pcl_param, &le_param);
	}
	// Boundary conditions
	MassForceBCParam_Constant mfp_param;
	/*// particle 1
	mfp_param.dof = DegreeOfFreedom::x;
	mfp_param.massForce = 0.01;
	mfp_param.index = 1;
	objects[0].addMassForceBC(&mfp_param);*/
	// Surface forces
	SurfaceForceBCParam_Constant sfp_param1;
	sfp_param1.dof = DegreeOfFreedom::x;
	sfp_param1.surfaceForce = 0.01;
	sfp_param1.index = 1;
	objects[0].addSurfaceForceBC(&sfp_param1);
	SurfaceForceBCParam_Linear sfp_param2;
	sfp_param2.dof = DegreeOfFreedom::x;
	sfp_param2.surfaceForce_begin = 0.0;
	sfp_param2.surfaceForce_end = 1.0;
	sfp_param2.t_length = time_step / 3.0 * 2.0;
	sfp_param2.index = 1;
	//objects[0].addSurfaceForceBC(&sfp_param1);
	// complete initialization
	objects[0].finish_init();

	OutputRequest output_request("test_result");
	Output *out1;
	std::vector<size_t> pcl_id;
	std::vector<unsigned long long> field_id;

	for (i = 0; i < 5; i++) pcl_id.push_back(i + 1);
	field_id.push_back(static_cast<unsigned long long>(OutputFieldType_1D_Hydromechanics::x));
	field_id.push_back(static_cast<unsigned long long>(OutputFieldType_1D_Hydromechanics::n));
	field_id.push_back(static_cast<unsigned long long>(OutputFieldType_1D_Hydromechanics::p));
	field_id.push_back(static_cast<unsigned long long>(OutputFieldType_1D_Hydromechanics::estress11));

	out1 = output_request.addOutput("test_output1", 10);
	out1->addObject(&objects[0], &field_id, &pcl_id);

	output_request.outputMeshGeometry(&mesh);

	Solver_1D_Hydromechanics_1D2_Explicit_FixedMem *solver;
	solver = new Solver_1D_Hydromechanics_1D2_Explicit_FixedMem(time_step, mesh, objects, output_request);
	solver->solve(0.1);
	delete solver;

	TmpDataToHdf5 trans;
	trans.init(output_request);
	trans.generate();
}


#endif