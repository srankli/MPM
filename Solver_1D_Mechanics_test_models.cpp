#ifdef _DEBUG
#include <iostream>

#include "Solver_1D_Mechanics_1D2_Explicit_FixedMem.h"
#include "TmpDataToHdf5.h"

void test_Solver_1D_Mechanics_1D2_Explicit_FixedMem(void)
{
	double time_step = 10;

	Mesh_1D2 mesh;
	double coords[] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
	mesh.initMesh(coords, sizeof(coords) / sizeof(coords[0]));
	// velocity boundary conditions
	VelocityBCParam_Fixed vbc;
	vbc.index = 1;
	vbc.dof = DegreeOfFreedom::x;
	mesh.addVelocityBC(&vbc);
	mesh.finish_init();

	std::vector<ObjectByParticle_1D_Mechanics> objects;
	ParticleParam_1D_Mechanics pcl_param;
	LinearElasticityModelParam le_param;

	objects.resize(1);

	objects[0].setName("testObject1");
	// particle 1
	pcl_param.x = 0.5;
	pcl_param.mass = 1.0;
	pcl_param.density = 2.0;
	pcl_param.momentum1 = 0.0;
	pcl_param.stress11 = 0.0;
	pcl_param.strain11 = 0.0;
	pcl_param.estrain11 = 0.0;
	pcl_param.pstrain11 = 0.0;
	le_param.E = 100.0;
	le_param.nu = 0.0;
	objects[0].addParticle(&pcl_param, &le_param);
	// particle 2
	pcl_param.x = 1.5;
	objects[0].addParticle(&pcl_param, &le_param);
	// particle 3
	pcl_param.x = 2.5;
	objects[0].addParticle(&pcl_param, &le_param);
	// particle 4
	pcl_param.x = 3.5;
	objects[0].addParticle(&pcl_param, &le_param);
	// particle 5
	pcl_param.x = 4.5;
	objects[0].addParticle(&pcl_param, &le_param);
	/*
	// Boundary conditions
	MassForceBCParam_Constant mfp_param;
	// particle 1
	mfp_param.dof = DegreeOfFreedom::x;
	mfp_param.massForce = 0.01;
	mfp_param.index = 1;
	objects[0].addMassForceBC(&mfp_param);
	// particle 2
	mfp_param.index = 2;
	objects[0].addMassForceBC(&mfp_param);
	// particle 3
	mfp_param.index = 3;
	objects[0].addMassForceBC(&mfp_param);
	// particle 4
	mfp_param.index = 4;
	objects[0].addMassForceBC(&mfp_param);
	// particle 5
	mfp_param.index = 5;
	objects[0].addMassForceBC(&mfp_param);
	*/
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
	// complete initialization
	objects[0].finish_init();

	OutputRequest output_request;
	Output *out1;
	std::vector<size_t> pcl_id;
	std::vector<unsigned long long> field_id;
	
	output_request.setFileName("test_result");
	output_request.setStepTime(time_step);
	pcl_id.push_back(1);
	pcl_id.push_back(2);
	pcl_id.push_back(3);
	pcl_id.push_back(4);
	pcl_id.push_back(5);
	//pcl_id.push_back(6);
	//pcl_id.push_back(7);
	field_id.push_back(static_cast<unsigned int>(OutputFieldType_1D_Mechanics::x));
	//field_id.push_back(2);
	//field_id.push_back(static_cast<unsigned int>(OutputFieldType_1D_Mechanics::mass));
	field_id.push_back(static_cast<unsigned int>(OutputFieldType_1D_Mechanics::isInMesh));
	field_id.push_back(static_cast<unsigned int>(OutputFieldType_1D_Mechanics::velocity1));
	field_id.push_back(static_cast<unsigned int>(OutputFieldType_1D_Mechanics::stress11));
	field_id.push_back(static_cast<unsigned int>(OutputFieldType_1D_Mechanics::strain11));
	
	out1 = output_request.addOutput("test_output1", 30);
	out1->addObject(&objects[0], &field_id, &pcl_id);
	
	output_request.finish_init();
	output_request.outputMeshGeometry(&mesh);
	
	Solver_1D_Mechanics_1D2_Explicit_FixedMem *solver;
	solver = new Solver_1D_Mechanics_1D2_Explicit_FixedMem(time_step, mesh, objects, output_request);
	solver->solve(0.01);

	delete solver;

	TmpDataToHdf5 trans;
	trans.init(output_request);
	trans.generate();
}

#endif