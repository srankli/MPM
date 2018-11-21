#ifdef _DEBUG
#include <iostream>

#include "Solver_1D_Hydromechanics_1D2_Explicit_FixedMem.h"
#include "TmpDataToHdf5.h"

// 1 point
void test_Solver_1D_Hydromechanics_1D2_Explicit_FixedMem_1(void)
{
	double time_step = 2.5;

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
	vbc2.index = 2;
	vbc2.v = -0.1;
	vbc2.dof = DegreeOfFreedom::x;
	//mesh.addVelocityBC(&vbc2);
	vbc2.v = -0.1;
	vbc2.dof = DegreeOfFreedom::x_f;
	mesh.addVelocityBC(&vbc2);
	VelocityBCParam_CubicSmooth vbc3;
	vbc3.index = 2;
	vbc3.dof = DegreeOfFreedom::x_f;
	vbc3.v_begin = 0.0;
	vbc3.v_end = -1.0;
	vbc3.t_length = time_step * 3.0 / 4.0;
	//mesh.addVelocityBC(&vbc3);
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
	pcl_param.Kf = 33.0;
	pcl_param.k = 0.0;
	le_param.E = 1.0;
	le_param.nu = 0.0;
	objects[0].addParticle(&pcl_param, &le_param);
	// Boundary conditions
	MassForceBCParam_Constant mfp_param;
	// particle 1
	mfp_param.dof = DegreeOfFreedom::x;
	mfp_param.massForce = 0.01;
	mfp_param.index = 1;
	//objects[0].addMassForceBC(&mfp_param);
	// Surface forces
	SurfaceForceBCParam_Constant sfp_param1;
	sfp_param1.dof = DegreeOfFreedom::x;
	sfp_param1.surfaceForce = -0.1;
	sfp_param1.index = 1;
	//objects[0].addSurfaceForceBC(&sfp_param1);
	SurfaceForceBCParam_Linear sfp_param2;
	sfp_param2.dof = DegreeOfFreedom::x;
	sfp_param2.surfaceForce_begin = 0.0;
	sfp_param2.surfaceForce_end = 1.0;
	sfp_param2.t_length = time_step / 3.0 * 2.0;
	sfp_param2.index = 1;
	//objects[0].addSurfaceForceBC(&sfp_param2);
	// complete initialization
	objects[0].finish_init();

	OutputRequest output_request("test_result");
	output_request.outputMeshGeometry(&mesh);

	Output *out1;
	std::vector<size_t> pcl_id;
	std::vector<unsigned long long> field_id;

	pcl_id.push_back(1);
	field_id.push_back(static_cast<unsigned long long>(OutputFieldType_1D_Hydromechanics::x));
	field_id.push_back(static_cast<unsigned long long>(OutputFieldType_1D_Hydromechanics::n));
	field_id.push_back(static_cast<unsigned long long>(OutputFieldType_1D_Hydromechanics::p));
	field_id.push_back(static_cast<unsigned long long>(OutputFieldType_1D_Hydromechanics::estress11));
	out1 = output_request.addOutput("test_output1", 10);
	out1->addObject(&objects[0], &field_id, &pcl_id);

	Solver_1D_Hydromechanics_1D2_Explicit_FixedMem *step1;
	step1 = new Solver_1D_Hydromechanics_1D2_Explicit_FixedMem(time_step, mesh, objects, output_request);
	step1->solve(0.1);
	delete step1;

	TmpDataToHdf5 trans;
	trans.init(output_request);
	trans.generate();
}


// 4 point, seepage, no deformation
void test_Solver_1D_Hydromechanics_1D2_Explicit_FixedMem_2(void)
{
	double time_step = 10.0;

	Mesh_1D2 mesh;
	double coords[] = { 0.0, 1.0, 2.0 };
	mesh.initMesh(coords, sizeof(coords) / sizeof(coords[0]));
	// Velocity boundary conditions
	VelocityBCParam_Fixed vbc;
	vbc.index = 1;
	vbc.dof = DegreeOfFreedom::x;
	mesh.addVelocityBC(&vbc);
	vbc.index = 2;
	//mesh.addVelocityBC(&vbc);
	vbc.index = 3;
	mesh.addVelocityBC(&vbc);
	mesh.finish_init();

	std::vector<ObjectByParticle_1D_Hydromechanics> objects;
	ParticleParam_1D_Hydromechanics pcl_param;
	LinearElasticityModelParam le_param;

	objects.resize(1);

	objects[0].setName("testObject1");
	// particle 1
	pcl_param.x = 0.25;
	pcl_param.mass_s = 1.0;
	pcl_param.density_s = 4.0;
	pcl_param.mass_f = 1.0;
	pcl_param.density_f = 4.0;
	pcl_param.Kf = 33.0; // the higher the bulk modulus of fluid, the smaller the time step
	pcl_param.k = 0.01; // the higher the permeability, the smaller the time step
	le_param.E = 1.0;
	le_param.nu = 0.0;
	objects[0].addParticle(&pcl_param, &le_param);
	// particle 2
	pcl_param.x = 0.75;
	objects[0].addParticle(&pcl_param, &le_param);
	// particle 3
	pcl_param.x = 1.25;
	objects[0].addParticle(&pcl_param, &le_param);
	// particle 4
	pcl_param.x = 1.75;
	objects[0].addParticle(&pcl_param, &le_param);
	// Boundary conditions
	// Surface forces
	SurfaceForceBCParam_Constant sfp_param1;
	sfp_param1.dof = DegreeOfFreedom::x_f;
	sfp_param1.surfaceForce = 0.1;
	sfp_param1.index = 4;
	objects[0].addSurfaceForceBC(&sfp_param1);
	SurfaceForceBCParam_Linear sfp_param2;
	sfp_param2.dof = DegreeOfFreedom::x_f;
	sfp_param2.surfaceForce_begin = 0.0;
	sfp_param2.surfaceForce_end = 1.0;
	sfp_param2.t_length = time_step / 3.0 * 2.0;
	sfp_param2.index = 1;
	//objects[0].addSurfaceForceBC(&sfp_param2);
	// complete initialization
	objects[0].finish_init();

	OutputRequest output_request("test_result");
	output_request.outputMeshGeometry(&mesh);

	Output *out1;
	std::vector<size_t> pcl_id;
	std::vector<unsigned long long> field_id;

	out1 = output_request.addOutput("test_output1", 10);
	pcl_id.push_back(1);
	pcl_id.push_back(2);
	pcl_id.push_back(3);
	pcl_id.push_back(4);
	field_id.push_back(static_cast<unsigned long long>(OutputFieldType_1D_Hydromechanics::x));
	field_id.push_back(static_cast<unsigned long long>(OutputFieldType_1D_Hydromechanics::n));
	field_id.push_back(static_cast<unsigned long long>(OutputFieldType_1D_Hydromechanics::p));
	field_id.push_back(static_cast<unsigned long long>(OutputFieldType_1D_Hydromechanics::estress11));
	field_id.push_back(static_cast<unsigned long long>(OutputFieldType_1D_Hydromechanics::velocity1_f));
	out1->addObject(&objects[0], &field_id, &pcl_id);

	Solver_1D_Hydromechanics_1D2_Explicit_FixedMem *step1;
	step1 = new Solver_1D_Hydromechanics_1D2_Explicit_FixedMem(time_step, mesh, objects, output_request);
	step1->solve(0.1);
	delete step1;

	TmpDataToHdf5 trans;
	trans.init(output_request);
	trans.generate();
}


// 1D consolidation, one side impermeable
void test_Solver_1D_Hydromechanics_1D2_Explicit_FixedMem_3(void)
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
	pcl_param.k = 0.0;
	le_param.E = 1.0;
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
	output_request.outputMeshGeometry(&mesh);

	Output *out1;
	std::vector<size_t> pcl_id;
	std::vector<unsigned long long> field_id;

	out1 = output_request.addOutput("test_output1", 10);
	for (i = 0; i < 5; i++) pcl_id.push_back(i + 1);
	field_id.push_back(static_cast<unsigned long long>(OutputFieldType_1D_Hydromechanics::x));
	field_id.push_back(static_cast<unsigned long long>(OutputFieldType_1D_Hydromechanics::n));
	field_id.push_back(static_cast<unsigned long long>(OutputFieldType_1D_Hydromechanics::p));
	field_id.push_back(static_cast<unsigned long long>(OutputFieldType_1D_Hydromechanics::estress11));
	out1->addObject(&objects[0], &field_id, &pcl_id);

	Solver_1D_Hydromechanics_1D2_Explicit_FixedMem *step1;
	step1 = new Solver_1D_Hydromechanics_1D2_Explicit_FixedMem(time_step, mesh, objects, output_request);
	step1->solve(0.1);
	delete step1;

	TmpDataToHdf5 trans;
	trans.init(output_request);
	trans.generate();
}


#endif