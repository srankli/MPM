#ifdef _DEBUG
#include <iostream>

#include "Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem.h"
#include "TmpDataToHdf5.h"

// only one element
void test_Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem_1(void)
{
	//size_t i, j;
	double time_step = 1.0;
	Mesh_R2D4 mesh;
	double x_coords[] = { 0.0, 1.0 };
	double y_coords[] = { 0.0, 1.0 };
	mesh.initMesh(x_coords, sizeof(x_coords) / sizeof(x_coords[0]),
		y_coords, sizeof(y_coords) / sizeof(y_coords[0]));

	// Velocity boundary conditions on mesh
	// Fixed the two nodes at the bottom
	VelocityBCParam_Fixed vbc_x;
	VelocityBCParam_Fixed vbc_y;
	vbc_x.dof = DegreeOfFreedom::x;
	vbc_x.index = 1;
	mesh.addVelocityBC(&vbc_x);
	vbc_y.dof = DegreeOfFreedom::y;
	vbc_y.index = 1;
	mesh.addVelocityBC(&vbc_y);
	vbc_x.index = 2;
	mesh.addVelocityBC(&vbc_x);
	vbc_y.index = 2;
	mesh.addVelocityBC(&vbc_y);
	// Assign constant velocity to nodes at the top 
	VelocityBCParam_Constant vbc_x2;
	VelocityBCParam_Constant vbc_y2;
	vbc_x2.v = 0.1;
	vbc_x2.dof = DegreeOfFreedom::x;
	vbc_x2.index = 3;
	mesh.addVelocityBC(&vbc_x2);
	vbc_y2.v = 0.0;
	vbc_y2.dof = DegreeOfFreedom::y;
	vbc_y2.index = 3;
	mesh.addVelocityBC(&vbc_y2);
	vbc_x2.index = 4;
	mesh.addVelocityBC(&vbc_x2);
	vbc_y2.index = 4;
	mesh.addVelocityBC(&vbc_y2);
	mesh.finish_init();

	// Initialize objects
	std::vector<ObjectByParticle_2D_Mechanics> objects;
	ParticleParam_2D_Mechanics pcl_param;
	LinearElasticityModelParam le_param;

	objects.resize(1);

	objects[0].setName("testObject1");
	// Add particles
	pcl_param.mass = 1.0;
	pcl_param.density = 4.0; // 0.5 by 0.5 per particle
	le_param.E = 100.0;
	le_param.nu = 0.3;
	pcl_param.x = 0.25;
	pcl_param.y = 0.25;
	objects[0].addParticle(&pcl_param, &le_param);
	//le_param.E = 200.0;
	pcl_param.x = 0.75;
	pcl_param.y = 0.25;
	objects[0].addParticle(&pcl_param, &le_param);
	//le_param.E = 300.0;
	pcl_param.x = 0.25;
	pcl_param.y = 0.75;
	objects[0].addParticle(&pcl_param, &le_param);
	//le_param.E = 400.0;
	pcl_param.x = 0.75;
	pcl_param.y = 0.75;
	objects[0].addParticle(&pcl_param, &le_param);
	objects[0].finish_init();

	// Output request
	OutputRequest output_request("test_result");
	Output *out1;
	std::vector<size_t> pcl_id;
	std::vector<unsigned long long> field_id;

	for (size_t i = 0; i < 4; i++) pcl_id.push_back(i + 1);
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::x));
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::y));
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::velocity1));
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::velocity2));
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::stress11));
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::stress12));
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::stress22));
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::strain11));
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::strain12));
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::strain22));
	out1 = output_request.addOutput("test_output1", 10);
	out1->addObject(&objects[0], &field_id, &pcl_id);

	output_request.outputMeshGeometry(&mesh);

	Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem *solver;
	solver = new Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem(
		time_step, mesh, objects, output_request);
	solver->solve(0.01);
	delete solver;

	TmpDataToHdf5 trans;
	trans.init(output_request);
	trans.generate();
}


// 2 by 2
void test_Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem_2(void)
{
	size_t i, j;
	double time_step = 1.0;
	Mesh_R2D4 mesh;
	double x_coords[] = { 0.0, 1.0, 2.0 };
	double y_coords[] = { 0.0, 1.0, 2.0 };
	mesh.initMesh(x_coords, sizeof(x_coords) / sizeof(x_coords[0]),
		y_coords, sizeof(y_coords) / sizeof(y_coords[0]));

	// velocity boundary conditions
	VelocityBCParam_Fixed vbc_x1;
	VelocityBCParam_Fixed vbc_y1;
	vbc_x1.dof = DegreeOfFreedom::x;
	vbc_x1.index = 1;
	mesh.addVelocityBC(&vbc_x1);
	vbc_y1.dof = DegreeOfFreedom::y;
	vbc_y1.index = 1;
	mesh.addVelocityBC(&vbc_y1);
	vbc_x1.index = 2;
	mesh.addVelocityBC(&vbc_x1);
	vbc_y1.index = 2;
	mesh.addVelocityBC(&vbc_y1);
	vbc_x1.index = 3;
	mesh.addVelocityBC(&vbc_x1);
	vbc_y1.index = 3;
	mesh.addVelocityBC(&vbc_y1);
	VelocityBCParam_Constant vbc_x2;
	VelocityBCParam_Constant vbc_y2;
	vbc_x2.dof = DegreeOfFreedom::x;
	vbc_x2.v = 0.1;
	vbc_x2.index = 7;
	mesh.addVelocityBC(&vbc_x2);
	vbc_y2.dof = DegreeOfFreedom::y;
	vbc_y2.v = 0;
	vbc_y2.index = 7;
	mesh.addVelocityBC(&vbc_y2);
	vbc_x2.index = 8;
	mesh.addVelocityBC(&vbc_x2);
	vbc_y2.index = 8;
	mesh.addVelocityBC(&vbc_y2);
	vbc_x2.index = 9;
	mesh.addVelocityBC(&vbc_x2);
	vbc_y2.index = 9;
	mesh.addVelocityBC(&vbc_y2);

	mesh.finish_init();

	std::vector<ObjectByParticle_2D_Mechanics> objects;
	ParticleParam_2D_Mechanics pcl_param;
	LinearElasticityModelParam le_param;

	objects.resize(1);

	objects[0].setName("testObject1");
	// Add particles
	pcl_param.mass = 1.0;
	pcl_param.density = 4.0; // 0.5 by 0.5 per particle
	le_param.E = 100.0;
	le_param.nu = 0.1;
	for (j = 0; j < 4; j++)
	{
		for (i = 0; i < 4; i++)
		{
			pcl_param.x = 0.25 + i * 0.5;
			pcl_param.y = 0.25 + j * 0.5;
			objects[0].addParticle(&pcl_param, &le_param);
		}
	}
	// Boundary conditions
	MassForceBCParam_Constant mfp_param;
	// particle 1
	mfp_param.dof = DegreeOfFreedom::x;
	mfp_param.massForce = 1.0;
	mfp_param.index = 1;
	for (j = 0; j < 4; j++)
	{
		for (i = 0; i < 4; i++)
		{
			//objects[0].addMassForceBC(&mfp_param);
			mfp_param.index++;
		}
	}
	
	// complete initialization
	objects[0].finish_init();

	// Output request
	OutputRequest output_request("test_result");
	Output *out1;
	std::vector<size_t> pcl_id;
	std::vector<unsigned long long> field_id;

	for (size_t i = 0; i < 16; i++) pcl_id.push_back(i + 1);
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::x));
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::y));
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::velocity1));
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::velocity2));
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::stress11));
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::stress12));
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::stress22));
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::strain11));
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::strain12));
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::strain22));
	out1 = output_request.addOutput("Output1", 20);
	out1->addObject(&objects[0], &field_id, &pcl_id);

	output_request.outputMeshGeometry(&mesh);

	Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem *solver;
	solver = new Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem(
		time_step, mesh, objects, output_request);
	solver->solve(0.1);
	delete solver;

	TmpDataToHdf5 trans;
	trans.init(output_request);
	trans.generate();
}


// Unixial compression
void test_Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem_3(void)
{
	size_t i, j;
	double time_step = 30.0;

	Mesh_R2D4 mesh;
	double x_coords[] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };
	double y_coords[] = { 0.0, 1.0, 2.0 };
	mesh.initMesh(x_coords, sizeof(x_coords) / sizeof(x_coords[0]),
		y_coords, sizeof(y_coords) / sizeof(y_coords[0]));

	// velocity boundary conditions
	VelocityBCParam_Fixed vbc_x1;
	VelocityBCParam_Fixed vbc_y1;
	vbc_x1.dof = DegreeOfFreedom::x;
	vbc_x1.index = 1;
	mesh.addVelocityBC(&vbc_x1);
	vbc_y1.dof = DegreeOfFreedom::y;
	vbc_y1.index = 1;
	mesh.addVelocityBC(&vbc_y1);
	vbc_x1.index = 11;
	mesh.addVelocityBC(&vbc_x1);
	vbc_y1.index = 11;
	mesh.addVelocityBC(&vbc_y1);
	vbc_x1.index = 21;
	mesh.addVelocityBC(&vbc_x1);
	vbc_y1.index = 21;
	mesh.addVelocityBC(&vbc_y1);
	mesh.finish_init();

	std::vector<ObjectByParticle_2D_Mechanics> objects;
	ParticleParam_2D_Mechanics pcl_param;
	LinearElasticityModelParam le_param;

	objects.resize(1);

	objects[0].setName("testObject1");
	// Add particles
	pcl_param.mass = 1.0;
	pcl_param.density = 4.0; // 0.5 by 0.5 per particle
	le_param.E = 100.0;
	le_param.nu = 0.3;
	for (j = 0; j < 4; j++)
	{
		for (i = 0; i < 6; i++)
		{
			pcl_param.x = 0.25 + i * 0.5;
			pcl_param.y = 0.25 + j * 0.5;
			objects[0].addParticle(&pcl_param, &le_param);
		}
	}
	//objects[0].addParticleFromTriElement(2.0, 3.464, 0.0, 0.0, 4.0, 0.0, &pcl_param, &le_param, 0.7);
	// Boundary conditions
	MassForceBCParam_Constant mfp_param;
	// particle 1
	mfp_param.dof = DegreeOfFreedom::x;
	mfp_param.index = 1;
	mfp_param.massForce = 0.2;
	for (i = 0; i < 6; i++)
	{
		for (j = 0; j < 4; j++)
		{
			//objects[0].addMassForceBC(&mfp_param);
			mfp_param.index++;
		}
	}
	SurfaceForceBCParam_Constant sfp_param1;
	sfp_param1.dof = DegreeOfFreedom::x;
	sfp_param1.surfaceForce = 5.0;
	sfp_param1.index = 6;
	objects[0].addSurfaceForceBC(&sfp_param1);
	sfp_param1.index = 12;
	objects[0].addSurfaceForceBC(&sfp_param1);
	sfp_param1.index = 18;
	objects[0].addSurfaceForceBC(&sfp_param1);
	sfp_param1.index = 24;
	objects[0].addSurfaceForceBC(&sfp_param1);
	// complete initialization
	objects[0].finish_init();

	// Output request
	OutputRequest output_request("test_result");
	Output *out;
	std::vector<size_t> pcl_id;
	std::vector<unsigned long long> field_id;


	for (size_t i = 0; i < 24; i++)	pcl_id.push_back(i + 1);
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::x));
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::y));
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::isInMesh));
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::velocity1));
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::velocity2));
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::stress11));
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::stress12));
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::stress22));
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::strain11));
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::strain12));
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::strain22));
	out = output_request.addOutput("Output1", 20);
	out->addObject(&objects[0], &field_id, &pcl_id);

	output_request.outputMeshGeometry(&mesh);

	Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem *solver;
	solver = new Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem(
		time_step, mesh, objects, output_request);
	solver->solve(0.01);
	delete solver;

	TmpDataToHdf5 trans;
	trans.init(output_request);
	trans.generate();
}

// Cantilever beam
void test_Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem_4(void)
{
	size_t i, j;
	double time_step = 100.0;

	Mesh_R2D4 mesh;
	double x_coords[] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
	double y_coords[] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
	mesh.initMesh(x_coords, sizeof(x_coords) / sizeof(x_coords[0]),
		y_coords, sizeof(y_coords) / sizeof(y_coords[0]));

	// velocity boundary conditions
	VelocityBCParam_Fixed vbc_x1;
	VelocityBCParam_Fixed vbc_y1;
	vbc_x1.dof = DegreeOfFreedom::x;
	vbc_x1.index = 45;
	mesh.addVelocityBC(&vbc_x1);
	vbc_y1.dof = DegreeOfFreedom::y;
	vbc_y1.index = 45;
	mesh.addVelocityBC(&vbc_y1);
	vbc_x1.index = 56;
	mesh.addVelocityBC(&vbc_x1);
	vbc_y1.index = 56;
	mesh.addVelocityBC(&vbc_y1);
	vbc_x1.index = 67;
	mesh.addVelocityBC(&vbc_x1);
	vbc_y1.index = 67;
	mesh.addVelocityBC(&vbc_y1);
	mesh.finish_init();

	std::vector<ObjectByParticle_2D_Mechanics> objects;
	ParticleParam_2D_Mechanics pcl_param;
	LinearElasticityModelParam le_param;

	objects.resize(1);

	objects[0].setName("testObject1");
	// Add particles
	pcl_param.mass = 1.0;
	pcl_param.density = 4.0; // 0.5 by 0.5 per particle
	le_param.E = 100.0;
	le_param.nu = 0.3;
	for (j = 0; j < 4; j++)
	{
		for (i = 0; i < 16; i++)
		{
			pcl_param.x = 0.25 + i * 0.5;
			pcl_param.y = 4.0 + 0.25 + j * 0.5;
			objects[0].addParticle(&pcl_param, &le_param);
		}
	}
	//objects[0].addParticleFromTriElement(2.0, 3.464, 0.0, 0.0, 4.0, 0.0, &pcl_param, &le_param, 0.7);
	// Boundary conditions
	MassForceBCParam_Constant mfp_param;
	// particle 1
	mfp_param.dof = DegreeOfFreedom::y;
	mfp_param.massForce = -0.2;
	mfp_param.index = 1;
	for (j = 0; j < 4; j++)
	{
		for (i = 0; i < 16; i++)
		{
			//objects[0].addMassForceBC(&mfp_param);
			mfp_param.index++;
		}
	}
	SurfaceForceBCParam_Constant sfp_param1;
	sfp_param1.dof = DegreeOfFreedom::y;
	sfp_param1.surfaceForce = -0.1;
	sfp_param1.index = 16;
	objects[0].addSurfaceForceBC(&sfp_param1);
	sfp_param1.index = 32;
	objects[0].addSurfaceForceBC(&sfp_param1);
	sfp_param1.index = 48;
	objects[0].addSurfaceForceBC(&sfp_param1);
	sfp_param1.index = 64;
	objects[0].addSurfaceForceBC(&sfp_param1);
	// complete initialization
	objects[0].finish_init();

	// Output request
	OutputRequest output_request("test_result");
	Output *out;
	std::vector<size_t> pcl_id;
	std::vector<unsigned long long> field_id;

	for (size_t i = 0; i < 4*16; i++)	pcl_id.push_back(i + 1);
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::x));
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::y));
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::isInMesh));
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::velocity1));
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::velocity2));
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::stress11));
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::stress12));
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::stress22));
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::strain11));
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::strain12));
	field_id.push_back((unsigned long long)(OutputFieldType_2D_Mechanics::strain22));
	out = output_request.addOutput("Output1", 100);
	out->addObject(&objects[0], &field_id, &pcl_id);

	output_request.outputMeshGeometry(&mesh);

	Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem *solver;
	solver = new Solver_2DPlaneStrain_Mechanics_R2D4_Explicit_FixedMem(
		time_step, mesh, objects, output_request);
	solver->solve(0.01);
	delete solver;

	TmpDataToHdf5 trans;
	trans.init(output_request);
	trans.generate();
}


#endif // _DEBUG