#ifndef __GAUSSINTEGRATIONMACRO_H__
#define __GAUSSINTEGRATIONMACRO_H__

// 1 dimensional integration
#define GUASS_INTEGRATION_1D_1ORDER(integration_function) (2.0 * integration_function##(0.0))

#define GUASS_INTEGRATION_1D_2ORDER(integration_function) \
	(integration_function##(-0.577350269189626) + integration_function##(0.577350269189626))

#define GUASS_INTEGRATION_1D_3ORDER(integration_function) \
	 (0.555555555555555 * integration_function##(-0.774596669241483) \
	+ 0.888888888888888 * integration_function##( 0.000000000000000) \
	+ 0.555555555555555 * integration_function##( 0.774596669241483))

// 2 dimensional integration
#define GUASS_INTEGRATION_2D_1ORDER(integration_function) (2.0 * 2.0 * integration_function##(0.0, 0.0))

#define GUASS_INTEGRATION_2D_2ORDER(integration_function) \
	 (integration_function##(-0.577350269189626, -0.577350269189626) \
	+ integration_function##(-0.577350269189626,  0.577350269189626) \
	+ integration_function##( 0.577350269189626, -0.577350269189626) \
	+ integration_function##( 0.577350269189626,  0.577350269189626))

#define GUASS_INTEGRATION_2D_3ORDER(integration_function) \
	 (0.555555555555555 * 0.555555555555555 * integration_function##(-0.774596669241483, -0.774596669241483) \
	+ 0.555555555555555 * 0.888888888888888 * integration_function##(-0.774596669241483,  0.000000000000000) \
	+ 0.555555555555555 * 0.555555555555555 * integration_function##(-0.774596669241483,  0.774596669241483) \
	+ 0.888888888888888 * 0.555555555555555 * integration_function##( 0.000000000000000, -0.774596669241483) \
	+ 0.888888888888888 * 0.888888888888888 * integration_function##( 0.000000000000000,  0.000000000000000) \
	+ 0.888888888888888 * 0.555555555555555 * integration_function##( 0.000000000000000,  0.774596669241483) \
	+ 0.555555555555555 * 0.555555555555555 * integration_function##( 0.774596669241483, -0.774596669241483) \
	+ 0.555555555555555 * 0.888888888888888 * integration_function##( 0.774596669241483,  0.000000000000000) \
	+ 0.555555555555555 * 0.555555555555555 * integration_function##( 0.774596669241483,  0.774596669241483))

#endif