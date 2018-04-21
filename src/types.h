#ifndef MCUBES_TYPES_H
#define MCUBES_TYPES_H

/* definerer attributter for et punkt */
typedef struct 
{
	double x_pos;
	double y_pos;
	double z_pos;

	double normal_x;
	double normal_y;
	double normal_z;
	double flux;
	bool inside;
} vertex;

/* definerer attributter for kraftkildene til metaballene */
typedef struct
{
	double x_pos;
	double y_pos;
	double z_pos;

	double power;
} metapoint;

#endif