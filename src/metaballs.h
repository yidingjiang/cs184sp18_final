#ifndef MCUBES_METABALLS_H
#define MCUBES_METABALLS_H
#include "types.h"
#include <stdio.h>
#include <math.h>

class mcubes_metaballs
{
public:
	/* konstrukt�r som initialiserer to baller til predefinerte verdier */
	mcubes_metaballs()
	{
		num_metapoints = 2;
		metapoints[0].power = 1.5;
		metapoints[0].x_pos = 0.0;
		metapoints[0].y_pos = -2.0;
		metapoints[0].z_pos = 0.0;

		metapoints[1].power = 1.5;
		metapoints[1].x_pos = -4.0;
		metapoints[1].y_pos = 2.0;
		metapoints[1].z_pos = 4.0;
	}

	~mcubes_metaballs();

	double iso_value;

	/* setter verdien for hvor grensen mellom ball og ikke ball skal g� */
	void set_iso_value(double iv)
	{
		this->iso_value = iv;
	}

	/* flytter kilden til et gitt metapunkt den angitte avstanden */
	inline void move_ball(int idx, double x, double y, double z)
	{
		metapoints[idx].x_pos += x;
		metapoints[idx].y_pos += y;
		metapoints[idx].z_pos += z;
	}

	/*	returnerer et punkt som er interpolert mellom to andre punkter, 
		b�de for normaler og for posisjon.
	*/
	inline vertex interpolate(vertex v1, vertex v2)
	{
		vertex v;
		double diff;

		diff = (this->iso_value - v1.flux) / (v2.flux - v1.flux);

		/* finner hvor p� linjestykket punktet ligger */
		v.x_pos = v1.x_pos + (v2.x_pos - v1.x_pos) * diff;
		v.y_pos = v1.y_pos + (v2.y_pos - v1.y_pos) * diff;
		v.z_pos = v1.z_pos + (v2.z_pos - v1.z_pos) * diff;
		v.flux = (v1.flux + v2.flux) * 0.5;

		/* regner ut en gjennomsnittsnormal ut fra normalene p� de enkelte punktene */
		v.normal_x = v1.normal_x + (v2.normal_x - v1.normal_x) * diff;
		v.normal_y = v1.normal_y + (v2.normal_y - v1.normal_y) * diff;
		v.normal_z = v1.normal_z + (v2.normal_z - v1.normal_z) * diff;

		return v;
	}

	/* henter ut flux-verdi for et gitt punkt i gridden */
	inline double get_vertex_value(vertex v)
	{
		double flux = 0.0;

		for (i = 0; i < num_metapoints; i++)
		{
			/* regner ut hvor langt det er til dette punktet fra det gitte metapunktet */
			length_x = metapoints[i].x_pos - v.x_pos;
			length_y = metapoints[i].y_pos - v.y_pos;
			length_z = metapoints[i].z_pos - v.z_pos;
		
			/* regner ut styrken som dette punktet p�virkes med */
			flux += fabs(metapoints[i].power) * metapoints[i].power / (
						length_x * length_x + 
						length_y * length_y +
						length_z * length_z + 1);
		}

		return flux;
	}

private:
	int i;
	/* antallet metapunkter */
	int num_metapoints;
	/* angir lengden mellom punkter */
	double length_x, length_y, length_z;
	/* selve metapunktene */
	metapoint metapoints[3];
};
#endif