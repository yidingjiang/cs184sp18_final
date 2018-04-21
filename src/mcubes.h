#include "metaballs.h"

class mcubes
{
public:
	/* konstruktør for objektet, tar initialiseringsverdier som beskrevet lenger ned */
	mcubes(	double start_x, double start_y, double start_z, double end_x, double end_y, double end_z,
			double step_x, double step_y, double step_z);
	virtual ~mcubes();

	/* ikke i bruk, tiltenkt for å kunne resize gridden dynamisk */
	void resize();
	/* tegner metaballene til opengl */
	void draw();
	/* regner ut hvorvidt punktene er innenfor eller utenfor */
	void computeMetaBalls();

	/* setter peker til metaballs-objekt og isoverdien man ønsker å benytte */
	void setMetaBalls(mcubes_metaballs *mb, double iso_value)
	{
		this->mb = mb;
		this->mb->set_iso_value(iso_value);
		this->metaballs_iso_value = iso_value;
	}

	void setWireframe(bool s)
	{
		this->wireframe = s;
	}

private:
	/* egenskaper for gridden, hvilke x,y,z-koordinater (i objekt-rommet) som gridden skal starte fra */
	double start_x;
	double start_y;
	double start_z;

	/* egenskaper for gridden, hvilke x,y,z-koordinater (i objekt-rommet) som gridden skal ende på */
	double end_x;
	double end_y;
	double end_z;

	/* hvor lange steg man skal ta når man beveger seg fra start_* til end_* */
	double step_x;
	double step_y;
	double step_z;

	/* størrelsen på gridden, regnes ut i konstruktøren */
	int size_x;
	int size_y;
	int size_z;

	/* isoverdien for metaballer */
	double metaballs_iso_value;

	/* hvorvidt vi skal tegne ting som en wireframe-struktur */
	bool wireframe;

	/* en peker til et metaball-objekt */
	mcubes_metaballs *mb;

	/* tabellene for oppslag av kanter og triangler */
	const static int edgeTable[256];
	const static int triTable[256][16];

	/* lagrer punktene som er generert */
	vertex *vertices;

	/* lagrer punktene som genereres fra en enkelt kube */
	vertex verts[12];
};