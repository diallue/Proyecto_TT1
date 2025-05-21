#include "..\include\Position.hpp"

/**
 * Calcula el vector de posición en coordenadas cartesianas a partir de coordenadas geodésicas.
 * @param lon Longitud geodésica (en radianes).
 * @param lat Latitud geodésica (en radianes).
 * @param h Altura sobre el elipsoide (en metros).
 * @return Vector de posición cartesiana (3x1, en metros).
 */
Matrix Position(double lon, double lat, double h) {
	double R_equ = R_EARTH;
	double f = F_EARTH;
	
	double e2     = f*(2.0-f);
	double CosLat = cos(lat);
	double SinLat = sin(lat);

	double N = R_equ / sqrt(1.0-e2*SinLat*SinLat);

	Matrix r(3,1);
	r(1,1) =  (         N+h)*CosLat*cos(lon);
	r(2,1) =  (         N+h)*CosLat*sin(lon);
	r(3,1) =  ((1.0-e2)*N+h)*SinLat;
	
	return r;
}