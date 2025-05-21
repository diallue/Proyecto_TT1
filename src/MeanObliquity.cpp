#include "..\include\MeanObliquity.hpp"

/**
 * Calcula la oblicuidad media de la eclíptica para una fecha dada.
 * @param Mjd_TT Fecha en días julianos modificados (Tiempo Terrestre).
 * @return Oblicuidad media en radianes.
 */
double MeanObliquity(double Mjd_TT) {
	
	double T = (Mjd_TT-MJD_J2000)/36525;

	double MOblq = RAD *( 84381.448/3600-(46.8150+(0.00059-0.001813*T)*T)*T/3600 );
	
	return MOblq;
}