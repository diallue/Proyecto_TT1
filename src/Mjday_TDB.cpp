#include "..\include\Mjday_TDB.hpp"

/**
 * Convierte el día juliano modificado en Tiempo Terrestre (Mjd_TT) a Tiempo Dinámico Baricéntrico (Mjd_TDB).
 * @param Mjd_TT Día juliano modificado en Tiempo Terrestre (MJD TT).
 * @return Día juliano modificado en Tiempo Dinámico Baricéntrico (MJD TDB).
 */
double Mjday_TDB(double Mjd_TT) {
	double T_TT = (Mjd_TT - MJD_J2000)/36525;
	
	double Mjd_TDB = Mjd_TT + ( 0.001658*sin(628.3076*T_TT + 6.2401) 
                 +   0.000022*sin(575.3385*T_TT+4.2970) 
                 +   0.000014*sin(1256.6152*T_TT + 6.1969) 
                 +   0.000005*sin(606.9777*T_TT+4.0212) 
                 +   0.000005*sin(52.9691*T_TT+0.4444) 
                 +   0.000002*sin(21.3299*T_TT+5.5431)
                 +   0.000010*sin(628.3076*T_TT+4.2490) )/86400;
				 
	return Mjd_TDB;
}