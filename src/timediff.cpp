#include "..\include\timediff.hpp"

/**
 * Calcula diferencias entre escalas de tiempo astronómicas y de navegación.
 * @param UT1_UTC Diferencia UT1-UTC [s] (DUT1)
 * @param TAI_UTC Diferencia TAI-UTC [s] (número de segundos intercalares + 32.184)
 * @return Tupla con diferencias de tiempo:
 *         [0] UT1-TAI [s]
 *         [1] UTC-GPS [s]
 *         [2] UT1-GPS [s]
 *         [3] TT-UTC [s]
 *         [4] GPS-UTC [s]
 */
tuple<double, double, double, double, double> timediff(double UT1_UTC, double TAI_UTC) {
	const double TT_TAI  = 32.184;          // TT-TAI time difference [s]
    const double GPS_TAI = -19.0;           // GPS-TAI time difference [s]
    const double TT_GPS = TT_TAI - GPS_TAI; // TT-GPS time difference [s]
    const double TAI_GPS = -GPS_TAI;        // TAI-GPS time difference [s]

    double UT1_TAI = UT1_UTC - TAI_UTC;     // UT1-TAI time difference [s]
    double UTC_TAI = -TAI_UTC;              // UTC-TAI time difference [s]
    double UTC_GPS = UTC_TAI - GPS_TAI;     // UTC_GPS time difference [s]
    double UT1_GPS = UT1_TAI - GPS_TAI;     // UT1-GPS time difference [s]
    double TT_UTC  = TT_TAI - UTC_TAI;      // TT-UTC time difference [s]
    double GPS_UTC = GPS_TAI - UTC_TAI;     // GPS-UTC time difference [s]

    return make_tuple(UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
}