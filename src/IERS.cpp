#include "..\include\IERS.hpp"

void IERS(double eop, double Mjd_UTC, double interp) {
	
	if (nargin == 2) {
	   interp = 'n';
	}

	if (interp =='l') {
		// linear interpolation
		mjd = (floor(Mjd_UTC));
		i = find(mjd==eop(4,:),1,'first');
		preeop = eop(:,i);
		nexteop = eop(:,i+1);
		mfme = 1440*(Mjd_UTC-floor(Mjd_UTC));
		fixf = mfme/1440;
		// Setting of IERS Earth rotation parameters
		// (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
		x_pole  = preeop(5)+(nexteop(5)-preeop(5))*fixf;
		y_pole  = preeop(6)+(nexteop(6)-preeop(6))*fixf;
		UT1_UTC = preeop(7)+(nexteop(7)-preeop(7))*fixf;
		LOD     = preeop(8)+(nexteop(8)-preeop(8))*fixf;
		dpsi    = preeop(9)+(nexteop(9)-preeop(9))*fixf;
		deps    = preeop(10)+(nexteop(10)-preeop(10))*fixf;
		dx_pole = preeop(11)+(nexteop(11)-preeop(11))*fixf;
		dy_pole = preeop(12)+(nexteop(12)-preeop(12))*fixf;
		TAI_UTC = preeop(13);
		
		x_pole  = x_pole/ARCS;  // Pole coordinate [rad]
		y_pole  = y_pole/ARCS;  // Pole coordinate [rad]
		dpsi    = dpsi/ARCS;
		deps    = deps/ARCS;
		dx_pole = dx_pole/ARCS; // Pole coordinate [rad]
		dy_pole = dy_pole/ARCS; // Pole coordinate [rad]
	} elseif (interp =='n') {   
		mjd = (floor(Mjd_UTC));
		i = find(mjd==eop(4,:),1,'first');
		eop = eop(:,i);
		// Setting of IERS Earth rotation parameters
		// (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
		x_pole  = eop(5)/ARCS;  // Pole coordinate [rad]
		y_pole  = eop(6)/ARCS;  // Pole coordinate [rad]
		UT1_UTC = eop(7);       // UT1-UTC time difference [s]
		LOD     = eop(8);       // Length of day [s]
		dpsi    = eop(9)/ARCS;
		deps    = eop(10)/ARCS;
		dx_pole = eop(11)/ARCS; // Pole coordinate [rad]
		dy_pole = eop(12)/ARCS; // Pole coordinate [rad]
		TAI_UTC = eop(13);      // TAI-UTC time difference [s]
	}
}