#ifndef _SAT_CONST_
#define _SAT_CONST_

// Mathematical constants
#define PI2          (2 * 3.141592653589793)       // 2pi
#define RAD          (3.141592653589793 / 180)     // Radians per degree
#define DEG          (180 / 3.141592653589793)     // Degrees per radian
#define ARCS         (3600 * 180 / 3.141592653589793)  // Arcseconds per radian

// General constants
#define MJD_J2000    51544.5                      // Modified Julian Date of J2000
#define T_B1950      -0.500002108                 // Epoch B1950
#define C_LIGHT      299792458.000000000          // Speed of light [m/s]; DE430
#define AU           149597870700.000000          // Astronomical unit [m]; DE430

// Physical parameters of the Earth, Sun, and Moon

// Equatorial radius and flattening
#define R_EARTH      6378.1363e3                  // Earth's radius [m]; DE430
#define F_EARTH      (1.0 / 298.257223563)        // Flattening; WGS-84
#define R_SUN        696000e3                     // Sun's radius [m]; DE430
#define R_MOON       1738e3                       // Moon's radius [m]; DE430

// Earth rotation (derivative of GMST at J2000; differs from inertial period by precession)
#define OMEGA_EARTH  (15.04106717866910 / 3600 * RAD)  // [rad/s]; WGS-84

// Gravitational coefficients
#define GM_EARTH     398600.435436e9             // [m^3/s^2]; DE430
#define GM_SUN       132712440041.939400e9       // [m^3/s^2]; DE430
#define GM_MOON      (GM_EARTH / 81.30056907419062) // [m^3/s^2]; DE430
#define GM_MERCURY   22031.780000e9              // [m^3/s^2]; DE430
#define GM_VENUS     324858.592000e9             // [m^3/s^2]; DE430
#define GM_MARS      42828.375214e9              // [m^3/s^2]; DE430
#define GM_JUPITER   126712764.800000e9          // [m^3/s^2]; DE430
#define GM_SATURN    37940585.200000e9           // [m^3/s^2]; DE430
#define GM_URANUS    5794548.600000e9            // [m^3/s^2]; DE430
#define GM_NEPTUNE   6836527.100580e9            // [m^3/s^2]; DE430
#define GM_PLUTO     977.0000000000009e9         // [m^3/s^2]; DE430

// Solar radiation pressure at 1 AU 
#define P_SOL        (1367 / C_LIGHT)             // [N/m^2] (~1367 W/m^2); IERS 96

#endif
