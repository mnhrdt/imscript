// 1. functions for converting between geographic coordinate systems
//
// the "basic" system is longitude, latitude, height over the ellipsoid
// other systems are expressed to and from this one
//
// Note: angles are always in degrees, distances are always in meters

// 1.1. map projections (coordinate conversions in 2d)

// UTMZ: Universal traverse mercator with a signed zone (sign=hemisphere)
int utmz_from_lonlat(double out_utm[2], double in_lonlat[2]);
void lonlat_from_utmz(double out_lonlat[2], double in_utm[2], int zone);

// Normalized latitude and longitude (plate carrée)
void nll_from_lonlat(double out_nll[2], double in_lonlat[2]);
void lonlat_from_nll(double out_lonlat[2], double in_nll[2]);

// Traditional mercator projection
void mercator_from_lonlat(double out_mercator[2], double in_lonlat[2]);
void lonlat_from_mercator(double out_lonlat[2], double in_mercator[2]);

// Traditional Gall-Peters (just for trolling)
void peters_from_lonlat(double out_peters[2], double in_lonlat[2]);
void lonlat_from_peters(double out_lonlat[2], double in_peters[2]);

// 1.2. three-dimensional reference frames

// Earth Centered Earth Fixed Cartesian reference frame
void xyz_from_llh(double out_xyz[3], double in_llh[3]);
void llh_from_xyz(double out_llh[3], double in_xyz[3]);

// Tangential reference system: Earth, North, Up
double zyz_from_enu(double out_xyz[3], double ref_xyz[3], double in_enu[3]);
double enu_from_xyz(double out_enu[3], double ref_xyz[3], double in_xyz[3]);



// 2. functions for extracting data from outer sources

// geoid: sea level over the WGS ellipsoid
double egm96(double lon, double lat);

// shuttle radar topography mission
double srtm4(double lon, double lat);
