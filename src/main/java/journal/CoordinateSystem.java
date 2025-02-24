package journal;

/**
 * A class to perform coordinates transformations between equatorial, ecliptic, galactic, and supergalatic systems. 
 * The orientation of the galactic plane follows Jia-Cheng Liu et al. 2010 (see http://arxiv.org/abs/1010.3773).
 * The orientation of the supergalactic pole follows Lahav et al 2000 (see https://arxiv.org/abs/astro-ph/9809343).
 */
public class CoordinateSystem {

    // J2000 position of the galactic pole in FK5 system. See Jia-Cheng Liu et al. 2010,
    // http://arxiv.org/abs/1010.3773. Distance to galactic center was set to 25830 +/- 500
    // light years (VERA collaboration, 2020)
    // FK5
    private static final double GALACTIC_POLE_RA_J2000 = 192.85948120833334 * Constant.DEG_TO_RAD;
    private static final double GALACTIC_POLE_DEC_J2000 = 27.128251194444445 * Constant.DEG_TO_RAD;
    private static final double GALACTIC_NODE_J2000 = 122.93191857 * Constant.DEG_TO_RAD;
    // Redefinition by Jia-Cheng Liu et al.
    //private static final double GALACTIC_POLE_RA_J2000 = 192.90297999208332 * Constant.DEG_TO_RAD;
    //private static final double GALACTIC_POLE_DEC_J2000 = 27.103109214444444 * Constant.DEG_TO_RAD;
    //private static final double GALACTIC_NODE_J2000 = 123.0075021536 * Constant.DEG_TO_RAD;
    // Galactic position of the supergalactic pole. See Lahav et al 2000, MNRAS 312 166-176, 
    // https://arxiv.org/abs/astro-ph/9809343. Definition by G. de Vaucouleurs 1991.
    private static final double SUPER_GALACTIC_POLE_RA = 47.37 * Constant.DEG_TO_RAD;
    private static final double SUPER_GALACTIC_POLE_DEC = 6.32 * Constant.DEG_TO_RAD;
    
    /**
     * Returns a 3x3 pure rotation matrix along axis X.
     * @param angle The angle to rotate.
     * @return The matrix.
     */
    public static double[][] getRotX(double angle) {
	return new double[][] {
	    new double[] {1.0, 0.0, 0.0},
	    new double[] {0.0, Math.cos(angle), Math.sin(angle)},
	    new double[] {0.0, -Math.sin(angle), Math.cos(angle)}};
    }

    /**
     * Returns a 3x3 pure rotation matrix along axis Y.
     * @param angle The angle to rotate.
     * @return The matrix.
     */
    public static double[][] getRotY(double angle) {
	return new double[][] {
	    new double[] {Math.cos(angle), 0.0, -Math.sin(angle)},
	    new double[] {0.0, 1.0, 0.0},
	    new double[] {Math.sin(angle), 0.0, Math.cos(angle)}};
    }

    /**
     * Returns a 3x3 pure rotation matrix along axis Z.
     * @param angle The angle to rotate.
     * @return The matrix.
     */
    public static double[][] getRotZ(double angle) {
	return new double[][] {
	    new double[] {Math.cos(angle), Math.sin(angle), 0.0},
	    new double[] {-Math.sin(angle), Math.cos(angle), 0.0},
	    new double[] {0.0, 0.0, 1.0}};
    }

    /**
     * Multiplication of a vector with a matrix
     * @param p The vector, with 3 (position) or 6 (position and velocity) components
     * @param m The 3x3 rotation matrix
     * @return The result of the rotation
     */
    public static double[] rotate(double[] p, double[][] m) {
	double[] out = new double[p.length];
	for (int i=0; i<p.length; i++) {
	    int ip = i;
	    if (i > 2) ip = i - 3;
	    out[i] = 0;
	    for (int j=0; j<3; j++) {
		out[i] += m[ip][j] * p[j];
	    }
	}
	return out;
    }
    
    /**
     * Transform coordinates from cartesian, x y z, to spherical, lon lat r
     * @param p x, y, z
     * @return lon, lat (radians), r
     */
    public static double[] cartesianToSpherical(double[] p) {
	double lon, lat;

	double x = p[0], y = p[1], z = p[2];
	if (y != 0.0 || x!= 0.0) {
	    double h = Math.hypot(x, y);
	    lon = Math.atan2(y, x);
	    lat = Math.atan2(z, h);
	} else {
	    lon = 0.0;
	    lat = Constant.PI_OVER_TWO;
	    if (z < 0.0) lat = -lat;
	}
	double rad = Math.sqrt(x * x + y * y + z * z);

	return new double[] {lon, lat, rad};
    }

    /**
     * Transforms coordinates from spherical, lon, lat, to cartesian, x y z
     * @param lon Longitude in radians
     * @param lat Latitude in radians
     * @return x y z coordinates, with unity as the norm of the vector
     */
    public static double[] sphericalToCartesian(double lon, double lat) {
	double rad = 1;
	double cl = Math.cos(lat);
	double x = rad * Math.cos(lon) * cl;
	double y = rad * Math.sin(lon) * cl;
	double z = rad * Math.sin(lat);

	return new double[] { x, y, z };
    }

    /**
     * Transforms coordinates from equatorial to ecliptic
     * @param p Rectangular coordinates
     * @param jd Julian day
     * @return Ecliptic rectangular coordinates
     */
    public static double[] equatorialToEcliptic(double[] p, double jd) {
	return rotate(p, getRotX(EarthAngles.trueObliquity(jd)));
    }

    /**
     * Transforms coordinates from ecliptic to equatorial
     * @param p Rectangular coordinates
     * @param jd Julian day
     * @return Equatorial rectangular coordinates
     */
    public static double[] eclipticToEquatorial(double[] p, double jd) {
	return rotate(p, getRotX(-EarthAngles.trueObliquity(jd)));
    }
    
    /**
     * Transforms coordinates from equatorial to horizontal
     * @param p Rectangular coordinates
     * @param jd Julian day
     * @param obsLon Observer's longitude in radians
     * @param obsLat Observer's latitude in radians
     * @return Horizontal rectangular coordinates
     */
    public static double[] equatorialToHorizontal(double[] p, double jd, double obsLon, double obsLat) {
	double[] out = rotate(rotate(p, getRotZ(EarthAngles.localApparentSiderealTime(jd, obsLon))), getRotY((Constant.PI_OVER_TWO - obsLat)));
	out = cartesianToSpherical(out);
	return sphericalToCartesian(Math.PI - out[0], out[1]);
    }

    /**
     * Transforms coordinates from horizontal to equatorial
     * @param p Rectangular coordinates
     * @param jd Julian day
     * @param obsLon Observer's longitude in radians
     * @param obsLat Observer's latitude in radians
     * @return Equatorial rectangular coordinates
     */
    public static double[] horizontalToEquatorial(double[] p, double jd, double obsLon, double obsLat) {
	double[] in = cartesianToSpherical(p);
	in = sphericalToCartesian(Math.PI - in[0], in[1]);
	return rotate(rotate(in, getRotY(-(Constant.PI_OVER_TWO - obsLat))), getRotZ(-EarthAngles.localApparentSiderealTime(jd, obsLon)));
    }
    
    /**
     * Transforms coordinates from J2000 equatorial to galactic
     * @param p Rectangular coordinates
     * @return Galactic rectangular coordinates
     */
    public static double[] equatorialJ2000ToGalactic(double[] p) {
	return rotate(rotate(rotate(p, getRotZ(GALACTIC_POLE_RA_J2000)), 
		getRotY(Constant.PI_OVER_TWO - GALACTIC_POLE_DEC_J2000)), getRotZ(Math.PI-GALACTIC_NODE_J2000));
    }

    /**
     * Transforms coordinates from galactic to J2000 equatorial
     * @param p Rectangular coordinates
     * @return Equatorial rectangular coordinates
     */
    public static double[] galacticToEquatorialJ2000(double[] p) {
	return rotate(rotate(rotate(p, getRotZ(-(Math.PI - GALACTIC_NODE_J2000))), 
		getRotY(-(Constant.PI_OVER_TWO - GALACTIC_POLE_DEC_J2000))), getRotZ(-GALACTIC_POLE_RA_J2000));
    }

    /**
     * Transforms coordinates from J2000 equatorial to supergalactic
     * @param p Rectangular coordinates
     * @return Supergalactic rectangular coordinates
     */
    public static double[] equatorialJ2000ToSupergalactic(double[] p) {
	return rotate(rotate(rotate(equatorialJ2000ToGalactic(p), getRotZ(SUPER_GALACTIC_POLE_RA)), 
		getRotY(Constant.PI_OVER_TWO - SUPER_GALACTIC_POLE_DEC)), getRotZ(Constant.PI_OVER_TWO));
    }

    /**
     * Transforms coordinates from supergalactic to J2000 equatorial
     * @param p Rectangular coordinates
     * @return J2000 equatorial rectangular coordinates
     */
    public static double[] supergalacticToEquatorialJ2000(double[] p) {
	return galacticToEquatorialJ2000(rotate(rotate(rotate(p, getRotZ(-Constant.PI_OVER_TWO)), 
		getRotY(-(Constant.PI_OVER_TWO - SUPER_GALACTIC_POLE_DEC))), getRotZ(-SUPER_GALACTIC_POLE_RA)));
    }

    /**
     * Transforms coordinates from J2000 equatorial to J1900
     * @param p Rectangular coordinates
     * @return J1900 rectangular coordinates
     */
    public static double[] equatorialJ2000ToEquatorialJ1900(double[] p) {
	double[] pa = EarthAngles.precessionAnglesFromJ2000(Constant.J2000 - Constant.JULIAN_DAYS_PER_CENTURY);
	return rotate(rotate(rotate(p, getRotZ(-pa[1])), getRotY(pa[2])), getRotZ(-pa[0]));
    }

    /**
     * Transforms coordinates from J1900 equatorial to J2000
     * @param p Rectangular coordinates
     * @return J2000 rectangular coordinates
     */
    public static double[] equatorialJ1900ToEquatorialJ2000(double[] p) {
	double[] pa = EarthAngles.precessionAnglesFromJ2000(Constant.J2000 - Constant.JULIAN_DAYS_PER_CENTURY);
	return rotate(rotate(rotate(p, getRotZ(pa[0])), getRotY(-pa[2])), getRotZ(pa[1]));
    }
    
    /**
     * Test program
     * @param s Not used
     */
    public static void main(String[] s) {
	try {
	    /*
	    Example fata from Horizons, for galactic better test from https://ned.ipac.caltech.edu/forms/calculator.html
	    Date__(UT)__HR:MN     R.A._(a-appar)_DEC.  Azi____(a-app)___Elev  L_Ap_Sid_Time        TDB-UT     ObsEcLon    ObsEcLat      GlxLon     GlxLat  L_Ap_Hour_Ang
 	    2016-Jan-20 00:00  m  205.73877 -35.68542  118.369371 -22.625219   7.6753665728     68.184474  217.2137137 -23.2621310  314.366606  26.140141   -6.040551182
 	    2000-Jan-20 00:00  m   31.19459 -13.82278  255.493146  -4.310376   7.6670225044     64.184432   23.8036275 -24.7875597  178.724053 -68.320060    5.587383271
	    */
	    double jd = 2457407.5;
	    double obsLon = -3.6879 * Constant.DEG_TO_RAD;
	    double obsLat = 40.408414 * Constant.DEG_TO_RAD;
	    double[] eq = sphericalToCartesian(205.73877 * Constant.DEG_TO_RAD, -35.68542 * Constant.DEG_TO_RAD);
	    // Uncomment for a second test
	    //jd = 2451563.5;
	    //eq = sphericalToCartesian(31.19459 * Constant.DEG_TO_RAD, -13.82278 * Constant.DEG_TO_RAD);
	    
	    // Equatorial to horizontal
	    double[] hor = equatorialToHorizontal(eq, jd, obsLon, obsLat);
	    double[] shor = cartesianToSpherical(hor);
	    System.out.println("Acimut:    " + Util.formatValue(shor[0] * Constant.RAD_TO_DEG, 3));
	    System.out.println("Elevation: " + Util.formatValue(shor[1] * Constant.RAD_TO_DEG, 3));
	    // Invert back to equatorial
	    eq = horizontalToEquatorial(hor, jd, obsLon, obsLat);
	    double[] seq = cartesianToSpherical(eq);
	    System.out.println("RA:  " + Util.formatValue(seq[0] * Constant.RAD_TO_DEG, 3));
	    System.out.println("DEC: " + Util.formatValue(seq[1] * Constant.RAD_TO_DEG, 3));
	    // To ecliptic, using the ecliptic of date (true obliquity), as in the above results from Horizons
	    double[] ecl = equatorialToEcliptic(eq, jd);
	    double[] secl = cartesianToSpherical(ecl);
	    System.out.println("Ecliptic longitude: " + Util.formatValue(secl[0] * Constant.RAD_TO_DEG, 3));
	    System.out.println("Ecliptic latitude:  " + Util.formatValue(secl[1] * Constant.RAD_TO_DEG, 3));
	    eq = horizontalToEquatorial(hor, jd, obsLon, obsLat);
	    seq = cartesianToSpherical(eq);
	    System.out.println("RA:  " + Util.formatValue(seq[0] * Constant.RAD_TO_DEG, 3));
	    System.out.println("DEC: " + Util.formatValue(seq[1] * Constant.RAD_TO_DEG, 3));
	    // To galactic
	    //eq = sphericalToCartesian(192.85948120833334 * Constant.DEG_TO_RAD, 27.128251194444445 * Constant.DEG_TO_RAD);
	    double[] gal = equatorialJ2000ToGalactic(eq);
	    double[] sgal = cartesianToSpherical(gal);
	    System.out.println("Galactic longitude: " + Util.formatValue(sgal[0] * Constant.RAD_TO_DEG, 3));
	    System.out.println("Galactic latitude:  " + Util.formatValue(sgal[1] * Constant.RAD_TO_DEG, 3));
	    eq = galacticToEquatorialJ2000(gal);
	    seq = cartesianToSpherical(eq);
	    System.out.println("RA:  " + Util.formatValue(seq[0] * Constant.RAD_TO_DEG, 3));
	    System.out.println("DEC: " + Util.formatValue(seq[1] * Constant.RAD_TO_DEG, 3));
	    // To supergalactic
	    double[] supergal = equatorialJ2000ToSupergalactic(eq);
	    double[] ssupergal = cartesianToSpherical(supergal);
	    System.out.println("Supergalactic longitude: " + Util.formatValue(ssupergal[0] * Constant.RAD_TO_DEG, 3));
	    System.out.println("Supergalactic latitude:  " + Util.formatValue(ssupergal[1] * Constant.RAD_TO_DEG, 3));
	    eq = supergalacticToEquatorialJ2000(supergal);
	    seq = cartesianToSpherical(eq);
	    System.out.println("RA:  " + Util.formatValue(seq[0] * Constant.RAD_TO_DEG, 3));
	    System.out.println("DEC: " + Util.formatValue(seq[1] * Constant.RAD_TO_DEG, 3));
	    
	    System.out.println();
	    double[] eq1900 = equatorialJ2000ToEquatorialJ1900(eq);
	    seq = cartesianToSpherical(eq1900);
	    System.out.println("RA:  " + Util.formatValue(seq[0] * Constant.RAD_TO_DEG, 3));
	    System.out.println("DEC: " + Util.formatValue(seq[1] * Constant.RAD_TO_DEG, 3));
	    eq = equatorialJ1900ToEquatorialJ2000(eq1900);
	    seq = cartesianToSpherical(eq);
	    System.out.println("RA:  " + Util.formatValue(seq[0] * Constant.RAD_TO_DEG, 3));
	    System.out.println("DEC: " + Util.formatValue(seq[1] * Constant.RAD_TO_DEG, 3));
	} catch (Exception exc) {
	    exc.printStackTrace();
	}
    }
}
