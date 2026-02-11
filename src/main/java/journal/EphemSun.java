package journal;

/**
 * A class to compute the ephemerides of the Sun. This class uses the code inside {@linkplain EphemReduction} by inheritance.
 */
public class EphemSun extends EphemReduction {

    public EphemSun(double jd_utc, double lon, double lat, double alt, TWILIGHT tw, TWILIGHT_MODE twm, int tz) {
	super(jd_utc, lon, lat, alt, tw, twm, tz);
    }

    // Sun data from "Planetary Programs and Tables" by Pierre Bretagnon and Jean-Louis Simon, Willman-Bell, 1986
    private static double[][] sunElements = { 
	new double[] { 403406.0, 0.0, 4.721964, 1.621043 }, new double[] { 195207.0, -97597.0, 5.937458, 62830.348067 }, 
	new double[] { 119433.0, -59715.0, 1.115589, 62830.821524 }, new double[] { 112392.0, -56188.0, 5.781616, 62829.634302 }, 
	new double[] { 3891.0, -1556.0, 5.5474, 125660.5691 }, new double[] { 2819.0, -1126.0, 1.512, 125660.9845 }, 
	new double[] { 1721.0, -861.0, 4.1897, 62832.4766 }, new double[] { 0.0, 941.0, 1.163, .813 }, 
	new double[] { 660.0, -264.0, 5.415, 125659.31 }, new double[] { 350.0, -163.0, 4.315, 57533.85 }, 
	new double[] { 334.0, 0.0, 4.553, -33.931 }, new double[] { 314.0, 309.0, 5.198, 777137.715 }, 
	new double[] { 268.0, -158.0, 5.989, 78604.191 }, new double[] { 242.0, 0.0, 2.911, 5.412 }, 
	new double[] { 234.0, -54.0, 1.423, 39302.098 }, new double[] { 158.0, 0.0, .061, -34.861 }, 
	new double[] { 132.0, -93.0, 2.317, 115067.698 }, new double[] { 129.0, -20.0, 3.193, 15774.337 }, 
	new double[] { 114.0, 0.0, 2.828, 5296.67 }, new double[] { 99.0, -47.0, .52, 58849.27 }, 
	new double[] { 93.0, 0.0, 4.65, 5296.11 }, new double[] { 86.0, 0.0, 4.35, -3980.7 }, 
	new double[] { 78.0, -33.0, 2.75, 52237.69 }, new double[] { 72.0, -32.0, 4.5, 55076.47 }, 
	new double[] { 68.0, 0.0, 3.23, 261.08 }, new double[] { 64.0, -10.0, 1.22, 15773.85 }, 
	new double[] { 46.0, -16.0, .14, 188491.03 }, new double[] { 38.0, 0.0, 3.44, -7756.55 }, 
	new double[] { 37.0, 0.0, 4.37, 264.89 }, new double[] { 32.0, -24.0, 1.14, 117906.27 }, 
	new double[] { 29.0, -13.0, 2.84, 55075.75 }, new double[] { 28.0, 0.0, 5.96, -7961.39 }, 
	new double[] { 27.0, -9.0, 5.09, 188489.81 }, new double[] { 27.0, 0.0, 1.72, 2132.19 }, 
	new double[] { 25.0, -17.0, 2.56, 109771.03 }, new double[] { 24.0, -11.0, 1.92, 54868.56 }, 
	new double[] { 21.0, 0.0, .09, 25443.93 }, new double[] { 21.0, 31.0, 5.98, -55731.43 }, 
	new double[] { 20.0, -10.0, 4.03, 60697.74 }, new double[] { 18.0, 0.0, 4.27, 2132.79 }, 
	new double[] { 17.0, -12.0, .79, 109771.63 }, new double[] { 14.0, 0.0, 4.24, -7752.82 }, 
	new double[] { 13.0, -5.0, 2.01, 188491.91 }, new double[] { 13.0, 0.0, 2.65, 207.81 }, 
	new double[] { 13.0, 0.0, 4.98, 29424.63 }, new double[] { 12.0, 0.0, .93, -7.99 }, 
	new double[] { 10.0, 0.0, 2.21, 46941.14 }, new double[] { 10.0, 0.0, 3.59, -68.29 }, 
	new double[] { 10.0, 0.0, 1.5, 21463.25 }, new double[] { 10.0, -9.0, 2.55, 157208.4 } 
    };
	
    /**
     * Geocentric position of the Sun. Mean equinox and ecliptic of date. Expansion is from "Planetary Programs and Tables", 
     * by Pierre Bretagnon and Jean-Louis Simon, Willman-Bell, 1986. The expansion is valid from 4000 B.C. to 8000 A.D, but 
     * has been modified to reproduce the DE440 integration. The ecliptic latitude is not supposed to be 0 as in the previous 
     * reference, but computed approximately from the first term of the VSOP solution in Jean Meeus's Astronomical Algorithms. 
     * Peak errors are below 2".
     * @return Array with ecliptic longitude, latitude (0), distance, and angular radius of the Sun.
     */
    @Override
    public double[] getBodyPosition() {
	double t = EarthAngles.toCenturiesRespectJ2000(jd_UT, true) / 100.0;
	double L = 0.0, R = 0.0;

	for (int i = 0; i < sunElements.length; i++) {
	    double variable = sunElements[i][2] + sunElements[i][3] * t;
	    double u = Util.normalizeRadians(variable);
	    double sU = Math.sin(u);
	    double cU = Math.cos(u);
	    L = L + sunElements[i][0] * sU;
	    R = R + sunElements[i][1] * cU;
	}

	double tmp = Util.normalizeRadians(62833.196168 * t);
	L = Util.normalizeRadians(4.9353929 + tmp + L * 1E-7);
	R = 1.0001026 + R * 1E-7;
	L -= 0.09 * Constant.ARCSEC_TO_RAD; // Basic transformation from mean dynamical equinox to FK5 J2000, as explained by Meeus
	L += getCorrectionToEclipticLongitude(jd_UT + EarthAngles.TTminusUT1(jd_UT) / Constant.SECONDS_PER_DAY);

	// Compute aberration, substracted later
	double aberration = (993 - 17 * Math.cos(3.10 + 62830.14 * t)) * 1E-7;
	
	return new double[] { L - aberration, getSolarGeocenticEclipticLatitude(t * 10), R, Math.atan(696250.0 / (R * Constant.AU)) };
    }
    
    /**
     * Geometric rectangular position of the Sun, mean ecliptic of date.
     * @param jd Julian day (TT).
     * @return x, y, z, vx, vy, vz, in AU and AU/d.
     */
    public static double[] getGeometricEclipticPositionEquinoxOfDate(double jd) {
	double t = (jd - Constant.J2000) / 3652500.0;
	double L = 0.0, R = 0.0, DL = 0.0, DR = 0.0;

	for (int i = 0; i < sunElements.length; i++) {
	    double variable = sunElements[i][2] + sunElements[i][3] * t;
	    double u = Util.normalizeRadians(variable);
	    double sU = Math.sin(u);
	    double cU = Math.cos(u);
	    L = L + sunElements[i][0] * sU;
	    R = R + sunElements[i][1] * cU;
	    DL = DL + sunElements[i][0] * sunElements[i][3] * cU;
	    DR = DR - sunElements[i][1] * sunElements[i][3] * sU;
	}

	double tmp = Util.normalizeRadians(62833.196168 * t);
	L = Util.normalizeRadians(4.9353929 + tmp + L * 1E-7);
	R = 1.0001026 + R * 1E-7;
	DL = (62833.196168 + DL * 1E-7) / 3652500.0;
	DR = (DR * 1E-7) / 3652500.0;	
	L -= 0.09 * Constant.ARCSEC_TO_RAD; // Basic transformation from mean dynamical equinox to FK5 J2000, as explained by Meeus
	double b = getSolarGeocenticEclipticLatitude(t * 10);
	L += getCorrectionToEclipticLongitude(jd);
	
	// For rectangular coordinates
	double x = R * Math.cos(L) * Math.cos(b);
	double y = R * Math.sin(L) * Math.cos(b);
	double z = R * Math.sin(b);
	double vx = DR * Math.cos(L) - DL * y;
	double vy = DR * Math.sin(L) + DL * x;
	double vz = 0.0;

	return new double[] { x, y, z, vx, vy, vz};
    }
    
    private static double getSolarGeocenticEclipticLatitude(double t) {
	// First 0-order B term from Meeus's Astronomical Algorithms, Appendix II, due to the Earth position respect E-M barycenter when the Moon has 
	// a high geocentric ecliptic latitude. Computed with the mean distance of moon from its ascending node. At most 0.6"
	return -280E-8 * Math.cos(3.199 + 84334.662 * t); // t in Julian millenia
    }
    
    private static double getCorrectionToEclipticLongitude(double jd) {
	return -(9.40 - 0.0000090 * jd + 291E-14 * jd * jd - 321E-21 * jd * jd * jd) * Constant.ARCSEC_TO_RAD; // To fit DE441 (Horizons)
    }
    
    /**
     * Returns the dates of the official (geocentric) equinoxes and solstices.
     * @return Dates of equinoxes and solstices, error always below 50s.
     */
    public double[] getEquinoxesAndSolstices() {
	double jdOld = jd_UT;
	double[] out = new double[4];

	double prec = 0.1 / Constant.SECONDS_PER_DAY; // Output precision 0.1s, accuracy around 1 minute
	JulianDay julDay = new JulianDay(jd_UT);
	int year = julDay.year;

	int[] months = new int[] {3, 9, 6, 12};
	double[] refs = new double[] {0, Math.PI, Constant.PI_OVER_TWO, 3 * Constant.PI_OVER_TWO}; // Will search for these RAs
	for (int i=0; i<4; i++) {
	    julDay = new JulianDay(year, months[i], 18);
	    double jd = julDay.getJulianDay();
	    setUTDate(jd);
	    double min = -1, minT = -1;
	    double stepDays = 0.25, lastRaDif = -1;
	    while (true) {
		EphemData data = doCalc(getBodyPosition(), true);
		double lon = data.rightAscension; // Using RA, replace with next line to use ecliptic longitude
		//double lon = Util.normalizeRadians(CoordinateSystem.cartesianToSpherical(CoordinateSystem.equatorialToEcliptic(CoordinateSystem.sphericalToCartesian(data.rightAscension, data.declination), jd))[0]);
		double raDif = Math.abs(refs[i] - lon);
		if (raDif > Math.PI) raDif = Constant.TWO_PI - raDif;
		if (raDif < min || min == -1) {
		    min = raDif;
		    minT = jd;
		}
		if (raDif > lastRaDif && lastRaDif >= 0) {
		    if (Math.abs(stepDays) < prec) {
			out[i] = minT;
			break;						
		    }
		    stepDays = -stepDays / 2;
		}
		lastRaDif = raDif;
		jd += stepDays;
		setUTDate(jd);
	    }
	}
	setUTDate(jdOld);
	return out;
    }

    /**
     * Test program
     * @param args Not used
     */
    public static void main(String[] args) {
	// Prepare input data
	int year = 2020, month = 6, day = 9, h = 18, m = 0, s = 0;
	JulianDay jday = new JulianDay(year, month, day);
	jday.setDayFraction((h + m / 60.0 + s / 3600.0) / 24.0);
	
	double jd_utc = jday.getJulianDay();
	double lon = -4; // degrees
	double lat = 40;
	double alt = 0; // m
	int tz = 3; // h
	TWILIGHT tw = TWILIGHT.HORIZON_34arcmin;
	TWILIGHT_MODE twm = TWILIGHT_MODE.TODAY_UT;
	
	// Compute the ephemerides data
	EphemSun sunEph = new EphemSun(jd_utc, lon, lat, alt, tw, twm, tz);
	EphemData sunData = EphemReduction.getEphemeris(sunEph);

	// Report
	System.out.println("Sun");
	System.out.println(sunData.toString());
	
	double[] eqxSols = sunEph.getEquinoxesAndSolstices();
	System.out.println("Sprint equinox:  " + EphemData.getDateAsString(eqxSols[0]));
	System.out.println("Autumn equinox:  " + EphemData.getDateAsString(eqxSols[1]));
	System.out.println("Summer solstice: " + EphemData.getDateAsString(eqxSols[2]));
	System.out.println("Winter solstice: " + EphemData.getDateAsString(eqxSols[3]));
	
/*
Sun
 Az:       285.78952°
 El:       17.423506°
 Dist:     1.0152142 au
 RA:       78.40923°
 DEC:      23.007484°
 Ill:      100.0%
 ang.R:    0.26245394°
 Rise:     2020/06/09 04:46:56 UT
 Set:      2020/06/09 19:43:59 UT
 Transit:  2020/06/09 12:15:21 UT (elev. 72.99472°)

Sprint equinox:  2020/03/20 03:49:37 UT
Autumn equinox:  2020/09/22 13:30:48 UT
Summer solstice: 2020/06/20 21:43:35 UT
Winter solstice: 2020/12/21 10:02:34 UT

Horizons:  2020-Jun-09 18:01:09.200 *    78.40923  23.00749  285.788988  17.372201  1.01521568029101   69.184687
 */
    }
}
