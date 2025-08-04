package journal;

/**
 * Computes multiple parameters related to the orientation of the Earth: TT minus UT1, 
 * nutation angles, mean obliquity, precession angles, and local apparent sidereal time
 * @version 1.1 Summer 2025: fixed an error in the nutation function, and output of TTminusUT1 set to 69.2s between 2018 and 2030
 */
public class EarthAngles {

    /**
     * Computes the difference between Terrestrial Time and Universal Time UT1.
     * A Fix is applied to the year interval 2018-2030, keeping TT-UT1 = 69.2s. A revision will be needed.
     * @param jd The Julian day
     * @return TT-UT1 in seconds
     */
    public static double TTminusUT1(double jd) {
	JulianDay julDay = new JulianDay(jd);
	int year = julDay.year;	
	int month = julDay.month;
	double day = julDay.day + julDay.getDayFraction();

	double TTminusUT1 = 0;
	double ndot = -25.858, c0 = 0.91072 * (ndot + 26.0) ;
	if (year < -500 || year >= 2200) {
		double u = (jd - 2385800.5) / 36525.0; // centuries since J1820
		TTminusUT1 = -20 + 32.0 * u * u;
	} else {
		double x = year + (month - 1 + (day - 1) / 30.0) / 12.0;
		double x2 = x * x, x3 = x2 * x, x4 = x3 * x, x8 = x4 * x4;
		if (year < 1600) {
			TTminusUT1 = 10535.328003 - 9.9952386275 * x + 0.00306730763 * x2 - 7.7634069836E-6 * x3 + 3.1331045394E-9 * x4 + 
				8.2255308544E-12 * x2 * x3 - 7.4861647156E-15 * x4 * x2 + 1.936246155E-18 * x4 * x3 - 8.4892249378E-23 * x8;
		} else {
			TTminusUT1 = -1027175.34776 + 2523.2566254 * x - 1.8856868491 * x2 + 5.8692462279E-5 * x3 + 3.3379295816E-7 * x4 + 
				1.7758961671E-10 * x2 * x3 - 2.7889902806E-13 * x2 * x4 + 1.0224295822E-16 * x3 * x4 - 1.2528102371E-20 * x8;
			// Apply a fix to TT-UT1 between 2018 and 2030 (for now) due to TT-UT1 no longer increasing at the previous rate
			if (year >= 2018 && year <= 2030) TTminusUT1 = 69.2;
		}
		c0 = 0.91072 * (ndot + 25.858) ;
	}
	double c = -c0 * Math.pow((jd - 2435109.0) / 36525.0, 2);
	if (year < 1955 || year > 2005) TTminusUT1 += c;
	return TTminusUT1;
    }
    
    /**
     * Transforms the input Julian day into centuries from J2000 in Terrestrial Time, to compute ephemerides
     * @param jd Julian day
     * @param UT True if input Julian day is in UT, otherwise TT is assumed
     * @return Centuries from J2000 in TT
     */
    public static double toCenturiesRespectJ2000(double jd, boolean UT) {
	if (UT) jd = jd + TTminusUT1(jd) / Constant.SECONDS_PER_DAY;
	return (jd - Constant.J2000) / Constant.JULIAN_DAYS_PER_CENTURY;
	
    }
    
    /**
     * Computes nutation in longitude and obliquity
     * @param jd Julian day in UT
     * @return Nutation angles in radians
     */
    public static double[] nutation(double jd) {
	double t = toCenturiesRespectJ2000(jd, true);
	
	// Mean longitude of the ascending node of the Moon
	double M1 = (124.90 - 1934.134 * t + 0.002063 * t * t) * Constant.DEG_TO_RAD;
	// 2 * Mean longitude of Sun
	double M2 = (201.11 + 72001.5377 * t + 0.00057 * t * t) * Constant.DEG_TO_RAD;
	
	// Compute approximate nutation
	double nutLon = (-(17.2026 + 0.01737 * t) * Math.sin(M1) + (-1.32012 + 0.00013 * t) * Math.sin(M2) + .2088 * Math.sin(2 * M1));
	double nutObl = ((9.2088 + .00091 * t) * Math.cos(M1) + (0.552204 - 0.00029 * t) * Math.cos(M2) - .0904 * Math.cos(2 * M1));
	
	return new double[] {nutLon * Constant.ARCSEC_TO_RAD, nutObl * Constant.ARCSEC_TO_RAD};
    }
    
    /**
     * Returns the mean obliquity in radians
     * @param jd Julian day in UT
     * @return Mean obliquity
     */
    public static double meanObliquity(double jd) {
	double t = toCenturiesRespectJ2000(jd, true);

	// IAU 1976 formulation, still used by Horizons
	//double eps0 = 84381.448;
	//double[] pol = {-468150., -590., 181320.};
	// J. Laskar's expansion comes from "Secular terms of classical planetary theories using the results of general theory," Astronomy and Astrophysics 157, 59070 (1986)
	double eps0 = 84381.448;
	double[] pol = {-468093., -155., 199925., -5138., -24967., -3905., 712., 2787., 579., 245. };
	// Capitaine et al. Astronomy and Astrophysics 412, 567-586, (2003), Hilton et al. 2006, 
	//double eps0 = 84381.406;
	//double[] pol = {-468367.69, -183.1, 200340., -5760., -43400.};

	double meanObliquity = 0;
	for (int i=0; i<pol.length; i++) {
	    meanObliquity += pol[i] * 0.01 * Math.pow(t * 0.01, i + 1);
	}
	return (meanObliquity + eps0) * Constant.ARCSEC_TO_RAD;
    }
    
    /**
     * True obliquity
     * @param jd Julian day
     * @return True obliquity in radians
     */
    public static double trueObliquity(double jd) {
	return meanObliquity(jd) + nutation(jd)[1];
    }
    
    /**
     * Computes the angles to correct for precession between J2000 and another date,
     * using the method by Lieske (IAU 1976)
     * @param jd Julian day in UT
     * @return The output precession angles
     */
    public static double[] precessionAnglesFromJ2000(double jd) {
	double t = toCenturiesRespectJ2000(jd, true);

	double zeta = 2306.2181 * t + 0.30188 * t * t + 0.017998 * t * t * t;
	double z = 2306.2181 * t + 1.09468 * t * t + 0.018203 * t * t * t;
	double theta = 2004.3109 * t - 0.42665 * t * t - 0.041833 * t * t * t;

	zeta = Util.normalizeRadians(zeta * Constant.ARCSEC_TO_RAD);
	z = Util.normalizeRadians(z * Constant.ARCSEC_TO_RAD);
	theta = Util.normalizeRadians(theta * Constant.ARCSEC_TO_RAD);

	return new double[] {zeta, z, theta};
    }

    /**
     * Returns the local apparent sidereal time
     * @param jd Julian day of calculations, in UTC (UT1 if possible)
     * @param obsLon Longitude of the observer in radians
     * @return Apparent sidereal time in radians
     */
    public static double localApparentSiderealTime(double jd, double obsLon) {
	// Obtain local apparent sidereal time
	double jd0 = Math.floor(jd - 0.5) + 0.5; // previous midnight
	double t0 = toCenturiesRespectJ2000(jd0, false); // centuries from previous midnight
	double secs = (jd - jd0) * Constant.SECONDS_PER_DAY;
	double gmst = (((((-6.2e-6 * t0) + 9.3104e-2) * t0) + 8640184.812866) * t0) + 24110.54841;
	double msday = 1.0 + (((((-1.86e-5 * t0) + 0.186208) * t0) + 8640184.812866) / (Constant.SECONDS_PER_DAY * 
		Constant.JULIAN_DAYS_PER_CENTURY));
	gmst = (gmst + msday * secs) * 15.0 * Constant.ARCSEC_TO_RAD;

	// IAU 1994 resolution C7 added two terms (dependent on the mean ascending node of the lunar orbit omega) 
	// to the equation of equinoxes, taking effect since 1997-02-27
	double dt = toCenturiesRespectJ2000(jd, true);
	double omega = 125.04452 - 1934.136261 * dt + 0.0020708 * dt * dt + (dt * dt * dt) / 450000;
	omega = Util.normalizeDegrees(omega);

	double nutLon = nutation(jd)[0]; // First element of the array returned by nutation
	double last = Util.normalizeRadians(
		gmst + obsLon + nutLon * Math.cos(meanObliquity(jd))
		+ 0.00264 * Math.sin(omega) * Constant.ARCSEC_TO_RAD 
		+ 0.000063 * Math.sin(2 * omega) * Constant.ARCSEC_TO_RAD
		);
	
	return last;
    }
    
    /**
     * Test program
     * @param args Not used
     */
    public static void main(String[] args) {
	try {
	    // Input Julian day
	    JulianDay julDay = new JulianDay(2016, 1, 20);
	    double jd = julDay.getJulianDay(); // 2457407.5
	    // Input longitude for the observer in radians. 0 = Greenwich
	    double lon = -(3 + 41 / 60.0 + 18.1 / 3600.0) * Constant.DEG_TO_RAD; // -3.68836 deg

	    // Compute values
	    double ttMinusUt1 = TTminusUT1(jd);
	    double[] nut = nutation(jd);
	    double obl = meanObliquity(jd);
	    double last = localApparentSiderealTime(jd, lon);

	    // Report values
	    // To check TT-UT1:
	    // 	https://maia.usno.navy.mil/ser7/deltat.data
	    // 	https://web.archive.org/web/20220918033245/http://asa.hmnao.com/SecK/DeltaT.html
	    System.out.println("TT-UT1:    "+(float) ttMinusUt1+" s");// 66.6s, 68.18 in Horizons
	    System.out.println("Nutation:  "+Util.formatDEC(nut[0], 2)+", "+
		    Util.formatDEC(nut[1], 2));                       // -00° 00' 00.72", -00° 00' 09.62"
	    System.out.println("Obliquity: "+Util.formatDEC(obl, 2)); // 23° 26' 13.93"
	    System.out.println("LAST:      "+Util.formatRA(last, 3)); // 07h 40m 31.144s, 07:40:31.1426 in USNO, 07 40 31.3197 in Horizons
	    System.out.println();
	    
	    // Check current sidereal time with the web page (for longitude = 0): https://www.localsiderealtime.com/
	    jd = new JulianDay(System.currentTimeMillis()).getJulianDay(); // Julian day now in UTC from the computer
	    double lastNow = localApparentSiderealTime(jd, 0); // For Greenwich
	    
	    System.out.println("GAST now:  "+Util.formatRA(lastNow, 3));
	} catch (Exception exc) {
	    exc.printStackTrace();
	}
    }
}
