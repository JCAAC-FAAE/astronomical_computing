package journal;

/**
 * A class to compute the ephemerides of the Moon. This class uses the code inside {@linkplain EphemReduction} by inheritance.
 */
public class EphemMoon extends EphemReduction {

    /** The set of phases to compute the moon phases */
    public enum MOONPHASE {
	/** New Moon phase */
	NEW_MOON ("New Moon:        ", 0),
	/** Crescent quarter phase */
	CRESCENT_QUARTER ("Crescent quarter:", 0.25),
	/** Full Moon phase */
	FULL_MOON ("Full Moon:       ", 0.5),
	/** Descent quarter phase */
	DESCENT_QUARTER ("Descent quarter: ", 0.75);

	/** Phase name */
	public String phaseName;
	/** Phase value */
	public double phase;

	private MOONPHASE(String name, double ph) {
	    phaseName = name;
	    phase = ph;
	}
    }
	
    public EphemMoon(double jd_utc, double lon, double lat, double alt, TWILIGHT tw, TWILIGHT_MODE twm, int tz) {
	super(jd_utc, lon, lat, alt, tw, twm, tz);
    }
	
    /**
     * Geocentric position of the Moon. Mean equinox and ecliptic of date.
     * Based on the Petter Duffet program (and mean elements from S. L. Moshier), 
     * modified to fit DE441 during millenia.
     * @return Array with ecliptic longitude, latitude (0), distance, and angular radius of the Sun
     */
    @Override
    protected double[] getBodyPosition() {
	double t = EarthAngles.toCenturiesRespectJ2000(jd_UT, true);
	
	// These expansions up to t^7 for the mean elements are taken from S. L. Moshier
	/* Mean elongation of moon = D */
	double x = (1.6029616009939659e+09 * t + 1.0722612202445078e+06);
	x += (((((-3.207663637426e-013 * t + 2.555243317839e-011) * t + 2.560078201452e-009) * t - 3.702060118571e-005) * t + 6.9492746836058421e-03) * t /* D, t^3 */
			- 6.7352202374457519e+00) * t * t; /* D, t^2 */
	double phase = Util.normalizeRadians(Constant.ARCSEC_TO_RAD * x);

	/* Mean distance of moon from its ascending node = F */
	x = (1.7395272628437717e+09 * t + 3.3577951412884740e+05);
	x += (((((4.474984866301e-013 * t + 4.189032191814e-011) * t - 2.790392351314e-009) * t - 2.165750777942e-006) * t - 7.5311878482337989e-04) * t /* F, t^3 */
			- 1.3117809789650071e+01) * t * t; /* F, t^2 */
	double node = Util.normalizeRadians(Constant.ARCSEC_TO_RAD * x);

	/* Mean anomaly of sun = l' (J. Laskar) */
	x = (1.2959658102304320e+08 * t + 1.2871027407441526e+06);
	x += ((((((((1.62e-20 * t - 1.0390e-17) * t - 3.83508e-15) * t + 4.237343e-13) * t + 8.8555011e-11) * t - 4.77258489e-8) * t - 1.1297037031e-5) * t + 8.7473717367324703e-05) * t - 5.5281306421783094e-01) * t * t;
	double sanomaly = Util.normalizeRadians(Constant.ARCSEC_TO_RAD * x);

	/* Mean anomaly of moon = l */
	x = (1.7179159228846793e+09 * t + 4.8586817465825332e+05);
	x += (((((-1.755312760154e-012 * t + 3.452144225877e-011) * t - 2.506365935364e-008) * t - 2.536291235258e-004) * t + 5.2099641302735818e-02) * t /* l, t^3 */
			+ 3.1501359071894147e+01) * t * t; /* l, t^2 */
	double anomaly = Util.normalizeRadians(Constant.ARCSEC_TO_RAD * x);

	/* Mean longitude of moon, re mean ecliptic and equinox of date = L */
	x = (1.7325643720442266e+09 * t + 7.8593980921052420e+05);
	x += (((((7.200592540556e-014 * t + 2.235210987108e-010) * t - 1.024222633731e-008) * t - 6.073960534117e-005) * t + 6.9017248528380490e-03) * t /* L, t^3 */
			- 5.6550460027471399e+00) * t * t; /* L, t^2 */
	double l = Util.normalizeRadians(Constant.ARCSEC_TO_RAD * x) * Constant.RAD_TO_DEG;
	
	// Now longitude, with the three main correcting terms of evection,
	// variation, and equation of year, plus other terms (error<0.01 deg)
	// P. Duffet's MOON program taken as reference for the periodic terms
	double E = 1.0 - (.002495 + 7.52E-06 * (t + 1.0)) * (t + 1.0), E2 = E * E;
	double td = t + 1, td2 = t * t;
	double M6 = td * Constant.JULIAN_DAYS_PER_CENTURY * 360.0 / 6.798363307E3;
	double NA = Util.normalizeRadians((2.59183275E2 - M6 + (2.078E-3 + 2.2E-6 * td) * td2) * Constant.DEG_TO_RAD);
	double C = NA + Constant.DEG_TO_RAD * (275.05 - 2.3 * td);
	
	l += 6.28875 * Math.sin(anomaly) + 1.274018 * Math.sin(2 * phase - anomaly) + .658309 * Math.sin(2 * phase);
	l +=  0.213616 * Math.sin(2 * anomaly) - E * .185596 * Math.sin(sanomaly) - 0.114336 * Math.sin(2 * node);
	l += .058793 * Math.sin(2 * phase - 2 * anomaly) + .057212 * E * Math.sin(2 * phase - anomaly - sanomaly) + .05332 * Math.sin(2 * phase + anomaly);
	l += .045874 * E * Math.sin(2 * phase - sanomaly) + .041024 * E * Math.sin(anomaly - sanomaly) - .034718 * Math.sin(phase) - E * .030465 * Math.sin(sanomaly + anomaly);
	l += .015326 * Math.sin(2 * (phase - node)) - .012528 * Math.sin(2 * node + anomaly) - .01098 * Math.sin(2 * node - anomaly) + .010674 * Math.sin(4 * phase - anomaly);
	l += .010034 * Math.sin(3 * anomaly) + .008548 * Math.sin(4 * phase - 2 * anomaly);
	l += -E * .00791 * Math.sin(sanomaly - anomaly + 2 * phase) - E * .006783 * Math.sin(2 * phase + sanomaly) + .005162 * Math.sin(anomaly - phase) + E * .005 * Math.sin(sanomaly + phase);
	l += .003862 * Math.sin(4 * phase) + E * .004049 * Math.sin(anomaly - sanomaly + 2 * phase) + .003996 * Math.sin(2 * (anomaly + phase)) + .003665 * Math.sin(2 * phase - 3 * anomaly);
	l += E * 2.695E-3 * Math.sin(2 * anomaly - sanomaly) + 2.602E-3 * Math.sin(anomaly - 2*(node+phase));
	l += E * 2.396E-3 * Math.sin(2*(phase - anomaly) - sanomaly) - 2.349E-3 * Math.sin(anomaly+phase);
	l += E * E * 2.249E-3 * Math.sin(2*(phase-sanomaly)) - E * 2.125E-3 * Math.sin(2*anomaly+sanomaly);
	l += -E * E * 2.079E-3 * Math.sin(2*sanomaly) + E * E * 2.059E-3 * Math.sin(2*(phase-sanomaly)-anomaly);
	l += -1.773E-3 * Math.sin(anomaly+2*(phase-node)) - 1.595E-3 * Math.sin(2*(node+phase));
	l += E * 1.22E-3 * Math.sin(4*phase-sanomaly-anomaly) - 1.11E-3 * Math.sin(2*(anomaly+node));
	l += 8.92E-4 * Math.sin(anomaly - 3 * phase) - E * 8.11E-4 * Math.sin(sanomaly + anomaly + 2 * phase);
	l += E * 7.61E-4 * Math.sin(4 * phase - sanomaly - 2 * anomaly);
	l += E2 * 7.04E-4 * Math.sin(anomaly - 2 * (sanomaly + phase));
	l += E * 6.93E-4 * Math.sin(sanomaly - 2 * (anomaly - phase));
	l += E * 5.98E-4 * Math.sin(2 * (phase - node) - sanomaly);
	l += 5.5E-4 * Math.sin(anomaly + 4 * phase) + 5.38E-4 * Math.sin(4 * anomaly);
	l += E * 5.21E-4 * Math.sin(4 * phase - sanomaly) + 4.86E-4 * Math.sin(2 * anomaly - phase);
	l += E2 * 7.17E-4 * Math.sin(anomaly - 2 * sanomaly);

	// Ecliptic longitude: additional terms to fit DE404
	l += 14.1983 * Math.cos(2.3232 * t + 0.4964) / 3600.0;
	double ph = 0.5582 + 0.11 * t;
	l += 7.2167 * Math.cos(33.8624 * t - ph) / 3600.0;
	// Additional terms to fit DE431
	l += -0.0759777 * Math.pow(t - 1.541336, 2.0) / 3600.0;
	// Previous term can be replaced by this one for a better fit between years -2000 to 6000 (only for DE431)
	//l += (-0.0001617 * t * t * t - 0.075925 * t * t + 0.3932 * t - 0.4375) / 3600.0;
	l += -0.01 * t * t * Math.cos(anomaly) / 3600.0;
	// Fit DE440
	l -= (0.88 - 0.156 * t - 0.00584 * t * t - 0.003128 * t * t * t) / 3600.0;
	l += 0.25 * t * Math.cos(anomaly) / 3600.0;
	
	double longitude = l * Constant.DEG_TO_RAD;
			
	// Now Moon parallax
	double p = .950724 + .051818 * Math.cos(anomaly) + .009531 * Math.cos(2 * phase - anomaly);
	p += .007843 * Math.cos(2 * phase) + .002824 * Math.cos(2 * anomaly);
	p += 0.000857 * Math.cos(2 * phase + anomaly) + E * .000533 * Math.cos(2 * phase - sanomaly);
	p += E * .000401 * Math.cos(2 * phase - anomaly - sanomaly) + E * .00032 * Math.cos(anomaly - sanomaly) - .000271 * Math.cos(phase);
	p += -E * .000264 * Math.cos(sanomaly + anomaly) - .000198 * Math.cos(2 * node - anomaly);
	p += 1.73E-4 * Math.cos(3 * anomaly) + 1.67E-4 * Math.cos(4*phase-anomaly);
	p += -E * 1.11E-4 * Math.cos(sanomaly) + 1.03E-4 * Math.cos(4 * phase - 2 * anomaly);
	p += -8.4E-5 * Math.cos(2 * anomaly - 2 * phase) - E * 8.3E-5 * Math.cos(2 * phase + sanomaly);
	p += 7.9E-5 * Math.cos(2 * phase + 2 * anomaly) + 7.2E-5 * Math.cos(4 * phase);
	p += E * 6.4E-5 * Math.cos(2 * phase - sanomaly + anomaly)
			- E * 6.3E-5 * Math.cos(2 * phase + sanomaly - anomaly);
	p += E * 4.1E-5 * Math.cos(sanomaly + phase) + E * 3.5E-5 * Math.cos(2 * anomaly - sanomaly);
	p += -3.3E-5 * Math.cos(3 * anomaly - 2 * phase) - 3E-5 * Math.cos(anomaly + phase);
	p += -2.9E-5 * Math.cos(2 * (node - phase)) - E * 2.9E-5 * Math.cos(2 * anomaly + sanomaly);
	p += E2 * 2.6E-5 * Math.cos(2 * (phase - sanomaly)) - 2.3E-5 * Math.cos(2 * (node - phase) + anomaly);
	p += E * 1.9E-5 * Math.cos(4 * phase - sanomaly - anomaly);
	
	// Parallax: additional terms to fit DE431+
	double pc = (20 - 1799 * t) * Constant.DEG_TO_RAD;
	p += 0.0000925 * Math.cos(pc);
	
	// So Moon distance in Earth radii is, more or less,
	double distance = 1.0 / Math.sin(p * Constant.DEG_TO_RAD);

	// Ecliptic latitude with nodal phase (error<0.01 deg)
	l = 5.128189 * Math.sin(node) + 0.280606 * Math.sin(node + anomaly) + 0.277693 * Math.sin(anomaly - node);
	l += .173238 * Math.sin(2 * phase - node) + .055413 * Math.sin(2 * phase + node - anomaly);
	l += .046272 * Math.sin(2 * phase - node - anomaly) + .032573 * Math.sin(2 * phase + node);
	l += .017198 * Math.sin(2 * anomaly + node) + .009267 * Math.sin(2 * phase + anomaly - node);
	l += .008823 * Math.sin(2 * anomaly - node) + E * .008247 * Math.sin(2 * phase - sanomaly - node) + .004323 * Math.sin(2 * (phase - anomaly) - node);
	l += .0042 * Math.sin(2 * phase + node + anomaly) + E * .003372 * Math.sin(node - sanomaly - 2 * phase);
	l += E * 2.472E-3 * Math.sin(2 * phase + node - sanomaly - anomaly);
	l += E * 2.222E-3 * Math.sin(2 * phase + node - sanomaly);
	l += E * 2.072E-3 * Math.sin(2 * phase - node - sanomaly - anomaly);
	l += E * 1.877E-3 * Math.sin(node - sanomaly + anomaly) + 1.828E-3 * Math.sin(4 * phase - node - anomaly);
	l += -E * 1.803E-3 * Math.sin(node + sanomaly) - 1.75E-3 * Math.sin(3 * node);
	l += E * 1.57E-3 * Math.sin(anomaly - sanomaly - node) - 1.487E-3 * Math.sin(node + phase);
	l += -E * 1.481E-3 * Math.sin(node + sanomaly + anomaly) + E * 1.417E-3 * Math.sin(node - sanomaly - anomaly);
	l += E * 1.35E-3 * Math.sin(node - sanomaly) + 1.33E-3 * Math.sin(node - phase);
	l += 1.106E-3 * Math.sin(node + 3 * anomaly) + 1.02E-3 * Math.sin(4 * phase - node);
	l += 8.33E-4 * Math.sin(node + 4 * phase - anomaly) + 7.81E-4 * Math.sin(anomaly - 3 * node);
	l += 6.7E-4 * Math.sin(node + 4 * phase - 2 * anomaly) + 6.06E-4 * Math.sin(2 * phase - 3 * node);
	l += 5.97E-4 * Math.sin(2 * (phase + anomaly) - node);
	l += E * 4.92E-4 * Math.sin(2 * phase + anomaly - sanomaly - node)
			+ 4.5E-4 * Math.sin(2 * (anomaly - phase) - node);
	l += 4.39E-4 * Math.sin(3 * anomaly - node) + 4.23E-4 * Math.sin(node + 2 * (phase + anomaly));
	l += 4.22E-4 * Math.sin(2 * phase - node - 3 * anomaly)
			- E * 3.67E-4 * Math.sin(sanomaly + node + 2 * phase - anomaly);
	l += -E * 3.53E-4 * Math.sin(sanomaly + node + 2 * phase) + 3.31E-4 * Math.sin(node + 4 * phase);
	l += E * 3.17E-4 * Math.sin(2 * phase + node - sanomaly + anomaly);
	l += E2 * 3.06E-4 * Math.sin(2 * (phase - sanomaly) - node) - 2.83E-4 * Math.sin(anomaly + 3 * node);
	double W1 = 4.664E-4 * Math.cos(NA);
	double W2 = 7.54E-5 * Math.cos(C);

	// Ecliptic latitude: additional terms to fit DE431
	double M1 = (124.90 - 1934.134 * t + 0.002063 * t * t) * Constant.DEG_TO_RAD;        
	l +=  -8 * Math.sin(M1) * Math.cos(node) / 3600.0; // For DE431
	l +=  -0.007 * t * t * Math.cos(node) / 3600.0; // For DE431
	l +=  0.48 * t * Math.cos(node) / 3600.0; // For DE440
	
	double latitude = l * Constant.DEG_TO_RAD * (1.0 - W1 - W2);
	
	double earthRadius = 6378.1366;
	double moonRadius = 1737.25;
	return new double[] {longitude, latitude, distance * earthRadius / Constant.AU, 
		Math.atan(moonRadius / (distance * earthRadius)), phase};
    }
    
    /**
     * Returns the moon age
     * @param sun Ephemerides of the Sun, or null for a moon age with less precision
     * @return Age in days
     */
    public double getMoonAge(EphemData sun) {
	double[] moon = getBodyPosition();
	double Psin = 29.530588853;
	if (sun != null) {
	    // Get Moon age, more accurate than 'phase' but we need the Sun position
	    return Util.normalizeRadians(moon[0] - sun.eclipticLongitude) * Psin / Constant.TWO_PI;
	} else {
	    // Use the phase variable as estimate, less accurate but this is used only when we don't need an accurate value
	    return moon[4] * Psin / Constant.TWO_PI;			
	}	
    }
    
    /**
     * Returns the instant of a given moon phase.
     * @param phase The phase.
     * @return The instant of that phase, accuracy around 1 minute or better.
     */
    public double getMoonPhaseTime(MOONPHASE phase) {
	double accuracy = 10 / (30 * Constant.SECONDS_PER_DAY); // 10s / lunar cycle length in s => 10s accuracy
	double refPhase = phase.phase;
	double oldJD = jd_UT;
	EphemSun sunEph = new EphemSun(jd_UT, obsLon * Constant.RAD_TO_DEG, obsLat * Constant.RAD_TO_DEG, obsAlt, 
		twilight, twilightMode, timeZone);
	while (true) {
	    double[] sun = sunEph.getBodyPosition();
	    double[] moon = getBodyPosition();
	    double age = Util.normalizeRadians((moon[0] - sun[0])) / Constant.TWO_PI - refPhase;
	    if (age < 0) age += 1;
	    if (age < accuracy || age > 1 - accuracy) break;
	    if (age < 0.5) {
		jd_UT -= age;
	    } else {
		jd_UT += 1-age;
	    }
	    sunEph.setUTDate(jd_UT);
	    setUTDate(jd_UT);
	}
	double out = jd_UT;
	setUTDate(oldJD);
	return out;
    }

    /**
     * Returns the orientation angles of the lunar disk figure. Illumination fraction is 
     * returned in the main program. Simplification of the method presented by Eckhardt, 
     * D. H., "Theory of the Libration of the Moon", Moon and planets 25, 3 (1981), without 
     * the physical librations of the Moon. Accuracy around 0.5 deg for each value. 
     * Moon and Sun positions must be computed before calling this method.
     * @return Optical libration in longitude, latitude, position angle of 
     * axis, bright limb angle, and paralactic angle.
     */
    public double[] getMoonDiskOrientationAngles(EphemData sun) {
	EphemData moon = getEphemeris(this);
	double moonLon = moon.eclipticLongitude, moonLat = moon.eclipticLatitude, 
		moonRA = moon.rightAscension, moonDEC = moon.declination;
	double sunRA = sun.rightAscension, sunDEC = sun.declination;

	// Obliquity of ecliptic
	double t = EarthAngles.toCenturiesRespectJ2000(jd_UT, true);
	double eps = EarthAngles.trueObliquity(jd_UT);
	// Moon's argument of latitude
	double F = (93.2720993 + 483202.0175273 * t - 0.0034029 * t * t - t * t * t / 3526000.0 + t * t * t * t / 863310000.0) * Constant.DEG_TO_RAD;
	// Moon's inclination
	double I = 1.54242 * Constant.DEG_TO_RAD;
	// Moon's mean ascending node longitude
	double omega = (125.0445550 - 1934.1361849 * t + 0.0020762 * t * t + t * t * t / 467410.0 - t * t * t * t / 18999000.0) * Constant.DEG_TO_RAD;

	double cosI = Math.cos(I), sinI = Math.sin(I);
	double cosMoonLat = Math.cos(moonLat), sinMoonLat = Math.sin(moonLat);
	double cosMoonDec = Math.cos(moonDEC), sinMoonDec = Math.sin(moonDEC);

	// Obtain optical librations lp and bp
	double W = moonLon - omega;
	double sinA = Math.sin(W) * cosMoonLat * cosI - sinMoonLat * sinI;
	double cosA = Math.cos(W) * cosMoonLat;
	double A = Math.atan2(sinA, cosA);
	double lp = Util.normalizeRadians(A - F);
	double sinbp = -Math.sin(W) * cosMoonLat * sinI - sinMoonLat * cosI;
	double bp = Math.asin(sinbp);

	// Obtain position angle of axis p
	double x = sinI * Math.sin(omega);
	double y = sinI * Math.cos(omega) * Math.cos(eps) - cosI * Math.sin(eps);
	double w = Math.atan2(x, y);
	double sinp = Math.hypot(x, y) * Math.cos(moonRA - w) / Math.cos(bp);
	double p = Math.asin(sinp);

	// Compute bright limb angle bl
	double bl = (Math.PI + Math.atan2(Math.cos(sunDEC) * Math.sin(moonRA - sunRA), Math.cos(sunDEC) * 
		sinMoonDec * Math.cos(moonRA - sunRA) - Math.sin(sunDEC) * cosMoonDec));

	// Paralactic angle par
	y = Math.sin(lst - moonRA);
	x = Math.tan(obsLat) * cosMoonDec - sinMoonDec * Math.cos(lst - moonRA);
	double par = 0.0;
	if (x != 0.0) {
	    par = Math.atan2(y, x);
	} else {
	    par = (y / Math.abs(y)) * Constant.PI_OVER_TWO;
	}
	return new double[] {lp, bp, p, bl, par};
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
	EphemMoon moonEph = new EphemMoon(jd_utc, lon, lat, alt, tw, twm, tz);
	EphemData moonData = EphemReduction.getEphemeris(moonEph);

	EphemSun sunEph = new EphemSun(jd_utc, lon, lat, alt, tw, twm, tz);
	EphemData sunData = EphemReduction.getEphemeris(sunEph);
	moonData.setIlluminationPhase(sunData);
	
	// Report
	System.out.println("Moon");
	System.out.println(moonData.toString());
	System.out.println(" Moon age: " + (float) moonEph.getMoonAge(sunData) + " days");
	System.out.println();
	
	System.out.println("Closest Moon phases:");
	System.out.println(" New moon:         " + EphemData.getDateAsString(moonEph.getMoonPhaseTime(MOONPHASE.NEW_MOON)));
	System.out.println(" Crescent quarter: " + EphemData.getDateAsString(moonEph.getMoonPhaseTime(MOONPHASE.CRESCENT_QUARTER)));
	System.out.println(" Full moon:        " + EphemData.getDateAsString(moonEph.getMoonPhaseTime(MOONPHASE.FULL_MOON)));
	System.out.println(" Descent quarter:  " + EphemData.getDateAsString(moonEph.getMoonPhaseTime(MOONPHASE.DESCENT_QUARTER)));
/*
Moon
 Az:       64.83035°
 El:       -57.527813°
 Dist:     0.0026317919 au
 RA:       313.0637°
 DEC:      -21.553576°
 Ill:      82.01859%
 ang.R:    0.25632995°
 Rise:     2020/06/09 23:17:41 UT
 Set:      2020/06/09 08:10:47 UT
 Transit:  2020/06/09 03:21:08 UT (elev. 26.579252°)

 Moon age: 18.86005 days

Closest Moon phases:
 New moon:         2020/06/21 06:41:23 UT
 Crescent quarter: 2020/05/30 03:30:03 UT
 Full moon:        2020/06/05 19:12:33 UT
 Descent quarter:  2020/06/13 06:23:07 UT
 */
    }
}
