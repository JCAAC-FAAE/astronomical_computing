package journal;

/**
 * A class to reduce the position of a body to obtain the ephemerides as visible from a given observer, include rise/set times.
 */
public class EphemReduction {

    /** The set of twilights to calculate (types of rise/set events) */
    public enum TWILIGHT {
	/** Identifier to compute rise/set times for astronomical twilight (center of the body at -18 degrees of geometrical elevation) */
	ASTRONOMICAL(-18),
	/** Identifier to compute rise/set times for nautical twilight (center of the body at -12 degrees of geometrical elevation) */
	NAUTICAL(-12),
	/** Identifier to compute rise/set times for civil twilight (center of the body at -6 degrees of geometrical elevation) */
	CIVIL(-6),
	/** The standard value of 34' for the refraction at the local horizon */
	HORIZON_34arcmin(-34.0 / 60.0);
	
	/** Target elevation of the center of the body in degrees */
	private double elevation;
	
	TWILIGHT(double elev) {
	    elevation = elev;
	}
    }

    /** Possible options to return the rise/set/transit times */
    public enum TWILIGHT_MODE {
	/** Closest events */
	CLOSEST, 
	/** Compute events for the current date in UT */
	TODAY_UT, 
	/** Compute events for the current date in LT */
	TODAY_LT;
    }

    /** The set of events to calculate (rise/set/transit events) */
    public enum EVENT { RISE, SET, TRANSIT }
	
    protected double jd_UT;
    private double nutLon;
    protected double obsLon, obsLat, obsAlt;
    protected double lst;
    protected TWILIGHT twilight;
    protected TWILIGHT_MODE twilightMode;
    protected int timeZone = 0; /** Time zone for option {@linkplain TWILIGHT_MODE#TODAY_LT}, LT-UT, hours. */

    /**
     * The constructor with the data for the ephemerides reduction process
     * @param jd_utc
     * @param lon
     * @param lat
     * @param alt
     * @param tw
     * @param twm
     * @param tz
     */
    public EphemReduction(double jd_utc, double lon, double lat, double alt, TWILIGHT tw, TWILIGHT_MODE twm, int tz) {
	obsLon = lon * Constant.DEG_TO_RAD;
	obsLat = lat * Constant.DEG_TO_RAD;
	obsAlt = alt;
	twilight = tw;
	twilightMode = twm;
	timeZone = tz;
	setUTDate(jd_utc);
    }

    /**
     * Sets the UT date from the provided Julian day and computes the nutation and sidereal time
     * @param jd The new Julian day in UT
     */
    public void setUTDate(double jd) {
	this.jd_UT = jd;
	double[] nut = EarthAngles.nutation(jd);
	nutLon = nut[0];
	lst = EarthAngles.localApparentSiderealTime(jd, obsLon);
    }
	
    /**
     * Compute the position of the body
     * @param pos Values for the ecliptic longitude, latitude, distance and so on from previous methods for the specific body
     * @param geocentric True to return geocentric position. Set this to false generally
     * @return The Ephem object with the output position. The rise/set/transit times returned here are only valid for non-moving bodies
     */
    public EphemData doCalc(double[] pos, boolean geocentric) {
	// Correct for nutation in longitude
	pos[0] = pos[0] + nutLon;

	// Ecliptic to equatorial coordinates using true obliquity
	double[] xyz = CoordinateSystem.eclipticToEquatorial(CoordinateSystem.sphericalToCartesian(pos[0], pos[1]), jd_UT);

	// Obtain topocentric rectangular coordinates
	double geocLat = (obsLat - .1925 * Math.sin(2 * obsLat) * Constant.DEG_TO_RAD);
	double geocR = 1.0 - Math.pow(Math.sin(obsLat), 2) / 298.257;
	double eradius = (geocR * Constant.EARTH_RADIUS + obsAlt * 0.001);
	double radiusAU = eradius / (pos[2] * Constant.AU);
	if (!geocentric) {
	    double cosLat = Math.cos(geocLat); 
	    double[] correction = new double[] {
		    radiusAU * cosLat * Math.cos(lst),
		    radiusAU * cosLat * Math.sin(lst),
		    radiusAU * Math.sin(geocLat)};
	    xyz[0] -= correction[0];
	    xyz[1] -= correction[1];
	    xyz[2] -= correction[2];
	}

	// Obtain spherical topocentric equatorial coordinates
	double[] sph = CoordinateSystem.cartesianToSpherical(xyz);
	double ra = sph[0], dec = sph[1], dist = pos[2] * sph[2];

	// Hour angle
	double angh = lst - ra;
	
	// Correct the equatorial position by diurnal aberration (< 0.3")
	if (!geocentric) {
	    double rotRate = Constant.SIDEREAL_DAY_LENGTH * Constant.TWO_PI / Constant.SECONDS_PER_DAY; // rad/s
	    double factor = rotRate * (eradius * 1000.0) / Constant.SPEED_OF_LIGHT; // v/c
	    double ddec = factor * Math.cos(geocLat) * Math.sin(dec) * Math.sin(angh);
	    if (Math.cos(dec) != 0.0) ra += factor * Math.cos(angh) * Math.cos(geocLat) / Math.cos(dec);
	    dec += ddec;
	}
	
	// Obtain azimuth and geometric alt
	double sinDec = Math.sin(dec), cosDec = Math.cos(dec);
	double sinLat = Math.sin(obsLat); 
	double cosLat = Math.cos(obsLat); 
	double h = sinLat * sinDec + cosLat * cosDec * Math.cos(angh);
	double alt = Math.asin(h);
	double azx = Math.cos(angh) * sinLat - sinDec * cosLat / cosDec;
	double azi = Math.PI + Math.atan2(Math.sin(angh), azx); // 0 = north

	if (geocentric) return new EphemData(azi, alt, -1, -1, -1, -1, Util.normalizeRadians(ra), dec, dist, pos[0], pos[1], pos[3]);

	// Get apparent elevation
	alt = getApparentElevation(alt);

	double tmp = twilight.elevation * Constant.DEG_TO_RAD;
	// Consider the angular radius (pos[3]) for rise, set, transit times when using the HORIZON_34arcmin twilight. 
	// Removing angular radius here would do calculations for the center of the disk instead of the lower/upper limb.
	if (twilight == TWILIGHT.HORIZON_34arcmin) tmp = tmp - pos[3];

	// Compute cosine of hour angle
	tmp = (Math.sin(tmp) - sinLat * sinDec) / (cosLat * cosDec);

	// Make calculations for the meridian
	double transit_alt = Math.asin(sinDec * sinLat + cosDec * cosLat);
	transit_alt = getApparentElevation(transit_alt);

	// Obtain the current transit event in time
	double transit = getTwilightEvent(ra, 0);

	// Make calculations for rise and set
	double rise = -1, set = -1;
	if (Math.abs(tmp) <= 1.0) {
	    double ang_hor = Math.abs(Math.acos(tmp));
	    rise = getTwilightEvent(ra, -ang_hor);
	    set = getTwilightEvent(ra, ang_hor);
	}

	EphemData out = new EphemData(azi, alt, rise, set, transit, transit_alt, 
		Util.normalizeRadians(ra), dec, dist, pos[0], pos[1], pos[3]);
	return out;
    }

    private double getTwilightEvent(double ra, double angh) {
	double celestialHoursToEarthTime = 1.0 / (Constant.SIDEREAL_DAY_LENGTH * Constant.TWO_PI);
	double jdToday_UT = Math.floor(jd_UT - 0.5) + 0.5;

	double eventTime = celestialHoursToEarthTime * Util.normalizeRadians(ra + angh - lst);
	double eventTimePrev = celestialHoursToEarthTime * (Util.normalizeRadians(ra + angh - lst) - Constant.TWO_PI);
	double eventDatePrev_UT = Math.floor(jd_UT + eventTimePrev - 0.5) + 0.5;

	if (Math.abs(eventTimePrev) < Math.abs(eventTime) && twilightMode == TWILIGHT_MODE.CLOSEST) eventTime = eventTimePrev;
	if (twilightMode == TWILIGHT_MODE.TODAY_UT) {
	    double eventDate_UT = Math.floor(jd_UT + eventTime - 0.5) + 0.5;
	    if (jdToday_UT != eventDate_UT) eventTime = -jd_UT - 1;
	    if (jdToday_UT == eventDatePrev_UT) eventTime = eventTimePrev;
	}
	if (twilightMode == TWILIGHT_MODE.TODAY_LT) {
	    double tz = timeZone / 24.0, jdToday_LT = Math.floor(jd_UT + tz - 0.5) + 0.5;
	    double eventDate_LT = Math.floor(jd_UT + tz + eventTime - 0.5) + 0.5;
	    if (jdToday_LT != eventDate_LT) eventTime = -jd_UT - 1;

	    double eventDatePrev_LT = Math.floor(jd_UT + tz + eventTimePrev - 0.5) + 0.5;
	    if (jdToday_LT == eventDatePrev_LT) eventTime = eventTimePrev;

	    double eventTimeNext = celestialHoursToEarthTime * (Util.normalizeRadians(ra + angh - lst) + Constant.TWO_PI);
	    double eventDateNext_LT = Math.floor(jd_UT + tz + eventTimeNext - 0.5) + 0.5;
	    if (jdToday_LT == eventDateNext_LT) eventTime = eventTimeNext;
	}

	return jd_UT + eventTime;
    }

    /** Corrects input geometric elevation for refraction if it is greater than -3 degrees, returning the apparent elevation */
    private double getApparentElevation(double alt) {
	if (alt <= -3 * Constant.DEG_TO_RAD) return alt;

	double altIn = alt, prevAlt = alt;
	int niter = 0;
	do {
	    double altOut = getGeometricElevation(alt);
	    alt = altIn - (altOut-alt);
	    niter ++;
	    if (Math.abs(prevAlt-alt) < 0.001 * Constant.DEG_TO_RAD) break;
	    prevAlt = alt;
	} while (niter < 8);

	return alt;
    }

    /** Compute geometric elevation from apparent elevation. Note ephemerides calculates geometric elevation, so an inversion is 
     * required, something achieved in method {@linkplain #getApparentElevation(double)} by iteration */
    private double getGeometricElevation(double alt) {
	double ps = 1010; // Pressure in mb
	double ts = 10 + 273.15; // Temperature in K
	double altDeg = alt * Constant.RAD_TO_DEG;

	// Bennet 1982 formulae for optical wavelengths, do the job but not accurate close to horizon
	double r = Constant.DEG_TO_RAD * Math.abs(Math.tan(Constant.PI_OVER_TWO - (altDeg + 7.31 / (altDeg + 4.4)) * Constant.DEG_TO_RAD)) / 60.0;
	double refr = r * (0.28 * ps / ts);
	return Math.min(alt - refr, Constant.PI_OVER_TWO);

	/*
	// Bennet formulae adapted to radio wavelenths. Use this for position in radio wavelengths
	// Reference for some values: http://icts-yebes.oan.es/reports/doc/IT-OAN-2003-2.pdf (Yebes 40m radiotelescope)
	double hs = 20; // Humidity %
	// Water vapor saturation pressure following Crane (1976), as in the ALMA memorandum
	double esat = 6.105 * Math.exp(25.22 * (ts - 273.15) / ts) + Math.pow(ts / 273.15, -5.31);
	double Pw = hs * esat / 100.0;

	double R0 = (16.01 / ts) * (ps - 0.072 * Pw + 4831 * Pw / ts) * Constant.ARCSEC_TO_RAD;
	double refr2 = R0 * Math.abs(Math.tan(Constant.PI_OVER_TWO - (altDeg + 5.9 / (altDeg + 2.5)) * Constant.DEG_TO_RAD));
	return Math.min(alt - refr, Constant.PI_OVER_TWO);
	*/
    }
    
    /**
     * Computes an accurate rise/set/transit time for a moving object
     * @param riseSetJD Start date for the event
     * @param index Event identifier
     * @param niter Maximum number of iterations
     * @return The Julian day in UT for the event, 1s accuracy
     */
    protected double obtainAccurateRiseSetTransit(double riseSetJD, EVENT index, int niter) {
	double step = -1;
	for (int i = 0; i< niter; i++) {
	    if (riseSetJD == -1) return riseSetJD; // -1 means no rise/set from that location
	    setUTDate(riseSetJD);
	    EphemData out = doCalc(getBodyPosition(), false);

	    double val = out.rise;
	    if (index == EVENT.SET) val = out.set;
	    if (index == EVENT.TRANSIT) val = out.transit;
	    step = Math.abs(riseSetJD - val);
	    riseSetJD = val;
	    if (step <= 1.0 / Constant.SECONDS_PER_DAY) break; // convergency reached
	}
	if (step > 1.0 / Constant.SECONDS_PER_DAY) return -1; // did not converge => without rise/set/transit in this date
	return riseSetJD;
    }
    
    protected double[] getBodyPosition() {
	return null;
    }
    
    /**
     * Computes the ephemerides for a body, including accurate rise/set/transit times
     * @param bodyRed The body reduction data
     * @return The ephemerides data
     */
    public static EphemData getEphemeris(EphemReduction bodyRed) {
	double jd_UT = bodyRed.jd_UT;
	EphemData bodyData = bodyRed.doCalc(bodyRed.getBodyPosition(), false);
	int niter = 15; // Number of iterations to get accurate rise/set/transit times
	bodyData.rise = bodyRed.obtainAccurateRiseSetTransit(bodyData.rise, EVENT.RISE, niter);
	bodyData.set = bodyRed.obtainAccurateRiseSetTransit(bodyData.set, EVENT.SET, niter);
	bodyData.transit = bodyRed.obtainAccurateRiseSetTransit(bodyData.transit, EVENT.TRANSIT, niter);
	if (bodyData.transit == -1) {
	    bodyData.transitElevation = 0;
	} else {
	    // Update Sun's maximum elevation
	    bodyRed.setUTDate(bodyData.transit);
	    bodyData.transitElevation = bodyRed.doCalc(bodyRed.getBodyPosition(), false).transitElevation;
	}
	bodyRed.setUTDate(jd_UT);
	return bodyData;
    }
}
