package journal;

/**
 * Obtain times of solar and lunar eclipses using a geometric 'brute-force' method independent from the ephemerides theory. 
 * An advantage is the possibility of calculating eclipses produced by other satellites (not the Moon).
 */
public class Eclipses {

    /** The set of eclipse types for both solar and lunar eclipses. */
    public enum ECLIPSE_TYPE {
	/** Constant ID for a total solar/lunar eclipse. */
	TOTAL,
	/** Constant ID for a partial solar/lunar eclipse. */
	PARTIAL,
	/** Constant ID for an annular solar eclipse. */
	ANNULAR,
	/** Identifier for a penumbral lunar eclipse. */
	PENUMBRAL,
	/** Constant ID for an inexistent solar eclipse, used internally. */
	NO_ECLIPSE;	
    }

    /** Instance variables */
    private ECLIPSE_TYPE type = ECLIPSE_TYPE.NO_ECLIPSE;
    private EphemMoon moonEph;
    private EphemSun sunEph;
    private EphemData moonData, sunData;
    private double step = 0.01 / Constant.SECONDS_PER_DAY;
    private double jdMax;
    private double[] events;
    private double maxMag = 0;
    private boolean isSolar;
    private boolean earthMoon = true;
    private boolean computeBackEvent = false;
    private boolean considerRiseSet = false;

    /**
     * Sets the desired accuracy of the iterative method in seconds for eclipses.
     * @param s Accuracy in seconds, must be greater than 0 and lower than 10. Default is 0.01s
     */
    public void setAccuracy(double s)  {
	if (s > 0 && s < 10) step = s / Constant.SECONDS_PER_DAY;
    }

    /**
     * Selects if the times returned consider if the Sun/Moon is above the horizon or not.
     * @param b Selection. Default is false (geometric events)
     */
    public void setConsiderRiseSet(boolean b) {
	considerRiseSet = b;
    }

    private boolean[] checkSolarEclipse(double jd_utc) {
	double mr = moonData.angularRadius;
	double sr = sunData.angularRadius;
	double dist = EphemData.getAngularDistance(moonData, sunData);
	if (dist > (sr + mr)) return new boolean[] { false, false }; // Not eclipsed
	if (considerRiseSet && sunData.elevation < -sr * 0.8) return new boolean[] { false, false }; // Not visible
	if (moonData.distance > sunData.distance) {
	    if (!computeBackEvent) return new boolean[] {false, false}; // No eclipse since the Moon is farther
	} else {
	    if (computeBackEvent) return new boolean[] {false, false};
	}

	boolean totality = false;
	double mag = 0;
	if ((dist + mr) <= sr || (dist + sr) <= mr) {
	    totality = true;
	    if (mr > sr) {
		type = ECLIPSE_TYPE.TOTAL;
	    } else {
		type = ECLIPSE_TYPE.ANNULAR;
	    }
	    double add = mr - (dist + sr);
	    mag = 1.0 + add / (sr * 2);
	} else {
	    if (type == ECLIPSE_TYPE.NO_ECLIPSE) type = ECLIPSE_TYPE.PARTIAL;
	    double add = mr - (dist - sr);
	    mag = add / (sr * 2);
	}
	if (mag > maxMag) {
	    maxMag = mag;
	    jdMax = jd_utc;
	}

	return new boolean[] { true, totality };
    }

    private boolean[] checkLunarEclipse(double jd_utc) {
	boolean insideShadow = false, totality = false;
	boolean insidePenumbra = false, totalityPenumbra = false;
	if (considerRiseSet && moonData.elevation < -moonData.angularRadius * 0.8) return new boolean[] { insidePenumbra, totalityPenumbra, insideShadow, totality }; // Not visible
	
	// Get distance between shadow cone direction and the Moon. Note ephemerides should be geocentric in this case
	double dist = EphemData.getAngularDistance(moonData, sunData);
	double equRadius = 6378.14;
	double polRadius = 6356.75;

	// The main calculation is to position the center of the Earth shadow cone. We consider this cone to be indeed an oval with a size
	// slightly larger than the Earth's equatorial and polar radius. This excess can be understood taking into account Earth surface
	// elevation and opacity of the atmosphere. Values fitted to lunar eclipses in 2007. Note AA supplement uses 1.02 in page 429, but more
	// precision is required (and a separation between polar and equatorial axis) to get an accuracy of few seconds.
	// See http://eclipse.gsfc.nasa.gov/LEcat5/shadow.html for a discussion
	double opaqueAtmFactorEqu = 1.0131;
	double opaqueAtmFactorPole = 1.015;
	double earthShadowConeSize = equRadius / (Constant.AU * Math.tan(sunData.angularRadius));
	double angRadiusMax = Math.atan2(equRadius / Constant.AU, moonData.distance) * (opaqueAtmFactorEqu - (moonData.distance / earthShadowConeSize));
	double angRadiusMin = Math.atan2(polRadius / Constant.AU, moonData.distance) * (opaqueAtmFactorPole - (moonData.distance / earthShadowConeSize));
	double penumbraAngRadius = 2.0 * sunData.angularRadius;
	double penumbraScaleMax = angRadiusMax + penumbraAngRadius;
	double penumbraScaleMin = angRadiusMin + penumbraAngRadius;

	double pa = 3.0 * Constant.PI_OVER_TWO - EphemData.getPositionAngle(moonData, sunData);
	double sx = Math.sin(pa) / angRadiusMax;
	double sy = Math.cos(pa) / angRadiusMin;
	double sr = 1.0 / Math.hypot(sx, sy);
	double srUmbra = sr;

	if (dist <= (sr + moonData.angularRadius)) insideShadow = true;
	if (dist <= (sr - moonData.angularRadius)) totality = true;

	sx = Math.sin(pa) / penumbraScaleMax;
	sy = Math.cos(pa) / penumbraScaleMin;
	sr = 1.0 / Math.hypot(sx, sy);

	if (dist <= (sr + moonData.angularRadius)) insidePenumbra = true;
	if (dist <= (sr - moonData.angularRadius)) totalityPenumbra = true;

	if ((insidePenumbra || totalityPenumbra) && type.ordinal() > ECLIPSE_TYPE.PENUMBRAL.ordinal()) type = ECLIPSE_TYPE.PENUMBRAL;
	if (insideShadow && type.ordinal() > ECLIPSE_TYPE.PARTIAL.ordinal()) {
	    type = ECLIPSE_TYPE.PARTIAL;
	    maxMag = 0; // Reset the magnitude, that was penumbral until this point
	}
	if (totality && type.ordinal() > ECLIPSE_TYPE.TOTAL.ordinal()) type = ECLIPSE_TYPE.TOTAL;

	double mag;
	if (type == ECLIPSE_TYPE.PENUMBRAL) {
	    double add = sr + moonData.angularRadius - dist;
	    mag = add / (moonData.angularRadius * 2);
	} else {
	    double add = srUmbra - moonData.angularRadius - dist;
	    mag = 1.0 + add / (moonData.angularRadius * 2);
	}
	if (mag > maxMag) {
	    maxMag = mag;
	    jdMax = jd_utc;
	}
	return new boolean[] { insidePenumbra, totalityPenumbra, insideShadow, totality };
    }

    /**
     * Obtain events for the next solar or lunar eclipses in TDB. Eclipse type (total, partial, annular) 
     * and events are set in static variables. Moon/Sun elevation above local horizon is not considered 
     * in this method, so output events could be not visible from the position of the observer.
     * @param jd_utc Input time in UTC before the eclipse
     * @param lon Longitude of the observer in degrees
     * @param lat Latitude of the observer
     * @param alt Altitude of the observer in m
     */
    public Eclipses(double jd_utc, double lon, double lat, double alt) {
	moonEph = new EphemMoon(jd_utc, lon, lat, alt, EphemReduction.TWILIGHT.HORIZON_34arcmin, EphemReduction.TWILIGHT_MODE.TODAY_LT, 0);
	sunEph = new EphemSun(jd_utc, lon, lat, alt, EphemReduction.TWILIGHT.HORIZON_34arcmin, EphemReduction.TWILIGHT_MODE.TODAY_LT, 0);
    }

    /**
     * Returns the next solar eclipse
     * @param jd_utc Initial date for the search
     */
    public void nextSolarEclipse(double jd_utc) {
	compute(jd_utc, true);
    }

    /**
     * Returns the next lunar eclipse
     * @param jd_utc Initial date for the search
     */
    public void nextLunarEclipse(double jd_utc) {
	compute(jd_utc, false);
    }

    private double stepForward(double jd_utc, boolean solar) {
	jd_utc += step; // Time step of the search
	moonEph.setUTDate(jd_utc);
	sunEph.setUTDate(jd_utc);
	moonData = moonEph.doCalc(moonEph.getBodyPosition(), !solar); 
	sunData = sunEph.doCalc(sunEph.getBodyPosition(), !solar);
	if (solar) {
	    // Correction by Xavier Jubier in the solar radius for solar eclipses
	    if (earthMoon) sunData.angularRadius *= 696250.0 / 696000.; 
	} else {
	    // For lunar eclipses, use as sun position the Earth's cone shadow (opposite of geocentric Sun)
	    sunData.rightAscension += Math.PI;
	    sunData.declination *= -1;
	}
	return jd_utc;
    }

    private double getTimeIncrementForNextPossibleEclipse(double jd_utc, boolean solar) {
	if (!earthMoon) return 0;
	// Accelerate calculations from the input date to one close to the next possible solar/lunar eclipse
	jd_utc = stepForward(jd_utc, solar);
	double dist = moonData.eclipticLongitude - sunData.eclipticLongitude;
	if (!solar) dist += Math.PI;
	return (1.0 - Util.normalizeRadians(dist) / Constant.TWO_PI) * 29.53 - 1;	
    }
    
    private void compute(double jd_utc, boolean solar) {
	// 4 events for solar eclipses: start/end of partial/total
	// 8 events for lunar eclipses: penumbral partial/penumbral total/partial/total
	int outSize = solar ? 4 : 8; 
	double[] out = new double[outSize];

	// Prepare and execute the iteration
	jd_utc += getTimeIncrementForNextPossibleEclipse(jd_utc, solar);	
	type = ECLIPSE_TYPE.NO_ECLIPSE;
	double jdPrev = -1, statusPrev = -1, outPrev[] = null;
	boolean detailedMode = false;
	double jd0 = jd_utc;
	do {
	    jd_utc = stepForward(jd_utc, solar);
	    boolean[] event = solar ? checkSolarEclipse(jd_utc) : checkLunarEclipse(jd_utc);

	    for (int i=0; i<outSize;i++) {
		if (i < outSize / 2) { // Start of events
		    if (event[i] && out[i] == 0.0) {
			out[i] = jd_utc;
			if (i == 0) break; // if the eclipse starts, do not let a second event to happen at this point
		    }
		} else { // End of events
		    int iStart = outSize - 1 - i;
		    if (!event[iStart] && out[i] == 0.0 && out[iStart] != 0.0) out[i] = jd_utc; // End of events
		}
	    }

	    // Generate an status id based on current events
	    double status = 0;
	    for (int i=0; i<out.length; i++) {
		if (out[i] != 0) status += Math.pow(10, i);
	    }

	    // If status changed then go to previous iteration in detailed mode, to avoid losing precision
	    if (statusPrev >= 0 && status != statusPrev && !detailedMode) {
		out = outPrev;
		jd_utc = jdPrev;
		detailedMode = true;
		continue;
	    }
	    // Still in detailed mode ? Then continue iterating
	    if (detailedMode && status == statusPrev) continue;

	    // Disable detailed mode and go forward in time in longer steps to search faster
	    detailedMode = false;
	    jdPrev = jd_utc;
	    statusPrev = status;
	    outPrev = out.clone();
	    jd_utc += 10 / Constant.SECONDS_PER_DAY; // 10s to search faster for the next event
	    if (earthMoon && out[0] == 0 && jd_utc - jd0 > 2) { // If no event found in 2 days, try in the next lunation
		jd_utc += getTimeIncrementForNextPossibleEclipse(jd_utc, solar);	
		jd0 = jd_utc;
	    }
	} while (out[outSize-1] == 0.0);

	// Re-estimate eclipse maximum accurately for the case of partial solar eclipses (*)
	if (solar && earthMoon && type == ECLIPSE_TYPE.PARTIAL) {
	    jd_utc = jdMax - 10 / Constant.SECONDS_PER_DAY;
	    double jdEnd = jdMax + 10 / Constant.SECONDS_PER_DAY;
	    while (jd_utc < jdEnd) {
		jd_utc = stepForward(jd_utc, solar);
		boolean[] event = checkSolarEclipse(jd_utc);

		// In critical cases the eclipse type may change from partial to total here, for a few seconds. 
		// Ex. Lezama, for the eclipse on August 2026 
		if (event[0] && out[0] == 0.0)
		    out[0] = jd_utc;
		if (event[1] && out[1] == 0.0)
		    out[1] = jd_utc;

		if (!event[1] && out[2] == 0.0 && out[1] != 0.0)
		    out[2] = jd_utc;
		if (!event[0] && out[3] == 0.0 && out[0] != 0.0)
		    out[3] = jd_utc;
	    }
	}

	events = out;
	isSolar = solar;
    }

    /**
     * Returns type of the current calculated eclipse.
     * @return Eclipse type, such as annular, partial, or total.
     */
    public ECLIPSE_TYPE getEclipseType() {
	return type;
    }

    /**
     * Returns the set of events. Times are given in UTC.
     * @return Events as an array, 4 for solar eclipses, 8 for lunar ones.
     */
    public double[] getEvents() {
	return events;
    }

    /**
     * Returns the name of a given event.
     * @param events Events computed.
     * @param index Index of the event, 0-3 for solar eclipses, 0-7 for lunar ones.
     * @return Event name.
     */
    public String getEventName(double[] events, int index) {
	// start/end of partial/total for solar eclipses, or penumbral partial/penumbral total/partial/total for lunar ones	
	if (events.length == 4) { // Solar eclipses
	    return "Partial phase start,Total phase start,Total phase end,Partial phase end".split(",")[index];
	} else { // Lunar eclipses
	    return ("Penumbral phase start,Total penumbral phase start,Partial phase start,Total phase start," + 
		    "Total phase end,Partial phase end,Total penumbral phase end,Penumbral phase end").split(",")[index];
	}
    }

    /**
     * Return eclipse maximum as Julian day.
     * @return Eclipse maximum in TDB.
     */
    public double getEclipseMaximum() {
	return jdMax;
    }

    /**
     * Returns the eclipse magnitude.
     * @return Eclipse magnitude.
     */
    public float getMagnitude() {
	return (float) maxMag;
    }

    /**
     * Returns if the computed eclipse is solar or not (lunar)
     * @return If it is a solar eclipse
     */
    public boolean isSolarEclipse() {
	return isSolar;
    }
    
    /**
     * Test program
     * @param args Not used
     */
    public static void main(String[] args) {
	// Prepare input data
	int year = 2026, month = 1, day = 1, h = 0, m = 0, s = 0;
	JulianDay jday = new JulianDay(year, month, day);
	jday.setDayFraction((h + m / 60.0 + s / 3600.0) / 24.0);

	double jd_utc = jday.getJulianDay();
	String locName = "Vitoria-Gasteiz";
	double lon = -2.67324; // degrees
	double lat = 42.85319;
	double alt = 516; // m

	// Compute the next eclipses during 3 years
	double jdEnd = jd_utc + 3 * 365.25;
	double eclTime = jd_utc;
	while (eclTime < jdEnd) {
	    Eclipses ecl = new Eclipses(eclTime, lon, lat, alt);
	    ecl.setConsiderRiseSet(true);
	    ecl.nextSolarEclipse(eclTime);
	    eclTime = ecl.getEclipseMaximum();
	    formatEclipse(ecl, locName);
	    eclTime += 20;
	}
	eclTime = jd_utc;
	while (eclTime < jdEnd) {
	    Eclipses ecl = new Eclipses(eclTime, lon, lat, alt);
	    ecl.setConsiderRiseSet(true);
	    ecl.nextLunarEclipse(eclTime);
	    eclTime = ecl.getEclipseMaximum();
	    formatEclipse(ecl, locName);
	    eclTime += 20;
	}
	System.out.println("TT-UT1 assumed: "+EarthAngles.TTminusUT1(eclTime)+" s");
    }

    /**
     * Formats an eclipse to show the information in the console.
     * @param ecl The eclipse instance with an eclipse computed.
     * @param locName Location where the calculations where done for.
     */
    public static void formatEclipse(Eclipses ecl, String locName) {
	ECLIPSE_TYPE type = ecl.getEclipseType();
	double[] times = ecl.getEvents();

	JulianDay jd = new JulianDay(ecl.getEclipseMaximum());
	String eclType = ecl.isSolarEclipse() ? "solar" : "lunar";
	System.out.println("Found a "+type.name()+" "+eclType+" eclipse on "+jd.toString()+" from "+locName);
	System.out.println("Eclipse magnitude: " + Util.formatValue(ecl.getMagnitude(), 4));
	for (int i=0; i<times.length; i++) {
	    if (times[i] <= 0) continue;
	    System.out.println(ecl.getEventName(times, i) + ": " + new JulianDay(times[i]).toString());
	}
	System.out.println();
    }
}
