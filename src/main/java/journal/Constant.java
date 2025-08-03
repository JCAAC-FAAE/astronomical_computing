package journal;

/**
 * Useful constants
 * @version 1.1 Summer 2025: very little change to AU contant to make it consistent with DE440, no practical effects
 */
public class Constant {
    
	/** Radians to degrees. */
	public static final double RAD_TO_DEG = 180.0 / Math.PI;

	/** Degrees to radians. */
	public static final double DEG_TO_RAD = 1.0 / RAD_TO_DEG;

	/* Arcseconds to radians */
	public static final double ARCSEC_TO_RAD = (DEG_TO_RAD / 3600.0);

	/** Astronomical Unit in km. As defined by JPL in DE440. */
	public static final double AU = 149597870.7;

	/** Earth equatorial radius in km. IERS 2003 Conventions. */
	public static final double EARTH_RADIUS = 6378.1366;

	/** Two times Pi. */
	public static final double TWO_PI = 2.0 * Math.PI;

	/** Pi divided by two. */
	public static final double PI_OVER_TWO = Math.PI / 2.0;

	/** Julian century conversion constant = 100 * days per year. */
	public static final double JULIAN_DAYS_PER_CENTURY = 36525.0;

	/** Seconds in one day. */
	public static final double SECONDS_PER_DAY = 86400;

	/** Light time in days for 1 AU. DE405 definition. */
	public static final double LIGHT_TIME_DAYS_PER_AU = 0.00577551833109;

	/** Our default epoch. The Julian Day which represents noon on 2000-01-01. */
	public static final double J2000 = 2451545.0;
	
	/** Speed of light in m/s, exact as it is defined. */
	public static final double SPEED_OF_LIGHT = 299792458.0;
	
	/** Length of a sidereal day in days according to IERS Conventions. */
	public static final double SIDEREAL_DAY_LENGTH = 1.00273781191135448;
}
