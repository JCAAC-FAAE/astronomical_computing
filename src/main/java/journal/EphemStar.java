package journal;

/**
 * A class to compute the ephemerides of a star, or any other catalog body. 
 * This class uses the code inside {@linkplain EphemReduction} by inheritance.
 */
public class EphemStar extends EphemReduction {

    private double ra;
    private double dec;
    
    public EphemStar(double jd_utc, double lon, double lat, double alt, TWILIGHT tw, TWILIGHT_MODE twm, int tz) {
	super(jd_utc, lon, lat, alt, tw, twm, tz);
    }

    /**
     * Sets the J2000 position of the body from catalog coordinates
     * @param ra Right ascension, hours
     * @param dec Declination, degrees
     */
    public void setJ2000Position(double ra, double dec) {
	this.ra = ra * 15.0 * Constant.DEG_TO_RAD;
	this.dec = dec * Constant.DEG_TO_RAD;
    }
    
    /**
     * Transforms coordinates from J2000 equatorial to mean equinox of date
     * @param p Rectangular coordinates
     * @return Mean equatorial rectangular coordinates
     */
    public double[] equatorialJ2000ToMeanEquatorialOfDate(double[] p) {
	double t = EarthAngles.toCenturiesRespectJ2000(jd_UT, true);
	double[] pa = EarthAngles.precessionAnglesFromJ2000(Constant.J2000 + t * Constant.JULIAN_DAYS_PER_CENTURY);
	return CoordinateSystem.rotate(CoordinateSystem.rotate(CoordinateSystem.rotate(p, CoordinateSystem.getRotZ(-pa[1])), CoordinateSystem.getRotY(pa[2])), CoordinateSystem.getRotZ(-pa[0]));
    }
    
    @Override
    public double[] getBodyPosition() {
	double[] cartesianJ2000 = CoordinateSystem.sphericalToCartesian(ra, dec);
	double[] meanEq = equatorialJ2000ToMeanEquatorialOfDate(cartesianJ2000);
	double[] meanEcl = CoordinateSystem.equatorialToEcliptic(meanEq, jd_UT);
	double[] meanEclSph = CoordinateSystem.cartesianToSpherical(meanEcl);
	return new double[] {meanEclSph[0], meanEclSph[1], 1E100, 0}; // Assume infinite distance and angular size 0
    }

    /**
     * Test program
     * @param args Not used
     */
    public static void main(String[] args) {
	// Julian day for current instant
	JulianDay jday = new JulianDay(System.currentTimeMillis());

	// Observer position (Madrid, as in theskylive.com)
	double lon = -(3 + 42 / 60.0 + 9.2 / 3600.0); // degrees
	double lat = 40 + 24 / 60.0 + 59.4 / 3600.0;
	double alt = 0; // m
	int tz = 2; // h

	// Catalog body, J2000 coordinates. Check Vega with https://theskylive.com/sky/stars/vega-alpha-lyrae-star
	String name = "Vega";
	double ra = 18.0 + 36 / 60.0 + 56 / 3600.0;
	double dec = 38 + 47 / 60.0 + 0 / 3600.0;
	
	// Compute the ephemerides data
	EphemStar bodyEph = new EphemStar(jday.getJulianDay(), lon, lat, alt, TWILIGHT.HORIZON_34arcmin, TWILIGHT_MODE.TODAY_UT, tz);
	bodyEph.setJ2000Position(ra, dec);
	EphemData bodyData = EphemReduction.getEphemeris(bodyEph);

	// Report
	System.out.println(name);
	System.out.println(bodyData.toString());
	
	System.out.println(" RA:  " + Util.formatRA(bodyData.rightAscension, 1));
	System.out.println(" DEC: " + Util.formatDEC(bodyData.declination, 0));
    }
}
