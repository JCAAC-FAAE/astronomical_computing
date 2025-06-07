package journal;

/**
 * A class to compute the ephemerides of the Sun. This class uses the code inside {@linkplain EphemReduction} by inheritance.
 */
public class EphemSunSimple extends EphemReduction {

    public EphemSunSimple(double jd_utc, double lon, double lat, double alt, TWILIGHT tw, TWILIGHT_MODE twm, int tz) {
	super(jd_utc, lon, lat, alt, tw, twm, tz);
    }

    @Override
    public double[] getBodyPosition() {
	double t = EarthAngles.toCenturiesRespectJ2000(jd_UT, true);

	// Mean longitude of Sun corrected for precession, taken from Moshier. 
	// Final accuracy up to +/- 27" or better over many millenia compared to VSOP87
	double xe = (129597742.283429 * t + 361679.198) + (-5.23e-6 * t - 2.04411e-2) * t * t;
	xe += (((((((((-8.66e-20 * t - 4.759e-17) * t + 2.424e-15) * t + 1.3095e-12) * t + 1.7451e-10) * t - 1.8055e-8) * t - 0.0000235316) * t + 0.000076) * t + 1.105414) * t + 5028.791959) * t;
	double lon = Util.normalizeRadians(Math.PI + Constant.ARCSEC_TO_RAD * xe) * Constant.RAD_TO_DEG;

	/* Mean anomaly of sun = l' (J. Laskar) */
	double x = (1.2959658102304320e+08 * t + 1.2871027407441526e+06);
	x += ((((((((1.62e-20 * t - 1.0390e-17) * t - 3.83508e-15) * t + 4.237343e-13) * t + 8.8555011e-11) * t - 4.77258489e-8) * t - 1.1297037031e-5) * t + 8.7473717367324703e-05) * t - 5.5281306421783094e-01) * t * t;
	double sanomaly = Util.normalizeRadians(Constant.ARCSEC_TO_RAD * x);

	double c = (1.9146 - .004817 * t - .000014 * t * t) * Math.sin(sanomaly);
	c = c + (.019993 - .000101 * t) * Math.sin(2 * sanomaly);
	c = c + .00029 * Math.sin(3.0 * sanomaly); // Correction to the mean ecliptic longitude

	// Now compute approximate aberration
	double d = -.00569;

	double slongitude = lon + c + d; // apparent longitude (error<0.003 deg)
	double slatitude = 0; // Sun's ecliptic latitude is always negligible
	double ecc = .016708617 - 4.2037E-05 * t - 1.236E-07 * t * t; // Eccentricity
	double v = sanomaly + c * Constant.DEG_TO_RAD; // True anomaly
	double sdistance = 1.000001018 * (1.0 - ecc * ecc) / (1.0 + ecc * Math.cos(v)); // In UA

	return new double[] {slongitude * Constant.DEG_TO_RAD, slatitude * Constant.DEG_TO_RAD, sdistance, Math.atan(696000.0 / (sdistance * Constant.AU))};
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
	EphemSunSimple sunEph = new EphemSunSimple(jd_utc, lon, lat, alt, tw, twm, tz);
	EphemData sunData = EphemReduction.getEphemeris(sunEph);

	// Report
	System.out.println("Sun");
	System.out.println(sunData.toString());
	
/*
Sun
 Az:       285.78873°
 El:       17.424834°
 Dist:     1.015206 au
 RA:       78.41086°
 DEC:      23.007696°
 Ill:      100.0%
 ang.R:    0.26256925°
 Rise:     2020/06/09 04:46:57 UT
 Set:      2020/06/09 19:43:59 UT
 Transit:  2020/06/09 12:15:22 UT (elev. 72.99493°)
 */
    }
}
