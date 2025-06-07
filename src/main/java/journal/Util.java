package journal;

import java.text.DecimalFormat;

/**
 * Utility functions
 * @version April, 2025 (added method getLineSeparator(), to return the line break string of the system)
 * @version September, 2024 (first version)
 */
public class Util {

    /**
     * Formats an integer with two digits with an optional suffix
     * @param i The integer
     * @param suffix An optional suffix
     * @return The formatted string
     */
    public static String fmt02(int i, String suffix) {
	String out = "" + i;
	if (i >= 0 && i < 10) out = "0" + out;
	if (suffix != null) out = out + suffix;
	return out;
    }
    
    /**
     * Reduce an angle in degrees to the range (0 &lt;= &lt; 360).
     *
     * @param d Value in degrees.
     * @return The reduced degree value.
     */
    public static double normalizeDegrees(double d) {
	d -= 360.0 * Math.floor(d / 360.0);
	if (d < 0.) d += 360.0;
	return d;
    }

    /**
     * Reduce an angle in radians to the range (0 - 2 Pi).
     *
     * @param r Value in radians.
     * @return The reduced radian value.
     */
    public static double normalizeRadians(double r) {
	double d = r - Constant.TWO_PI * Math.floor(r / Constant.TWO_PI);
	if (d < 0.) d += Constant.TWO_PI;
	return d;
    }

    /**
     * Format a real number.
     *
     * @param val Numerical value
     * @param decimals Number of decimal places.
     * @return String with the adequate format.
     */
    public static String formatValue(double val, int decimals) {
	DecimalFormat formatter;
	if (decimals > 0) {
	    String out = repeatString("0", decimals);
	    formatter = new DecimalFormat("##0." + out);
	} else {
	    formatter = new DecimalFormat("##0");
	}

	String out = formatter.format(val);
	out = out.replace(",", ".");

	return out;
    }
    
    /**
     * Repeats a given string s n times.
     * @param s The string.
     * @param n Number of times to repeat.
     * @return n times the string s.
     */
    public static String repeatString(String s, int n) {
        final StringBuilder sb = new StringBuilder(s);
        for(int i = 1; i < n; i++) {
            sb.append(s);
        }
        return sb.toString();

    }
    
    /**
     * Format declination.
     *
     * @param dec Declination in radians. Must be in the range -Pi/2 to +Pi/2.
     * @param nsec Number of decimal places in arcseconds.
     * @return String with the format $##&deg; ##' ##.##...'' ($ is the sign).
     */
    public static String formatDEC(double dec, int nsec) {
	DecimalFormat formatter, formatter0;
	if (nsec > 0) {
	    String decimal = repeatString("0", nsec);
	    formatter = new DecimalFormat("00." + decimal);
	    formatter0 = new DecimalFormat("00");			
	} else {
	    formatter = new DecimalFormat("00");
	    formatter0 = formatter;
	}
	double dec_d = Math.abs(dec) * Constant.RAD_TO_DEG;
	double dec_m = (dec_d - Math.floor(dec_d)) * 60.0;
	double dec_s = (dec_m - Math.floor(dec_m)) * 60.0;
	dec_d = Math.floor(dec_d);
	dec_m = Math.floor(dec_m);

	String out = "" + formatter0.format(dec_d) + "\u00b0 " + formatter0.format(dec_m) + "' " + formatter.format(dec_s) + "\"";
	out = out.replaceAll(",", ".");
	if (dec < 0.0) out = "-" + out;

	return out;
    }

    /**
     * Format right ascension. Significant digits are adapted to common
     * ephemeris precision.
     *
     * @param ra Right ascension in radians.
     * @param nsec Number of decimal places in seconds of time.
     * @return String with the format ##h ##m ##.##...s.
     */
    public static String formatRA(double ra, int nsec) {
	DecimalFormat formatter;
	if (nsec > 0) {
	    String dec = repeatString("0", nsec);
	    formatter = new DecimalFormat("00." + dec);
	} else {
	    formatter = new DecimalFormat("00");
	}
	DecimalFormat formatter0 = new DecimalFormat("00");
	ra = normalizeRadians(ra);
	double ra_h = ra * Constant.RAD_TO_DEG / 15.0;
	double ra_m = (ra_h - Math.floor(ra_h)) * 60.0;
	double ra_s = (ra_m - Math.floor(ra_m)) * 60.0;
	ra_h = Math.floor(ra_h);
	ra_m = Math.floor(ra_m);

	String out = "" + formatter0.format(ra_h) + "h " + formatter0.format(ra_m) + "m " + formatter.format(ra_s) + "s";
	out = out.replaceAll(",", ".");

	return out;
    }
    
    /**
     * Get the system line separator.
     * @return Line separator.
     */
    public static String getLineSeparator() {
	return System.getProperty("line.separator");
    }

    /**
     * Returns declination in radians given degrees, minutes, and arcseconds. A
     * minus sign can be set in degrees for southern positions.
     *
     * @param decg Degrees.
     * @param min Arcminutes.
     * @param sec Arcseconds
     * @return Declination in radians
     */
    public static double parseDeclination(String decg, double min, double sec) {
	double g = Double.parseDouble(decg);
	double dec = Math.abs(g) + min / 60.0 + sec / 3600.0;
	dec = dec * Constant.DEG_TO_RAD;
	if (decg.indexOf("-") >= 0) dec = -dec;
	return dec;
    }

    /**
     * Returns right ascension in radians given hours, minutes, and seconds of
     * time.
     *
     * @param hour Hours.
     * @param min Minutes.
     * @param sec Seconds.
     * @return Right ascension value in radians.
     */
    public static double parseRightAscension(double hour, double min, double sec) {
	double ra = Math.abs(hour) + min / 60.0 + sec / 3600.0;
	ra = ra * 15 * Constant.DEG_TO_RAD;
	return ra;

    }

    /**
     * Separates a RA value into hours, minutes, and seconds of time.
     * @param ra The input right ascension in radians.
     * @return H, M, and S of time.
     */
    public static double[] getHMS(double ra) {
	double[] out = new double[3];
	double rah = ra * Constant.RAD_TO_DEG / 15.0;
	out[0] = (int) rah;
	double ram = (rah - out[0]) * 60.0;
	out[1] = (int) ram;
	out[2] = (ram - out[1]) * 60.0;
	return out;
    }

    /**
     * Separates a DEC value into degrees, minutes, and seconds of arc.
     * @param dec The input declination in radians.
     * @return D, M, and S of arc, and sign in the fourth component. The
     * sign is 1 if DEC is positive, and -1 if it is negative.
     */
    public static double[] getDMSs(double dec) {
	double[] out = new double[4];
	double ded = Math.abs(dec * Constant.RAD_TO_DEG);
	out[0] = (int) ded;
	double dem = (ded - out[0]) * 60.0;
	out[1] = (int) dem;
	out[2] = (dem - out[1]) * 60.0;
	out[3] = 1;
	if (dec < 0.0) out[3] = -1;
	return out;
    }
}
