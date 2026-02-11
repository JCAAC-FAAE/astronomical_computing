package journal;

/**
 * A class to hold the results of the ephemerides calculations.
 */
public class EphemData {

    /** Values for azimuth, elevation, rise, set, and transit for the Sun. Angles in radians, rise ... 
     * as Julian days in UT. Distance in AU */
    public double azimuth, elevation, rise, set, transit, transitElevation, distance, rightAscension, 
    	declination, illuminationPhase, eclipticLongitude, eclipticLatitude, angularRadius;
    
    /**
     * Main constructor
     * @param azi Azimuth
     * @param ele Elevation
     * @param rise2 Rise
     * @param set2 Set
     * @param transit2 Transit
     * @param transit_alt Transit elevation
     * @param ra Right ascension
     * @param dec Declination
     * @param dist Distance
     * @param eclLon Ecliptic longitude
     * @param eclLat Ecliptic latitude
     * @param angR Angular radius
     */
    public EphemData(double azi, double ele, double rise2, double set2,
	    double transit2, double transit_alt, double ra, double dec,
	    double dist, double eclLon, double eclLat, double angR) {
	azimuth = azi;
	elevation = ele;
	rise = rise2;
	set = set2;
	transit = transit2;
	transitElevation = transit_alt;
	rightAscension = ra;
	declination = dec;
	distance = dist;
	illuminationPhase = 100;
	eclipticLongitude = eclLon;
	eclipticLatitude = eclLat;
	angularRadius = angR;
    }
    
    @Override
    public String toString() {
	String degSymbol = "\u00b0";
	String lsep = Util.getLineSeparator();
	StringBuilder sb = new StringBuilder();
	sb.append(" Az:       "+(float) (azimuth * Constant.RAD_TO_DEG)+degSymbol + lsep);
	sb.append(" El:       "+(float) (elevation * Constant.RAD_TO_DEG)+degSymbol + lsep);
	sb.append(" Dist:     "+(float) (distance)+" au" + lsep);
	sb.append(" RA:       "+(float) (rightAscension * Constant.RAD_TO_DEG)+degSymbol + lsep);
	sb.append(" DEC:      "+(float) (declination * Constant.RAD_TO_DEG)+degSymbol + lsep);
	sb.append(" Ill:      "+(float) (illuminationPhase)+"%" + lsep);
	sb.append(" ang.R:    "+(float) (angularRadius * Constant.RAD_TO_DEG)+degSymbol + lsep);
	sb.append(" Rise:     "+EphemData.getDateAsString(rise) + lsep);
	sb.append(" Set:      "+EphemData.getDateAsString(set) + lsep);
	sb.append(" Transit:  "+EphemData.getDateAsString(transit)+" (elev. "+(float) (transitElevation * Constant.RAD_TO_DEG)+degSymbol+")" + lsep);
	return sb.toString();
    }
    
    /**
     * Sets the illumination phase field for the body
     * @param sun The Ephem object for the Sun
     */
    public void setIlluminationPhase(EphemData sun) {
	double dlon = rightAscension - sun.rightAscension;
	double cosElong = (Math.sin(sun.declination) * Math.sin(declination) + 
		Math.cos(sun.declination) * Math.cos(declination) * Math.cos(dlon));

	double rsun = sun.distance;
	double rbody = distance;
	// Use elongation cosine as trick to solve the rectangle and get rp (distance body - sun)
	double rp = Math.sqrt(-(cosElong * 2.0 * rsun * rbody - rsun * rsun - rbody * rbody));

	double dph = ((rp * rp + rbody * rbody - rsun * rsun) / (2.0 * rp * rbody));
	illuminationPhase = 100 * (1.0 + dph) * 0.5;
    }
    
    /**
     * Returns the apparent magnitude assuming the body is an asteroid. See Meeus' Astronomical Algorithms, pages 216 and 217
     * @param sun The ephemeris object for the Sun
     * @param magAbs Absolute magnitude of the body
     * @param magSlope Magnitude slope of the body
     */
    public double getAsteroidApparentMagnitude(EphemData sun, double magAbs, double magSlope) {
	double dlon = rightAscension - sun.rightAscension;
	double cosElong = (Math.sin(sun.declination) * Math.sin(declination) + 
		Math.cos(sun.declination) * Math.cos(declination) * Math.cos(dlon));

	double rsun = sun.distance;
	double rbody = distance;
	// Use elongation cosine as trick to solve the rectangle and get rp (distance body - sun)
	double rp = Math.sqrt(-(cosElong * 2.0 * rsun * rbody - rsun * rsun - rbody * rbody));

	double dph = ((rp * rp + rbody * rbody - rsun * rsun) / (2.0 * rp * rbody));
	illuminationPhase = 100 * (1.0 + dph) * 0.5;
	
	double phaseAngle = Math.acos(dph);
	double tmp = Math.tan(phaseAngle * 0.5);
	double phi1 = Math.exp(-3.33 * Math.pow(tmp, 0.63));
	double phi2 = Math.exp(-1.87 * Math.pow(tmp, 1.22));
	double apMag = magAbs + (float) (5.0 * Math.log10(distance * rp)) - 2.5 * Math.log10(phi1 * (1.0 - magSlope) + phi2 * magSlope);	
	return apMag;
    }
    
    /**
     * Returns a date as a string yyyy/mm/dd hh:mm:ss UT
     * @param jd The Julian day
     * @return The String
     */
    public static String getDateAsString(double jd) {
	if (jd == -1) return "NO RISE/SET/TRANSIT FOR THIS OBSERVER/DATE";

	JulianDay julDay = new JulianDay(jd);
	double[] hms = Util.getHMS(julDay.getDayFraction() * Constant.TWO_PI);
	String out = Util.fmt02(julDay.year, "/");
	out += Util.fmt02(julDay.month, "/");
	out += Util.fmt02(julDay.day, " ");
	out += Util.fmt02((int) hms[0], ":");
	out += Util.fmt02((int) hms[1], ":");
	out += Util.fmt02((int) (hms[2] + 0.5), " UT");
	return out;
    }
}
