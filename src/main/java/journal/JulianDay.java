package journal;

/**
 * A class to store date/time values with a precision of one second or better. Support to
 * retrieve the Julian day number is given for the Gregorian and Julian calendars.
 * <P>
 * Selectable dates are not limited in years, but there are some invalid dates (non
 * existent) in the civil calendar, between October 5, 1582 and October 14, 1582.
 * <P>
 * In the constructors and fields the year is entered without considering that year 0 does 
 * not exist (astronomical year). Year 0 is 1 B.C.
 * @version September, 2024 (first version)
 */
public class JulianDay {

    /** Instance field */
    public int year, month, day, hour, minute;
    /** Instance field */
    double second;

    /**
     * Constructor for a Julian day at 0h
     * @param y Year
     * @param m Month
     * @param d Day
     */
    public JulianDay(int y, int m, int d) {
	year = y;
	month = m;
	day = d;
	hour = 0;
	minute = 0;
	second = 0;
	if (isInvalid()) throw new IllegalArgumentException("Date is invalid");
    }

    /**
     * Constructor for a date given as a long value, representing the number of
     * milliseconds elapsed from 1970, January 1, at 0h UTC
     * @param t Time from 1970-1-1 in milliseconds
     */
    public JulianDay(long t) {
	setFromJd(JulianDay.dateToJulianDay(1970, 1, 1, false) + t / 86400000.0);
    }

    /**
     * Constructor for a Julian day
     * @param jd Julian day number
     */
    public JulianDay(double jd) {
	setFromJd(jd);
    }

    private void setFromJd(double jd) {
	// The conversion formulas are from Meeus, Chapter 7
	double z = (int) (Math.abs(jd) + 0.5);
	if (jd < 0) z = -z;
	double a = z;
	if (z >= 2299161.0) { // Gregorian
	    int a2 = (int) ((z - 1867216.25) / 36524.25);
	    a += 1 + a2 - (int) (a2 / 4.0);
	}
	double b = a + 1524;
	int c = (int) ((b - 122.1) / 365.25);
	int d = (int) (c * 365.25);
	int e = (int) ((b - d) / 30.6001);

	double f = jd + 0.5 - z;
	double exactDay = f + b - d - (int) (30.6001 * e);
	day = (int) exactDay;
	month = e - 1;
	if (month > 12) month = month - 12;
	year = c - 4715;
	if (month > 2) year = year - 1;
	setDayFraction(exactDay - day);	
    }

    /**
     * Returns the day fraction, from 0 to 1
     * @return
     */
    public double getDayFraction() {
	return (hour + minute / 60.0 + second / 3600.0) / 24.0;
    }

    /**
     * Sets the day fraction (to set hour, minute, second)
     * @param f Day fraction, from 0 to 1 (exclusive)
     */
    public void setDayFraction(double f) {
	double frac = f * 24.0;
	hour = (int) frac;
	frac = (frac - hour) * 60.0;
	minute = (int) frac;
	second = (frac - minute) * 60.0;
    }    

    /**
     * Convert this instance of {@linkplain JulianDay} to a Julian day number. Dates 
     * for Gregorian or Julian calendars are handled automatically
     * @return The Julian day that corresponds to this {@linkplain JulianDay} instance
     */
    public double getJulianDay() {
	return dateToJulianDay(year, month, day, isJulian()) + getDayFraction();
    }

    /**
     * Check if the {@linkplain JulianDay} instance contains an invalid date. A date is
     * invalid between October 5, 1582 and October 14, 1582.
     * @return true if the date is invalid, false otherwise.
     */
    public boolean isInvalid() {
	return (year == 1582 && month == 10 && (day >= 5 && day < 15));
    }

    /**
     * True if the instance represents a date in the Julian calendar, before October, 15, 1582
     * @return
     */
    public boolean isJulian() {
	if (year < 1582) return true;
	if (year == 1582 && month < 10) return true;
	if (year == 1582 && month == 10 && day < 15) return true;
	return false;
    }

    /**
     * Convert this {@linkplain JulianDay} instance to a String as YYYY-MM-DD hh:mm:ss.sss
     * @return A date/time String as YYYY-MM-DD hh:mm:ss.sss
     */
    public String toString() {
	return ""+year + "-" + Util.fmt02(month, "-") + Util.fmt02(day, " ") + 
		Util.fmt02(hour, ":") + Util.fmt02(minute, ":") + Util.formatValue(second, 3);
    }

    /**
     * Convert a date (0h) to a Julian Day. See Meeus, Astronomical Algorithms, chapter 7
     * @param year The year
     * @param month The month
     * @param day The day
     * @param julian true = Julian calendar, false for Gregorian. If not sure, enter false
     * @return The Julian day for the date and calendar specified
     */
    public static double dateToJulianDay(int year, int month, int day, boolean julian) {
	if (month < 3) {
	    year = year - 1;
	    month = month + 12;
	}
	int a = year / 100;
	int b = 0;
	if (!julian) b = 2 - a + a / 4;
	return (int) (365.25 * (year + 4716)) + (int) (30.6001 * (month + 1)) + day + b - 1524.5;
    }    

    /**
     * Test program
     * @param s Not used
     */
    public static void main(String[] s) {
	try {
	    System.out.println("JD TEST 1: GREGORIAN (NOW)");
	    JulianDay jd = new JulianDay(System.currentTimeMillis()); // UTC
	    System.out.println("JD:      " + jd.getJulianDay());
	    System.out.println("STR:     " + jd.toString());
	    System.out.println("Julian?  " + jd.isJulian());
	    System.out.println("Invalid? " + jd.isInvalid());
	    jd.setFromJd(jd.getJulianDay());
	    jd.setDayFraction(jd.getDayFraction());
	    System.out.println("NOW:     " + jd.toString());

	    System.out.println();
	    System.out.println("JD TEST 2: JULIAN");
	    JulianDay jdJulian = new JulianDay(1582, 9, 15);
	    jdJulian.setDayFraction(0.25);
	    System.out.println("JD:      " + jdJulian.getJulianDay());
	    System.out.println("STR:     " + jdJulian.toString());
	    System.out.println("Julian?  " + jdJulian.isJulian());
	    System.out.println("Invalid? " + jdJulian.isInvalid());
	    jdJulian.setFromJd(jdJulian.getJulianDay());
	    jdJulian.setDayFraction(jdJulian.getDayFraction());
	    System.out.println("NOW:     " + jdJulian.toString());

	    System.out.println();
	    System.out.println("JD TEST 3: INVALID DATE");
	    new JulianDay(1582, 10, 10);

	    /*
            JD TEST 1: GREGORIAN (NOW)
            JD:      2460550.127662118
            STR:     2024-08-27 15:03:50.007
            Julian?  false
            Invalid? false
            NOW:     2024-08-27 15:03:50.007

            JD TEST 2: JULIAN
            JD:      2299140.75
            STR:     1582-09-15 06:00:0.000
            Julian?  true
            Invalid? false
            NOW:     1582-09-15 06:00:0.000

            JD TEST 3: INVALID DATE
            java.lang.IllegalArgumentException: Date is invalid
            	at journal.JulianDay.<init>(JulianDay.java:33)
            	at journal.JulianDay.main(JulianDay.java:181)	    
	     */
	} catch (Exception exc) {
	    exc.printStackTrace();
	}
    }
}
