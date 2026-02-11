package journal;

import java.util.Arrays;

/** A class to compute the ephemerides of the planets and minor bodies with numerical integration, using {@linkplain EphemReduction} by inheritance */
public class EphemPlanets extends EphemReduction {

    private static final double C2 = Math.pow(Constant.SPEED_OF_LIGHT * 0.001 * Constant.SECONDS_PER_DAY / Constant.AU, 2); // AU^2 / d^2
    private static final double GMS = 2.9591220828411956E-4; // GMS in DE441, AU^3/d^2, assuming Msun = 1
    private static final double SOLAR_SYSTEM_DATA_EPOCH = Constant.J2000;
    private static final double[][] SOLAR_SYSTEM_DATA = new double[][] {
	// x, y, z, vx, vy, vz (AU, AU/d, respect Sun, JD 2451545.0 TDB), mass (respect Sun), radius (km), index in array bodyName. From DE440. Do NOT change order
	new double[] {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 696000, 0}, // Sun
	new double[] {-0.1300936053754522, -0.40059372164232543, -0.20048930201672596, 0.021366395668016163, -0.004926299692875428, -0.004847433077772866,  1.660120825489089E-7, 2440.53, 1}, // Mercury
	new double[] {-0.718302296345389, -0.04627424670211335, 0.02464063845542861, 7.981175157753219E-4, -0.018491837481062413, -0.008369735338020125,  2.447838287796944E-6, 6051.8, 2}, // Venus
	new double[] {1.390715921746287, 0.001401217626814569, -0.036960167196011424, 6.71499521033585E-4, 0.013814037515614361, 0.006317900433310847,  3.2271560829138995E-7, 3396.19, 3}, // Mars
	new double[] {4.001177161126057, 2.7365787240216024, 1.0755122808242419, -0.004568313526752718, 0.005881462129979568, 0.0026323030159255195,  9.547919099414246E-4, 71492.0, 4}, // Jupiter
	new double[] {6.406408859532808, 6.174657792651239, 2.274770783705428, -0.00429235187384325, 0.0035283445659715523, 0.0016419315191857945,  2.8588567002459455E-4, 60268.0, 5}, // Saturn
	new double[] {14.431856614381783, -12.506266259007408, -5.681690059289029, 0.00267810559142015, 0.002462004302669807, 0.0010404094481152937,  4.3662496132221186E-5, 25559.0, 6}, // Uranus
	new double[] {16.812046968052883, -22.980100505749114, -9.824427653612803, 0.002579274259047737, 0.001668425282316438, 6.188152032295604E-4,  5.151383772654574E-5, 24764.0, 7}, // Neptune
	new double[] {-9.87535222992358, -27.978868119163504, -5.753691421762491, 0.0030287508460142658, -0.001127593278936313, -0.001265129364676525,  7.350478973158631E-9, 1188.3, 8}, // Pluto
	new double[] {-1.771350992727098E-01, 8.874285223255191E-01, 3.847428990882070E-01, -1.720762506872895E-02, -2.898167717572411E-03, -1.256395052182784E-03,  3.003489615465139E-6, 6378.1366, 10}, // Earth
	new double[] {-1.790843809223965E-01, 8.856456304126460E-01, 3.842341853815847E-01, -1.683595459141215E-02, -3.282865544741707E-03, -1.430425208901575E-03,  3.6943033501098785E-8, 1737.4, 11}, // Moon
	//new double[] {-2.379327705915647, 5.456711318631370E-01, 7.412254807065680E-01, -3.584228273182221E-03, -9.845217307745823E-03, -3.904543826033526E-03,  4.719142276709673E-10, 469.7, 12}, // Ceres
	new double[] {-1.037925696095382, -0.1268092611155217, -0.07404282940865300, 4.227374301983514E-03, -1.412107207092049E-02, -5.145341154826345E-03, 0, 2, 13}, // Apophis, JPL elements (Horizons)
	//new double[] {-2.594982968100301E-01, -2.604198475982948, -1.161962651205037, 7.551958208072309E-03, 5.004810090124190E-03, 2.674350852369668E-03, 0, 0.03, 14}, // 2024 YR4, JPL elements
    };
    private static final String[] BODY_NAME = "Sun,Mercury,Venus,Mars,Jupiter,Saturn,Uranus,Neptune,Pluto,EMB,Earth,Moon,Ceres,Apophis,2024 YR4,Pallas,Juno,Vesta,Iris,Hygiea".split(","); // Keep consistent with SOLAR_SYSTEM_DATA !

    private int body = 0;
    private double lastPos[][], lastJD;
    private int timeStepSeconds = 8640;

    public EphemPlanets(double jdUT, double lon, double lat, double alt, TWILIGHT tw, TWILIGHT_MODE twm, int tz) {
	super(jdUT, lon, lat, alt, tw, twm, tz);
    }

    // Clone a 2d double array, some work is needed for this in Java ...
    private double[][] clone(double[][] in) {
	double[][] out = new double[in.length][];
	for (int i = 0; i < in.length; i ++) {
	    out[i] = in[i].clone();
	}
	return out;
    }
    
    private double[][] integrateTo(double[][] data, double jdStart, double jdEnd) {
	if (data == null) {
	    data = clone(SOLAR_SYSTEM_DATA);
	    jdStart = SOLAR_SYSTEM_DATA_EPOCH;
	}
	int nBodies = data.length;
	double[][] pos = clone(data);
	double jd = jdStart;
	if (jdEnd == jd) return pos;
	double speed = Math.abs(timeStepSeconds) / Constant.SECONDS_PER_DAY;
	if (jdEnd < jd) speed = -speed;
	if (Math.abs(jdEnd - jd) < Math.abs(speed)) speed = jdEnd - jd;

	// The integration schema use an explicit Runge-Kutta (RK). See the Butcher tables of some at https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
	// Bogacki-Shampine 5-4 method, by P. Bogacki, L.F. Shampine, An efficient Runge-Kutta (4,5) pair, Computers and Mathematics with Applications,
	// Volume 32, Issue 6, September 1996, Pages 15-28
	double[] ak = new double[] { 0, 1 / 6., 4 / 27., 1053 / 1372., 1960 / 3861., 4617 / 20480., 482048 / 414219., 7267 / 94080.,
		    2 / 27., -162 / 343., 42 / 143., 58653 / 366080., -29421 / 29068., 2152 / 5985.,
		    183 / 1372., -4 / 11., 63099 / 585728., 666106 / 1039181., 387 / 44800., 68 / 297., 81 / 352., 8152137 / 19744439., 24353 / 124800.,
		    597 / 22528., -30942 / 79937., 4440339 / 15491840., 174197 / 959244., 0, 587 / 8064. };
	double[] bk = new double[] { 587 / 8064.0, 0, 4440339 / 15491840., 24353 / 124800., 387 / 44800., 2152 / 5985., 7267 / 94080., 0 };

	// Here we substract J2000 date to integration start/end dates. This is EXTREMELY IMPORTANT to save digits in double precision
	jd -= Constant.J2000;
	jdEnd -= Constant.J2000;

	// Arrays, define once here to improve memory handling
	double[][] poorNewGravity = new double[nBodies][3];
	double[][] poorNewPosVel = new double[nBodies][6];
	double[][] futurePosVel = new double[nBodies][6];
	double[][][] fpv = new double[bk.length][nBodies][6];
	while (true) {
	    // Reset previous integration step
	    for (int i = 0; i < nBodies; i++) {
		Arrays.fill(futurePosVel[i], 0.0);
		Arrays.fill(poorNewPosVel[i], 0.0);
	    }

	    evaluateRK(ak, bk, pos, poorNewPosVel, poorNewGravity, futurePosVel, fpv, speed);

	    // Make integration effective on vectors and date
	    for (int i = 0; i < nBodies; i++) {
		for (int j = 0; j < 6; j++) {
		    pos[i][j] += futurePosVel[i][j];
		}
	    }
	    jd += speed;

	    // Control end of integration
	    double difEnd = Math.abs(jdEnd - jd);
	    if (difEnd < 1E-6) break;
	    if (difEnd < Math.abs(speed)) speed = (jdEnd - jd);
	}
	return pos;
    }

    private void evaluateRK(double[] ak, double[] bk, double[][] pos, double[][] poorNewPosVel, double[][] poorNewGravity,
	    double[][] futurePosVel, double[][][] fpv, double speed) {
	int nBodies = pos.length;
	int akl = ak.length;
	int bkl = bk.length;
	boolean fsal = bk[bkl - 1] == 0.0; // FSAL = First Same At Last in the Runge-Kutta, to skip last evaluation step
	boolean triangular = akl > bkl;
	int delta0 = bk.length - 2;
	for (int k = 0; k < bkl; k++) {
	    if (k > 0) {
		double a = ak[k];
		int km1 = k - 1;
		int ai0 = k + delta0;
		for (int i = 0; i < nBodies; i++) {
		    for (int j = 0; j < 6; j++) {
			poorNewPosVel[i][j] = fpv[km1][i][j] * a;
		    }
		    if (!triangular || k == 1) continue;

		    // Add the rest of terms if they are present
		    int ai = ai0;
		    int delta = delta0;
		    for (int an = 1; an < k; an++) {
			if (akl <= ai) break;
			double akn = ak[ai];
			if (akn != 0) {
			    int kn = km1 - an;
			    for (int j = 0; j < 6; j++) {
				poorNewPosVel[i][j] += fpv[kn][i][j] * akn;
			    }
			}
			if (km1 == an) break;
			delta--;
			ai += delta;
		    }
		}
	    }
	    if (k == bkl - 1 && fsal) continue;

	    // Compute accelerations
	    computeAccelerations(pos, poorNewPosVel, poorNewGravity);

	    // Update fpv for next RK step, and the output pv vector
	    for (int i = 0; i < nBodies; i++) {
		for (int j = 0; j < 3; j++) {
		    int jp3 = j + 3;
		    fpv[k][i][j] = (pos[i][jp3] + poorNewPosVel[i][jp3]) * speed;
		    fpv[k][i][jp3] = poorNewGravity[i][j] * speed;
		}
		if (bk[k] == 0) continue;
		for (int j = 0; j < 6; j++) {
		    futurePosVel[i][j] += fpv[k][i][j] * bk[k];
		}
	    }
	}
    }

    private void computeAccelerations(double[][] pos, double[][] newPos, double[][] gravity) {
	int nBodies = pos.length;

	// Compute Newtonian gravity and the distance dependent term of the relativistic correction
	Arrays.fill(gravity[0], 0.0); // Reset gravity of Sun, that should be the first body
	double[] dPos = new double[6];
	for (int i = 0; i < nBodies - 1; i++) {
	    boolean isun = (int) pos[i][8] == 0;
	    boolean iearth = (int) pos[i][8] == 10;
	    for (int j = i + 1; j < nBodies; j++) {
		boolean jmoonOrAster = (int) pos[j][8] >= 11;
		double d2 = 0;
		for (int k=0; k<3; k ++) {
		    dPos[k] = pos[i][k] - pos[j][k];
		    if (newPos != null) dPos[k] += newPos[i][k] - newPos[j][k];
		    d2 += dPos[k] * dPos[k];
		}
		double d = Math.sqrt(d2);
		double d3 = d2 * d;

		// Add relativistic terms from the Sun using Damour-Deruelle (1985) GR formulation (or Brumberg 1972, Sitarsky 1983)
		if (isun) {
		    double v2 = 0;
		    for (int k=3; k<6; k ++) {
			dPos[k] = pos[i][k] - pos[j][k];
			if (newPos != null) dPos[k] += newPos[i][k] - newPos[j][k];
			v2 += dPos[k] * dPos[k];
		    }
		    double pv = dPos[0] * dPos[3] + dPos[1] * dPos[4] + dPos[2] * dPos[5];

		    double a = GMS * (4.0 / d) - v2;
		    double b = 4.0 * pv;
		    double prefac = -GMS / (d3 * C2); // AU^-2
		    for (int k=0; k<3; k ++) {
			gravity[j][k] = prefac * (a * dPos[k] + b * dPos[k + 3]); // Set gravity for all (j) bodies except the Sun
		    }
		}

		double gmsrj = -GMS * pos[j][6] / d3;
		double gmsri = -GMS * pos[i][6] / d3;
		if (isun && (int) pos[j][8] == 9) gmsri *= 1.000000061783; // Approximately account for the effects of Earth tides on the barycenter
		for (int k=0; k < 3; k ++) {
		    gravity[i][k] += gmsrj * dPos[k];
		    gravity[j][k] -= gmsri * dPos[k];
		}

		// Other forces: oblateness of Earth and NGF
		if (iearth && jmoonOrAster) oblatenessJ2(dPos, gravity, d, j, i, pos, -1.08262539E-3);
		if (isun && (int) pos[j][8] == 13) apophisNGF(dPos, d, gravity[j], -1);
	    }
	}
    }

    private void oblatenessJ2(double[] pv, double[][] gravity, double r, int affectedBody, int mainBody, double[][] pos, double j2) {
	// https://en.wikipedia.org/wiki/Geopotential_model#The_deviations_of_Earth.27s_gravitational_field_from_that_of_a_homogeneous_sphere
	double r2 = r * r;
	double costheta2 = (pv[2] * pv[2]) / r2;
	double j2e0 = j2 * Math.pow(pos[mainBody][7] / Constant.AU, 2) * pos[mainBody][6];
	double j2eMain = 1.5 * j2e0 * GMS / (r2 * r2 * r);
	double j2eFac = 5.0 * costheta2 - 1.0;
	double[] accj2 = new double[] {j2eMain * j2eFac * pv[0], j2eMain * j2eFac * pv[1], j2eMain * (j2eFac - 2.0) * pv[2]};
	double sc = -pos[affectedBody][6] / pos[mainBody][6];
	for (int k = 0; k < 3; k++) {
	    gravity[affectedBody][k] += accj2[k];
	    gravity[mainBody][k] += accj2[k] * sc;
	}
    }

    // Calculates accelerations on Apophis due to non-gravitational forces, see Marsden et al. (1973), Astron. J. 78, 211-225.
    private void apophisNGF(double[] pv, double r, double[] gravity, double f) {
	// r dot v vector
	double rv = 0;
	for (int k = 0; k < 3; k++) {
	    rv += pv[k] * pv[k + 3];
	}

	// Within-orbital-plane transverse vector components
	double r2 = r * r;
	double tx = r2 * pv[3] - rv * pv[0];
	double ty = r2 * pv[4] - rv * pv[1];
	double tz = r2 * pv[5] - rv * pv[2];

	// Multiplication factors. NGF (A) values are read in AU/s^2. a3 = 0 for Apophis
	double a1 = 4.999999873689E-13 / r;
	double a2 = -2.901766720242E-14 / Math.sqrt(tx * tx + ty * ty + tz * tz);

	// X, Y and Z components of non-gravitational acceleration
	gravity[0] += f * (a1 * pv[0] + a2 * tx) / r2;
	gravity[1] += f * (a1 * pv[1] + a2 * ty) / r2;
	gravity[2] += f * (a1 * pv[2] + a2 * tz) / r2;
    }

    /** Sets the body to compute, from 0 (Sun) to the number of bodies supported minus 1.
     * @param i Body index. Default is 0 (Sun). */
    public void setBodyIndex(int i) {
	body = i;
    }

    /** Selects the integration time step in seconds.
     * @param seconds Integration step. Default is 8640 seconds or 0.1 days. */
    public void setTimeStep(int seconds) {
	timeStepSeconds = seconds;
    }

    /** Returns the number of bodies supported.
     * @return Number of bodies. */
    public static int getNumberOfBodies() {
	return SOLAR_SYSTEM_DATA.length;
    }

    /** Returns the name of a body.
     * @param index Index of the body, from 0 to the number of bodies supported minus 1.
     * @return Body name. */
    public String getBodyName(int index) {
	return BODY_NAME[(int) SOLAR_SYSTEM_DATA[index][8]];
    }

    /** Correction for precession using the algorithms by Laskar, 1986, consistent with the work by Bretagnon (Sun position) and JPL
     * @param in Input equatorial vector, positions and velocities
     * @param jd Julian day in UT
     * @param toJ2000 True to transform to J2000, false to transform from J2000 to this date
     * @return New vector */
    public static double[] precessionLaskarToOrFromJ2000(double[] in, double jd, boolean toJ2000) {
	/* Evaluation of Laskar precession angles */
	double[] pApol = new double[] { 0, -8.66e-10, -4.759e-8, 2.424e-7, 1.3095e-5, 1.7451e-4, -1.8055e-3, -0.235316, 0.07732, 111.1971, 50290.966 };
	double[] Wpol = new double[] { 6.6402e-16, -2.69151e-15, -1.547021e-12, 7.521313e-12, 6.3190131e-10, -3.48388152e-9, -1.813065896e-7,
		2.75036225e-8, 7.4394531426e-5, -0.042078604317, 3.052112654975 };
	double[] zpol = new double[] { 1.2147e-16, 7.3759e-17, -8.26287e-14, 2.503410e-13, 2.4650839e-11, -5.4000441e-11, 1.32115526e-9,
		-5.998737027e-7, -1.6242797091e-5, 0.002278495537, 0.0 };
	double t = EarthAngles.toCenturiesRespectJ2000(jd, true) * 0.1;
	double pA = pApol[0], W = Wpol[0], z = zpol[0];
	for (int i = 1; i < pApol.length; i ++) {
	    pA = pA * t + pApol[i];
	    W = W * t + Wpol[i];
	    z = z * t + zpol[i];
	}
	pA *= Constant.ARCSEC_TO_RAD * t;

	// Rotation angles to or from J2000, note input and output are equatorial
	double j2000_ut = Constant.J2000 - EarthAngles.TTminusUT1(Constant.J2000);
	double[] rotAngles = null;
	if (!toJ2000) rotAngles = new double[] {EarthAngles.meanObliquity(j2000_ut), W, z, -W - pA, -EarthAngles.meanObliquity(jd)};
	if (toJ2000) rotAngles = new double[] {EarthAngles.meanObliquity(jd), W + pA, -z, -W, -EarthAngles.meanObliquity(j2000_ut)};

	/* Implementation by elementary rotations using expansions. First rotate about the x axis from the initial equator to the ecliptic */
	double[] out = CoordinateSystem.rotate(in, CoordinateSystem.getRotX(rotAngles[0]));
	/* Rotate about z axis to the node */
	out = CoordinateSystem.rotate(out, CoordinateSystem.getRotZ(rotAngles[1]));
	/* Rotate about new x axis by the inclination of the moving ecliptic on the ecliptic for the initial time */
	out = CoordinateSystem.rotate(out, CoordinateSystem.getRotX(rotAngles[2]));
	/* Rotate about new z axis back from the node */
	out = CoordinateSystem.rotate(out, CoordinateSystem.getRotZ(rotAngles[3]));
	/* Rotate about x axis to final equator */
	return CoordinateSystem.rotate(out, CoordinateSystem.getRotX(rotAngles[4]));
    }

    private double[] equatorialJ2000ToEclipticOfDate(double[] pv, double jd_UT) {
	return CoordinateSystem.rotate(precessionLaskarToOrFromJ2000(pv, jd_UT, false), CoordinateSystem.getRotX(EarthAngles.meanObliquity(jd_UT)));
    }

    /** Geocentric position of the planets.
     * @return Array with ecliptic longitude, latitude, distance, and angular radius */
    @Override
    public double[] getBodyPosition() {
	double t = EarthAngles.toCenturiesRespectJ2000(jd_UT, true);
	double jdTT = Constant.J2000 + t * Constant.JULIAN_DAYS_PER_CENTURY;
	double[][] pos = lastPos;

	if (lastJD != jdTT) { // Avoid a full re-integration when re-evaluating positions to compute rise/set times
	    pos = integrateTo(lastPos, lastJD, jdTT);
	    if (Math.abs(lastJD - jdTT) > 1) { 
		// But do not save the results of these little integrations during rise/set computations, since they will affect Apophis
		lastPos = pos;
		lastJD = jdTT;
	    }
	}

	// Get the geocentric ecliptic position of the Sun for the equinox of date. For years between 1200 and 2800
	// the numerical integration is used in case the Earth body is present. Outside this interval, or when the
	// Earth is not included in the integration, the algorithm from Planetary Programs and Tables (PPT) is used
	// instead, since it is more accurate than the numerical integration model
	double[] sunPos = new double[6];
	boolean usePPT = Math.abs(t) > 8 || (int) pos[9][8] != 10;
	if (usePPT) sunPos = EphemSun.getGeometricEclipticPositionEquinoxOfDate(jdTT);

	// Get the heliocentric ecliptic position of the body for the equinox of date
	double[] pv = new double[6];
	for (int i=0; i<6; i++) {
	    pv[i] = pos[body][i] - pos[0][i]; // Relative to the Sun, equatorial J2000
	    if (!usePPT) sunPos[i] = pos[0][i] - pos[9][i];
	}
	double[] meanEcl = equatorialJ2000ToEclipticOfDate(pv, jd_UT);
	if (!usePPT) sunPos = equatorialJ2000ToEclipticOfDate(sunPos, jd_UT);

	// Get the geocentric ecliptic position of the body for the equinox of date
	for (int i=0; i<6; i++) {
	    pv[i] = meanEcl[i] + sunPos[i];
	}
	double[] meanEclSph = CoordinateSystem.cartesianToSpherical(pv);

	// First order light-time correction
	for (int i=0; i<3; i++) {
	    pv[i] -= pv[i+3] * meanEclSph[2] * Constant.LIGHT_TIME_DAYS_PER_AU;
	}
	meanEclSph = CoordinateSystem.cartesianToSpherical(pv);

	return new double[] {meanEclSph[0], meanEclSph[1], meanEclSph[2], Math.atan(pos[body][7] / (meanEclSph[2] * Constant.AU))};
    }

    /** Test program
     * @param args Not used */
    public static void main(String[] args) {
	// Prepare input data
	int year = 2029, month = 4, day = 13, h = 21, m = 38, s = 0; // For the close encounter of Apophis
	JulianDay jday = new JulianDay(year, month, day);
	jday.setDayFraction((h + m / 60.0 + s / 3600.0) / 24.0);

	double jdUTC = jday.getJulianDay();
	double lon = -(3 + 42 / 60.0); // degrees, for Madrid
	double lat = 40 + 26 / 60.0;
	double alt = 0; // m
	int tz = 2; // h
	TWILIGHT tw = TWILIGHT.HORIZON_34arcmin;
	TWILIGHT_MODE twm = TWILIGHT_MODE.TODAY_UT;

	// Report calculation time
	double jdTT = Constant.J2000 + EarthAngles.toCenturiesRespectJ2000(jdUTC, true) * Constant.JULIAN_DAYS_PER_CENTURY;
	System.out.println("Date (UTC): "+jday.toString());
	System.out.println("JD (UTC):   "+jdUTC);
	System.out.println("JD (TT):    "+jdTT);
	System.out.println();

	// Compute the ephemerides data
	EphemPlanets planEph = new EphemPlanets(jdUTC, lon, lat, alt, tw, twm, tz);
	EphemSun sunEph = new EphemSun(jdUTC, lon, lat, alt, tw, twm, tz);
	EphemData sunData = EphemReduction.getEphemeris(sunEph);
	for (int body = 0; body < EphemPlanets.getNumberOfBodies(); body ++) {
	    planEph.setBodyIndex(body);
	    EphemData planetData = EphemReduction.getEphemeris(planEph);
	    if (body != 0) planetData.setIlluminationPhase(sunData);

	    // Report
	    System.out.println(planEph.getBodyName(body));
	    System.out.println(planetData.toString());
	    if (planEph.getBodyName(body).equals("Apophis")) {
		double mag = planetData.getAsteroidApparentMagnitude(sunData, 19.09, 0.24);
		System.out.println(" App. Mag: " + (float) mag);
	    }
	    if (planEph.getBodyName(body).equals("2024 YR4")) {
		double mag = planetData.getAsteroidApparentMagnitude(sunData, 23.92, 0.15);
		System.out.println(" App. Mag: " + (float) mag);
	    }
	}
/*
Mars.      Horizons: RA 179.238743, DEC 3.211346, Dist 0.6712987 (distance in Horizons follows a different criteria)
 Az:       154.49199°
 El:       50.060345°
 Dist:     0.6712703 au
 RA:       179.23877°
 DEC:      3.2113345°
 Ill:      98.10554%
 ang.R:    0.0019377233°
 Rise:     2029/04/13 16:29:30 UT
 Set:      2029/04/13 04:59:35 UT
 Transit:  2029/04/13 22:42:06 UT (elev. 52.7938°)
 
Apophis.   Horizons: RA 125.12850, DEC 23.17325, Dist 0.00021921426064
 Az:       253.27751°
 El:       53.756172°
 Dist:     2.1911896E-4 au
 RA:       125.12876°
 DEC:      23.17038°
 Ill:      56.757015%
 ang.R:    0.0034958057°
 Rise:     NO RISE/SET/TRANSIT FOR THIS OBSERVER/DATE
 Set:      2029/04/13 04:48:18 UT
 Transit:  2029/04/13 00:42:40 UT (elev. 20.622501°)
 */
    }    
}
