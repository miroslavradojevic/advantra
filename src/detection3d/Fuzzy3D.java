package detection3d;

import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ColorProcessor;

import java.io.File;
import java.util.Arrays;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 12/13/13
 * Time: 8:54 AM
 */
public class Fuzzy3D {

	public float iDiff;

	static int 		N = 101; 				// nr points covering out membership functions
	static float 	TRIANGLE_HALF_W = 0.5f; // half width of the membership function triangle
	static float[] 	x;      				// serve as x axis for membership functions (range 0 to 1)
	static float[] 	agg;
	static float[] 	v_YES;
	static float[] 	v_MAYBE;
	static float[] 	v_NO;

	// membership fuzzy sets (rules will call these categories "YES", "MAYBE", "NO")
	private static float[] q_YES;
	private static float[] q_MAYBE;
	private static float[] q_NO;

	String exportFile = System.getProperty("user.home")+ File.separator+"fuzzy3d.dat";

	public Fuzzy3D(float _iDiff) {

		iDiff = _iDiff;
		x = new float[N];
		for (int i=0; i<N; i++) x[i] = i * (1f / (N-1));

		/*
		 "YES" - define membership function: from TRIANGLE_HALF_W to 1.0
		 */
		q_YES = new float[N];
		float A = 1 - TRIANGLE_HALF_W;
		float B = 1;
		float C = 1 + TRIANGLE_HALF_W;
		for (int i=0; i<N; i++) {
			if (x[i]<=A) q_YES[i] = 0;
			else if (x[i]>A && x[i]<=B) q_YES[i] = (x[i]-A) / TRIANGLE_HALF_W;
				//else if (x[i]>B && x[i]<=C) q_YES[i] = 1 - (x[i]-B) / TRIANGLE_HALF_W;
			else q_YES[i] = 0;
		}

		/*
			"MAYBE"
		 */
		q_MAYBE = new float[N];
		A = .5f - TRIANGLE_HALF_W;
		B = .5f;
		C = .5f + TRIANGLE_HALF_W;
		boolean flag = false;
		for (int i=0; i<N; i++) {
			if (x[i]<=A) q_MAYBE[i] = 0;
			else if (x[i]>A && x[i]<=B) q_MAYBE[i] = (x[i]-A) / TRIANGLE_HALF_W;
				//else if (x[i]> B && !flag) { flag = true; q_MAYBE[i] = 1; }
			else if (x[i]>B && x[i]<=C) q_MAYBE[i] = 1 - (x[i]-B) / TRIANGLE_HALF_W;
			else q_MAYBE[i] = 0;
		}

		/*
			"NO"
		 */
		q_NO = new float[N];
		A = .0f - TRIANGLE_HALF_W;
		B = .0f;
		C = .0f + TRIANGLE_HALF_W;
		for (int i=0; i<N; i++) {
			//if (x[i]<=A) q_NO[i] = 0;
			//else if (x[i]>A && x[i]<=B) q_NO[i] = (x[i]-A) / TRIANGLE_HALF_W;
			if (x[i]<=B) q_NO[i] = 1;
			else if (x[i]>B && x[i]<=C) q_NO[i] = 1 - (x[i]-B) / TRIANGLE_HALF_W;
			else q_NO[i] = 0;
		}

		agg = new float[N];
		v_YES = new float[N];
		v_MAYBE = new float[N];
		v_NO = new float[N];

	}

	/*
        fuzzification - input linguistic variable
	 */

	private float h_low(float theta)
	{
		if (theta<=0)
			return 1;
		else if(theta>0 && theta<=iDiff/2)
			return 1 - theta / (iDiff/2);
		else
			return 0;
	}

	private float h_mid(float theta)
	{
		if (theta<=0)
			return 0;
		else if (theta>0 && theta<=iDiff/2)
			return (theta) / (iDiff/2);
			//else if (theta>(dDiff/2) && theta<=iDiff)
			//	return 1;
		else if (theta>iDiff/2 && theta<=iDiff)
			return 1 - (theta-iDiff/2) / (iDiff/2);
		else
			return 0;
	}

	private float h_hgh(float theta)
	{
		if (theta<=iDiff/2)
			return 0;
		else if (theta>iDiff/2 && theta<=iDiff)
			return (theta-iDiff/2) / (iDiff/2);
		else
			return 1;

	}

	/*
		defuzzification
	 */
	private float[] fi_YES(float z)
	{
		//float[] out = new float[N];

		Arrays.fill(v_YES, 0);

		for (int i=0; i<N; i++)
			if (q_YES[i]<=z)
				v_YES[i] = q_YES[i];
			else
				v_YES[i] = z;

		return v_YES;
	}

	private float q_YES(float fi) // check paper for notation - sample out function fi ranges 0-1 at index corresponding to fi value
	{
		fi = (fi>1f)? 1f : fi ;
		fi = (fi<0f)? 0f : fi ;
		return q_YES[Math.round((fi/1f)*N)];
	}

	private float[] fi_MAYBE(float z)
	{

		Arrays.fill(v_MAYBE, 0);

		for (int i=0; i<N; i++)
			if (q_MAYBE[i]<=z)
				v_MAYBE[i] = q_MAYBE[i];
			else
				v_MAYBE[i] = z;

		return v_MAYBE;
	}

	private float q_MAYBE(float fi)
	{
		fi = (fi>1f)? 1f : fi ;
		fi = (fi<0f)? 0f : fi ;
		return q_MAYBE[Math.round((fi/1f)*N)];
	}

	private float[] fi_NO(float z)
	{

		Arrays.fill(v_NO, 0);

		for (int i=0; i<N; i++)
			if (q_NO[i]<=z)
				v_NO[i] = q_NO[i];
			else
				v_NO[i] = z;

		return v_NO;
	}

	private float q_NO(float fi)
	{
		fi = (fi>1f)? 1f : fi ;
		fi = (fi<0f)? 0f : fi ;
		return q_NO[Math.round((fi/1f)*N)];
	}

	/*
		FLS rules
	 */
	public float bifurcationess(float theta0, float theta1, float theta2, float theta3) {

		Arrays.fill(agg, 0);

		float[] cur;
		float mu;

		// Rule 1: all are high 				(theta0=="HIGH")&&(theta1=="HIGH")&&(theta2=="HIGH")&&(theta3=="HIGH") => J=="YES"
		mu = min (h_hgh(theta0), h_hgh(theta1), h_hgh(theta2), h_hgh(theta3)); cur = fi_YES(mu); accumulate(cur, agg);

		// Rule 2: central is fair 				(theta0=="FAIR")&&(theta1=="HIGH")&&(theta2=="HIGH")&&(theta3=="HIGH") => J=="YES"
		mu = min (h_mid(theta0), h_hgh(theta1), h_hgh(theta2), h_hgh(theta3)); cur = fi_YES(mu); accumulate(cur, agg);

		// Rule 3: one of the branches is fair	(theta0=="HIGH")&&(theta1=="FAIR")&&(theta2=="HIGH")&&(theta3=="HIGH") => J=="MAYBE"
		mu = max(
						min (h_hgh(theta0), h_mid(theta1), h_hgh(theta2), h_hgh(theta3)),
						min (h_hgh(theta0), h_hgh(theta1), h_mid(theta2), h_hgh(theta3)),
						min (h_hgh(theta0), h_hgh(theta1), h_hgh(theta2), h_mid(theta3))
		);
		cur = fi_MAYBE(mu); accumulate(cur, agg);

		// Rule 4: two of the branches are fair (at least) (theta0=="HIGH")&&(theta1=="FAIR")&&(theta2=="FAIR")&&(theta3=="HIGH") => J=="NO"
		// should cover the case when three are FAIR as well, third is HIGH
		mu = max(
						min (h_hgh(theta0), h_mid(theta1), h_mid(theta2), h_hgh(theta3)), // 12
						min (h_hgh(theta0), h_mid(theta2), h_mid(theta3), h_hgh(theta1)), // 23
						min (h_hgh(theta0), h_mid(theta1), h_mid(theta3), h_hgh(theta2))  // 13
		);
		cur = fi_NO(mu); accumulate(cur, agg);

		// Rule 5: at least one branch is low 			(theta0=="HIGH")&&(theta1=="LOW")&&(theta2=="")&&(theta3=="") => J=="NO"
		mu = max(
						h_low(theta1), // 12 // , h_hgh(theta3)  , h_mid(theta1), h_mid(theta2)
						h_low(theta2),//min (h_hgh(theta0), h_mid(theta2), h_mid(theta3)), // 23 // , h_mid(theta3)
						h_low(theta3) //min (h_hgh(theta0), h_mid(theta1), h_mid(theta3))  // 13 // , h_mid(theta3)
		);
		cur = fi_NO(mu); accumulate(cur, agg);

		// find centroid of agg (defuzzification)
		float cx = 0;
		float cx_norm = 0;
		for (int i=1; i<N; i++) {
			cx += .5f * (agg[i]+agg[i-1]) * (x[i]-x[i-1]) * .5f * (x[i]+x[i-1]);
			cx_norm += .5f * (agg[i]+agg[i-1]) * (x[i]-x[i-1]);
		}

		if (cx_norm>Float.MIN_VALUE)
			cx /= cx_norm;

		return cx;

	}

//	// visualize decision making (3 inputs) with the RGB cube
//	public float[] decisionOutcomeTh123(float theta1, float theta2, float theta3){
//		return out_NMY;
//	}

	public ImagePlus decisionCube(){

		int W = 128;
		int H = 128;
		int L = 128;

		ImageStack is_cube = new ImageStack(W, H);

		for (int l=0; l<L; l++) {
			is_cube.addSlice(new ColorProcessor(W, H));
		}

		// fill the values in
		for (int x=0; x<W; x++) {
			for (int y=0; y<H; y++) {
				for (int z=0; z<L; z++) {

					// x will represent theta1
					float th_1 = ((float) x / W) * iDiff;
					float th_2 = ((float) y / H) * iDiff;
					float th_3 = ((float) z / L) * iDiff;

					float[] out_NMY = new float[3]; // out_NMY[0] ~ NO, out_NMY[1] ~ MAYBE, out_NMY[2] ~ YES

					float centroid = bifurcationess(iDiff, th_1, th_2, th_3);

					out_NMY[0] = q_NO(centroid);
					out_NMY[1] = q_MAYBE(centroid);
					out_NMY[2] = q_YES(centroid);

					// convert to colors
					int red = Math.round(out_NMY[2] * 255);
					red = (red>255)? 255 : red;
					red = (red<0  )? 0   : red;

					int green = Math.round(out_NMY[1] * 255);
					green = (green>255)? 255 : green;
					green = (green<0  )? 0   : green;

					int blue = Math.round(out_NMY[0] * 255);
					blue = (blue>255)? 255 : blue;
					blue = (blue<0  )? 0   : blue;

					int vox = ((red & 0xff)<<16)+((green & 0xff)<<8) + (blue & 0xff);

					is_cube.setVoxel(x, y, z, vox);

				}
			}
		}

		ImagePlus rgb_cube = new ImagePlus("decision_cube_rgb", is_cube);
		return rgb_cube;

	}


	/*
		some tools used for fuzzy set operations, to operate on the membership values
	 */
	public static final void accumulate(float[] values, float[] accumulator)
	{
		for (int i=0; i<N; i++) if (values[i]>accumulator[i]) accumulator[i] = values[i];
	}

	/*
	OR operation in rules
	 */
	public static final float max(float in1, float in2)
	{
		return Math.max(in1, in2);
	}

	public static final float max(float in1, float in2, float in3)
	{
		return Math.max(in1, max(in2, in3));
	}

	public static final float max(float in1, float in2, float in3, float in4)
	{
		return Math.max(in1, max(in2, in3, in4));
	}

	public static final float max(float in1, float in2, float in3, float in4, float in5)
	{
		return Math.max(in1, max(in2, in3, in4, in5));
	}

	public static final float max(float in1, float in2, float in3, float in4, float in5, float in6)
	{
		return Math.max(in1, max(in2, in3, in4, in5, in6));
	}

	public static final float max(float in1, float in2, float in3, float in4, float in5, float in6, float in7)
	{
		return Math.max(in1, max(in2, in3, in4, in5, in6, in7));
	}

	/*
	AND operation in rules
	 */

	public static final float min(float in1, float in2)
	{
		return Math.min(in1, in2);
	}

	public static final float min(float in1, float in2, float in3)
	{
		return Math.min(in1, min(in2, in3));
	}

	public static final float min(float in1, float in2, float in3, float in4)
	{
		return Math.min(in1, min(in2, in3, in4));
	}

	public static final float min(float in1, float in2, float in3, float in4, float in5)
	{
		return Math.min(in1, min(in2, in3, in4, in5));
	}

	public static final float min(float in1, float in2, float in3, float in4, float in5, float in6)
	{
		return Math.min(in1, min(in2, in3, in4, in5, in6));
	}

	public static final float min(float in1, float in2, float in3, float in4, float in5, float in6, float in7)
	{
		return Math.min(in1, min(in2, in3, in4, in5, in6, in7));
	}

}
