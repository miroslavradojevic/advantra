package detection;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Plot;

import java.awt.*;
import java.io.*;
import java.util.Arrays;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 9/1/13
 * Time: 9:32 AM
 */
public class Fuzzy {

	// parameters
	public float iDiff;//, ratio_iDiff, dDiff;

	static int 		N = 101; 	// nr points covering out membership functions
	static float 	TRIANGLE_HALF_W = 0.5f;
	static float[] 	x;      // serve as x axis for membership functions
	static float[] 	agg;
	static float[] 	v_YES;
	static float[] 	v_MAYBE;
	static float[] 	v_NO;

	// membership fuzzy sets (rules will call these categories "YES", "MAYBE", "NO")
	private static float[] q_YES;
	private static float[] q_MAYBE;
	private static float[] q_NO;

	String exportFile = "fuzzy.log";

	/*
		there could be three categories "YES" "MAYBE" "NO" refering to input being a junction
	 */

	public Fuzzy(float iDiff1) // , float ratio_iDiff1
	{

		iDiff = iDiff1;
		//ratio_iDiff = ratio_iDiff1;
		//dDiff = ratio_iDiff*iDiff;
		//dDiff = (dDiff<2)? 2 : dDiff ;
		//dDiff = (dDiff>=iDiff)? iDiff-1 : dDiff ;

		x = new float[N];
		for (int i=0; i<N; i++) x[i] = i * (1f / (N-1)); // -TRIANGLE_HALF_W + i * ((1f+2*TRIANGLE_HALF_W)/(N-1));

		// "YES"     - TRIANGLE_HALF_W from 1.0
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

		// "MAYBE"     +/- TRIANGLE_HALF_W around 0.5
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

		// "NO" 		+ TRIANGLE_HALF_W after 0.0
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

		//for (int i=0; i<N; i++) agg[i] = Math.max(Math.max(q_NO[i], q_MAYBE[i]), q_YES[i]);
		//M = 0;
		//for (int i=1; i<N; i++) M += .5f*(Ma[i]+Ma[i-1]) * (x[i]-x[i-1]);

	}

	/*
		membership fuzzy sets, 3 categories "HIGH" "MID" and "LOW"
		refering to how much above the background level is the intensity at input
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

	private float[] fi_MAYBE(float z)
	{
//		float[] out = new float[N];
		Arrays.fill(v_MAYBE, 0);
		for (int i=0; i<N; i++)
			if (q_MAYBE[i]<=z)
				v_MAYBE[i] = q_MAYBE[i];
			else
				v_MAYBE[i] = z;
		return v_MAYBE;
	}

	private float[] fi_NO(float z)
	{
//		float[] out = new float[N];
		Arrays.fill(v_NO, 0);
		for (int i=0; i<N; i++)
			if (q_NO[i]<=z)
				v_NO[i] = q_NO[i];
			else
				v_NO[i] = z;
		return v_NO;
	}

//	private float min6(float in1, float in2, float in3, float in4, float in5, float in6)
//	{
//		return Math.min(in1, Math.min(in2, Math.min(in3, Math.min(in4, Math.min(in5, in6)))));
//	}

	/*
	set of rules... how much it is "VERY HIGH", "HIGH", "MODERATE", "LOW" according to each rule
	 */

	public float bifurcationess(
		float theta0,
		float theta1, float theta2,
		float theta3, float theta4,
		float theta5, float theta6,
		boolean plotIt
	)
	{

		Arrays.fill(agg, 0);

		float[] cur;
		float mu;

		// all are high  => YES (theta0=="HIGH")&&(theta1=="HIGH")&&(theta2=="HIGH")&&...&&(theta6=="HIGH") => J=="YES"
		mu = min (h_hgh(theta0), h_hgh(theta1), h_hgh(theta2), h_hgh(theta3), h_hgh(theta4), h_hgh(theta5), h_hgh(theta6)); cur = fi_YES(mu); accumulate(cur, agg);

		// allow one to be mid, rest are high => MAYBE
		mu = min (h_mid(theta0), h_hgh(theta1), h_hgh(theta2), h_hgh(theta3), h_hgh(theta4), h_hgh(theta5), h_hgh(theta6)); cur = fi_MAYBE(mu); accumulate(cur, agg);
		mu = min (h_hgh(theta0), h_mid(theta1), h_hgh(theta2), h_hgh(theta3), h_hgh(theta4), h_hgh(theta5), h_hgh(theta6)); cur = fi_MAYBE(mu); accumulate(cur, agg);
		mu = min (h_hgh(theta0), h_hgh(theta1), h_mid(theta2), h_hgh(theta3), h_hgh(theta4), h_hgh(theta5), h_hgh(theta6)); cur = fi_MAYBE(mu); accumulate(cur, agg);
		mu = min (h_hgh(theta0), h_hgh(theta1), h_hgh(theta2), h_mid(theta3), h_hgh(theta4), h_hgh(theta5), h_hgh(theta6)); cur = fi_MAYBE(mu); accumulate(cur, agg);
		mu = min (h_hgh(theta0), h_hgh(theta1), h_hgh(theta2), h_hgh(theta3), h_mid(theta4), h_hgh(theta5), h_hgh(theta6)); cur = fi_MAYBE(mu); accumulate(cur, agg);
		mu = min (h_hgh(theta0), h_hgh(theta1), h_hgh(theta2), h_hgh(theta3), h_hgh(theta4), h_mid(theta5), h_hgh(theta6)); cur = fi_MAYBE(mu); accumulate(cur, agg);
		mu = min (h_hgh(theta0), h_hgh(theta1), h_hgh(theta2), h_hgh(theta3), h_hgh(theta4), h_hgh(theta5), h_mid(theta6)); cur = fi_MAYBE(mu); accumulate(cur, agg);

		// if one is low => NO - that means there is dissconntinuity between
		mu = min (h_low(theta0), h_hgh(theta1), h_hgh(theta2), h_hgh(theta3), h_hgh(theta4), h_hgh(theta5), h_hgh(theta6)); cur = fi_NO(mu); accumulate(cur, agg);
		mu = min (h_hgh(theta0), h_low(theta1), h_hgh(theta2), h_hgh(theta3), h_hgh(theta4), h_hgh(theta5), h_hgh(theta6)); cur = fi_NO(mu); accumulate(cur, agg);
		mu = min (h_hgh(theta0), h_hgh(theta1), h_low(theta2), h_hgh(theta3), h_hgh(theta4), h_hgh(theta5), h_hgh(theta6)); cur = fi_NO(mu); accumulate(cur, agg);
		mu = min (h_hgh(theta0), h_hgh(theta1), h_hgh(theta2), h_low(theta3), h_hgh(theta4), h_hgh(theta5), h_hgh(theta6)); cur = fi_NO(mu); accumulate(cur, agg);
		mu = min (h_hgh(theta0), h_hgh(theta1), h_hgh(theta2), h_hgh(theta3), h_low(theta4), h_hgh(theta5), h_hgh(theta6)); cur = fi_NO(mu); accumulate(cur, agg);
		mu = min (h_hgh(theta0), h_hgh(theta1), h_hgh(theta2), h_hgh(theta3), h_hgh(theta4), h_low(theta5), h_hgh(theta6)); cur = fi_NO(mu); accumulate(cur, agg);
		mu = min (h_hgh(theta0), h_hgh(theta1), h_hgh(theta2), h_hgh(theta3), h_hgh(theta4), h_hgh(theta5), h_low(theta6)); cur = fi_NO(mu); accumulate(cur, agg);

		if(false) {
		// allow two from different clusters to be mid => YES
		mu = min (h_hgh(theta0), h_mid(theta1), h_hgh(theta2), h_mid(theta3), h_hgh(theta4), h_hgh(theta5), h_hgh(theta6)); cur = fi_YES(mu); accumulate(cur, agg);// 1,3
		mu = min (h_hgh(theta0), h_mid(theta1), h_hgh(theta2), h_hgh(theta3), h_mid(theta4), h_hgh(theta5), h_hgh(theta6)); cur = fi_YES(mu); accumulate(cur, agg);// 1,4
		mu = min (h_hgh(theta0), h_mid(theta1), h_hgh(theta2), h_hgh(theta3), h_hgh(theta4), h_mid(theta5), h_hgh(theta6)); cur = fi_YES(mu); accumulate(cur, agg);// 1,5
		mu = min (h_hgh(theta0), h_mid(theta1), h_hgh(theta2), h_hgh(theta3), h_hgh(theta4), h_hgh(theta5), h_mid(theta6)); cur = fi_YES(mu); accumulate(cur, agg);// 1,6

		mu = min (h_hgh(theta0), h_hgh(theta1), h_mid(theta2), h_mid(theta3), h_hgh(theta4), h_hgh(theta5), h_hgh(theta6)); cur = fi_YES(mu); accumulate(cur, agg);// 2,3
		mu = min (h_hgh(theta0), h_hgh(theta1), h_mid(theta2), h_hgh(theta3), h_mid(theta4), h_hgh(theta5), h_hgh(theta6)); cur = fi_YES(mu); accumulate(cur, agg);// 2,4
		mu = min (h_hgh(theta0), h_hgh(theta1), h_mid(theta2), h_hgh(theta3), h_hgh(theta4), h_mid(theta5), h_hgh(theta6)); cur = fi_YES(mu); accumulate(cur, agg);// 2,5
		mu = min (h_hgh(theta0), h_hgh(theta1), h_mid(theta2), h_hgh(theta3), h_hgh(theta4), h_hgh(theta5), h_mid(theta6)); cur = fi_YES(mu); accumulate(cur, agg);// 2,6

		mu = min (h_hgh(theta0), h_hgh(theta1), h_hgh(theta2), h_mid(theta3), h_hgh(theta4), h_mid(theta5), h_hgh(theta6)); cur = fi_YES(mu); accumulate(cur, agg);// 3,5
		mu = min (h_hgh(theta0), h_hgh(theta1), h_hgh(theta2), h_mid(theta3), h_hgh(theta4), h_hgh(theta5), h_mid(theta6)); cur = fi_YES(mu); accumulate(cur, agg);// 3,6
		mu = min (h_hgh(theta0), h_hgh(theta1), h_hgh(theta2), h_hgh(theta3), h_mid(theta4), h_mid(theta5), h_hgh(theta6)); cur = fi_YES(mu); accumulate(cur, agg);// 4,5
		mu = min (h_hgh(theta0), h_hgh(theta1), h_hgh(theta2), h_hgh(theta3), h_mid(theta4), h_hgh(theta5), h_mid(theta6)); cur = fi_YES(mu); accumulate(cur, agg);// 4,6
		}

		// two from the same cluster mid => NO
		mu = min (h_hgh(theta0), h_mid(theta1), h_mid(theta2), h_hgh(theta3), h_hgh(theta4), h_hgh(theta5), h_hgh(theta6)); cur = fi_NO(mu); accumulate(cur, agg);
		mu = min (h_hgh(theta0), h_hgh(theta1), h_hgh(theta2), h_mid(theta3), h_mid(theta4), h_hgh(theta5), h_hgh(theta6)); cur = fi_NO(mu); accumulate(cur, agg);
		mu = min (h_hgh(theta0), h_hgh(theta1), h_hgh(theta2), h_hgh(theta3), h_hgh(theta4), h_mid(theta5), h_mid(theta6)); cur = fi_NO(mu); accumulate(cur, agg);

		// two from the same cluster low => NO (this would cover the case when theyre all low as well)
		mu = min (h_hgh(theta0), h_low(theta1), h_low(theta2), h_hgh(theta3), h_hgh(theta4), h_hgh(theta5), h_hgh(theta6)); cur = fi_NO(mu); accumulate(cur, agg);
		mu = min (h_hgh(theta0), h_hgh(theta1), h_hgh(theta2), h_low(theta3), h_low(theta4), h_hgh(theta5), h_hgh(theta6)); cur = fi_NO(mu); accumulate(cur, agg);
		mu = min (h_hgh(theta0), h_hgh(theta1), h_hgh(theta2), h_hgh(theta3), h_hgh(theta4), h_low(theta5), h_low(theta6)); cur = fi_NO(mu); accumulate(cur, agg);

		// find acc centroid (defuzzification)
		float cx = 0;
		float cx_norm = 0;
		for (int i=1; i<N; i++) {
			cx += .5f * (agg[i]+agg[i-1]) * (x[i]-x[i-1]) * .5f * (x[i]+x[i-1]);
			cx_norm += .5f * (agg[i]+agg[i-1]) * (x[i]-x[i-1]);
		}

		if (cx_norm>Float.MIN_VALUE)
			cx /= cx_norm;

//		// Plot
		if (plotIt) {
			Plot p = new Plot("", "x", "acc", x, agg);
			p.draw();
			p.setLineWidth(3);
			p.setColor(Color.RED);
			p.addPoints(new float[]{cx, cx}, new float[]{0f, 1f}, Plot.LINE);
			p.show();
		}

		return cx;

	}

	private static void accumulate(float[] values, float[] accumulator)
	{
		for (int i=0; i<N; i++) if (values[i]>accumulator[i]) accumulator[i] = values[i];
	}

	private float min(float in1, float in2)
	{
		return Math.min(in1, in2);
	}

	private float min(float in1, float in2, float in3)
	{
		return Math.min(in1, min(in2, in3));
	}

	private float min(float in1, float in2, float in3, float in4)
	{
		return Math.min(in1, min(in2, in3, in4));
	}

	private float min(float in1, float in2, float in3, float in4, float in5)
	{
		return Math.min(in1, min(in2, in3, in4, in5));
	}

	private float min(float in1, float in2, float in3, float in4, float in5, float in6)
	{
		return Math.min(in1, min(in2, in3, in4, in5, in6));
	}

	private float min(float in1, float in2, float in3, float in4, float in5, float in6, float in7)
	{
		return Math.min(in1, min(in2, in3, in4, in5, in6, in7));
	}

	/*
		demo
	 */
	public void demo(){

		float bkg = 40;

		// will make plots and export variables for plotting in R

		int demo_theta_size = (int) Math.round(3 * iDiff);
		float[] demo_theta = new float[demo_theta_size];
		float[] demo_low_theta = new float[demo_theta_size];
		float[] demo_mid_theta = new float[demo_theta_size];
		float[] demo_hgh_theta = new float[demo_theta_size];

		for (int i=0; i<demo_theta_size; i++) {
			demo_theta[i] 		= -iDiff+i;
			demo_low_theta[i] 	= h_low(demo_theta[i]);
			demo_mid_theta[i] 	= h_mid(demo_theta[i]);
			demo_hgh_theta[i] 	= h_hgh(demo_theta[i]);
		}

		// plots
		ImageStack plot_is = new ImageStack(528, 255);
		Plot p = new Plot("", "theta", "h[theta]", demo_theta, demo_low_theta);
		p.setLineWidth(3);
		p.draw();

		p.setColor(Color.BLUE);
		p.addPoints(demo_theta, demo_mid_theta, Plot.LINE);
		p.draw();

		p.setColor(Color.RED);
		p.addPoints(demo_theta, demo_hgh_theta, Plot.LINE);
		p.draw();

		plot_is.addSlice(p.getProcessor());
		//plot_is.addSlice(new Plot("", "theta", "h_mig[theta]", demo_theta, demo_mid_theta).getProcessor());
		//plot_is.addSlice(new Plot("", "theta", "h_hgh[theta]", demo_theta, demo_hgh_theta).getProcessor());

		Plot p1 = new Plot("", "x", "Q[x]", x, q_YES);
		p1.setLineWidth(3);
		p1.draw();

		p1.setColor(Color.BLUE);
		p1.addPoints(x, q_MAYBE, Plot.LINE);
		p1.draw();

		p1.setColor(Color.RED);
		p1.addPoints(x, q_NO, Plot.LINE);
		p1.draw();

		plot_is.addSlice(p1.getProcessor());
		//plot_is.addSlice(new Plot("", "x", "Q_MAYBE[x]", ).getProcessor());
		//plot_is.addSlice(new Plot("", "x", "Q_NO[x]", ).getProcessor());
		new ImagePlus("fuzzy_membership_mapping", plot_is).show();

		// file export
		PrintWriter writer = null;
		try {
			writer = new PrintWriter(exportFile);  writer.print("");   writer.close();
		} catch (FileNotFoundException ex) {}

		try {
			PrintWriter out;

			out = new PrintWriter(new BufferedWriter(new FileWriter(exportFile, true)));
			// fuzzy membership sets
			out.print("theta");			for (int idx=0; idx<demo_theta.length; idx++)    	out.print(" " + demo_theta[idx]);      out.print("\n");
			out.print("h_low[theta]"); 	for (int idx=0; idx<demo_low_theta.length; idx++)   out.print(" " + demo_low_theta[idx]);  out.print("\n");
			out.print("h_mid[theta]"); 	for (int idx=0; idx<demo_mid_theta.length; idx++)   out.print(" " + demo_mid_theta[idx]);  out.print("\n");
			out.print("h_hgh[theta]"); 	for (int idx=0; idx<demo_hgh_theta.length; idx++)   out.print(" " + demo_hgh_theta[idx]);  out.print("\n");
			// fuzzy decision sets
			out.print("x");				for (int idx=0; idx<x.length; idx++)    	out.print(" " + x[idx]);      		out.print("\n");
			out.print("Q_YES"); 		for (int idx=0; idx<q_YES.length; idx++)    out.print(" " + q_YES[idx]);  	out.print("\n");
			out.print("Q_MAYBE"); 		for (int idx=0; idx<q_MAYBE.length; idx++)  out.print(" " + q_MAYBE[idx]);  		out.print("\n");
			out.print("Q_NO"); 			for (int idx=0; idx<q_NO.length; idx++)		out.print(" " + q_NO[idx]);  	out.print("\n");
			//
			out.close();

		} catch (IOException e1) {}

		IJ.log("exported to : \n" + new File(exportFile).getAbsolutePath()+ " \n ");

	}

}
