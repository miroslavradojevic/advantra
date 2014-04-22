package detection2d;

import aux.Stat;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Plot;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.awt.*;
import java.util.Arrays;

/**
 * Created by miroslav on 4/10/14.
 */
public class Fuzzy2D {

    // two inputs ncc (average ncc along the streamline), and streamline likelihood
	// output 1 - saying whether one streamline (that takes ncc and likelihood) is ON or OFF
	// output 2 - saying whether up to three streamlines taken are END, NONE, or BIF

	// define the average ncc range (values below 0 are cut off)
	float ncc_start = 0;
    float ncc_end 	= 1;
	float ncc_high; // trapezoid
	float ncc_low;

    float likelihood_start 	= 0;
    float likelihood_end 	= 1;
	float likelihood_high; // trapezoid
	float likelihood_low;

	float 	output_sigma;

    // output1 can be off or on
	float 	output1_range_start;
	float 	output1_range_ended;
	float 	output1_off  = 0f;
	float 	output1_none = 1f;
	float   output1_on 	 = 2f;

	float 	output2_range_start;
	float 	output2_range_ended;
	float 	output2_endpoint = 1f;
	float   output2_nonpoint = 2f;
	float   output2_bifpoint = 3f;

    int      N1 = 100;         		// number of points to aggregate, hardcoded, need to cover the range well
	int		 N2 = N1;

    // fuzzification parameters
    float[] 	x1;      	    	// serve as x axis for membership, aggregation functions for output1 (range output1_range_start to output1_range_ended)
    float[] 	x2;      	    	// serve as x axis for membership, aggregation functions for output2 (range output2_range_start to output2_range_ended)
    public float[] 	agg1;            // serve as aggregation for output1
    public float[] 	agg2;            // serve as aggregation for output2

    // output categories
    float[] 	v_off;
	float[]		v_none;
    float[] 	v_on;
	float[]		v_endpoint;
	float[]		v_nonpoint;
	float[]		v_bifpoint;

    // membership fuzzy sets for output
    float[] 	q_off;
    float[] 	q_none;
    float[] 	q_on;
    float[] 	q_endpoint;
    float[] 	q_nonpoint;
    float[] 	q_bifpoint;

    public Fuzzy2D(
								float _ncc_high,
								float _ncc_low,
								float _likelihood_high,
								float _likelihood_low,
								float _output_sigma
								)
    {
        ncc_high 	= _ncc_high;
        ncc_low 	= _ncc_low;
		likelihood_high = _likelihood_high;
		likelihood_low 	= _likelihood_low;
		output_sigma = _output_sigma;
		output1_range_start = output1_off-3*output_sigma;
		output1_range_ended = output1_on+3*output_sigma;
		output2_range_start = output2_endpoint-3*output_sigma;
		output2_range_ended = output2_bifpoint+3*output_sigma;

        x1 = new float[N1]; // output1 range
		x2 = new float[N2]; // output2 range

        for (int i=0; i<N1; i++) x1[i] = output1_range_start + i * ( (output1_range_ended-output1_range_start) / (N1-1));
        for (int i=0; i<N2; i++) x2[i] = output2_range_start + i * ( (output2_range_ended-output2_range_start) / (N2-1));

        // output membership functions
        q_off = new float[N1];
        for (int i=0; i<N1; i++)
			q_off[i] = (float) Math.exp(-Math.pow(x1[i] - output1_off, 2) / (2 * Math.pow(output_sigma, 2)));

		q_none = new float[N1];
		for (int i=0; i<N1; i++)
			q_none[i] = (float) Math.exp(-Math.pow(x1[i] - output1_none, 2) / (2 * Math.pow(output_sigma, 2)));

		q_on = new float[N1];
        for (int i=0; i<N1; i++)
			q_on[i] = (float) Math.exp(-Math.pow(x1[i] - output1_on, 2) / (2 * Math.pow(output_sigma, 2)));

		q_endpoint = new float[N2];
		for (int i=0; i<N2; i++)
			q_endpoint[i] = (float) Math.exp(-Math.pow(x2[i] - output2_endpoint, 2) / (2 * Math.pow(output_sigma, 2)));

		q_nonpoint = new float[N2];
		for (int i=0; i<N2; i++)
			q_nonpoint[i] = (float) Math.exp(-Math.pow(x2[i] - output2_nonpoint, 2) / (2 * Math.pow(output_sigma, 2)));

		q_bifpoint = new float[N2];
		for (int i=0; i<N2; i++)
			q_bifpoint[i] = (float) Math.exp(-Math.pow(x2[i] - output2_bifpoint, 2) / (2 * Math.pow(output_sigma, 2)));

        // aux arrays used for calculations
        agg1 = new float[N1];
        agg2 = new float[N2];

        v_off = new float[N1];
		v_none = new float[N1];
        v_on  = new float[N1];

		v_endpoint = new float[N2];
		v_nonpoint = new float[N2];
		v_bifpoint = new float[N2];

    }

    /*
    fuzzification of the ncc input
     */
    private float h_ncc_high(float ncc)
    {
        if (ncc>=ncc_high)
            return 1;
        else if (ncc<ncc_high && ncc>=ncc_low)
			return (ncc-ncc_low)/(ncc_high-ncc_low);
		else
            return 0;
    }

	private float h_ncc_low(float ncc)
	{
		if (ncc<ncc_low)
			return 1;
		else if (ncc>=ncc_low && ncc<ncc_high)
			return (ncc-ncc_high)/(ncc_low-ncc_high);
		else
			return 0;
	}

    /*
    fuzzification of the likelihood input
     */
    private float h_likelihood_high(float likelihood)
    {
        if (likelihood>=likelihood_high)
            return 1;
        else if (likelihood<likelihood_high && likelihood>=likelihood_low)
            return (likelihood-likelihood_low)/(likelihood_high-likelihood_low);
		else
			return 0;
	}

	private float h_likelihood_low(float likelihood)
	{
		if (likelihood<likelihood_low)
			return 1;
		else if (likelihood>=likelihood_low && likelihood<likelihood_high)
			return (likelihood-likelihood_high)/(likelihood_low-likelihood_high);
		else
			return 0;
	}

    public void showFuzzification()
    {
        int NN = 51;

        // ncc input fuzzification
        float[] xx = new float[NN];
        float[] yy_on = new float[NN];
        float[] yy_off = new float[NN];

        for (int ii=0; ii<xx.length; ii++) xx[ii] = ncc_start + ii * ( (ncc_end-ncc_start) / (NN-1));
        // on
        for (int ii=0; ii<xx.length; ii++) yy_on[ii] = h_ncc_high(xx[ii]);
        //off
        for (int ii=0; ii<xx.length; ii++) yy_off[ii] = h_ncc_low(xx[ii]);

        Plot p = new Plot("ncc", "", "", xx, yy_on, Plot.LINE);
        p.setColor(Color.RED);
        p.draw();
        p.setColor(Color.BLUE);
        p.addPoints(xx, yy_off, Plot.LINE);
        p.show();

        // likelihood input fuzzification from 0 to 1
        for (int ii=0; ii<xx.length; ii++) xx[ii] = likelihood_start + ii * ( (likelihood_end-likelihood_start) / (NN-1));
        // on
        for (int ii=0; ii<xx.length; ii++) yy_on[ii] = h_likelihood_high(xx[ii]);
        //off
        for (int ii=0; ii<xx.length; ii++) yy_off[ii] = h_likelihood_low(xx[ii]);

        p = new Plot("likelihood", "", "", xx, yy_on, Plot.LINE);
        p.setColor(Color.RED);
        p.draw();
        p.setColor(Color.BLUE);
        p.addPoints(xx, yy_off, Plot.LINE);
        p.show();

    }

    public void showDefuzzification()
    {
        Plot p = new Plot("output1", "", "", x1, q_on, Plot.LINE);
        p.setColor(Color.RED);
        p.draw();
		p.setColor(Color.GREEN);
		p.addPoints(x1, q_none, Plot.LINE);
		p.draw();
        p.setColor(Color.BLUE);
        p.addPoints(x1, q_off, Plot.LINE);
        p.show();


		p = new Plot("output2", "", "", x2, q_bifpoint, Plot.LINE);
		p.setColor(Color.RED);
		p.draw();
		p.setColor(Color.GREEN);
		p.addPoints(x2, q_nonpoint, Plot.LINE);
		p.setColor(Color.BLUE);
		p.addPoints(x2, q_endpoint, Plot.LINE);
		p.show();

    }

//    public void showAgg1()
//    {
//        Plot p = new Plot("agg", "", "", x1, agg1, Plot.LINE);
//        p.show();
//    }

//	public void showAgg2()
//	{
//		Plot p = new Plot("agg", "", "", x2, agg2, Plot.LINE);
//		p.show();
//	}

//	public void showIt(float ncc_1, float likelihood_1)
//	{
//		float[] tst = new float[3];
//		critpointScore(ncc_1, likelihood_1, tst, true);
//		showAgg2();
//	}

//	public void showIt(float ncc_1, float likelihood_1, float ncc_2, float likelihood_2)
//	{
//		float[] tst = new float[3];
//		critpointScore(ncc_1, likelihood_1, ncc_2, likelihood_2, tst);
//		showAgg2();
//	}

//	public void showIt(float ncc_1, float likelihood_1, float ncc_2, float likelihood_2, float ncc_3, float likelihood_3)
//	{
//		float[] tst = new float[3];
//		critpointScore(ncc_1, likelihood_1, ncc_2, likelihood_2, ncc_3, likelihood_3, tst);
//		showAgg2();
//	}
//
//	public void showIt(float ncc_1, float likelihood_1, float ncc_2, float likelihood_2, float ncc_3, float likelihood_3, float ncc_4, float likelihood_4)
//	{
//		float[] tst = new float[3];
//		critpointScore(ncc_1, likelihood_1, ncc_2, likelihood_2, ncc_3, likelihood_3, ncc_4, likelihood_4, tst);
//		showAgg2();
//	}

	public void showDefuzzificationSurface(){

		int w = 200; 	// 0 / 1
		int h = w; 		// 0 / 1

		ImageStack is_out = new ImageStack(w, h);

		float[][] t = new float[w][h];
		float[][] t_off = new float[w][h];
		float[][] t_none = new float[w][h];
		float[][] t_on = new float[w][h];
		float[] dummy = new float[3];

		for (int x=0; x<w; x++) {

			for (int y=0; y<h; y++) {

				float ncc = ncc_start + x * ((ncc_end-ncc_start)/(w-1));
				float likelihood = likelihood_start + y * ((likelihood_end-likelihood_start)/(h-1));
				branchStrength(ncc, likelihood, dummy);
//				outputFuzzified(t[x][y], dummy);
				t[x][y] = branchStrength(ncc,likelihood);
				t_off[x][y] 	= dummy[0];
				t_none[x][y] 	= dummy[1];
				t_on[x][y] 		= dummy[2];

			}

		}

		is_out.addSlice("branch strength defuzzified", new FloatProcessor(t));
		is_out.addSlice("off", new FloatProcessor(t_off));
		is_out.addSlice("none", new FloatProcessor(t_none));
		is_out.addSlice("on", new FloatProcessor(t_on));

		new ImagePlus("output", is_out).show();

	}

    private float[] fi_off(float z)
    {
        Arrays.fill(v_off, 0);

        for (int i=0; i<N1; i++)
            if (q_off[i]<=z)
                v_off[i] = q_off[i];
            else
                v_off[i] = z;

        return v_off;
    }

	private float[] fi_none(float z)
	{
		Arrays.fill(v_none, 0);

		for (int i=0; i<N1; i++)
			if (q_none[i]<=z)
				v_none[i] = q_none[i];
			else
				v_none[i] = z;

		return v_none;
	}

    private float[] fi_on(float z)
    {
        Arrays.fill(v_on, 0);

        for (int i=0; i<N1; i++)
            if (q_on[i]<=z)
                v_on[i] = q_on[i];
            else
                v_on[i] = z;

        return v_on;
    }

	private float[] fi_endpoint(float z)
	{
		Arrays.fill(v_endpoint, 0);

		for (int i=0; i<N2; i++)
			if (q_endpoint[i]<=z)
				v_endpoint[i] = q_endpoint[i];
			else
				v_endpoint[i] = z;

		return v_endpoint;
	}

	private float[] fi_nonpoint(float z)
	{
		Arrays.fill(v_nonpoint, 0);

		for (int i=0; i<N2; i++)
			if (q_nonpoint[i]<=z)
				v_nonpoint[i] = q_nonpoint[i];
			else
				v_nonpoint[i] = z;
		return v_nonpoint;
	}

	private float[] fi_bifpoint(float z)
	{
		Arrays.fill(v_bifpoint, 0);

		for (int i=0; i<N2; i++)
			if (q_bifpoint[i]<=z)
				v_bifpoint[i] = q_bifpoint[i];
			else
				v_bifpoint[i] = z;

		return v_bifpoint;
	}

    /*
    *
    *
    *
    *
    *
     */

    private float branchStrength(float ncc_in, float likelihood_in)
    {
        Arrays.fill(agg1, 0);
        float[] cur;
        float mu;

        // apply rules
		mu = min(h_ncc_low(ncc_in), 	h_likelihood_low(likelihood_in));	cur = fi_off(mu);	accumulate(cur, agg1);
		mu = min(h_ncc_high(ncc_in),	h_likelihood_high(likelihood_in)); 	cur = fi_on(mu);  	accumulate(cur, agg1);
		mu = min(h_ncc_low(ncc_in),		h_likelihood_high(likelihood_in)); 	cur = fi_none(mu);  accumulate(cur, agg1);
        // finished rules

        return get_agg1_centroid();
    }

	public void branchStrength(float ncc_in, float likelihood_in, float[] output_off_none_on) // output_off_none_on: [off_score, none_score, on_score]
	{
		float defuzz = branchStrength(ncc_in, likelihood_in);
		// use the same formulas that were used for q_on/q_off
		output_off_none_on[0] = (float) Math.exp(-Math.pow(defuzz - output1_off, 2) / (2 * Math.pow(output_sigma, 2)));
		output_off_none_on[1] = (float) Math.exp(-Math.pow(defuzz - output1_none, 2) / (2 * Math.pow(output_sigma, 2)));
		output_off_none_on[2] = (float) Math.exp(-Math.pow(defuzz - output1_on, 2) / (2 * Math.pow(output_sigma, 2)));
	}

	/*
    *
    *
    *
    *
    *
     */

	private float critpointScore(float ncc_1, float likelihood_1)
	{

		float[] t = new float[3];

		branchStrength(ncc_1, likelihood_1, t);
		float b1_is_off 	= t[0];
		float b1_is_none 	= t[1];
		float b1_is_on 		= t[2];

        Arrays.fill(agg2, 0);
		float mu;
		float[] cur;

        // rules
		mu = b1_is_off;     cur = fi_nonpoint(mu);   	accumulate(cur, agg2);
		mu = b1_is_none;    cur = fi_nonpoint(mu);   	accumulate(cur, agg2);
		mu = b1_is_on;      cur = fi_endpoint(mu);    	accumulate(cur, agg2);
        return get_agg2_centroid();

	}

    public ImageProcessor critpointScore(float ncc_1, float likelihood_1, float[] output_endpoint_nonpoint_bifpoint, boolean verbose)
    {
        float defuzz = critpointScore(ncc_1, likelihood_1);
		output_endpoint_nonpoint_bifpoint[0] = (float) Math.exp(-Math.pow(defuzz - output2_endpoint, 2) / (2 * Math.pow(output_sigma, 2)));
		output_endpoint_nonpoint_bifpoint[1] = (float) Math.exp(-Math.pow(defuzz - output2_nonpoint, 2) / (2 * Math.pow(output_sigma, 2)));
		output_endpoint_nonpoint_bifpoint[2] = (float) Math.exp(-Math.pow(defuzz - output2_bifpoint, 2) / (2 * Math.pow(output_sigma, 2)));
		if (verbose) {
			Plot p = new Plot("agg", "", "", x2, agg2, Plot.LINE);
			return p.getProcessor();
		}
		else return null;
    }

	private float critpointScore(
										float ncc_1, float likelihood_1,
									 	float ncc_2, float likelihood_2)
	{

        float[] t = new float[3];

        branchStrength(ncc_1, likelihood_1, t);
        float b1_is_off 	= t[0];
        float b1_is_none 	= t[1];
        float b1_is_on 		= t[2];

        branchStrength(ncc_2, likelihood_2, t);
        float b2_is_off 	= t[0];
        float b2_is_none 	= t[1];
        float b2_is_on 		= t[2];

        Arrays.fill(agg2, 0);
        float mu;
        float[] cur;

        // rules
        mu = min(b1_is_off, 	b2_is_off);     cur = fi_nonpoint(mu); 	accumulate(cur, agg2);
		mu = min(b1_is_off, 	b2_is_on);      cur = fi_endpoint(mu); 	accumulate(cur, agg2);
		mu = min(b1_is_on, 		b2_is_off);     cur = fi_endpoint(mu); 	accumulate(cur, agg2);
		mu = min(b1_is_on, 		b2_is_on);      cur = fi_nonpoint(mu); 	accumulate(cur, agg2);
		mu = max(b1_is_none, 	b2_is_none);	cur = fi_nonpoint(mu); 	accumulate(cur, agg2);

        return get_agg2_centroid();

    }

    public ImageProcessor critpointScore(
									float ncc_1, float likelihood_1,
                                  	float ncc_2, float likelihood_2, float[] output_endpoint_nonpoint_bifpoint, boolean verbose)
    {
        float defuzz = critpointScore(ncc_1, likelihood_1, ncc_2, likelihood_2);
		output_endpoint_nonpoint_bifpoint[0] = (float) Math.exp(-Math.pow(defuzz - output2_endpoint, 2) / (2 * Math.pow(output_sigma, 2)));
		output_endpoint_nonpoint_bifpoint[1] = (float) Math.exp(-Math.pow(defuzz - output2_nonpoint, 2) / (2 * Math.pow(output_sigma, 2)));
		output_endpoint_nonpoint_bifpoint[2] = (float) Math.exp(-Math.pow(defuzz - output2_bifpoint, 2) / (2 * Math.pow(output_sigma, 2)));
		if (verbose) {
			Plot p = new Plot("agg", "", "", x2, agg2, Plot.LINE);
			return p.getProcessor();
		}
		else return null;
    }

    private float critpointScore(	float ncc_1, float likelihood_1,
                                    float ncc_2, float likelihood_2,
                                    float ncc_3, float likelihood_3)
    {

        float[] t = new float[3];

        branchStrength(ncc_1, likelihood_1, t);
        float b1_is_off 	= t[0];
        float b1_is_none 	= t[1];
        float b1_is_on 		= t[2];

        branchStrength(ncc_2, likelihood_2, t);
        float b2_is_off 	= t[0];
        float b2_is_none 	= t[1];
        float b2_is_on 		= t[2];

        branchStrength(ncc_3, likelihood_3, t);
        float b3_is_off 	= t[0];
        float b3_is_none 	= t[1];
        float b3_is_on 		= t[2];


        Arrays.fill(agg2, 0);
        float mu;
        float[] cur;

        // rules
        mu = min(b1_is_off, b2_is_off, b3_is_off);     cur = fi_nonpoint(mu); accumulate(cur, agg2);
		mu = min(b1_is_off, b2_is_off, b3_is_on);      cur = fi_endpoint(mu); accumulate(cur, agg2);
		mu = min(b1_is_off, b2_is_on,  b3_is_off);     cur = fi_endpoint(mu); accumulate(cur, agg2);
		mu = min(b1_is_off, b2_is_on,  b3_is_on);      cur = fi_nonpoint(mu); accumulate(cur, agg2);

		mu = min(b1_is_on, b2_is_off, b3_is_off);      cur = fi_endpoint(mu); accumulate(cur, agg2);
        mu = min(b1_is_on, b2_is_off, b3_is_on);       cur = fi_nonpoint(mu); accumulate(cur, agg2);
        mu = min(b1_is_on, b2_is_on,  b3_is_off);      cur = fi_nonpoint(mu); accumulate(cur, agg2);
        mu = min(b1_is_on, b2_is_on,  b3_is_on);       cur = fi_bifpoint(mu); accumulate(cur, agg2);

		mu = max(b1_is_none, b2_is_none, b3_is_none);  cur = fi_nonpoint(mu); accumulate(cur, agg2);

        return get_agg2_centroid();

    }

    public ImageProcessor critpointScore(    float ncc_1, float likelihood_1,
                                   float ncc_2, float likelihood_2,
                                   float ncc_3, float likelihood_3, float[] output_endpoint_nonpoint_bifpoint, boolean verbose)
    {
        float defuzz = critpointScore(ncc_1, likelihood_1, ncc_2, likelihood_2, ncc_3, likelihood_3);
		output_endpoint_nonpoint_bifpoint[0] = (float) Math.exp(-Math.pow(defuzz - output2_endpoint, 2) / (2 * Math.pow(output_sigma, 2)));
		output_endpoint_nonpoint_bifpoint[1] = (float) Math.exp(-Math.pow(defuzz - output2_nonpoint, 2) / (2 * Math.pow(output_sigma, 2)));
		output_endpoint_nonpoint_bifpoint[2] = (float) Math.exp(-Math.pow(defuzz - output2_bifpoint, 2) / (2 * Math.pow(output_sigma, 2)));
		if (verbose) {
			Plot p = new Plot("agg", "", "", x2, agg2, Plot.LINE);
			return p.getProcessor();
		}
		else return null;
    }

    private float critpointScore(    float ncc_1, float likelihood_1,
                                     float ncc_2, float likelihood_2,
                                     float ncc_3, float likelihood_3,
                                     float ncc_4, float likelihood_4)
    {

        float[] t = new float[3];

        branchStrength(ncc_1, likelihood_1, t);
        float b1_is_off 	= t[0];
        float b1_is_none 	= t[1];
        float b1_is_on 		= t[2];

        branchStrength(ncc_2, likelihood_2, t);
        float b2_is_off 	= t[0];
        float b2_is_none 	= t[1];
        float b2_is_on 		= t[2];

		branchStrength(ncc_3, likelihood_3, t);
        float b3_is_off 	= t[0];
        float b3_is_none 	= t[1];
        float b3_is_on 		= t[2];

        branchStrength(ncc_4, likelihood_4, t);
        float b4_is_off 	= t[0];
        float b4_is_none 	= t[1];
        float b4_is_on 		= t[2];

        Arrays.fill(agg2, 0);
        float mu;
        float[] cur;

        // rules
        mu = min(b1_is_off, 	b2_is_off, 	b3_is_off, 	b4_is_off);      cur = fi_nonpoint(mu);  accumulate(cur, agg2); // 0000
		mu = min(b1_is_off, 	b2_is_off, 	b3_is_off, 	b4_is_on);       cur = fi_endpoint(mu);  accumulate(cur, agg2); // 0001 (end)
		mu = min(b1_is_off, 	b2_is_off, 	b3_is_on, 	b4_is_off);      cur = fi_endpoint(mu);  accumulate(cur, agg2); // 0010 (end)
		mu = min(b1_is_off, 	b2_is_off, 	b3_is_on,  	b4_is_on);       cur = fi_nonpoint(mu);  accumulate(cur, agg2); // 0011

		mu = min(b1_is_off, 	b2_is_on,   b3_is_off,  b4_is_off);      cur = fi_endpoint(mu);  accumulate(cur, agg2); // 0100 (end)
		mu = min(b1_is_off, 	b2_is_on, 	b3_is_off,  b4_is_on);       cur = fi_nonpoint(mu);  accumulate(cur, agg2); // 0101
		mu = min(b1_is_off, 	b2_is_on, 	b3_is_on,  	b4_is_off);      cur = fi_nonpoint(mu);  accumulate(cur, agg2); // 0110
		mu = min(b1_is_off, 	b2_is_on, 	b3_is_on,  	b4_is_on);       cur = fi_bifpoint(mu);  accumulate(cur, agg2); // 0111 (bif)

		mu = min(b1_is_on, 		b2_is_off, 	b3_is_off, 	b4_is_off);      cur = fi_endpoint(mu);  accumulate(cur, agg2); // 1000 (end)
		mu = min(b1_is_on,  	b2_is_off,  b3_is_off,  b4_is_on);       cur = fi_nonpoint(mu);  accumulate(cur, agg2); // 1001
        mu = min(b1_is_on,  	b2_is_off, 	b3_is_on,  	b4_is_off);      cur = fi_nonpoint(mu);  accumulate(cur, agg2); // 1010
        mu = min(b1_is_on,  	b2_is_off, 	b3_is_on,  	b4_is_on);       cur = fi_bifpoint(mu);  accumulate(cur, agg2); // 1011 (bif)

		mu = min(b1_is_on,  	b2_is_on, 	b3_is_off,  b4_is_off);      cur = fi_nonpoint(mu);  accumulate(cur, agg2); // 1100
        mu = min(b1_is_on,  	b2_is_on, 	b3_is_off,  b4_is_on);       cur = fi_bifpoint(mu);  accumulate(cur, agg2); // 1101 (bif)
        mu = min(b1_is_on,  	b2_is_on, 	b3_is_on,  	b4_is_off);      cur = fi_bifpoint(mu);  accumulate(cur, agg2); // 1110 (bif)
        mu = min(b1_is_on,  	b2_is_on, 	b3_is_on,  	b4_is_on);       cur = fi_bifpoint(mu);  accumulate(cur, agg2); // 1111 (bif)

		mu = min(b1_is_none,	b2_is_on, 	b3_is_on, 	b4_is_on);		 cur = fi_bifpoint(mu);  accumulate(cur, agg2); // x111 (bif)
		mu = min(b1_is_on,		b2_is_none, b3_is_on, 	b4_is_on);		 cur = fi_bifpoint(mu);  accumulate(cur, agg2); // 1x11 (bif)
		mu = min(b1_is_on,		b2_is_on, 	b3_is_none, b4_is_on);		 cur = fi_bifpoint(mu);  accumulate(cur, agg2); // 11x1 (bif)
		mu = min(b1_is_on,		b2_is_on, 	b3_is_on,   b4_is_none);     cur = fi_bifpoint(mu);  accumulate(cur, agg2); // 111x (bif)

		mu = min(b1_is_none,	1-b2_is_on, 1-b3_is_on,  1-b4_is_on);    cur = fi_nonpoint(mu);  accumulate(cur, agg2); // x1_1_1_
		mu = min(1-b1_is_on,	b2_is_none, 1-b3_is_on,  1-b4_is_on);    cur = fi_nonpoint(mu);  accumulate(cur, agg2); // 1_x1_1_
		mu = min(1-b1_is_on,	1-b2_is_on, b3_is_none,  1-b4_is_on);    cur = fi_nonpoint(mu);  accumulate(cur, agg2); // 1_1_x1_
		mu = min(1-b1_is_on,	1-b2_is_on, 1-b3_is_on,  b4_is_none);    cur = fi_nonpoint(mu);  accumulate(cur, agg2); // 1_1_1_x

        return get_agg2_centroid();

    }

    public ImageProcessor critpointScore(	   float ncc_1, float likelihood_1,
                                   float ncc_2, float likelihood_2,
                                   float ncc_3, float likelihood_3,
                                   float ncc_4, float likelihood_4, float[] output_endpoint_nonpoint_bifpoint, boolean verbose)
    {
        float defuzz = critpointScore(ncc_1, likelihood_1, ncc_2, likelihood_2, ncc_3, likelihood_3, ncc_4, likelihood_4);
		output_endpoint_nonpoint_bifpoint[0] = (float) Math.exp(-Math.pow(defuzz - output2_endpoint, 2) / (2 * Math.pow(output_sigma, 2)));
		output_endpoint_nonpoint_bifpoint[1] = (float) Math.exp(-Math.pow(defuzz - output2_nonpoint, 2) / (2 * Math.pow(output_sigma, 2)));
		output_endpoint_nonpoint_bifpoint[2] = (float) Math.exp(-Math.pow(defuzz - output2_bifpoint, 2) / (2 * Math.pow(output_sigma, 2)));
		if (verbose) {
			Plot p = new Plot("agg", "", "", x2, agg2, Plot.LINE);
			return p.getProcessor();
		}
		else return null;

    }

    private float get_agg1_centroid()
    {
        float centroid = 0;
        float norm = 0;
        for (int i=0; i<N1; i++) {
            centroid += x1[i]*agg1[i];
            norm += agg1[i];
        }

        if (norm>Float.MIN_VALUE) return centroid / norm;
        else return Stat.average(x1);

    }

	private float get_agg2_centroid()
	{
		float centroid = 0;
		float norm = 0;
		for (int i=0; i<N2; i++) {
			centroid += x2[i]*agg2[i];
			norm += agg2[i];
		}

		if (norm>Float.MIN_VALUE) return centroid / norm;
		else return Stat.average(x2);

	}

    private static final float min(float in1, float in2)
    {
        return Math.min(in1, in2);
    }

    private static final float min(float in1, float in2, float in3)
    {
        return Math.min(in1, min(in2, in3));
    }

    private static final float min(float in1, float in2, float in3, float in4)
    {
        return Math.min(in1, min(in2, in3, in4));
    }


    private static final void accumulate(float[] values, float[] accumulator)
    {
        for (int i=0; i<values.length; i++) if (values[i]>accumulator[i]) accumulator[i] = values[i];
    }

	private static final float max(float in1, float in2)
	{
		return Math.max(in1, in2);
	}

	private static final float max(float in1, float in2, float in3)
	{
		return Math.max(in1, max(in2, in3));
	}

	private static final float max(float in1, float in2, float in3, float in4)
	{
		return Math.max(in1, max(in2, in3, in4));
	}

}
