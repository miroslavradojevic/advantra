package detection2d;

import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Plot;
import ij.process.FloatProcessor;

import java.awt.*;
import java.util.Arrays;

/**
 * Created by miroslav on 4/10/14.
 */
public class FuzzyBranch2D {

    // two inputs ncc, likelihood, can be high or low
    float ncc_start = -1; // define the range
    float ncc_end = +1;
    float ncc_high_mean;
	float ncc_high_sigma;
	float ncc_low_mean;
	float ncc_low_sigma;

    float likelihood_start = 0;
    float likelihood_end = 1;
	float likelihood_high_mean;
	float likelihood_high_sigma;
	float likelihood_low_mean;
	float likelihood_low_sigma;

    // one output can be off or on
	float   output_on_mean = 1f;
	float 	output_on_sigma;
	float 	output_off_mean = 0f;
	float	output_off_sigma;
	float 	output_start;// = output_off_mean-3*output_off_sigma;
	float 	output_end;// = output_on_mean+3*output_on_sigma;

    int      N;         // number of points to aggregate

    // fuzzification parameters
    float[] 	x;      	    // serve as x axis for membership, aggregation functions (range output_start to output_end)
    public float[] 	agg;            // serve as aggregation

    // output categories
    float[] 	v_off;
    float[] 	v_on;

    // membership fuzzy sets for output (rules will call these categories "OFF", "ON")
    float[] q_off;
    float[] q_on;

    public FuzzyBranch2D(
								int _N,
								float _ncc_high_mean,
								float _ncc_high_sigma,
								float _ncc_low_mean,
								float _ncc_low_sigma,
								float _likelihood_high_mean,
								float _likelihood_high_sigma,
								float _likelihood_low_mean,
								float _likelihood_low_sigma,
								float _output_on_sigma
								)
    {
        N = _N;
        ncc_high_mean = _ncc_high_mean;
        ncc_high_sigma = _ncc_high_sigma;
        ncc_low_mean = _ncc_low_mean;
        ncc_low_sigma = _ncc_low_sigma;
		likelihood_high_mean = _likelihood_high_mean;
		likelihood_high_sigma = _likelihood_high_sigma;
		likelihood_low_mean = _likelihood_low_mean;
		likelihood_low_sigma = _likelihood_low_sigma;
		output_on_sigma = _output_on_sigma;
		output_off_sigma = output_on_sigma;
		output_start = output_off_mean-3*output_off_sigma;
		output_end = output_on_mean+3*output_on_sigma;

        x = new float[N]; // output range mu_off-3*sig_off, mu_on+3*sig_on , also used in centroid calculation
        for (int i=0; i<N; i++) x[i] = output_start + i * ( (output_end-output_start) / (N-1));

        // output membership functions
        q_off = new float[N];
        for (int i=0; i<N; i++)
			q_off[i] = (float) Math.exp(-Math.pow(x[i] - output_off_mean, 2) / (2 * Math.pow(output_off_sigma, 2)));

        q_on = new float[N];
        for (int i=0; i<N; i++)
			q_on[i] = (float) Math.exp(-Math.pow(x[i] - output_on_mean, 2) / (2 * Math.pow(output_on_sigma, 2)));

        // aux arrays used for calculations
        agg = new float[N];
        v_off = new float[N];
        v_on = new float[N];

    }

    /*
    fuzzification of the ncc input
     */
    private float h_ncc_high(float ncc)
    {
        if (ncc>=ncc_high_mean)
            return 1;
        else
            return (float) Math.exp(-Math.pow(ncc - ncc_high_mean, 2)/(2*Math.pow(ncc_high_sigma, 2)));
    }

	private float h_ncc_low(float ncc)
	{
		if (ncc<=ncc_low_mean)
			return 1;
		else
			return (float) Math.exp(-Math.pow(ncc - ncc_low_mean, 2)/(2*Math.pow(ncc_low_sigma, 2)));
	}

    /*
    fuzzification of the likelihood input
     */
    private float h_likelihood_high(float likelihood)
    {
        if (likelihood>=likelihood_high_mean)
            return 1;
        else
            return (float) Math.exp(-Math.pow(likelihood - likelihood_high_mean, 2)/(2*Math.pow(likelihood_high_sigma, 2)));
    }

	private float h_likelihood_low(float likelihood)
	{
		if (likelihood<=likelihood_low_mean)
			return 1;
		else
			return (float) Math.exp(-Math.pow(likelihood - likelihood_low_mean, 2)/(2*Math.pow(likelihood_low_sigma, 2)));
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

        Plot p = new Plot("ncc fuzzification", "", "", xx, yy_on, Plot.LINE);
        p.setColor(Color.RED);
        p.draw();
        p.setColor(Color.BLUE);
        p.addPoints(xx, yy_off, Plot.LINE);
        p.show();

        // likelihood input fuzzification from 0 to 1
        for (int ii=0; ii<xx.length; ii++) xx[ii] = 0 + ii * ( (+1f-(0)) / (NN-1));

        // on
        for (int ii=0; ii<xx.length; ii++) yy_on[ii] = h_likelihood_high(xx[ii]);
        //off
        for (int ii=0; ii<xx.length; ii++) yy_off[ii] = h_likelihood_low(xx[ii]);

        p = new Plot("likelihood fuzzification", "", "", xx, yy_on, Plot.LINE);
        p.setColor(Color.RED);
        p.draw();
        p.setColor(Color.BLUE);
        p.addPoints(xx, yy_off, Plot.LINE);
        p.show();

    }

    public void showDefuzzification()
    {
        Plot p = new Plot("defuzzification", "", "", x, q_on, Plot.LINE);
        p.setColor(Color.RED);
        p.draw();
        p.setColor(Color.BLUE);
        p.addPoints(x, q_off, Plot.LINE);
        p.show();
    }

    public void showAgg()
    {
        Plot p = new Plot("agg", "", "", x, agg, Plot.LINE);
        p.show();
    }

	public void showDefuzzificationSurface(){

		int w = 200; // -1 / 1
		int h = 100; //  0 / 1

		ImageStack is_out = new ImageStack(w, h);

		float[][] t = new float[w][h];
		float[][] t_off = new float[w][h];
		float[][] t_on = new float[w][h];
		float[] dummy = new float[2];

		for (int x=0; x<w; x++) {

			for (int y=0; y<h; y++) {

				float ncc = ncc_start + x * ((ncc_end-ncc_start)/(w-1));
				float likelihood = likelihood_start + y * ((likelihood_end-likelihood_start)/(h-1));
				t[x][y] = branchStrengthDefuzzified(ncc, likelihood);
				branchStrengthFuzzified(t[x][y], dummy);
				t_off[x][y] = dummy[0];
				t_on[x][y] = dummy[1];

			}

		}

		is_out.addSlice("branch strength defuzzified", new FloatProcessor(t));
		is_out.addSlice("branch strength off", new FloatProcessor(t_off));
		is_out.addSlice("branch strength on", new FloatProcessor(t_on));

		new ImagePlus("output", is_out).show();

	}

    private float[] fi_off(float z)
    {
        Arrays.fill(v_off, 0);

        for (int i=0; i<N; i++)
            if (q_off[i]<=z)
                v_off[i] = q_off[i];
            else
                v_off[i] = z;

        return v_off;
    }

    private float[] fi_on(float z)
    {
        Arrays.fill(v_on, 0);

        for (int i=0; i<N; i++)
            if (q_on[i]<=z)
                v_on[i] = q_on[i];
            else
                v_on[i] = z;

        return v_on;
    }

    public float branchStrengthDefuzzified(float ncc_in, float likelihood_in)
    {
        Arrays.fill(agg, 0);
        float[] cur;
        float mu;

        // apply rules
//        mu = min (h_ncc_low(ncc_in),   		h_likelihood_low(likelihood_in)  ) ; cur = fi_off(mu); accumulate(cur, agg);
//        mu = min (h_ncc_low(ncc_in),   		h_likelihood_high(likelihood_in) ) ; cur = fi_off(mu); accumulate(cur, agg);
//        mu = min (h_ncc_high(ncc_in),     	h_likelihood_low(likelihood_in)  ) ; cur = fi_off(mu); accumulate(cur, agg);

		mu = h_ncc_low(ncc_in) ; cur = fi_off(mu); accumulate(cur, agg);
		mu = h_likelihood_low(likelihood_in); cur = fi_off(mu); accumulate(cur, agg);
		mu = min (h_ncc_high(ncc_in),     	h_likelihood_high(likelihood_in) ) ; cur = fi_on(mu);  accumulate(cur, agg);  // and
        // finished rules

		float centroid = 0;
		float norm = 0;
		for (int i=0; i<N; i++) {
			centroid += x[i]*agg[i];
			norm += agg[i];
		}

		return centroid/norm;

    }

	public void branchStrengthFuzzified(float ncc_in, float likelihood_in, float[] output_off_on) // [off_score, on_score]
	{
		float defuzz = branchStrengthDefuzzified(ncc_in, likelihood_in);
		// use the same formulas that were used for q_on/q_off
		output_off_on[0] = (float) Math.exp(-Math.pow(defuzz - output_off_mean, 2) / (2 * Math.pow(output_off_sigma, 2)));
		output_off_on[1] = (float) Math.exp(-Math.pow(defuzz - output_on_mean, 2) / (2 * Math.pow(output_on_sigma, 2)));
	}

	public void branchStrengthFuzzified(float _defuzz, float[] output_off_on) // [off_score, on_score]
	{
//		float defuzz = branchStrengthDefuzzified(ncc_in, likelihood_in);
		// use the same formulas that were used for q_on/q_off
		output_off_on[0] = (float) Math.exp(-Math.pow(_defuzz - output_off_mean, 2) / (2 * Math.pow(output_off_sigma, 2)));
		output_off_on[1] = (float) Math.exp(-Math.pow(_defuzz - output_on_mean, 2) / (2 * Math.pow(output_on_sigma, 2)));
	}

	public float endpointDefuzzified(float ncc_1, float likelihood_1)
	{

		float[] t = new float[2];

		branchStrengthFuzzified(ncc_1, likelihood_1, t);
		float b1_is_off = t[0];
		float b1_is_on = t[1];

		// rules if there is one
		float mu;
		float[] cur;

		mu = b1_is_off; cur = fi_off(mu); accumulate(cur, agg);
		mu = b1_is_on; cur = fi_on(mu); accumulate(cur, agg);

		float centroid = 0;
		float norm = 0;
		for (int i=0; i<N; i++) {
			centroid += x[i]*agg[i];
			norm += agg[i];
		}

		return centroid/norm;

	}

	public float endpointDefuzzified(float ncc_1, float likelihood_1,
									 float ncc_2, float likelihood_2)
	{

	}



    public static final float min(float in1, float in2)
    {
        return Math.min(in1, in2);
    }

    private static final void accumulate(float[] values, float[] accumulator)
    {
        for (int i=0; i<values.length; i++) if (values[i]>accumulator[i]) accumulator[i] = values[i];
    }

}
