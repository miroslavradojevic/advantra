package detection2d;

import ij.gui.Plot;

import java.awt.*;
import java.util.Arrays;

/**
 * Created by miroslav on 4/10/14.
 */
public class FuzzyBranch2D {

    // two inputs of different nature, with different parameters: ncc, likelihood
    float ncc_start = -1; // define the range
    float ncc_end = +1;
    float    mu_ncc, sig_ncc;

    float start_likelihood = 0;
    float end_likelihood = 1;
    float    sig_likelihood, mu_likelihood;

    // one output: strength
    float    mu_on=1f, mu_off=0f, sig_on=0.3f, sig_off=0.3f;
    float    ;

    int      N;         // number of points to aggregate
//    int      L = 1;		// number of outputs - saying how much it is on or off
    int[]    out_idxs;  // indexes of out points

    // fuzzification parameters
         // high
    // todo add low separate parameters so that low has it's own function


    float[] 	x;      	    // serve as x axis for membership functions (range 0 to 4)
    public float[] 	agg;            // serve as aggregation

    // L outputs
    float[] 	v_off;
    float[] 	v_on;

    // membership fuzzy sets for L outputs (rules will call these categories "OFF", "ON")
    float[] q_off;
    float[] q_on;

    public FuzzyBranch2D(int _N, float _mu_ncc, float _sig_ncc, float _mu_likelihood, float _sig_likelihood)
    {
        N = _N;

        mu_likelihood = _mu_likelihood;
        sig_likelihood = _sig_likelihood;

        mu_ncc = _mu_ncc;
        sig_ncc = _sig_ncc;

        x = new float[N]; // output range mu_off-3*sig_off, mu_on+3*sig_on
        for (int i=0; i<N; i++) x[i] = (mu_off-3*sig_off) + i * ((float)(L-1) / (N-1));

        System.out.println(Arrays.toString(x));

        out_idxs = new int[L];
        int step = (N-1)/(L-1);
        for (int i=0; i<L; i++) {
            out_idxs[i] = i*step;
            System.out.println("idx " + out_idxs[i]);
        }

        // output membership functions
        q_off = new float[N];
        for (int i=0; i<N; i++) {
            q_off[i] = (float) Math.exp(-Math.pow(x[i] - mu_off, 2)/(2*Math.pow(sig_off, 2)));
        }

        q_on = new float[N];
        for (int i=0; i<N; i++) {
            q_on[i] = (float) Math.exp(-Math.pow(x[i] - mu_on, 2)/(2*Math.pow(sig_on, 2)));
        }

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
        if (ncc>=mu_ncc)
            return 1;
        else //if(theta>mu_ON) //  && theta<=mu_OFF
            return (float) Math.exp(-Math.pow(ncc - mu_ncc, 2)/(2*Math.pow(sig_ncc, 2)));
        //else
        //    return 0;
    }

    /*
    fuzzification of the likelihood input
     */
    private float h_likelihood_high(float likelihood)
    {
        if (likelihood>=mu_likelihood)
            return 1;
        else //if(theta>mu_ON) //  && theta<=mu_OFF
            return (float) Math.exp(-Math.pow(likelihood - mu_likelihood, 2)/(2*Math.pow(sig_likelihood, 2)));
        //else
        //    return 0;
    }

    public void showFuzzification()
    {
        int NN = 51;

        // ncc input fuzzification
        float[] xx = new float[NN];
        float[] yy_on = new float[NN];
        float[] yy_off = new float[NN];

        for (int ii=0; ii<xx.length; ii++) xx[ii] = -1 + ii * ( (+1f-(-1)) / (NN-1));

        // on
        for (int ii=0; ii<xx.length; ii++) yy_on[ii] = h_ncc_high(xx[ii]);
        //off
        for (int ii=0; ii<xx.length; ii++) yy_off[ii] = 1-h_ncc_high(xx[ii]);

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
        for (int ii=0; ii<xx.length; ii++) yy_off[ii] = 1-h_likelihood_high(xx[ii]);

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

    public void critpointScores(float ncc_in, float likelihood_in, float[] out)
    {
        Arrays.fill(agg, 0);
        float[] cur;
        float mu;

        // apply rules
        mu = min (1-h_ncc_high(ncc_in),   1-h_likelihood_high(likelihood_in)  ) ; cur = fi_off(mu); accumulate(cur, agg);
        System.out.println("00:"+mu); System.out.println("ncc :"+h_ncc_high(ncc_in)+" likelihood "+h_likelihood_high(likelihood_in));
        mu = min (1-h_ncc_high(ncc_in),   h_likelihood_high(likelihood_in)    ) ; cur = fi_off(mu); accumulate(cur, agg); System.out.println("01:" + mu);

        mu = min (h_ncc_high(ncc_in),     1-h_likelihood_high(likelihood_in)  ) ; cur = fi_off(mu); accumulate(cur, agg); System.out.println("10:" + mu);
        mu = min (h_ncc_high(ncc_in),     h_likelihood_high(likelihood_in)    ) ; cur = fi_on(mu);  accumulate(cur, agg); System.out.println("11:" + mu);
        // finished rules

        for (int i=0; i<L; i++) out[i] = agg[out_idxs[i]];

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
