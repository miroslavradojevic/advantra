package detection2d;

import conn.AmiraParameters;
import ij.gui.Plot;

import java.awt.*;
import java.util.Arrays;

/**
 * Created by miroslav on 3/7/14.
 * there can be 1,2,3 or 4 inputs, every input is one dimensional and graded as ON or OFF
 * parameters of gradation are given as arguments muON, sigmaON, muOFF, sigmaOFF (for one variable)
 *
 */
public class Fuzzy2D {

    static int      N;          // number of points
    static int[]    out_idxs;   //
    static int      L = 5;      // number of outputs

    // fuzzification parameters
    static float    mu_ON, sig_ON, mu_OFF, sig_OFF;     // gaussian standard deviation
    // defuzzification parameters
    static float    sig_OUT = 0.25f;
    static float    mu_NON = 0f;
    static float    mu_END = 1f;
    static float    mu_BDY = 2f;
    static float    mu_BIF = 3f;
    static float    mu_CRS = 4f;

    static float[] 	x;      	    // serve as x axis for membership functions (range 0 to 4)
    static float[] 	agg;            // serve as aggregation

    // L outputs
    static float[] 	v_NON;
    static float[] 	v_END;
    static float[] 	v_BDY;
    static float[] 	v_BIF;
    static float[] 	v_CRS;

    // membership fuzzy sets for L outputs (rules will call these categories "NON", "END", "BDY", "BIF", "CRS")
    private static float[] q_NON;
    private static float[] q_END;
    private static float[] q_BDY;
    private static float[] q_BIF;
    private static float[] q_CRS;

    public Fuzzy2D(int _N, float _mu_ON, float _sig_ON, float _mu_OFF, float _sig_OFF)
    {
        N = _N;

        mu_ON = _mu_ON;
        sig_ON = _sig_ON;

        mu_OFF = _mu_OFF;
        sig_OFF = _sig_OFF;

        x = new float[N];
        for (int i=0; i<N; i++) x[i] = i * ((float)(L-1) / (N-1));

        out_idxs = new int[L];
        int step = (N-1)/(L-1);
        for (int i=0; i<L; i++) {
            out_idxs[i] = i*step;
        }

        // output membership functions
        q_NON = new float[N];
        for (int i=0; i<N; i++) {
            q_NON[i] = (float) Math.exp(-Math.pow(x[i] - mu_NON, 2)/(2*Math.pow(sig_OUT, 2)));
        }

        q_END = new float[N];
        for (int i=0; i<N; i++) {
            q_END[i] = (float) Math.exp(-Math.pow(x[i] - mu_END, 2)/(2*Math.pow(sig_OUT, 2)));
        }

        q_BDY = new float[N];
        for (int i=0; i<N; i++) {
            q_BDY[i] = (float) Math.exp(-Math.pow(x[i] - mu_BDY, 2)/(2*Math.pow(sig_OUT, 2)));
        }

        q_BIF = new float[N];
        for (int i=0; i<N; i++) {
            q_BIF[i] = (float) Math.exp(-Math.pow(x[i] - mu_BIF, 2)/(2*Math.pow(sig_OUT, 2)));
        }

        q_CRS = new float[N];
        for (int i=0; i<N; i++) {
            q_CRS[i] = (float) Math.exp(-Math.pow(x[i] - mu_CRS, 2)/(2*Math.pow(sig_OUT, 2)));
        }

        // aux arrays used for calculations
        agg = new float[N];
        v_NON = new float[N];
        v_END = new float[N];
        v_BDY = new float[N];
        v_BIF = new float[N];
        v_CRS = new float[N];

    }

    /*
        fuzzification - input linguistic variable (1 dimensional input)
	 */

    private float h_off(float theta)
    {
        if (theta>=mu_OFF)
            return 1;
        else //if(theta>0 && theta<=1)
            return (float) Math.exp(-Math.pow(theta - mu_OFF, 2)/(2*Math.pow(sig_OFF, 2)));
//        else
//            return 0f;

    }

    private float h_on(float theta)
    {
        if (theta<=mu_ON)
            return 1;
        else //if(theta>mu_ON) //  && theta<=mu_OFF
            return (float) Math.exp(-Math.pow(theta - mu_ON, 2)/(2*Math.pow(sig_ON, 2)));
        //else
        //    return 0;
    }

    // todo h_OFF can be a small fuzzy logic system itself 2 inputs, 1 output

    public void showFuzzification()
    {
        int NN = 101;

        float[] xx = new float[NN];

        for (int ii=0; ii<xx.length; ii++) xx[ii] = mu_ON + ii * (mu_OFF / (NN-1));

        // on
        float[] yy_on = new float[NN];
        for (int ii=0; ii<xx.length; ii++) yy_on[ii] = h_on(xx[ii]);
        //off
        float[] yy_off = new float[NN];
        for (int ii=0; ii<xx.length; ii++) yy_off[ii] = h_off(xx[ii]);

        Plot p = new Plot("fuzzification", "", "", xx, yy_on, Plot.LINE);
        p.setColor(Color.RED);
        p.draw();
        p.setColor(Color.BLUE);
        p.addPoints(xx, yy_off, Plot.LINE);
        p.show();

    }

    public void showDefuzzification()
    {
        Plot p = new Plot("defuzzification", "", "", x, q_NON);
        p.addPoints(x, q_END, Plot.LINE);
        p.addPoints(x, q_BDY, Plot.LINE);
        p.addPoints(x, q_BIF, Plot.LINE);
        p.addPoints(x, q_CRS, Plot.LINE);
        p.show();
    }

    /*
		defuzzification
	 */
    private float[] fi_NON(float z)
    {
        Arrays.fill(v_NON, 0);

        for (int i=0; i<N; i++)
            if (q_NON[i]<=z)
                v_NON[i] = q_NON[i];
            else
                v_NON[i] = z;

        return v_NON;
    }

    private float[] fi_END(float z)
    {
        Arrays.fill(v_END, 0);

        for (int i=0; i<N; i++)
            if (q_END[i]<=z)
                v_END[i] = q_END[i];
            else
                v_END[i] = z;

        return v_END;
    }

    private float[] fi_BDY(float z)
    {
        Arrays.fill(v_BDY, 0);

        for (int i=0; i<N; i++)
            if (q_BDY[i]<=z)
                v_BDY[i] = q_BDY[i];
            else
                v_BDY[i] = z;

        return v_BDY;
    }

    private float[] fi_BIF(float z)
    {
        Arrays.fill(v_BIF, 0);

        for (int i=0; i<N; i++)
            if (q_BIF[i]<=z)
                v_BIF[i] = q_BIF[i];
            else
                v_BIF[i] = z;

        return v_BIF;
    }

    private float[] fi_CRS(float z)
    {
        Arrays.fill(v_CRS, 0);

        for (int i=0; i<N; i++)
            if (q_CRS[i]<=z)
                v_CRS[i] = q_CRS[i];
            else
                v_CRS[i] = z;

        return v_CRS;
    }

    /*
        apply rules
     */
    public void critpointScores(float theta1, float[] out)
    {
        Arrays.fill(agg, 0);
        float[] cur;
        float mu;

        // apply rules
        mu = h_off(theta1); cur = fi_NON(mu); accumulate(cur, agg);
        mu = h_on(theta1);  cur = fi_END(mu); accumulate(cur, agg);
        // finished rules

        for (int i=0; i<L; i++) out[i] = agg[out_idxs[i]];

    }

    public void critpointScores(float theta1, float theta2, float[] out)
    {
        Arrays.fill(agg, 0);
        float[] cur;
        float mu;

        // apply rules
        mu = min (h_off(theta1), h_off(theta2)                 ) ; cur = fi_NON(mu); accumulate(cur, agg);
        mu = min (h_off(theta1), h_on(theta2)                  ) ; cur = fi_END(mu); accumulate(cur, agg);

        mu = min (h_on(theta1),  h_off(theta2)                 ) ; cur = fi_END(mu); accumulate(cur, agg);
        mu = min (h_on(theta1),  h_on(theta2)                  ) ; cur = fi_BDY(mu); accumulate(cur, agg);
        // finished rules

        for (int i=0; i<L; i++) out[i] = agg[out_idxs[i]];

    }

    public void critpointScores(float theta1, float theta2, float theta3, float[] out)
    {
        Arrays.fill(agg, 0);
        float[] cur;
        float mu;

        // apply rules
        mu = min (h_off(theta1), h_off(theta2), h_off(theta3)                  ) ; cur = fi_NON(mu); accumulate(cur, agg);
        mu = min (h_off(theta1), h_off(theta2), h_on(theta3)                   ) ; cur = fi_END(mu); accumulate(cur, agg);
        mu = min (h_off(theta1), h_on(theta2),  h_off(theta3)                  ) ; cur = fi_END(mu); accumulate(cur, agg);
        mu = min (h_off(theta1), h_on(theta2),  h_on(theta3)                   ) ; cur = fi_BDY(mu); accumulate(cur, agg);

        mu = min (h_on(theta1),  h_off(theta2), h_off(theta3)                  ) ; cur = fi_END(mu); accumulate(cur, agg);
        mu = min (h_on(theta1),  h_off(theta2), h_on(theta3)                   ) ; cur = fi_BDY(mu); accumulate(cur, agg);
        mu = min (h_on(theta1),  h_on(theta2),  h_off(theta3)                  ) ; cur = fi_BDY(mu); accumulate(cur, agg);
        mu = min (h_on(theta1),  h_on(theta2),  h_on(theta3)                   ) ; cur = fi_BIF(mu); accumulate(cur, agg);
        // finished rules

        for (int i=0; i<L; i++) out[i] = agg[out_idxs[i]];

    }

    public void critpointScores(float theta1, float theta2, float theta3, float theta4, float[] out) {

        Arrays.fill(agg, 0);
        float[] cur;
        float mu;

        // apply rules
        mu = min (h_off(theta1), h_off(theta2), h_off(theta3), h_off(theta4)                  ) ; cur = fi_NON(mu); accumulate(cur, agg); // h_on(theta5)
        mu = min (h_off(theta1), h_off(theta2), h_off(theta3), h_on(theta4)                   ) ; cur = fi_END(mu); accumulate(cur, agg); //1-h_off(theta5)
        mu = min (h_off(theta1), h_off(theta2), h_on(theta3),  h_off(theta4)                  ) ; cur = fi_END(mu); accumulate(cur, agg); //    ,  1-h_off(theta5)
        mu = min (h_off(theta1), h_off(theta2), h_on(theta3),  h_on(theta4)                   ) ; cur = fi_BDY(mu); accumulate(cur, agg); // 1-h_off(theta5)
        mu = min (h_off(theta1), h_on(theta2),  h_off(theta3), h_off(theta4)                  ) ; cur = fi_END(mu); accumulate(cur, agg); //,  1-h_off(theta5)
        mu = min (h_off(theta1), h_on(theta2),  h_off(theta3), h_on(theta4)                   ) ; cur = fi_BDY(mu); accumulate(cur, agg); // h_on(theta5)
        mu = min (h_off(theta1), h_on(theta2),  h_on(theta3),  h_off(theta4)                  ) ; cur = fi_BDY(mu); accumulate(cur, agg); // h_on(theta5)
        mu = min (h_off(theta1), h_on(theta2),  h_on(theta3),  h_on(theta4)                   ) ; cur = fi_BIF(mu); accumulate(cur, agg); //  h_on(theta5)

        mu = min (h_on(theta1), h_off(theta2),  h_off(theta3), h_off(theta4)                   ) ; cur = fi_END(mu); accumulate(cur, agg); //,  1-h_off(theta5)
        mu = min (h_on(theta1), h_off(theta2),  h_off(theta3), h_on(theta4)                    ) ; cur = fi_BDY(mu); accumulate(cur, agg); // ,   h_on(theta5)
        mu = min (h_on(theta1), h_off(theta2),  h_on(theta3),  h_off(theta4)                   ) ; cur = fi_BDY(mu); accumulate(cur, agg); // ,  h_on(theta5)
        mu = min (h_on(theta1), h_off(theta2),  h_on(theta3),  h_on(theta4)                    ) ; cur = fi_BIF(mu); accumulate(cur, agg); //   h_on(theta5)
        mu = min (h_on(theta1), h_on(theta2),   h_off(theta3), h_off(theta4)                   ) ; cur = fi_BDY(mu); accumulate(cur, agg); // ,  h_on(theta5)
        mu = min (h_on(theta1), h_on(theta2),   h_off(theta3), h_on(theta4)                    ) ; cur = fi_BIF(mu); accumulate(cur, agg);  //,   h_on(theta5)
        mu = min (h_on(theta1), h_on(theta2),   h_on(theta3),  h_off(theta4)                   ) ; cur = fi_BIF(mu); accumulate(cur, agg); //,h_on(theta5)
        mu = min (h_on(theta1), h_on(theta2),   h_on(theta3),  h_on(theta4)                    ) ; cur = fi_CRS(mu); accumulate(cur, agg); // ,   h_on(theta5)
        // finished rules

        for (int i=0; i<L; i++) out[i] = agg[out_idxs[i]];

    }

//        if (true) {
//            Plot p = new Plot("", "", "", x, agg);
//            p.setLimits(0, (L-1), 0, 1);
//            p.show();
//            System.out.println(Arrays.toString(agg));
//        }

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

    public static final void accumulate(float[] values, float[] accumulator)
    {
        for (int i=0; i<N; i++) if (values[i]>accumulator[i]) accumulator[i] = values[i];
    }

}
