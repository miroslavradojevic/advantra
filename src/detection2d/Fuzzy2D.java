package detection2d;

import conn.AmiraParameters;
import ij.gui.Plot;

import java.util.Arrays;

/**
 * Created by miroslav on 3/7/14.
 */
public class Fuzzy2D {

    static int      N;          // number of points
    static int[]      out_idxs;
    static int      L = 5; // number of outputs
    static float    sig;     // gaussian standard deviation
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

    public Fuzzy2D(int _N, float _sig)
    {
        N = _N;
        sig = _sig;

        x = new float[N];
        for (int i=0; i<N; i++) x[i] = i * ((float)(L-1) / (N-1));

        out_idxs = new int[L];
        int step = (N-1)/(L-1);
//        System.out.println(""+step);
        for (int i=0; i<L; i++) {
            out_idxs[i] = i*step;
//            System.out.println(i+" "+x[out_idxs[i]]);
        }

        // output membership functions
        q_NON = new float[N];
        for (int i=0; i<N; i++) {
            q_NON[i] = (float) Math.exp(-(x[i]-0)*(x[i]-0)/(2*sig*sig));
        }

        q_END = new float[N];
        for (int i=0; i<N; i++) {
            q_END[i] = (float) Math.exp(-(x[i]-1)*(x[i]-1)/(2*sig*sig));
        }

        q_BDY = new float[N];
        for (int i=0; i<N; i++) {
            q_BDY[i] = (float) Math.exp(-(x[i]-2)*(x[i]-2)/(2*sig*sig));
        }

        q_BIF = new float[N];
        for (int i=0; i<N; i++) {
            q_BIF[i] = (float) Math.exp(-(x[i]-3)*(x[i]-3)/(2*sig*sig));
        }

        q_CRS = new float[N];
        for (int i=0; i<N; i++) {
            q_CRS[i] = (float) Math.exp(-(x[i]-4)*(x[i]-4)/(2*sig*sig));
        }

        // aux arrays used for calcualtoins
        agg = new float[N];
        v_NON = new float[N];
        v_END = new float[N];
        v_BDY = new float[N];
        v_BIF = new float[N];
        v_CRS = new float[N];

    }

    /*
        fuzzification - input linguistic variable
	 */

    private float h_off(float theta)
    {
        if (theta<=0)
            return 1;
        else if(theta>0 && theta<=1)
            return (float) Math.exp(-(theta-0)*(theta-0)/(2*sig*sig)); // mu = 0
        else
            return 0;
    }

    private float h_on(float theta)
    {
        if (theta<=0)
            return 0;
        else if(theta>0 && theta<=1)
            return (float) Math.exp(-(theta-1)*(theta-1)/(2*sig*sig)); // mu = 1
        else
            return 1;
    }


    public void showFuzzification()
    {
        int NN = 101;

        float[] xx = new float[NN];

        for (int ii=0; ii<xx.length; ii++) xx[ii] = ii * (1f / (NN-1));

        // on
        float[] yy_on = new float[NN];
        for (int ii=0; ii<xx.length; ii++) yy_on[ii] = h_on(xx[ii]);
        //off
        float[] yy_off = new float[NN];
        for (int ii=0; ii<xx.length; ii++) yy_off[ii] = h_off(xx[ii]);

        Plot p = new Plot("fuzzification", "", "", xx, yy_on);
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
    public void critpointScores(float theta1, float theta2, float theta3, float theta4, float theta5, float[] out) {

        Arrays.fill(agg, 0);

        float[] cur;
        float mu;

        // apply rules
        mu = min (h_off(theta1), h_off(theta2), h_off(theta3), h_off(theta4)                  ) ; cur = fi_NON(mu); accumulate(cur, agg); // h_on(theta5)
        mu = min (h_off(theta1), h_off(theta2), h_off(theta3), h_on(theta4),   1-h_off(theta5)) ; cur = fi_END(mu); accumulate(cur, agg);
        mu = min (h_off(theta1), h_off(theta2), h_on(theta3),  h_off(theta4),  1-h_off(theta5)) ; cur = fi_END(mu); accumulate(cur, agg);
        mu = min (h_off(theta1), h_off(theta2), h_on(theta3),  h_on(theta4)                   ) ; cur = fi_BDY(mu); accumulate(cur, agg); // 1-h_off(theta5)
        mu = min (h_off(theta1), h_on(theta2),  h_off(theta3), h_off(theta4),  1-h_off(theta5)) ; cur = fi_END(mu); accumulate(cur, agg);
        mu = min (h_off(theta1), h_on(theta2),  h_off(theta3), h_on(theta4)                   ) ; cur = fi_BDY(mu); accumulate(cur, agg); // h_on(theta5)
        mu = min (h_off(theta1), h_on(theta2),  h_on(theta3),  h_off(theta4)                  ) ; cur = fi_BDY(mu); accumulate(cur, agg); // h_on(theta5)
        mu = min (               h_on(theta2),  h_on(theta3),  h_on(theta4)                   ) ; cur = fi_BIF(mu); accumulate(cur, agg); // h_off(theta1), h_on(theta5)

        mu = min (h_on(theta1), h_off(theta2), h_off(theta3), h_off(theta4),  1-h_off(theta5) ) ; cur = fi_END(mu); accumulate(cur, agg);
        mu = min (h_on(theta1), h_off(theta2), h_off(theta3), h_on(theta4)                    ) ; cur = fi_BDY(mu); accumulate(cur, agg); // ,   h_on(theta5)
        mu = min (h_on(theta1), h_off(theta2), h_on(theta3),  h_off(theta4)                   ) ; cur = fi_BDY(mu); accumulate(cur, agg); // ,  h_on(theta5)
        mu = min (h_on(theta1),                h_on(theta3),  h_on(theta4)                    ) ; cur = fi_BIF(mu); accumulate(cur, agg);  //h_off(theta2), ,   h_on(theta5)
        mu = min (h_on(theta1), h_on(theta2),  h_off(theta3), h_off(theta4)                   ) ; cur = fi_BDY(mu); accumulate(cur, agg);  // ,  h_on(theta5)
        mu = min (h_on(theta1), h_on(theta2),                 h_on(theta4)                    ) ; cur = fi_BIF(mu); accumulate(cur, agg);  //h_off(theta3), ,   h_on(theta5)
        mu = min (h_on(theta1), h_on(theta2),  h_on(theta3)                                   ) ; cur = fi_BIF(mu); accumulate(cur, agg);  //h_off(theta4),h_on(theta5)
        mu = min (h_on(theta1), h_on(theta2),  h_on(theta3),  h_on(theta4)                    ) ; cur = fi_CRS(mu); accumulate(cur, agg);  // ,   h_on(theta5)

        mu = h_off(theta5); cur = fi_NON(mu); accumulate(cur, agg);
        // finished rules

//        if (true) {
//            Plot p = new Plot("", "", "", x, agg);
//            p.setLimits(0, (L-1), 0, 1);
//            p.show();
//            System.out.println(Arrays.toString(agg));
//        }

        float[] c_scores = new float[L];
        for (int i=0; i<L; i++) out[i] = agg[out_idxs[i]];
//        return c_scores;

    }

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
