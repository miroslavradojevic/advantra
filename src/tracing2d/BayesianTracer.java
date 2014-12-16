package tracing2d;

import detection2d.SemiCircle;
import ij.IJ;
import ij.ImageStack;
import ij.gui.Plot;

/**
 * Created by miroslav on 16-12-14.
 */
public class BayesianTracer {

    public static float[]       radiuses;
    public static int[]         nsamp2;
    public static float[][]     tt; // templates
    public static float[]       tta;

    private static SemiCircle scirc = new SemiCircle(1.2f);
    public static float cross_sampling = .5f;

    public static void init(int rmin, int rmax, float sigma_r)
    {

        float radstep = .5f;// rmax - rmin + 1;
        int nrads = 0;
        for (float rr = rmin; rr <= rmax; rr+=radstep) nrads++;

        radiuses = new float[nrads];
        nsamp2 = new int[nrads];
        tt = new float[nrads][];
        tta = new float[nrads];

        int idx=0;
        for (float rr = rmin; rr <= rmax; rr+=radstep) {

            radiuses[idx] = rr;
            nsamp2[idx] = Math.round((3*rr)/cross_sampling);

            tt[idx] = new float[2*nsamp2[idx]+1];
            tta[idx] = 0;

            for (int j = 0; j < tt[idx].length; j++) {
                tt[idx][j] = (float) Math.exp(-Math.pow((j-nsamp2[idx])*cross_sampling,2)/(2*radiuses[idx]*radiuses[idx]));
                tta[idx]+=tt[idx][j];
            }

            for (int j = 0; j < tt[idx].length; j++) { // normalize
                tt[idx][j] /= tta[idx];
            }

            tta[idx] = 0;
            for (int j = 0; j < tt[idx].length; j++) {
                tta[idx] += tt[idx][j];
            }
            tta[idx] = tta[idx]/tt[idx].length;

            idx++;

        }

    }

    public static ImageStack show_templates()
    {
        int wd = new Plot("", "", "", new float[1], new float[1]).getProcessor().getWidth();
        int ht = new Plot("", "", "").getProcessor().getHeight();
        ImageStack is = new ImageStack(wd, ht);

        String profiles = "templates<-c(";

        for (int ridx = 0; ridx < tt.length; ridx++) {

            int LL = tt[ridx].length;
            int RR = (LL-1)/2;

            float[] xx = new float[LL];
            xx[RR] = 0;
            for (int i = 1; i <= RR; i++) {
                xx[RR-i] = - i * cross_sampling;
                xx[RR+i] = + i * cross_sampling;
            }

            Plot plt = new Plot("", "", "", xx, tt[ridx], Plot.LINE);
            is.addSlice("r="+ IJ.d2s(radiuses[ridx], 2), plt.getProcessor());

                for (int k = 0; k < xx.length; k++)         profiles += xx[k] + ",";
                for (int k = 0; k < tt[ridx].length; k++)   profiles += tt[ridx][k] + ",";

        }

        if (profiles != null && profiles.length() > 1) profiles = profiles.substring(0, profiles.length() - 1);
        profiles += ")";
        IJ.log(profiles);

        return is;
    }

    public static float local_radius(float x,
                                     float y,
                                     int nang,
                                     float[][] img_xy
    )
    {

        // will give the radius estimation based on the correlations obtained at different directions
        float max_corr = -1;
        float max_corr_radius = Float.NaN;

        for (int iang = 0; iang < nang; iang++) {
            float ang = (float) (Math.PI/nang);
            float vx = ;
            float vy = ;
            for (int irad = 0; irad < radiuses.length; irad++) {

            }
        }

        return max_corr_radius;
    }

    public static void linear_tracer(float       x,
                                     float       y,
                                     float       r,
                                     float       vx,
                                     float       vy,
                                     float[][]   img_xy,
                                     float       sig_deg,
                                     float[][][] xt,
                                     float[][]   wt)
    {

        // state iteration xt[iter][sample]->[x,y,vx,vy,r]

        int Niters  = xt.length;
        int Nstates = xt[0].length; // allocated max number of states

        // initialize states for iter=0, non-existing states will be filled with NaNs
        // there is one initial state with linear_tracer
        xt[0][0][0] = x;
        xt[0][0][1] = y;

        float vnorm = (float) Math.sqrt(vx*vx+vy*vy);
        xt[0][0][2] = vx/vnorm;
        xt[0][0][3] = vy/vnorm;

        xt[0][0][4] = r;

        // fill the rest of the samples with NaN - starts from one point here
        for (int i = 1; i < xt[0].length; i++) {
            xt[0][i][0] = Float.NaN;
            xt[0][i][1] = Float.NaN;
            xt[0][i][2] = Float.NaN;
            xt[0][i][3] = Float.NaN;
            xt[0][i][4] = Float.NaN;
        }

        // loop the rest
        int iter = 1;

        while (iter<xt.length) {
//            bayesian_iteration();
            iter++;
        }

    }

}
