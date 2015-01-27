package tracing2d;

import aux.Interpolator;
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

    public static void zncc(         float x,
                                     float y,
                                     int   nang,
                                     float[][] img_xy,
                                     float[][] aux_vals,    // has to be allocated properly
                                     float[] out_corr_radius
    )
    {

        out_corr_radius[0] = -1;
        out_corr_radius[1] = 0;

        for (int iang = 0; iang < nang; iang++) {

            float ang = (float) (iang * (Math.PI/nang));
            float vx = (float) Math.cos(ang);
            float vy = (float) Math.sin(ang);

            // ncc calculation
            zncc(   x,y,vx,vy,
                    img_xy,
                    aux_vals,
                    out_corr_radius
            );

        }

    }

    public static void zncc( float      x,
                             float      y,
                             float      vx,
                             float      vy,
                             float[][]  img_xy,
                             float[][]  img_vals,
                             float[]    out_ncc_radius
                             ){

        out_ncc_radius[0] = -1;
        out_ncc_radius[1] = Float.NaN;

        int imgW = img_xy.length;
        int imgH = img_xy[0].length;

        for (int irad = 0; irad < tt.length; irad++) {

            // check if this sampling range falls into the image at all - borders are enough to check
            float p1x = x + nsamp2[irad]*cross_sampling*vy; // x+L*max_step*vy
            float p1y = y - nsamp2[irad]*cross_sampling*vx; // y-L*max_step*vx
            float p2x = x - nsamp2[irad]*cross_sampling*vy; // x-L*max_step*vy
            float p2y = y + nsamp2[irad]*cross_sampling*vx; // y+L*max_step*vx

            if ((p1x>=0 && p1x<=imgW-1 && p1y>=0 && p1y<=imgH-1) &&
                    (p2x>=0 && p2x<=imgW-1 && p2y>=0 && p2y<=imgH-1)) {

                img_vals[irad][nsamp2[irad]] = Interpolator.interpolateAt(x,y,img_xy);
                float avgVals = img_vals[irad][nsamp2[irad]];

                for (int shift = 1; shift <= nsamp2[irad]; shift++) {

                    img_vals[irad][nsamp2[irad]+shift] = Interpolator.interpolateAt(
                            x + shift * cross_sampling * vy,
                            y - shift * cross_sampling * vx,
                            img_xy
                    );
                    avgVals += img_vals[irad][nsamp2[irad]+shift];

                    img_vals[irad][nsamp2[irad]-shift] = Interpolator.interpolateAt(
                            x - shift * cross_sampling * vy,
                            y + shift * cross_sampling * vx,
                            img_xy
                    );
                    avgVals += img_vals[irad][nsamp2[irad]-shift];

                }

                avgVals /= (float) (2*nsamp2[irad]+1);

                float corra = 0;
                float corrb = 0;
                float corrc = 0;

                for (int i = 0; i < ((2*nsamp2[irad]+1)); i++) {
                    corra += (img_vals[irad][i]-avgVals) * (tt[irad][i]-tta[irad]);
                    corrb += Math.pow(img_vals[irad][i]-avgVals,2);
                    corrc += Math.pow(tt[irad][i]-tta[irad],2);
                }

                float ncc = (float) (corra / Math.sqrt(corrb*corrc));

                if (ncc>out_ncc_radius[0]) {
                    out_ncc_radius[0] = ncc;
                    out_ncc_radius[1] = radiuses[irad];
                }

            }

        }

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
