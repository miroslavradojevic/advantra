package tracing2d;

import aux.Interpolator;
import aux.Stat;
import aux.Tools;
import detection2d.MyColors;
import detection2d.SemiCircle;
import ij.IJ;
import ij.ImageStack;
import ij.gui.Line;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.gui.Plot;

import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by miroslav on 27-11-14.
 * this one is used to delineate the surrounding there are multiple hypotheses diverging
 * used to cover the critical point regions
 */
public class BayesianTracerMulti {

    private static float neuron_radius  = -1;
    private static float nbhood_scale   = 2;
    private static float nbhood_radius  = -1;
    private static Expansion expan = null;
    public static float[] sigmas;// = new float[];
//    public static float[]   sigratios   = new float[]{1/8f, 1/4f, 1/2f, 1/1f};  // different shapes
//    public static float[]  sstep        = new float[]{.5f, .75f, 1.0f};         // scale of the measurement (keep it ascending)

    public static int Ni                = -1;                                    // how many times to step with radius scaled steps to make sense EMPIRICAL
    public static int ns=-1, ns2=-1;

    private static SemiCircle   scirc       = new SemiCircle();

//    private static int         Rmax     = -1;           // needs initialisation
//    private static float[]     pred_steps = null;       // depends on Rmax
//    private static float       step_big = Float.NaN;    // used to change bayesian filtering prediction steps depending on the iteration index
//    private static float       step_sml = 1.2f;         // these are limit values
    private static float        cross_profile_step        = 0.5f;
    private static float        semi_circ_step = 1.2f;
    private  static float ring = .3f;
    private static float dd = 1f;

   private static float[][]   tt;//       = new float[][]; // templates to be matched when filtering 2*Rcnt+1
    private static float[]     tta;//      = new float[];

    public enum Expansion {OUTER, INNER}

    public static Overlay extractAt(int _x, int _y, float[][] _img_xy, boolean _show_tracks, int _Ns, float _sigma_deg) // make a version for threading
    {
        float[][][]     xt      = new float[Ni+1][_Ns][4];    // bayesian filtering states (x,y, vx, vy)
        float[][]       wt      = new float[Ni+1][_Ns];       // bayesian filtering weights
        byte[][]        pt      = new byte[Ni+1][_Ns];        // bayesian filtering parent pointer - index of the parent in previous iteration
        ArrayList[][]   ft      = new ArrayList[Ni+1][_Ns];   // bayesian filtering follow up point
        for (int i = 0; i < ft.length; i++) for (int j = 0; j < ft[i].length; j++) ft[i][j] = new ArrayList<Byte>();
        float[][]       mt      = new float[Ni+1][_Ns];       // measurements, likelihoods - correlations measured during sequential filtering todo expell this later if not necessary

        spherical_wavefront_2d(_x, _y, _img_xy, _sigma_deg, expan, xt, wt, pt, ft, mt);

//        ArrayList<float[][]>    delin_locs      = new ArrayList<float[][]>();
//        ArrayList<Boolean>      delin_complete  = new ArrayList<Boolean>();
//        ArrayList<float[]>      delin_scores    = new ArrayList<float[]>();
//        ArrayList<float[]>      delin_values    = new ArrayList<float[]>();
//        ArrayList<float[]>      delin_terminals = new ArrayList<float[]>();
//        ArrayList<Float>        delin_wascores  = new ArrayList<Float>();
//        ArrayList<Float>        delin_wavalues  = new ArrayList<Float>();

        boolean[][] et                          = new boolean[Ni+1][_Ns];     // extracting track clusters
        byte[][] trace_map = new byte[Ni+1][_Ns];

//        extract(_img_xy, _x, _y, _R, ratioC, dd, xt, wt, pt, mt, et,trace_map);

        // local min-max
        float mn = Float.POSITIVE_INFINITY;
        float mx = Float.NEGATIVE_INFINITY;

        int nrad = Math.round(nbhood_radius);
        for (int xx = _x-nrad; xx <= _x+nrad; xx++) {
            if (xx>=0 && xx<_img_xy.length) {
                for (int yy = _y-nrad; yy <= _y+nrad; yy++) {
                    if (yy>=0 && yy<_img_xy[0].length) {
                        float val = _img_xy[xx][yy];
                        if (val<mn) mn = val;
                        if (val>mx) mx = val;
                    }
                }
            }
        }

//        Tools.descendingFloat(delin_wavalues);
//        float score = (delin_wavalues.size()>1)? delin_wavalues.get(0)-delin_wavalues.get(1) : (delin_wavalues.size()==1)? delin_wavalues.get(0)-mn : 0 ;
//        String resume = "";
//        for (int i = 0; i < delin_wavalues.size(); i++) resume+= IJ.d2s(delin_wavalues.get(i),3) + "  ";
//        resume += IJ.d2s(mn,3) + "";
//        resume += " $" + IJ.d2s(score,3);
//        IJ.log(resume);

//        float rad = BayesianTracerMulti.pred_path;

        Overlay ov = new Overlay(); // output

//        TextRoi tr = new TextRoi(_x+rad, _y-rad, rad, rad, resume, new Font("TimesRoman", Font.PLAIN, 1));
//        tr.setStrokeColor(Color.RED);
//        ov.add(tr);

        Overlay ov_xt = viz_xt(xt, wt, pt, trace_map, _show_tracks);
        for (int i = 0; i < ov_xt.size(); i++) ov.add(ov_xt.get(i));
//        Overlay ov_delin = viz_delin(delin_locs, delin_values, delin_terminals);
//        for (int i = 0; i < ov_delin.size(); i++) ov.add(ov_delin.get(i));

        OvalRoi circ = new OvalRoi(_x-nbhood_radius+.5, _y-nbhood_radius+.5, 2*nbhood_radius, 2*nbhood_radius);
        circ.setStrokeColor(Color.RED);
        ov.add(circ);
        float ratio_in = (expan==Expansion.INNER)? 1-ring : (expan==Expansion.OUTER)? ring : Float.NaN;
        OvalRoi circ_center = new OvalRoi(_x-nbhood_radius*ratio_in+.5, _y-nbhood_radius*ratio_in+.5, 2*nbhood_radius*ratio_in, 2*nbhood_radius*ratio_in);
        circ_center.setFillColor(new Color(1,0,0,.1f));

        ov.add(circ_center);

        return ov;

    }

    public static void spherical_wavefront_2d(
            float       at_x,
            float       at_y,
            float[][]   img_xy,
            float       sigma_deg,      // tracing parameter
            Expansion   expan,          // OUTWARD
            float[][][] xt,             // Ni+1 x Ns X 4 (x,y,vx,vy)
            float[][]   wt,             // Ni+1 x Ns
            byte[][]    pt,             // Ni+1 x Ns
            ArrayList[][] ft,           // Ni+1 x Ns
            float[][]   mt              // Ni+1 x Ns
    )
    {

        int Nstates = xt[0].length;

        // initialise states for iter=0
        for (int i = 0; i < Nstates; i++) {

            float ang = (float) ((i*2*Math.PI)/Nstates);

            // initialize depending on the expansion type
            if (expan==Expansion.OUTER) {
                xt[0][i][0] = at_x; // x
                xt[0][i][1] = at_y; // y
                xt[0][i][2] = (float) Math.cos(ang); // vx
                xt[0][i][3] = (float) Math.sin(ang); // vy
            }
            else if (expan==Expansion.INNER) {
                float vx = (float) Math.cos(ang);
                float vy = (float) Math.sin(ang);
                xt[0][i][0] = at_x + nbhood_radius * vx; // x
                xt[0][i][1] = at_y + nbhood_radius * vy; // y
                xt[0][i][2] = -vx;
                xt[0][i][3] = -vy;
            }
            else {
                System.out.println("problem - wrong expansion!!");
                System.exit(0);
            }

            wt[0][i] = 1f/(float)Nstates;

            pt[0][i] = (byte) 255; // dummy parent - should not be accessed for iteration = 0, just to fill up the array

            mt[0][i] = 0;

        }

        int iter = 1;

        while (iter<=Ni) {

            bayesian_iteration(at_x, at_y, img_xy,iter,semi_circ_step, sigma_deg, xt,wt,pt,ft,mt);

            iter++;
        }

    }

    private static Overlay viz_xt(float[][][] xt, float[][] wt, byte[][] pt, byte[][] _et, boolean show_tracks)
    {

        Overlay ov = new Overlay();
        for (int i = 0; i < xt.length; i++) { // iterations

            Stat.min_max_normalize(wt[i]);
            for (int j = 0; j <xt[i].length; j++) {

                OvalRoi p = new OvalRoi(xt[i][j][0]-dd/2+.5, xt[i][j][1]-dd/2+.5, dd, dd);

//                p.setFillColor(new Color(1f, 1f, 1f, wt[i][j]));
                if (i==xt.length-1) { // _et[i][j]!=(byte)255 _et[i][j] & 0xff
                    p.setFillColor(new Color(MyColors.getR(0), MyColors.getG(0), MyColors.getB(0), wt[i][j]));
                }
                else {
                    p.setFillColor(new Color(1f, 1f, 1f, wt[i][j]));
                }
                ov.add(p);

//                PointRoi pt = new PointRoi(xt[i][j][0]+.5, xt[i][j][1]+.5);
//                pt.setStrokeColor(new Color(1f,1f,1f,wt[i][j]));
//                pt.setFillColor(new Color(1f, 1f, 1f, wt[i][j]));
//                ov.add(pt);

            }
        }

        if (show_tracks) {

            for (int i = 1; i < xt.length; i++) {
                for (int j = 0; j < xt[i].length; j++) {
                    float currx = xt[i][j][0]+.5f;
                    float curry = xt[i][j][1]+.5f;
                    float prevx = xt[i-1][pt[i][j]&0xff][0]+.5f;
                    float prevy = xt[i-1][pt[i][j]&0xff][1]+.5f;
                    ov.add(new Line(currx, curry, prevx, prevy));
                }
            }

        }

        return ov;
    }

    private static void bayesian_iteration(
            float xc, float yc, //float rc,
            float[][]               _likelihood_xy,
            int                     _iter,  // will mark the array index where the value will be filled in
            float                   _iter_step,
//            int                     _Ns,
            float                   _prior_sigma_deg,
            float[][][]             _xt,    // outputs
            float[][]               _wt,    // weights
            byte[][]                _pt,    // parent index
            ArrayList[][]           _ft,    // follow-up point
            float[][]               _mt     // measurement
    )
    {

        int _Ns = _xt[0].length;

        scirc.reset(_iter_step); // costly
        // some aux values for bayesian filtering (slows down a lot that they are allocated here - can be outside if the step would stay constant)
        float[][]   _transition_xy       = new float[_Ns*scirc.NN][4];   // Ns*Nscirc x 4 store all the predictions at current iteration
        float[]     _ptes                = new float[_Ns*scirc.NN];      // Ns*Nscirc        store the probabilities temporarily
        float[]     _lhoods              = new float[_Ns*scirc.NN];      // Ns*Nscirc        store the likelihoods
        byte[]      _parents             = new byte[_Ns*scirc.NN];       // Ns*Nscirc        store the index of the parent

//        int L = (tt[0].length - 1)/2; // templates are allocated as Nsigmas x (2L+1) and contain L in its dimensions, take it from the first shape
        float[] valI = new float[ns]; // profile values for the image samples
        int imgW = _likelihood_xy.length;
        int imgH = _likelihood_xy[0].length;

        // take states at _iter-1
        float[][] prev_x = _xt[_iter-1];
        float[]   prev_w = _wt[_iter-1];

        int Ns = prev_x.length;

        int count = 0;
        double sum_lhoods = 0;
        for (int i = 0; i < prev_x.length; i++) {

            float x  = prev_x[i][0];
            float y  = prev_x[i][1];
            float vx = prev_x[i][2];
            float vy = prev_x[i][3];

            // prior is calculated here, possible to change it so that it is averaged(vx,vy)
            scirc.set(x, y, vx, vy, _prior_sigma_deg);

            for (int j = 0; j < scirc.NN; j++) {

                _transition_xy[count][0] = scirc.p[j][0]; // start_xy[i][0] + _sph2d.locsXY.get(j)[0];
                _transition_xy[count][1] = scirc.p[j][1]; // start_xy[i][1] + _sph2d.locsXY.get(j)[1];
                _transition_xy[count][2] = scirc.v[j][0]; // vx
                _transition_xy[count][3] = scirc.v[j][1]; // vy

                _parents[count] = (byte) i; // distribution index of the parent stored in byte variable

                /*
                    likelihoood as a correlation score at different scales (sampling in charge of scaling) and different shapes (template def.)
                 */

                _lhoods[count] = 0; // will take the highest calculated correlation out of all possible correlations

                // check if this sampling range falls into the image at all - borders are enough to check
                float p1x = _transition_xy[count][0] + ns2*cross_profile_step*_transition_xy[count][3]; // x+L*max_step*vy
                float p1y = _transition_xy[count][1] - ns2*cross_profile_step*_transition_xy[count][2]; // y-L*max_step*vx

                float p2x = _transition_xy[count][0] - ns2*cross_profile_step*_transition_xy[count][3]; // x-L*max_step*vy
                float p2y = _transition_xy[count][1] + ns2*cross_profile_step*_transition_xy[count][2]; // y+L*max_step*vx

                if ((p1x>=0 && p1x<=imgW-1 && p1y>=0 && p1y<=imgH-1) && (p2x>=0 && p2x<=imgW-1 && p2y>=0 && p2y<=imgH-1)) // both limits fall within the image - can calculate likelihood
                {
                    for (int sig_idx = 0; sig_idx < sigmas.length; sig_idx++) {

                        Arrays.fill(valI,0);
                        valI[ns2] = Interpolator.interpolateAt(_transition_xy[count][0], _transition_xy[count][1], _likelihood_xy);
                        float avgI = valI[ns2];

                        for (int shift = 1; shift <= ns2; shift++) {

                            valI[ns2+shift] = Interpolator.interpolateAt(
                                    _transition_xy[count][0] + shift * cross_profile_step * _transition_xy[count][3],
                                    _transition_xy[count][1] - shift * cross_profile_step * _transition_xy[count][2],
                                    _likelihood_xy
                            );

                            avgI += valI[ns2+shift];

                            valI[ns2-shift] = Interpolator.interpolateAt(
                                    _transition_xy[count][0] - shift * cross_profile_step * _transition_xy[count][3],
                                    _transition_xy[count][1] + shift * cross_profile_step * _transition_xy[count][2],
                                    _likelihood_xy
                            );

                            avgI += valI[ns2-shift];

                        }
                        avgI = avgI / (float)(ns);

                        // use template that corresponds to this sigma
                        float corra = 0;
                        float corrb = 0;
                        float corrc = 0;

                        for (int profileidx = 0; profileidx <ns; profileidx++) {
                            corra += (valI[profileidx]-avgI) * (tt[sig_idx][profileidx]-tta[sig_idx]);
                            corrb += Math.pow(valI[profileidx]-avgI, 2);
                            corrc += Math.pow(tt[sig_idx][profileidx]-tta[sig_idx], 2);
                        }

                        float ncc = (float) (corra / Math.sqrt(corrb*corrc));

                        if (ncc>_lhoods[count])
                            _lhoods[count] = ncc;

                    }
                }

//                for (int k = 0; k < sstep.length; k++) {
//                        for (int sigmaidx = 0; sigmaidx < sigratios.length; sigmaidx++) {}
//                }

                sum_lhoods += _lhoods[count];

                _ptes[count] = prev_w[i] * scirc.w[j]; // w stores priors, possible to immediately multiply posterior here unless it is all zeros

                count++;
            }
        }

        // just in case check - if all the lhoods were zero then stay with the priors
        if (sum_lhoods>Float.MIN_VALUE) for (int i = 0; i < _ptes.length; i++) _ptes[i] *= _lhoods[i];

        Stat.probability_distribution(_ptes); // normalize

        // take best Ns (we'll always have enough to take)
        int[] sort_idx = Tools.descending(_ptes);      // ptes will be sorted as a side effect

        int cnt = 0;
        for (int i = 0; i < sort_idx.length; i++) { //

            float currx = _transition_xy[sort_idx[i]][0];
            float curry = _transition_xy[sort_idx[i]][1];

            if ((currx-xc)*(currx-xc)+(curry-yc)*(curry-yc)<=nbhood_radius*nbhood_radius) {

                _xt[_iter][cnt][0] = _transition_xy[sort_idx[i]][0];
                _xt[_iter][cnt][1] = _transition_xy[sort_idx[i]][1];
                _xt[_iter][cnt][2] = _transition_xy[sort_idx[i]][2];
                _xt[_iter][cnt][3] = _transition_xy[sort_idx[i]][3];
                _wt[_iter][cnt]    = _ptes[i]; // because they are already sorted
                _pt[_iter][cnt]    = _parents[sort_idx[i]];
                _ft[_iter-1][_pt[_iter][cnt]&0xff].add((byte)cnt);
                _mt[_iter][cnt]    = _lhoods[sort_idx[i]];
                cnt++;
                if (cnt==Ns) break;

            }

        }

        Stat.probability_distribution(_wt[_iter]); // _wt will be normalized as a side effect

    }

    public static void init(int _R, Expansion expan_type) // _R is neuron radius in pix
    {

        neuron_radius = _R;
        nbhood_radius = nbhood_scale *_R;
        expan = expan_type;

        sigmas = new float[_R]; // sigma will range from 1:_R
        for (int i = 0; i < sigmas.length; i++) sigmas[i] = i + 1;
        ns2 = (int) Math.ceil(sigmas[sigmas.length - 1] / cross_profile_step);
        ns = ns2*2+1;
        Ni = Math.round(nbhood_radius/semi_circ_step);
//        step_big = 0.75f*Rmax*sstep[sstep.length-1];
//        lbd = (float) (Math.log(step_sml/step_big)/(1-Ni));
//        pred_steps = new float[Ni];
//        pred_path = 0;
//        for (int i = 0; i < Ni; i++) pred_path += step;

        tt = new float[sigmas.length][ns];
        tta = new float[sigmas.length];

        for (int i = 0; i < sigmas.length; i++) {
            float ag = 0;
            for (int j = 0; j < ns; j++) {
                tt[i][j] = (float) ((1/(Math.sqrt(2*Math.PI)*sigmas[i])) * Math.exp(-Math.pow(j-ns2,2)/(2*Math.pow(sigmas[i],2))));
                ag += tt[i][j];
            }
            tta[i] = ag/ns;
        }

    }

    public static ImageStack show_templates()
    {
        int wd = new Plot("", "", "").getProcessor().getWidth();
        int ht = new Plot("", "", "").getProcessor().getHeight();
        ImageStack is = new ImageStack(wd, ht);

        int LL = tt[0].length;
        int RR = (LL-1)/2;

        float[][] xx = new float[sigmas.length][LL];
        for (int i = 0; i < sigmas.length; i++) {
            xx[i][RR] = 0;
            for (int j = 1; j <= RR; j++) {xx[i][RR-j] = - j * cross_profile_step; xx[i][RR+j] = + j * cross_profile_step;}
        }

        String profiles = "templates<-c(";

        for (int j = 0; j < tt.length; j++) {
                Plot plt = new Plot("", "", "", xx[j], tt[j], Plot.LINE);
                is.addSlice("sigma="+ IJ.d2s(sigmas[j],2), plt.getProcessor());

                for (int k = 0; k < LL; k++) profiles += xx[j][k] + ",";
                for (int k = 0; k < LL; k++) profiles += tt[j][k] + ",";

        }

        if (profiles != null && profiles.length() > 1) profiles = profiles.substring(0, profiles.length() - 1);
        profiles += ")";
        IJ.log(profiles);
        return is;
    }

}