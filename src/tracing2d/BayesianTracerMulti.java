package tracing2d;

import aux.Interpolator;
import aux.Stat;
import aux.Tools;
import detection2d.SemiCircle;
import ij.IJ;
import ij.ImageStack;
import ij.gui.Plot;

import java.util.Arrays;

/**
 * Created by miroslav on 27-11-14.
 * this one is used to delineate the surrounding there are multiple hypotheses diverging
 * used to cover the critical point regions
 */
public class BayesianTracerMulti {

    public static float[]   sigratios   = new float[]{1/8f, 1/4f, 1/2f, 1/1f};          // different shapes
    public static float[]  sstep       = new float[]{.5f, .75f, 1.0f};                 // scale of the measurement (keep it ascending)
    public static int Ni = 6;                                                           // how many times to step with radius scaled steps to make sense EMPIRICAL
    private static SemiCircle scirc = new SemiCircle();

    private static int         Rmax = -1; // needs initialisation
    private static float[][]   tt      = new float[sigratios.length][]; // templates to be matched when filtering 2*Rcnt+1
    private static float[]     tta     = new float[sigratios.length];

    public static void run2d(
            float       at_x,
            float       at_y,
            float[][]   img_xy,

//            SemiCircle  circ,           // will contain the radius

//            int         Rmax,
            float       sigma_deg,      // tracing parameter
//            int         Ni,             // number of the iterations excluding the initial one
            int         Ns,             // number of the states to maintain

//            float[][]   tt,             // list of 1d templates nr_templates x template size PRECOMUTED, given outside for efficiency reasons
//            float[]     tta,            // list of precomputed template averages

//            float[][]   transition_xy,
//            float[]     ptes,
//            float[]     lhoods,
//            byte[]      parents,

            float[][][] xt,             // Ni+1 x Ns X 4 (x,y,vx,vy)
            float[][]   wt,             // Ni+1 x Ns
            byte[][]    pt              // Ni+1 x Ns
    )
    {



        float big_iter_step = 0.75f*Rmax*sstep[sstep.length-1]; // stepping further in motion model when tracing is defined with scale
        float small_iter_step = 1.2f; // related to the resolution, just enought to capture
        float lbd = (float) (Math.log(small_iter_step/big_iter_step)/(1-Ni));

        // number of locations and allocated to cover the pi region, allocate the states for iter=0
        for (int i = 0; i < Ns; i++) {

            float ang = (float) ((i*2*Math.PI)/Ns);

            xt[0][i][0] = at_x; // x
            xt[0][i][1] = at_y; // y
            xt[0][i][2] = (float) Math.cos(ang); // vx
            xt[0][i][3] = (float) Math.sin(ang); // vy

            wt[0][i] = 1/(float)Ns;

            pt[0][i] = (byte) 255; // dummy parent - should not be accessed for iteration = 0, just to fill up the array

        }


        int iter = 1;

        while (iter<=Ni) {

            float iter_step = (float) (big_iter_step*Math.exp(-lbd*(iter-1)));

            bayesian_iteration(img_xy,iter,iter_step, Ns, sigma_deg, xt,wt,pt);// xt[iter], wt[iter], pt[iter]

            iter++;
        }

    }

    private static void bayesian_iteration(

            float[][]               _likelihood_xy,

            int                     _iter,  // will mark the array index where the value will be filled in

            float                   _iter_step,
            int                     _Ns,
//            SemiCircle _scirc,              // step is contained here
            float                   _prior_sigma_deg,

//            float[][]               _tt,    // templates used for the matching (min-max normalized gaussian profiles, 2L+1 length)
//            float[]                 _tta,   // precomputed averages

//            float[][]               _transition_xy, // auxilliary variables
//            float[]                 _ptes,
//            float[]                 _lhoods,
//            byte[]                  _parents,

            float[][][]             _xt,    // outputs
            float[][]               _wt,
            byte[][]                _pt

    )
    {

        scirc.reset(_iter_step);
        // some aux values for bayesian filtering (slows down a lot that they are allocated here - can be outside if the step would stay constant)
        float[][]   _transition_xy       = new float[_Ns*scirc.NN][4];   // Ns*Nscirc x 4 store all the predictions at current iteration
        float[]     _ptes                = new float[_Ns*scirc.NN];      // Ns*Nscirc        store the probabilities temporarily
        float[]     _lhoods              = new float[_Ns*scirc.NN];      // Ns*Nscirc        store the likelihoods
        byte[]      _parents             = new byte[_Ns*scirc.NN];       // Ns*Nscirc        store the index of the parent

        int L = (tt[0].length - 1)/2; // templates are allocated as Nsigmas x (2L+1) and contain L in its dimensions, take it from the first shape
        float[] valI = new float[2*L+1]; // profile values for the image samples
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

                for (int k = 0; k < sstep.length; k++) { // loop sampling schemes (scales)

                    // sampling around _transition_xy[count][0],_transition_xy[count][1]

                    // check if this sampling range falls into the image at all - borders are enough to check
                    float p1x = _transition_xy[count][0] + L*sstep[k]*_transition_xy[count][3]; // x+L*max_step*vy
                    float p1y = _transition_xy[count][1] - L*sstep[k]*_transition_xy[count][2]; // y-L*max_step*vx

                    float p2x = _transition_xy[count][0] - L*sstep[k]*_transition_xy[count][3]; // x-L*max_step*vy
                    float p2y = _transition_xy[count][1] + L*sstep[k]*_transition_xy[count][2]; // y+L*max_step*vx

                    if ((p1x>=0 && p1x<=imgW-1 && p1y>=0 && p1y<=imgH-1) && (p2x>=0 && p2x<=imgW-1 && p2y>=0 && p2y<=imgH-1)) // both limits fall within the image - can calculate likelihood
                    {

                        Arrays.fill(valI,0); // to store the samples

                        valI[L] = Interpolator.interpolateAt(_transition_xy[count][0], _transition_xy[count][1], _likelihood_xy);
                        float avgI = valI[L];

                        for (int shift = 1; shift <= L; shift++) {

                            valI[L+shift] = Interpolator.interpolateAt(
                                    _transition_xy[count][0] + shift * sstep[k] * _transition_xy[count][3],
                                    _transition_xy[count][1] - shift * sstep[k] * _transition_xy[count][2],
                                    _likelihood_xy
                            );

                            avgI += valI[L+shift];

                            valI[L-shift] = Interpolator.interpolateAt(
                                    _transition_xy[count][0] - shift * sstep[k] * _transition_xy[count][3],
                                    _transition_xy[count][1] + shift * sstep[k] * _transition_xy[count][2],
                                    _likelihood_xy
                            );

                            avgI += valI[L-shift];

                        }
                        avgI = avgI / (float)(2*L+1);

                        // take the highest likelihood for different sigmas
                        for (int sigmaidx = 0; sigmaidx < sigratios.length; sigmaidx++) {

                            // use template that corresponds to this sigma
                            float corra = 0;
                            float corrb = 0;
                            float corrc = 0;

                            for (int profileidx = 0; profileidx <2*L+1; profileidx++) {
                                corra += (valI[profileidx]-avgI) * (tt[sigmaidx][profileidx]-tta[sigmaidx]);
                                corrb += Math.pow(valI[profileidx]-avgI, 2);
                                corrc += Math.pow(tt[sigmaidx][profileidx]-tta[sigmaidx], 2);
                            }

                            float ncc = (float) (corra / Math.sqrt(corrb*corrc));

                            if (ncc>_lhoods[count])
                                _lhoods[count] = ncc;

                        }

                    }

                }

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

        for (int i = 0; i < Ns; i++) {

            _xt[_iter][i][0] = _transition_xy[sort_idx[i]][0];
            _xt[_iter][i][1] = _transition_xy[sort_idx[i]][1];
            _xt[_iter][i][2] = _transition_xy[sort_idx[i]][2];
            _xt[_iter][i][3] = _transition_xy[sort_idx[i]][3];
            _wt[_iter][i]    = _ptes[i]; // because they are already sorted
            _pt[_iter][i]    = _parents[sort_idx[i]];

        }

        Stat.probability_distribution(_wt[_iter]); // _wt will be normalized as a side effect

    }

    public static void init(int _L) // templates will be 2L+1 samples
    {
        Rmax = _L; // scale

        tt = new float[sigratios.length][2*_L+1];
        tta = new float[sigratios.length];
        for (int i = 0; i < sigratios.length; i++) {
            float ag = 0;
            for (int j = 0; j < 2*_L+1; j++) {
                tt[i][j] = (float) Math.exp(-Math.pow(j-_L,2)/(2*Math.pow(sigratios[i]*_L,2)));
                ag += tt[i][j];
            }
            tta[i] = ag/(2*_L+1);
        }

    }

    public static ImageStack show_templates(float[][] _templates)
    {
        int wd = new Plot("", "", "").getProcessor().getWidth();
        int ht = new Plot("", "", "").getProcessor().getHeight();
        ImageStack is = new ImageStack(wd, ht);

        float[] xx = new float[_templates[0].length];
        for (int i = 0; i < xx.length; i++) {
            xx[i] = i;
        }

        for (int i = 0; i < _templates.length; i++) {
            Plot plt = new Plot("", "", "", xx, _templates[i], Plot.LINE);
            is.addSlice("sig="+ IJ.d2s(sigratios[i],2)+"xL", plt.getProcessor());

        }
        return is;
    }

}
