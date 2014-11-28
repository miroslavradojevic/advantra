package tracing2d;

import aux.Interpolator;
import aux.Stat;
import detection2d.SemiCircle;
import ij.ImageStack;

/**
 * Created by miroslav on 27-11-14.
 * this one is used to delineate the surrounding there are multiple hypotheses diverging
 * used to cover the critical point regions
 */
public class BayesianTracerMulti {

    public static void run2d(
            float       at_x,
            float       at_y,
            float[][]   img_xy,
            SemiCircle  circ,           // will contain the radius
            float       sigma_deg,      // tracing parameter
            int         Ni,             // number of the iterations excluding the initial one
            int         Ns,             // number of the states to maintain
            float[][]   tt,             // list of 1d templates nr_templates x template size
            float[][][] xt,             // Ni+1 x Ns X 4 (x,y,vx,vy)
            float[][]   wt,             // Ni+1 x Ns
            byte[][]    pt              // Ni+1 x Ns
    )
    {

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
            bayesian_iteration(img_xy, iter, circ, sigma_deg, tt, xt, wt, pt); // xt[iter], wt[iter], pt[iter]
            iter++;
        }

    }

    private static void bayesian_iteration(

            float[][]               _likelihood_xy,

            int                     _iter,  // will mark the array index where the value will be filled in

            SemiCircle _scirc,              // step is contained here
            float                   _prior_sigma_deg,

            float[][]               _tt,    // templates used for the matching (min-max normalized gaussian profiles)
            float[]                 _tta,   // precomputed averages

            float[][][]             _xt,    // outputs
            float[][]               _wt,
            byte[][]                _pt

    )
    {

        // templates  - this can be given as an argument
//        float[] t = new float[]{0, .5f, 0};
//        float ta = .5f/3;
        float[] mes = new float[_tt[0].length]; // will be measurement from the image sampled with respect to the template size
        float   mesa = Float.NaN;

        // take states at _iter-1
        float[][] prev_x = _xt[_iter-1];
        float[]   prev_w = _wt[_iter-1];

        int Ns = prev_x.length;

        // TODO allocation in the loop!!! to speed up
        float[][]   transition_xy       = new float[prev_x.length * _scirc.NN][4]; // all the predictions are stored here
        float[]     ptes                = new float[prev_x.length * _scirc.NN];
        float[]     lhoods              = new float[prev_x.length * _scirc.NN];
        byte[]      parents             = new byte[prev_x.length * _scirc.NN];
        // TODO
        int count = 0;
        double sum_lhoods = 0;
        for (int i = 0; i < prev_x.length; i++) {

            float x  = prev_x[i][0];
            float y  = prev_x[i][1];
            float vx = prev_x[i][2];
            float vy = prev_x[i][3];

            // prior is calculated here, possible to change it so that it is averaged(vx,vy)
            _scirc.set(x, y, vx, vy, _prior_sigma_deg);

            for (int j = 0; j < _scirc.NN; j++) {

                transition_xy[count][0] = _scirc.p[j][0]; // start_xy[i][0] + _sph2d.locsXY.get(j)[0];
                transition_xy[count][1] = _scirc.p[j][1]; // start_xy[i][1] + _sph2d.locsXY.get(j)[1];
                transition_xy[count][2] = _scirc.v[j][0]; // vx
                transition_xy[count][3] = _scirc.v[j][1]; // vy

                parents[count] = (byte) i; // distribution index of the parent stored in byte variable

                // immediately calculate the likelihood - correlation score
                if ( // check if it is within the image before interpolating
                        transition_xy[count][0]<0                           ||
                                transition_xy[count][0]>_likelihood_xy.length-1     ||
                                transition_xy[count][1]<0                           ||
                                transition_xy[count][1]>_likelihood_xy[0].length-1
                        ){
                    lhoods[count] = 0;
                }
                else {


                    boolean use_pix = false;
                    if (use_pix) {
                        lhoods[count] = Interpolator.interpolateAt(transition_xy[count][0], transition_xy[count][1], _likelihood_xy);
                    }
                    else {

                        // todo additional check there are some that would go out!!! when using measurement model sampling
                        // calculate correlation as measurement
                        // do a multiscale gaussian, several sigmas, and widths
                        m[0] = Interpolator.interpolateAt(transition_xy[count][0]+transition_xy[count][3], transition_xy[count][1]-transition_xy[count][2], _likelihood_xy);
                        m[1] = Interpolator.interpolateAt(transition_xy[count][0], transition_xy[count][1], _likelihood_xy);
                        m[2] = Interpolator.interpolateAt(transition_xy[count][0]-transition_xy[count][3], transition_xy[count][1]+transition_xy[count][2], _likelihood_xy);
                        ma = (m[0] + m[1] + m[2]) / 3f;

                        lhoods[count] = (float) (((m[0]-ma)*(t[0]-ta) + (m[1]-ma)*(t[1]-ta) + (m[2]-ma)*(t[2]-ta)) /
                                Math.sqrt(
                                        (Math.pow(m[0]-ma,2)+Math.pow(m[1]-ma,2)+Math.pow(m[2]-ma,2)) *
                                                (Math.pow(t[0]-ta,2)+Math.pow(t[1]-ta,2)+Math.pow(t[2]-ta,2))
                                )+1);
                    }

                }

                sum_lhoods += lhoods[count];

                ptes[count] = prev_w[i] * _scirc.w[j]; // w stores priors, possible to immediately multiply posterior here unless it is all zeros

                count++;
            }
        }

        // just in case check
        if (sum_lhoods>Float.MIN_VALUE) {
            for (int i = 0; i < ptes.length; i++) {
                ptes[i] = ptes[i] * lhoods[i];
            }
        }

        Stat.probability_distribution(ptes); // normalize

        // now limit the number of estimates to best Ns
//        if (transition_xy.length>Ns) {

        // reduce, take best Ns
        int[] sort_idx = descending(ptes);      // ptes will be sorted as a side effect

//            float[][]   selection_xy    = new float[Ns][4];
//            float[]     selection_w     = new float[Ns];

        for (int i = 0; i < Ns; i++) {

            _xt[_iter][i][0] = transition_xy[sort_idx[i]][0];
            _xt[_iter][i][1] = transition_xy[sort_idx[i]][1];
            _xt[_iter][i][2] = transition_xy[sort_idx[i]][2];
            _xt[_iter][i][3] = transition_xy[sort_idx[i]][3];
            _wt[_iter][i]    = ptes[i]; // because they are already sorted
            _pt[_iter][i]    = parents[sort_idx[i]];
//                selection_xy[i][0] = transition_xy[sort_idx[i]][0];
//                selection_xy[i][1] = transition_xy[sort_idx[i]][1];
//                selection_w[i] = ptes[i];           //

        }


//            _Xt_xy.add(selection_xy);
        Stat.probability_distribution(_wt[_iter]);
//            _wt_xy.add(selection_w);

        // final estimate will be...
//        float[]     selection_e     = new float[4];
//        for (int i = 0; i < Ns; i++) {
//            selection_e[0] += selection_w[i] * selection_xy[i][0];
//            selection_e[1] += selection_w[i] * selection_xy[i][1];
//        }

//            _est_xy.add(selection_e);
//        }

    }

    public static float[][] generate_templates(float _radius)
    {

        int L = (int) Math.ceil(_radius/.5f);

        float[] sigs  = new float[]{_radius/4f, _radius/2f, _radius};

        float[][] templates = new float[sigs.length][2*L+1];

        for (int i = 0; i < sigs.length; i++) {
            for (int j = 0; j < 2*L+1; j++) {
                templates[i][j] = (float) Math.exp(-Math.pow(j-L,2)/(2*sigs[i]));
            }
        }

        return templates;

    }

    public static ImageStack show_templates(float[][] _templates)
    {
        ImageStack is = new ImageStack(600, 300);
        return is;
    }


}
