package tracing2d;

import aux.Interpolator;
import aux.Stat;
import ij.ImagePlus;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.gui.PointRoi;

import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by miroslav on 6-10-14.
 * take image and set of parameters ot initialize
 *
 */
public class BayesianTracer2D {
    // will trace in 2d from the seed point defined as (x, y, vx, vy)
    float           D;                      //      = 5f;
    float           prior_sigma_deg;        //      = 35f;
    int             Nt;                     //      = 100;
    int             MAX_ITER;               //      = 100;
    int             MIN_LEN;                //      = 5;
    int             MARGIN          = 20;

    int             W, H;
    float[][]       likelihood_xy;
    Sphere2D        sph2d;

    // output trace
    ArrayList<float[][]>    Xt_xy       = new ArrayList<float[][]>();   // this one is recursively updated, starts with 1 but keeps N elements distribution
    ArrayList<float[]>      wt_xy       = new ArrayList<float[]>();     // weights of the states that describe tube configuration
    ArrayList<float[]>      est_xy      = new ArrayList<float[]>();

    public BayesianTracer2D(
            String      _scales_list,               //  comma separated values
            ImagePlus   _inimg,                     //
            float       _D,                         //
            float       _prior_sigma_deg,           //
            int         _Nt,                        //
            int         _MAX_ITER,
            int         _MIN_LEN,
            boolean     _use_original
    )
    {

        // use _scales_list and _inimg to form neuriteness - will be used as the source for the measurement
        String[] scales1 = _scales_list.split(","); if (scales1.length==0) return;
        float[] scales = new float[scales1.length];
        for (int i=0; i<scales1.length; i++) scales[i] = Float.valueOf(scales1[i]);

        W = _inimg.getWidth();
        H = _inimg.getHeight();

        // likelihood load
        likelihood_xy = new float[W][H]; 	// x~column, y~row
        if(_use_original) {
            byte[] read = (byte[]) _inimg.getProcessor().getPixels(); // use 8bit original
            for (int idx=0; idx<read.length; idx++) {
                likelihood_xy[idx%W][idx/W] = read[idx] & 0xff;
            }
        }
        else {
            float[] read = (float[]) Calc.neuriteness(_inimg, scales).getProcessor().getPixels(); // calculate neuriteness
            for (int idx=0; idx<read.length; idx++) {
                likelihood_xy[idx%W][idx/W] = read[idx];
            }
        }

        D = _D;
        prior_sigma_deg = _prior_sigma_deg;
        Nt = _Nt;
        MAX_ITER = _MAX_ITER;
        MIN_LEN = _MIN_LEN;

        sph2d = new Sphere2D(D, 1.0f); // second argument is scale - here fixed to 1f, this is geometrical component for tracking

    }

    public int run(  // return the outcome label
            float _x,
            float _y,
            float _vx,
            float _vy,

            int _init_label,
            int[][] _region_map // this one is used to check when to end the trace
    )
    {
        // init
        float[][] start_xy = new float[1][2];
        start_xy[0][0] = _x;
        start_xy[0][1] = _y;
        Xt_xy.clear();
        Xt_xy.add(start_xy);

        float[] start_w = new float[1];
        start_w[0] = 1f;
        wt_xy.clear();
        wt_xy.add(start_w);

        float[] start_est = new float[2];
        start_est[0] = _x;
        start_est[1] = _y;
        est_xy.clear();
        est_xy.add(start_est);

        /*
            recursive tracking
         */

        int out_label = Integer.MAX_VALUE;

        int iter = 1;
        while (iter<=MAX_ITER) {

            Overlay track_iter = new Overlay();

            // append Xt_xy, wt_xy, est_xy

            if (iter==1) {

                bayesian_iteration(Nt, sph2d, prior_sigma_deg, likelihood_xy,
                        Xt_xy, wt_xy, est_xy,               // outputs
                        new float[] {_vx, _vy});            // priors with manually set direction (used at the start)

            }
            else {

                bayesian_iteration(Nt, sph2d, prior_sigma_deg, likelihood_xy,
                        Xt_xy, wt_xy, est_xy                // outputs
                );                                          // priors with manually set direction (used at the start)

            }

            // check if it is in/out by some margin or NaN
            float last_x = est_xy.get(est_xy.size()-1)[0];
            float last_y = est_xy.get(est_xy.size()-1)[1];

            if (Float.isNaN(last_x) || Float.isNaN(last_y)) {
                System.out.print("STOP, coords were: " + last_x + " , " + last_y + " at # estimates " + est_xy.size() + "    ");
                break; // out_label stays Integer.MAX_VALUE
            }

            float x1 = last_x - 0;
            float x2 = W - last_x;

            float y1 = last_y - 0;
            float y2 = H - last_y;

            if (x1<MARGIN || x2<MARGIN || y1<MARGIN || y2<MARGIN) {
                //System.out.print("STOP, reached the image margin.");
                break; // out_label stays Integer.MAX_VALUE
            }

            if (
                    iter>MIN_LEN &&
                    _region_map[Math.round(last_x)][Math.round(last_y)]!=0 &&
                    _region_map[Math.round(last_x)][Math.round(last_y)]!=_init_label
            ){
                //System.out.print("STOP, reached another CP region.");
                out_label = _region_map[Math.round(last_x)][Math.round(last_y)];
                break;
            }

            iter++;

        }

        return out_label;

    }

    private static void bayesian_iteration(

            int                     Nt,
            Sphere2D                _sph2d,
            float                   _prior_sigma_deg,
            float[][]               _likelihood_xy,

            ArrayList<float[][]>    _Xt_xy,
            ArrayList<float[]>      _wt_xy,
            ArrayList<float[]>      _est_xy,

            float[]                 _vxy     // priors with manually set direction (used at the start)
    )
    {

        float[][] start_xy = _Xt_xy.get(_Xt_xy.size()-1);

        int N = start_xy.length;

        float[][] transition_xy = new float[N * _sph2d.N][2]; // all the predictions are stored here

        int count = 0;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < _sph2d.N; j++) {
                transition_xy[count][0] = start_xy[i][0] + _sph2d.locsXY.get(j)[0];
                transition_xy[count][1] = start_xy[i][1] + _sph2d.locsXY.get(j)[1];
                count++;
            }
        }

        // assign p prev for each - inherit it from wt
        float[] start_w = _wt_xy.get(_wt_xy.size()-1);

        N = start_w.length;

        float[] ptes = new float[N * _sph2d.N];

        count = 0;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < _sph2d.N; j++) {
                ptes[count] = start_w[i];
                count++;
            }
        }

        // mul. priors (depending on the direction)
        float[] priors_per_sphere = new float[_sph2d.N];

        count=0;
        for (int i = 0; i < N; i++) {

            // direction for this point
            _sph2d.getPriors(_vxy, _prior_sigma_deg, priors_per_sphere);

            for (int j = 0; j < _sph2d.N; j++) {

                ptes[count] *= priors_per_sphere[j];
                count++;

            }
        }

//          REPLACED!!!
//        for (int i = 0; i < transition_xy.length; i++) {
//            ptes[i] *= Interpolator.interpolateAt(transition_xy[i][0], transition_xy[i][1], _likelihood_xy); // apply the likelihoods
//        }


        // fix: can happen that all the likelihoods are zero - that would cancel all the info in priors (all the posteriors reduced to zero)
        // if it is >0 and still equal - priors are propagated
        // if all are zero => skip the multiplication stage
        float check_sum = 0;
        boolean apply_likelihoods;
        for (int i = 0; i < transition_xy.length; i++)
            check_sum += Interpolator.interpolateAt(transition_xy[i][0], transition_xy[i][1], _likelihood_xy);

        if (check_sum<=Float.MIN_VALUE) apply_likelihoods = false;
        else apply_likelihoods = true;

        if (apply_likelihoods) for (int i = 0; i < transition_xy.length; i++)
            ptes[i] *= Interpolator.interpolateAt(transition_xy[i][0], transition_xy[i][1], _likelihood_xy); // likelihoods: vesselness, or the original values or anything else, interpolated at each location

        Stat.probability_distribution(ptes);

        float[]     selection_e     = new float[2];
        // now limit the number of estimates to best Nt
        if (transition_xy.length>Nt) {

            // reduce, take best Nt
            int[] sort_idx = descending(ptes);      // ptes will be sorted as a side effect

            float[][]   selection_xy    = new float[Nt][2];
            float[]     selection_w     = new float[Nt];

            for (int i = 0; i < Nt; i++) {
                selection_xy[i][0] = transition_xy[sort_idx[i]][0];
                selection_xy[i][1] = transition_xy[sort_idx[i]][1];
                selection_w[i] = ptes[i];           // because they are already sorted
            }

            _Xt_xy.add(selection_xy);
            Stat.probability_distribution(selection_w);
            _wt_xy.add(selection_w);

            // final estimate will be the mean
            for (int i = 0; i < Nt; i++) {
                selection_e[0] += selection_w[i] * selection_xy[i][0];
                selection_e[1] += selection_w[i] * selection_xy[i][1];
            }

            if (Float.isNaN(selection_e[0]) || Float.isNaN(selection_e[1])) {
                System.out.println("first est_xy! adding NaN, >Nt");
            }
            _est_xy.add(selection_e);

        }
        else {
            // keep all of them
            _Xt_xy.add(transition_xy);
            _wt_xy.add(ptes);

            // final estimate will be the mean
            for (int i = 0; i < transition_xy.length; i++) {
                selection_e[0] += ptes[i] * transition_xy[i][0];
                selection_e[1] += ptes[i] * transition_xy[i][1];
            }
//            if (Float.isNaN(selection_e[0]) || Float.isNaN(selection_e[1])) {
//                boolean test = true;
//                for (int i = 0; i < ptes.length; i++) {
//                    test = test & !Float.isNaN(ptes[i]);
//                }
//                System.out.println("first est_xy! adding NaN, <=Nt (keep all) ptes ok? " + test + "   " + transition_xy.length + "  " + ptes.length + "   " + Arrays.toString(ptes));
//            }
            _est_xy.add(selection_e);

        }



    }

    private static void bayesian_iteration(

            int                     Nt,
            Sphere2D                _sph2d,
            float                   _prior_sigma_deg,
            float[][]               _likelihood_xy,

            ArrayList<float[][]>    _Xt_xy,
            ArrayList<float[]>      _wt_xy,
            ArrayList<float[]>      _est_xy

    )
    {

        float[][] start_xy = _Xt_xy.get(_Xt_xy.size()-1); // start from the current distribution

        int N = start_xy.length;                        // there will be N spheres

        float[][] transition_xy = new float[N * _sph2d.N][2]; // all the predictions are stored here

        int count = 0;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < _sph2d.N; j++) {
                transition_xy[count][0] = start_xy[i][0] + _sph2d.locsXY.get(j)[0];
                transition_xy[count][1] = start_xy[i][1] + _sph2d.locsXY.get(j)[1];
                count++;
            }
        }

        // assign p prev for each - inherit it from wt
        float[] start_w = _wt_xy.get(_wt_xy.size()-1);

        N = start_w.length;

        float[] ptes = new float[N * _sph2d.N];

        count = 0;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < _sph2d.N; j++) {
                ptes[count] = start_w[i];
                count++;
            }
        }

        // mul. priors (depending on the direction)
        float[] priors_per_sphere = new float[_sph2d.N];

        float[] _vxy = new float[2];

        count=0;
        for (int i = 0; i < N; i++) {

            // direction for this point

            _vxy[0] = 0;
            _vxy[1] = 0;

            int nr_avg = 0;
            while (nr_avg<5) {



                if (nr_avg==0) {
                    _vxy[0] += start_xy[i][0] - _est_xy.get(_est_xy.size()-2)[0];
                    _vxy[1] += start_xy[i][1] - _est_xy.get(_est_xy.size()-2)[1];
                    nr_avg++;
                }
                else {
                    int last = _est_xy.size()-2-nr_avg;
                    if (last>=0) {
                        _vxy[0] += _est_xy.get(_est_xy.size()-2-(nr_avg-1) )[0] - _est_xy.get(_est_xy.size()-2-(nr_avg) )[0];
                        _vxy[1] += _est_xy.get(_est_xy.size()-2-(nr_avg-1) )[1] - _est_xy.get(_est_xy.size()-2-(nr_avg) )[1];
                        nr_avg++;
                    }
                    else {
                        break;
                    }

                }

            }

            if (nr_avg>0) {
                _vxy[0] = _vxy[0] / nr_avg;
                _vxy[1] = _vxy[1] / nr_avg;
            }

            _sph2d.getPriors(_vxy, _prior_sigma_deg, priors_per_sphere);

            for (int j = 0; j < _sph2d.N; j++) {

                ptes[count] *= priors_per_sphere[j];
                count++;

            }
        }

        // fix: can happen that all the likelihoods are zero - that would cancel all the info in priors (all the posteriors reduced to zero)
        // if it is >0 and still equal - priors are propagated
        // if all are zero => skip the multiplication stage
        float check_sum = 0;
        boolean apply_likelihoods;
        for (int i = 0; i < transition_xy.length; i++)
            check_sum += Interpolator.interpolateAt(transition_xy[i][0], transition_xy[i][1], _likelihood_xy);

        if (check_sum<=Float.MIN_VALUE) apply_likelihoods = false;
        else apply_likelihoods = true;

        if (apply_likelihoods) for (int i = 0; i < transition_xy.length; i++)
            ptes[i] *= Interpolator.interpolateAt(transition_xy[i][0], transition_xy[i][1], _likelihood_xy); // likelihoods: vesselness, or the original values or anything else, interpolated at each location

        Stat.probability_distribution(ptes);

        float[]     selection_e     = new float[2];
        // now limit the number of estimates to best Nt
        if (transition_xy.length>Nt) {

            // reduce, take best Nt
            int[] sort_idx = descending(ptes);      // ptes will be sorted as a side effect



            float[][]   selection_xy    = new float[Nt][2];
            float[]     selection_w     = new float[Nt];

            for (int i = 0; i < Nt; i++) {
                selection_xy[i][0] = transition_xy[sort_idx[i]][0];
                selection_xy[i][1] = transition_xy[sort_idx[i]][1];
                selection_w[i] = ptes[i];           // because they are already sorted
            }


            _Xt_xy.add(selection_xy);
            Stat.probability_distribution(selection_w);
            _wt_xy.add(selection_w);

            // final estimate will be the mean
            for (int i = 0; i < Nt; i++) {
                selection_e[0] += selection_w[i] * selection_xy[i][0];
                selection_e[1] += selection_w[i] * selection_xy[i][1];
            }

//            if (Float.isNaN(selection_e[0]) || Float.isNaN(selection_e[1])) {
//                System.out.println("adding NaN " + _est_xy.size());
//            }

            _est_xy.add(selection_e);

        }
        else {
            // keep all of them
            _Xt_xy.add(transition_xy);
            _wt_xy.add(ptes);

            // final estimate will be the mean
            for (int i = 0; i < transition_xy.length; i++) {
                selection_e[0] += ptes[i] * transition_xy[i][0];
                selection_e[1] += ptes[i] * transition_xy[i][1];
            }

//            if (Float.isNaN(selection_e[0]) || Float.isNaN(selection_e[1])) {
//                System.out.println("adding NaN " + _est_xy.size());
//            }

            _est_xy.add(selection_e);

        }

    }

    public static int[] descending(float[] a)
    {

        // prepare array with indexes first
        int[] idx = new int[a.length];
        for (int i=0; i<idx.length; i++) idx[i] = i;

        for (int i = 0; i < a.length-1; i++) {
            for (int j = i+1; j < a.length; j++) {
                if (a[j]>a[i]) { // desc.
                    float temp 	= a[i];
                    a[i]		= a[j];
                    a[j] 		= temp;

                    int temp_idx 	= idx[i];
                    idx[i] 			= idx[j];
                    idx[j]			= temp_idx;
                }
            }
        }

        return idx;

    }

}
