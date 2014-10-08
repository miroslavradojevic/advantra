package tracing2d;

import aux.Interpolator;
import aux.Stat;
import ij.ImagePlus;
import ij.gui.OvalRoi;
import ij.gui.Overlay;

import java.awt.*;
import java.util.ArrayList;

/**
 * Created by miroslav on 6-10-14.
 * take image and set of parameters ot initialize
 *
 */
public class Tracker {



    float           D;//                   = 5f;
    float           prior_sigma_deg;//     = 35f;
    int             Nt;// = 100;
//    boolean         save_midresults     = true;
    int             MAX_ITER = 500;
//    boolean         ADD_PARTICLES = false;
    int             MARGIN = 20;

    float[][]       likelihood_xy;
    Sphere2D        sph2d;

    // track components
    ArrayList<float[][]>    Xt_xy       = new ArrayList<float[][]>();   // this one is recursively updated, starts with 1 but keeps N elements distribution
    ArrayList<float[]>      wt_xy       = new ArrayList<float[]>();     // weights of the states that describe tube configuration
    ArrayList<float[]>      est_xy      = new ArrayList<float[]>();
    Overlay                 trace       = new Overlay();


    public Tracker(
            String      _scales_list,               //  comma separated values
            ImagePlus   _inimg,                     //
            float       _D,                         //
            float       _prior_sigma_deg,           //
            int         _Nt                         //
    )
    {

        // use _scales_list and _inimg to form neuriteness - will be used as the source for the measurement
        String[] scales1 = _scales_list.split(","); if (scales1.length==0) return;
        float[] scales = new float[scales1.length];
        for (int i=0; i<scales1.length; i++) scales[i] = Float.valueOf(scales1[i]);
        ImagePlus imout = Calc.neuriteness(_inimg, scales);
//        imout.show();
        likelihood_xy = new float[imout.getWidth()][imout.getHeight()]; 	// x~column, y~row
        float[] read = (float[]) imout.getProcessor().getPixels();
        for (int idx=0; idx<read.length; idx++) {
            likelihood_xy[idx%imout.getWidth()][idx/imout.getWidth()] = read[idx];
        }

        D = _D;
        prior_sigma_deg = _prior_sigma_deg;
        Nt = _Nt;

        sph2d = new Sphere2D(D, 1.0f); // second argument is scale - here fixed to 1f, this is geometrical component for tracking

    }

//    public int[] tracktest(int a){
//        int[] out = new int[5];
//        for (int i = 0; i < out.length; i++) {
//            out[i] = a+i;
//        }
//        return out;
//    }


    public int track(float _x, float _y, float _vx, float _vy) { // return the outcome label

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

        trace.clear();

        /*
            recursive tracking
         */
        int iter = 1;
        while (iter<=MAX_ITER) {

//            Overlay track_iter = new Overlay();

            if (iter==1) {

                trace = bayesian_iteration(Nt, sph2d, prior_sigma_deg, likelihood_xy,
                        Xt_xy, wt_xy, est_xy,
                        new float[] {_vx, _vy}); // priors with manually set direction (used at the start)

            }
            else {

                trace = bayesian_iteration(Nt, sph2d, prior_sigma_deg, likelihood_xy,
                        Xt_xy, wt_xy, est_xy
                ); // priors with manually set direction (used at the start)

            }

            // check if it is in/out by some margin
            float last_x = est_xy.get(est_xy.size()-1)[0];
            float last_y = est_xy.get(est_xy.size()-1)[1];

            float x1 = last_x - 0;
            float x2 = likelihood_xy.length - last_x;

            float y1 = last_y - 0;
            float y2 = likelihood_xy[0].length - last_y;

            if (x1<MARGIN || x2<MARGIN || y1<MARGIN || y2<MARGIN) {
//                System.out.println("stop, reached the margin.");
                break;
            }

//            if (ADD_PARTICLES) {
//                for (int i = 0; i < track_iter.size(); i++) curr_ov.add(track_iter.get(i));
//            }

//            System.out.println("iter " + iter + " ->   " + est_xy.get(est_xy.size()-1)[0] + " , " + est_xy.get(est_xy.size()-1)[1]);

            iter++;

        } // while()

        return iter;

    }

    private static Overlay bayesian_iteration(

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

        // likelihoods - vesselness interpolated at each location

        for (int i = 0; i < transition_xy.length; i++) {
            ptes[i] *= Interpolator.interpolateAt(transition_xy[i][0], transition_xy[i][1], _likelihood_xy);
        }

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

            _est_xy.add(selection_e);

        }



        // initialize output for this iteration
        Overlay ov = new Overlay();

        float[] radiuses = ptes.clone();
        Stat.normalize(radiuses);
        for (int i = 0; i < transition_xy.length; i++) {
            OvalRoi pt = new OvalRoi(transition_xy[i][0]+.5f-0.5f*radiuses[i], transition_xy[i][1]+.5-0.5f*radiuses[i], radiuses[i], radiuses[i]);
            pt.setStrokeColor(new Color(1f,1f,1f,radiuses[i]));
            pt.setFillColor(new Color(1f,1f,1f, radiuses[i]));
            ov.add(pt);
        }

        OvalRoi centroid = new OvalRoi(selection_e[0]-.5+.5f, selection_e[1]-.5+.5f,1,1);
        centroid.setFillColor(Color.YELLOW);
        ov.add(centroid);

        return ov;

    }

    private static Overlay bayesian_iteration(

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

        // likelihoods - vesselness interpolated at each location

        for (int i = 0; i < transition_xy.length; i++) {
            ptes[i] *= Interpolator.interpolateAt(transition_xy[i][0], transition_xy[i][1], _likelihood_xy);
        }

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

            _est_xy.add(selection_e);

        }

        // initialize output for this iteration
        Overlay ov = new Overlay();

        float[] radiuses = ptes.clone();
        Stat.normalize(radiuses);
        for (int i = 0; i < transition_xy.length; i++) {
            OvalRoi pt = new OvalRoi(transition_xy[i][0]+.5f-0.5f*radiuses[i], transition_xy[i][1]+.5-0.5f*radiuses[i], radiuses[i], radiuses[i]);
            pt.setStrokeColor(new Color(1f,1f,1f,radiuses[i]));
            pt.setFillColor(new Color(1f,1f,1f, radiuses[i]));
            ov.add(pt);
        }

        //OvalRoi centroid = new OvalRoi(selection_e[0]-.5+.5f, selection_e[1]-.5+.5f,1,1);
        //centroid.setFillColor(Color.YELLOW);
        //ov.add(centroid);

        return ov;

    }

    public static int[] descending(float[] a) {

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
