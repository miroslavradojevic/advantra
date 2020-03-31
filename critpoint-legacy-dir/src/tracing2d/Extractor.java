package tracing2d;

import aux.Interpolator;
import aux.Stat;
import aux.Tools;
import detection2d.MyColors;
import ij.IJ;
import ij.gui.*;
import ij.process.ByteProcessor;

import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by miroslav on 29-11-14.
 * extract local number of streamlines, their intensities based on the cloud of the locations covering, average intensities for each streamline
 * possible to multithread for selected number od locations
 * possible to call only for one particular location
 */
public class Extractor extends Thread {

    private int begN, endN;

    public  static int[][] 	 i2loc;
//    public  static int       R;
//    public  static int       Ns;
//    public  static float     sigma_deg;
//    public  static int       nstreams;
    private static int       dim;

//    private static int merging_margin_step = 1; // after this stream iteration index merging of the confluent paths is possible
//    private static float merging_margin_dist = 1f;

    public static ArrayList[] i2scores;     // i2scores[0]      is ArrayList<float[]>
    public static ArrayList[] i2locEnd;  // i2terminals[0]      is ArrayList<float[]>
    public static ArrayList[] i2range;// min-max Ni lowest / Ni highest in the neighbourhood defined with the range

    // number of iterations is maximum Ni+1 by definition
    private static float[][] wgts_outer; // set of weights assigned to each streamline length (Ni+1) high towards the end
//    private static float[][] wgts_inner; // high at the beginnining

    public Extractor (int n0, int n1)
    {
        begN = n0;
        endN = n1;
    }

    /* initializator for usage at one particular location */
    public static void loadTemplate(
            int     _dim,
            int     _radius,
//            float   _sigma_deg,
            BayesianTracerMulti.CircExpansion _expan_type
    )
    {

        i2loc        = null;
        i2scores     = null;
        i2locEnd     = null;
        i2range      = null;
        dim         = _dim;
//        R           = _radius;
//        Ns          = _nsamples;
//        sigma_deg   = _sigma_deg;
        BayesianTracerMulti.init(_radius);

        weights();

    }

    private  static void weights()
    {
        wgts_outer = new float[BayesianTracerMulti.Ni+1][];

        for (int i = 0; i < wgts_outer.length; i++) {

            wgts_outer[i] = new float[i+1];
            float sum = 0;

            for (int j = 0; j <= i; j++) {
                float ratio = (i>0)? (float)j/i : 0;
                wgts_outer[i][j] = (float) (1 - Math.exp(-2*ratio));
                sum += wgts_outer[i][j];
            }
            for (int j = 0; j <= i; j++) {
                wgts_outer[i][j] /= sum;
            }

        }
    }

    /* initializator for usage at selected set of locations (threaded) */
    public static void loadTemplate(
            int[][] _i2loc,
            int     _radius,
//            float   _sigma_deg,
            BayesianTracerMulti.CircExpansion _expan_type
    )
    {
        i2loc   = _i2loc;
        i2scores = new ArrayList[i2loc.length];
        i2locEnd = new ArrayList[i2loc.length];
        i2range  = new ArrayList[i2loc.length];

        for (int i = 0; i < i2loc.length; i++) {
            i2scores[i] = new ArrayList<Float>();
            i2locEnd[i] = new ArrayList<float[]>();
            i2range[i]  = new ArrayList<float[]>();
        }
        
        dim     = i2loc[0].length;
//        R       = _radius;
//        Ns      = _nsamples;
//        sigma_deg   = _sigma_deg;
        BayesianTracerMulti.init(_radius);

        weights();

    }

    public void run()
    {
        //allocate xt, etc...
    }

    private static int pickmax(float[] _wti, boolean[] _eti) // will pick max from the current set of weights that are not examined
    {
        float max_score = Float.NEGATIVE_INFINITY;
        int max_index = -1;

        for (int j = 0; j < _wti.length; j++) { // loop the last iteration - loop all the states/weights from the distribution
            if (!_eti[j]) {
                if (_wti[j] > max_score) {
                    max_score = _wti[j];
                    max_index = j;
                }
            }
        }
        return max_index;
    }

    private static int  extract(
            float[][]               _in_xy,
            float                   _xc,
            float                   _yc,
            float                   _r,
            float                   _ratio_center,
            float                   _nearness,
            float[][][]             _xt,                // states
            float[][]               _wt,                // weights (posteriors)
            byte[][]                _pt,                // parent (needed for top-bottom)
            float[][]               _mt,                // measurements (likelihoods, correlations)
            boolean[][]                _et,             // examined map (auxiliary map for extraction methods)
            byte[][]                _trace_map
//            ArrayList<float[][]>    _delin_locs,
//            ArrayList<float[]>      _delin_scores,
//            ArrayList<float[]>      _delin_values,
//            ArrayList<Float>        _region_averages,
//            ArrayList<float[]>      _region_terminals
//            ArrayList<Float>        _delin_wascores,
//            ArrayList<Float>        _delin_wavalues
    )
    {

        for (int i = 0; i < _et.length; i++) // reset examined map
            for (int j = 0; j < _et[i].length; j++)
            {
                _et[i][j] = false;
                _trace_map[i][j] = (byte) 255;
            }

//        _delin_locs.clear();
//        _delin_complete.clear();
//        _delin_scores.clear();
//        _delin_values.clear();
//        _region_terminals.clear();

        ArrayList<Integer> path_i = new ArrayList<Integer>();
        ArrayList<Integer> path_s = new ArrayList<Integer>();

        boolean all_examined;
        int curr_streams = 0; // stream counter

        do {

            int max_idx_samp = -1;
            int max_idx_iter = -1;

            for (int i = _wt.length-1; i >= 0; i--) { // loop bayesian iterations backwards

                max_idx_iter = i;

                max_idx_samp = pickmax(_wt[i], _et[i]); // pick next max among those that were not examined in iteration i

                if (max_idx_samp!=-1) break;
            }

            all_examined = false;

            if (max_idx_samp!=-1 && max_idx_iter!=-1) { // max was found at iteration max_idx_iter and sample index max_idx_samp

                float startx = _xt[max_idx_iter][max_idx_samp][0];
                float starty = _xt[max_idx_iter][max_idx_samp][1];

                if (Math.sqrt(Math.pow(startx - _xc, 2) + Math.pow(starty - _yc, 2))/_r<=_ratio_center) {

                    int sidx = max_idx_samp;
                    boolean is_overlapping = false;
                    boolean is_close = false;
//                    boolean passes_center = false;

                    byte label_found = (byte)255;

                    path_i.clear();
                    path_s.clear();

                    for (int iidx = max_idx_iter; iidx >=0 ; iidx--) { // top-bottom

                        _et[iidx][sidx] = true;
                        path_i.add(iidx);
                        path_s.add(sidx);

                        float currx = _xt[iidx][sidx][0];
                        float curry = _xt[iidx][sidx][1];

                        float relative_r = (float) (Math.sqrt(Math.pow(currx - _xc, 2) + Math.pow(curry - _yc, 2)) / _r);
//                        if (!passes_center) {  passes_center = relative_r <= _ratio_center; }

                        if (relative_r > _ratio_center && !is_close && !is_overlapping) { // allow merging below merging_margin_step of top-bottom
                            if (_trace_map[iidx][sidx]!=(byte)255) {
//                            System.out.println("^^there was an overlap!");
                                label_found = _trace_map[iidx][sidx];
                                is_overlapping = true;
                            }
                            else {
                                for (int titer = 0; titer < _xt.length; titer++) {
                                    for (int t = 0; t < _xt[titer].length; t++) {
                                        if (_trace_map[titer][t]!=(byte)255) {

                                            if (!(titer==iidx && t == sidx)) {

                                                float checkx = _xt[titer][t][0];
                                                float checky = _xt[titer][t][1];

                                                if (Math.sqrt(Math.pow(checkx - _xc, 2) + Math.pow(checky - _yc, 2)) / _r > _ratio_center) {
                                                    if (Math.pow(checkx - currx, 2) + Math.pow(checky - curry, 2) <= Math.pow(_nearness, 2)) {
                                                        label_found = _trace_map[titer][t];
                                                        is_close = true;
                                                    }
                                                }

                                            }
                                            else IJ.log("wrong!!");

                                        }
                                    }
                                }
                            }
                        }

//                        if (is_overlapping) break; // it is guaranteed that the all the iterations backwards are marked as examined
//                        if (!is_close) { // no need to store them further - won't be added anyway, just to fill _et examined map
//                            temp_locs[0][iidx] = _xt[iidx][sidx][0];
//                            temp_locs[1][iidx] = _xt[iidx][sidx][1];
//                            temp_scores[iidx] = _mt[iidx][sidx];
//                            temp_values[iidx] = Interpolator.interpolateAt(_xt[iidx][sidx][0], _xt[iidx][sidx][1], _in_xy);
//                        }

                        sidx = _pt[iidx][sidx]&0xff; // recursion backtrack

                    }// top-bottom

                    if (is_overlapping || is_close) {
                        for (int i = 0; i < path_i.size(); i++) {
                            _trace_map[path_i.get(i)][path_s.get(i)] = label_found;
                        }
                    }
                    else {
                        for (int i = 0; i < path_i.size(); i++) {
                            _trace_map[path_i.get(i)][path_s.get(i)] = (byte)curr_streams;
                        }
                        curr_streams++;
                    }

//                    if (!is_overlapping && !is_close) { // will occur up to 4 times
//
//                        // locs xy * length
//
//                        float[][]   add_locs    = new float[temp_locs.length][path_len];
//                        float[]     add_scores  = new float[path_len];
//                        float[]     add_values  = new float[path_len];
//
//                        float add_wascores = 0;
//                        float add_wavalues = 0;
//
//                        for (int i = 0; i < path_len; i++) {
//
//                            add_locs[0][i] = temp_locs[0][i];
//                            add_locs[1][i] = temp_locs[1][i];
//                            add_scores[i] = temp_scores[i];
//                            add_values[i] = temp_values[i];
//
//                            add_wascores += wgts_outer[path_len-1][i] * add_scores[i]; // weights for L long vectors are at wgts[L-1] index
//                            add_wavalues += wgts_outer[path_len-1][i] * add_values[i];
//
//                        }
//
//                        _delin_locs.add(add_locs);
//
//                        // complete
////                    _delin_complete.add(path_len == _xt.length);
//
//                        // scores
////                        _delin_scores.add(add_scores);
////                        _delin_wascores.add(add_wascores);
//
//
//                        // values
//                        _delin_values.add(add_values);
//                        _delin_wavalues.add(add_wavalues);
//
//                        // terminal locations // if (path_len == _xt.length)
//                        _delin_terminals.add(new float[]{temp_locs[0][0],temp_locs[1][0]}); // {_xt[max_idx_iter][max_idx_samp][0], _xt[max_idx_iter][max_idx_samp][1]}
//
//                        curr_streams++;
//
//                    }




                }
                else _et[max_idx_iter][max_idx_samp] = true;

//                int path_len = max_idx_iter + 1;
//
//                // auxilliary array to fill up during the back track
//                float[][]   temp_locs   = new float[2][path_len]; // x=0,y=1
//                float[]     temp_scores = new float[path_len];
//                float[]     temp_values = new float[path_len];

            }
            else all_examined = true;

        }
        while (!all_examined); // curr_streams<1 &&

//        IJ.log(curr_streams + " streams");
return curr_streams;
    }

    private static void arrange(ArrayList<Float> _vals, ArrayList<Boolean> _compl, ArrayList<Float> _A, ArrayList<Float> _B)
    {
        _A.clear();
        _B.clear();
        int idx = 0;
        for (int i = 0; i < _vals.size(); i++) {

            if (_compl.get(i).booleanValue()==true) {
                if (_A.size()==0) _A.add(_vals.get(i));
                else { // append so that the highest one stray on the top
                    idx = 0;
                    while (idx<_A.size() && _A.get(idx)>_vals.get(i)) {
                        idx++;
                    }
                    _A.add(idx, _vals.get(i));
                }
            }
            else {
                if (_B.size()==0) _B.add(_vals.get(i));
                else {
                    idx = 0;
                    while (idx<_B.size() && _B.get(idx)>_vals.get(i)) {
                        idx++;
                    }
                    _B.add(idx, _vals.get(i));
                }
            }

        }
    }

    private static Overlay viz_delin(ArrayList<float[][]> _delin_locs, ArrayList<float[]> _delin_values, ArrayList<float[]> _delin_terminals)
    {
        Overlay ov = new Overlay();

//        System.out.println(_delin_locs.size() + "delineations");

        for (int i = 0; i < _delin_locs.size(); i++) {

            float[] plgx = _delin_locs.get(i)[0];
            for (int j = 0; j < plgx.length; j++) plgx[j] += .5;

            float[] plgy = _delin_locs.get(i)[1];
            for (int j = 0; j < plgy.length; j++) plgy[j] += .5;

            for (int j = 0; j < _delin_locs.get(i)[0].length; j++) {

//                OvalRoi or = new OvalRoi(delin.get(i).get(j)[0]-1f+.5f, delin.get(i).get(j)[1]-1f+.5f, 2f, 2f);
//                or.setFillColor(Color.RED);
//                or.setStrokeColor(Color.RED);

//                PointRoi pt = new PointRoi(delin.get(i).get(j)[0]+.5f, delin.get(i).get(j)[1]+.5f);
//                pt.setFillColor(Color.RED);
//                pt.setStrokeColor(Color.RED);

                float atx = plgx[j]; //_delin_locs.get(i)[0][j]+ .5f;
                float aty = plgy[j]; //_delin_locs.get(i)[1][j]+ .5f;
                float scr = .5f;//_delin_values.get(i)[j];

//                OvalRoi ovr = new OvalRoi(atx-.5*scr, aty-.5*scr, scr, scr);  ovr.setFillColor(Color.RED);

//                if (i==0) ovr.setFillColor(Color.RED);
//                if (i==1) ovr.setFillColor(Color.ORANGE);
//                if (i==2) ovr.setFillColor(Color.GREEN);
//                if (i==3) ovr.setFillColor(Color.BLUE);
//                ov.add(ovr);

            }

            PolygonRoi plg = new PolygonRoi(plgx, plgy, Roi.FREELINE);
            plg.setStrokeColor(Color.RED);
            plg.setFillColor(Color.RED);
            plg.setStrokeWidth(.25f);
            ov.add(plg);

        }

        for (int i = 0; i < _delin_terminals.size(); i++) {
            OvalRoi ovr = new OvalRoi(_delin_terminals.get(i)[0]-.5+.5, _delin_terminals.get(i)[1]-.5+.5, 1, 1);
            ovr.setFillColor(Color.RED);
            ov.add(ovr);
        }

//        int xc = Math.round(_delin_locs.get(0)[0][0]);
//        int yc = Math.round(_delin_locs.get(0)[1][0]);

//        float[] vals = new float[];
//        for (int xx = xc-Math.round(rad); xx <= xc+Math.round(rad); xx++) {
//            for (int yy = yc-Math.round(rad); yy <= yc+Math.round(rad); yy++) {
//
//            }
//        }

        return ov;

    }

}


