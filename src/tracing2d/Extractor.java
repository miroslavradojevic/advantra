package tracing2d;

import aux.Interpolator;
import aux.Stat;
import ij.IJ;
import ij.gui.*;

import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by miroslav on 29-11-14.
 * extract local number of streamlines, their intensities, average intensities for each streamline
 * possible to multithread for selected number od locations
 */
public class Extractor extends Thread {

    private int begN, endN;

    public  static int[][] 	 i2loc;
    public  static int       R;
    public  static int       Ns;
    public  static float     sigma_deg;
    public  static int       nstreams;
    private static int       dim;

    private static int merging_margin_step = 1; // after this stream iteration index merging of the confluent paths is possible
    private static float merging_margin_dist = 1f;

    public static ArrayList[] i2scores;     // i2scores[0]      is ArrayList<float[]>
    public static ArrayList[] i2locEnd;  // i2terminals[0]      is ArrayList<float[]>
    public static ArrayList[] i2range;// min-max Ni lowest / Ni highest in the neighbourhood defined with the range

    // number of iterations is maximum Ni+1 by definition
    private static float[][] wgts_outer; // set of weights assigned to each streamline length (Ni+1) high towards the end
    private static float[][] wgts_inner; // high at the beginnining

    public Extractor (int n0, int n1)
    {
        begN = n0;
        endN = n1;
    }

    /* initializator for usage at one particular location */
    public static void loadTemplate(
            int     _dim,
            int     _radius,
            int     _nsamples,
            float   _sigma_deg,
            int     _nstreams
    )
    {

        i2loc        = null;
        i2scores     = null;
        i2locEnd     = null;
        i2range      = null;
        dim         = _dim;
        R           = _radius;
        Ns          = _nsamples;
        sigma_deg   = _sigma_deg;
        nstreams    = _nstreams;
        BayesianTracerMulti.init(R);

        weights();

    }

    private  static void weights()
    {
        wgts_outer = new float[BayesianTracerMulti.Ni+1][];
        wgts_inner = new float[BayesianTracerMulti.Ni+1][];

        for (int i = 0; i < wgts_outer.length; i++) {

            wgts_outer[i] = new float[i+1];
            wgts_inner[i] = new float[i+1];
            float sum = 0;

            for (int j = 0; j <= i; j++) {
                float ratio = (i>0)? (float)j/i : 0;
                wgts_outer[i][j] = (float) (1 - Math.exp(-2*ratio));
                sum += wgts_outer[i][j];
            }
            for (int j = 0; j <= i; j++) {
                wgts_outer[i][j] /= sum;
                wgts_inner[i][i-j] = wgts_outer[i][j];
            }

        }
    }

    /* initializator for usage at selected set of locations (threaded) */
    public static void loadTemplate(
            int[][] _i2loc,
            int     _radius,
            int     _nsamples,
            float   _sigma_deg,
            int     _nstreams
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
        R       = _radius;
        Ns      = _nsamples;
        sigma_deg   = _sigma_deg;
        nstreams    = _nstreams;
        BayesianTracerMulti.init(R);

        weights();

    }

    public void run()
    {

        float[][][] xt      = new float[BayesianTracerMulti.Ni+1][Ns][4];    // bayesian filtering states (x,y, vx, vy)
        float[][]   wt      = new float[BayesianTracerMulti.Ni+1][Ns];       // bayesian filtering states
        byte[][]    pt      = new byte[BayesianTracerMulti.Ni+1][Ns];        // bayesian filtering parent pointer - distribution index of the parent state in previous iteration
        boolean[][] et      = new boolean[BayesianTracerMulti.Ni+1][Ns];     // extracting streamlines record



    }

    public static Overlay extractAt(int _x, int _y, float[][] _img_xy, boolean _show_tracks, BayesianTracerMulti.Expansion expan) // make a version for threading
    {

        float[][][] xt      = new float[BayesianTracerMulti.Ni+1][Ns][4];    // bayesian filtering states (x,y, vx, vy)
        float[][]   wt      = new float[BayesianTracerMulti.Ni+1][Ns];       // bayesian filtering states
        byte[][]    pt      = new byte[BayesianTracerMulti.Ni+1][Ns];        // bayesian filtering parent pointer - distribution index of the parent state in previous iteration
        float[][]   mt      = new float[BayesianTracerMulti.Ni+1][Ns];       // measurements, likelihoods - correlations measured during sequential filtering

        BayesianTracerMulti.spherical_wavefront_2d(_x, _y, _img_xy, sigma_deg, Ns, expan, xt, wt, pt, mt);

        ArrayList<float[][]>    delin_locs      = new ArrayList<float[][]>();
        ArrayList<Boolean>      delin_complete  = new ArrayList<Boolean>();
        ArrayList<float[]>      delin_scores    = new ArrayList<float[]>();
        ArrayList<float[]>      delin_values    = new ArrayList<float[]>();
        ArrayList<float[]>      delin_terminals = new ArrayList<float[]>();
        ArrayList<Float>        delin_wascores  = new ArrayList<Float>();
        ArrayList<Float>        delin_wavalues  = new ArrayList<Float>();

        boolean[][] et                          = new boolean[BayesianTracerMulti.Ni+1][Ns];     // extracting streamlines record

        int nrad = extract(_img_xy, nstreams, expan, xt, wt, pt, mt, et, delin_locs, delin_complete, delin_scores, delin_values, delin_terminals, delin_wascores, delin_wavalues);
        float mn = Float.POSITIVE_INFINITY;
        float mx = Float.NEGATIVE_INFINITY;

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

        System.out.println("\n#\tCPLT\tE(V)");
        for (int i = 0; i < delin_complete.size(); i++)
            System.out.println((i+1) + "\t" + delin_complete.get(i) + "\t" + IJ.d2s(delin_wavalues.get(i),3)); //  + Arrays.toString(delin_values.get(i))

        System.out.println("------------------");


        ArrayList<Float> A = new ArrayList<Float>();
        ArrayList<Float> B = new ArrayList<Float>();
        arrange(delin_wavalues, delin_complete, A, B);

        System.out.println("MIN->" + mn);
        System.out.print("A->"); for (int i = 0; i < A.size(); i++) System.out.print(A.get(i) + "  "); System.out.println("");
        System.out.print("B->"); for (int i = 0; i < B.size(); i++) System.out.print(B.get(i) + "  "); System.out.println("");
        System.out.println("MAX->" + mx);

        Overlay ov = new Overlay();
        Overlay ov_xt = viz_xt(xt, wt, pt, et, _show_tracks);
        for (int i = 0; i < ov_xt.size(); i++) ov.add(ov_xt.get(i));
        Overlay ov_delin = viz_delin(delin_locs, delin_scores, delin_values, delin_terminals);
        for (int i = 0; i < ov_delin.size(); i++) ov.add(ov_delin.get(i));

        return ov;

    }

    private static int  extract(
            float[][]               _in_xy,
            int                     _nr_streams,
            BayesianTracerMulti.Expansion _expan,
            float[][][]             _xt,                // states
            float[][]               _wt,                // weights (posteriors)
            byte[][]                _pt,                // parent
            float[][]               _mt,                // measurements (likelihoods, correlations)
            boolean[][]             _et,                // existing in the back track (auxillliary map)
            ArrayList<float[][]>    _delin_locs,
            ArrayList<Boolean>      _delin_complete,
            ArrayList<float[]>      _delin_scores,
            ArrayList<float[]>      _delin_values,
            ArrayList<float[]>      _delin_terminals,
            ArrayList<Float>        _delin_wascores,
            ArrayList<Float>        _delin_wavalues
    )
    {

        float[][] wgts = (_expan== BayesianTracerMulti.Expansion.OUTER)? wgts_outer : (_expan== BayesianTracerMulti.Expansion.INNER)? wgts_inner : null;

        int merging_margin_idx = (_expan== BayesianTracerMulti.Expansion.OUTER)? merging_margin_step : (_expan== BayesianTracerMulti.Expansion.INNER)? _xt.length-1-merging_margin_step : null;

        for (int i = 0; i < _et.length; i++) // reset examined map
            for (int j = 0; j < _et[i].length; j++)
                if (_et[i][j]) _et[i][j] = false;

        _delin_locs.clear();
        _delin_complete.clear();
        _delin_scores.clear();
        _delin_values.clear();
        _delin_terminals.clear();

        int curr_streams = 0;

        boolean all_examined;

        do {

            float max_score;
            int max_idx = -1;
            int max_iter = -1;

            for (int i = _wt.length-1; i > merging_margin_idx; i--) { // loop bayesian iterations backwards

                max_score = Float.NEGATIVE_INFINITY;
                max_idx = -1;
                max_iter = -1;

                for (int j = 0; j < _wt[i].length; j++) { // loop the last iteration - loop all the states/weights from the distribution
                    if (!_et[i][j]) {
                        if (_wt[i][j] > max_score) {
                            max_score = _wt[_wt.length-1][j];
                            max_idx = j;
                            max_iter = i;
                        }
                    }
                }
                if (max_idx!=-1) {
                    break;
                }
            }

            all_examined = false;


            if (max_idx!=-1 && max_iter!=-1) { // max was found in the loop at one index

                boolean is_overlapping = false;
                boolean is_close = false;
                int sidx = max_idx;

                // auxilliary array to fill up during the back track
                float[][]   temp_locs = new float[2][max_iter+1]; // x=0,y=1
                float[]     temp_scores = new float[max_iter+1];
                float[]     temp_values = new float[max_iter+1];

                for (int iidx = max_iter; iidx >=0 ; iidx--) {

                    if (!is_close && iidx>merging_margin_idx) {
                        if (_et[iidx][sidx]) {
                            is_overlapping = true;
                        }
                        else {
                            for (int t = 0; t < _xt[iidx].length; t++) {
                                if (_et[iidx][t] && t!=sidx && Math.pow(_xt[iidx][t][0]-_xt[iidx][sidx][0],2)+Math.pow(_xt[iidx][t][1]-_xt[iidx][sidx][1],2)<=Math.pow(merging_margin_dist,2)) {
                                    is_close = true;
                                    break;
                                }
                            }
                        }
                    }

                    _et[iidx][sidx] = true;

                    if (is_overlapping) break;

                    if (!is_close) { // no need to store them  won't be added anyway, just to fill _et

                        temp_locs[0][iidx] = _xt[iidx][sidx][0];
                        temp_locs[1][iidx] = _xt[iidx][sidx][1];

                        temp_scores[iidx] = _mt[iidx][sidx];

                        temp_values[iidx] = Interpolator.interpolateAt(_xt[iidx][sidx][0], _xt[iidx][sidx][1], _in_xy);

                    }

                    sidx = _pt[iidx][sidx]&0xff;

                }

                if (!is_overlapping && !is_close) { // will occur up to 4 times

                    // locs xy * length
                    float[][]   add_locs    = new float[temp_locs.length][max_iter+1];
                    float[]     add_scores  = new float[max_iter+1];
                    float[]     add_values  = new float[max_iter+1];

                    float add_wascores = 0;
                    float add_wavalues = 0;

                    for (int i = 0; i < max_iter+1; i++) {

                        add_locs[0][i] = temp_locs[0][i];
                        add_locs[1][i] = temp_locs[1][i];
                        add_scores[i] = temp_scores[i];
                        add_values[i] = temp_values[i];

                        add_wascores += wgts[max_iter][i] * add_scores[i];
                        add_wavalues += wgts[max_iter][i] * add_values[i];

                    }

                    _delin_locs.add(add_locs);

                    // complete
                    _delin_complete.add(temp_locs[0].length == _xt.length);

                    // scores
                    _delin_scores.add(add_scores);
                    _delin_wavalues.add(add_wavalues);

                    // values
                    _delin_values.add(add_values);
                    _delin_wascores.add(add_wascores);

                    // terminal locations
                    if (temp_locs[0].length == _xt.length) _delin_terminals.add(new float[]{_xt[max_iter][max_idx][0], _xt[max_iter][max_idx][1]});

                    curr_streams++;

                }

            }
            else all_examined = true;

        }
        while (curr_streams<_nr_streams && !all_examined);

        // return neighbourhood radius
        int nradius = (int) Math.round(Math.sqrt(Math.pow(_delin_terminals.get(0)[0] - _delin_locs.get(0)[0][0], 2) + Math.pow(_delin_terminals.get(0)[1] - _delin_locs.get(0)[1][0], 2)));

        return nradius;

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

    private static Overlay viz_delin(ArrayList<float[][]> _delin_locs, ArrayList<float[]> _delin_scores, ArrayList<float[]> _delin_values, ArrayList<float[]> _delin_terminals)
    {
        Overlay ov = new Overlay();

        for (int i = 0; i < _delin_locs.size(); i++) {

            for (int j = 0; j < _delin_locs.get(i)[0].length; j++) {

//                OvalRoi or = new OvalRoi(delin.get(i).get(j)[0]-1f+.5f, delin.get(i).get(j)[1]-1f+.5f, 2f, 2f);
//                or.setFillColor(Color.RED);
//                or.setStrokeColor(Color.RED);

//                PointRoi pt = new PointRoi(delin.get(i).get(j)[0]+.5f, delin.get(i).get(j)[1]+.5f);
//                pt.setFillColor(Color.RED);
//                pt.setStrokeColor(Color.RED);

                float atx = _delin_locs.get(i)[0][j]+ .5f;
                float aty = _delin_locs.get(i)[1][j]+ .5f;
                float scr = 1;//_delin_values.get(i)[j];

                OvalRoi ovr = new OvalRoi(atx-.5*scr, aty-.5*scr, scr, scr);
                ovr.setFillColor(Color.YELLOW);
                ov.add(ovr);

//                TextRoi tr = new TextRoi(atx-.5*scr, aty-.5*scr, scr, scr, IJ.d2s(scr, 1), new Font("TimesRoman", Font.PLAIN, 1));
//                tr.setStrokeColor(Color.BLUE);
//                ov.add(tr);

            }

//            PolygonRoi plg = new PolygonRoi(plgx, plgy, Roi.FREELINE);
//            plg.setStrokeColor(Color.RED);
//            plg.setFillColor(Color.RED);
//            plg.setStrokeWidth(.25f);
//            ov.add(plg);

        }

        for (int i = 0; i < _delin_terminals.size(); i++) {
            OvalRoi ovr = new OvalRoi(_delin_terminals.get(i)[0]-.25, _delin_terminals.get(i)[1]-.25, .5, .5);
            ovr.setFillColor(Color.RED);
            ov.add(ovr);
        }

        float rad = (float) Math.sqrt(Math.pow(_delin_terminals.get(0)[0]-_delin_locs.get(0)[0][0],2)+Math.pow(_delin_terminals.get(0)[1]-_delin_locs.get(0)[1][0],2));
        OvalRoi circ = new OvalRoi(_delin_locs.get(0)[0][0]-rad, _delin_locs.get(0)[1][0]-rad, 2*rad, 2*rad);
        circ.setStrokeColor(Color.RED);
        ov.add(circ);

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

    private static Overlay viz_xt(float[][][] xt, float[][] wt, byte[][] pt, boolean[][] _et, boolean show_tracks)
    {

        float rad = .5f;

        Overlay ov = new Overlay();
        for (int i = 0; i < xt.length; i++) { // iterations

            Stat.min_max_normalize(wt[i]);
            for (int j = 0; j <xt[i].length; j++) {

                OvalRoi p = new OvalRoi(xt[i][j][0]-rad+.5, xt[i][j][1]-rad+.5, 2*rad, 2*rad);

//                if (i==xt.length-1) // last one in red
//                    p.setFillColor(new Color(1f,0f,0f,wt[i][j]));
//                else
                    p.setFillColor(new Color(1f,1f,1f,wt[i][j]));

//                if (!_et[i][j]) p.setFillColor(new Color(1f, 1f, 0f, wt[i][j]));

                ov.add(p);

//                PointRoi pt = new PointRoi(xt[i][j][0]+.5, xt[i][j][1]+.5);
//                pt.setStrokeColor(new Color(1f,1f,1f,wt[i][j]));
//                pt.setFillColor(new Color(1f, 1f, 1f, wt[i][j]));
//                ov.add(pt);

            }
        }

        if (show_tracks) {

            float[] trackx = new float[xt.length];
            float[] tracky = new float[xt.length];

            for (int i = 0; i < Ns; i++) {
                int curr = i;
                for (int t = xt.length-1; t >= 0; t--) {
                    trackx[t] = xt[t][curr][0] + .5f;
                    tracky[t] = xt[t][curr][1] + .5f;
                    curr = pt[t][curr]&0xff;
                }
                PolygonRoi proi = new PolygonRoi(trackx, tracky, Roi.FREELINE);
                proi.setFillColor(Color.WHITE);
                proi.setStrokeColor(Color.WHITE);
                ov.add(proi);
            }
        }

        return ov;
    }

}
