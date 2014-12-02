package tracing2d;

import aux.Stat;
import ij.gui.*;

import java.awt.*;
import java.util.ArrayList;

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

    private static int merging_margin_idx = 1; // after this stream iteration index merging of the confluent paths is possible
    private static float merging_margin = 1f;

    public static ArrayList[] i2scores;     // i2scores[0]      is ArrayList<Float>
    public static ArrayList[] i2locEnd;  // i2terminals[0]   is ArrayList<float[]>
    public static ArrayList[] i2range;// min-max Ni lowest / Ni highest in the neighbourhood defined with the range

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
    }

    public void run()
    {

        float[][]   tt  = new float[BayesianTracerMulti.Ni+1][2];            // aux to store the traces
        float[][][] xt      = new float[BayesianTracerMulti.Ni+1][Ns][4];    // bayesian filtering states (x,y, vx, vy)
        float[][]   wt      = new float[BayesianTracerMulti.Ni+1][Ns];       // bayesian filtering states
        byte[][]    pt      = new byte[BayesianTracerMulti.Ni+1][Ns];        // bayesian filtering parent pointer - distribution index of the parent state in previous iteration
        boolean[][] et      = new boolean[BayesianTracerMulti.Ni+1][Ns];     // extracting streamlines record



    }

    public static Overlay get_streamlines(int _x, int _y, float[][] _img_xy, boolean _show_tracks)
    {

        float[][][] xt      = new float[BayesianTracerMulti.Ni+1][Ns][4];    // bayesian filtering states (x,y, vx, vy)
        float[][]   wt      = new float[BayesianTracerMulti.Ni+1][Ns];       // bayesian filtering states
        byte[][]    pt      = new byte[BayesianTracerMulti.Ni+1][Ns];        // bayesian filtering parent pointer - distribution index of the parent state in previous iteration
        boolean[][] et      = new boolean[BayesianTracerMulti.Ni+1][Ns];     // extracting streamlines record

        BayesianTracerMulti.run2d(_x, _y, _img_xy,  sigma_deg, Ns,  xt, wt, pt);

        ArrayList<float[][]> delin_locs = new ArrayList<float[][]>();

        Overlay ov = new Overlay();

        extract_delin_locs(nstreams, xt, wt, pt, et, delin_locs);

        System.out.print(delin_locs.size() + "x ");

        Overlay ov_delin = viz_delin(delin_locs);
        for (int i = 0; i < ov_delin.size(); i++) ov.add(ov_delin.get(i));
        Overlay ov_xt = viz_xt(xt, wt, pt, _show_tracks);
        for (int i = 0; i < ov_xt.size(); i++) ov.add(ov_xt.get(i));

        return ov;

    }

    private static void  extract_delin_locs(
            int                     _nr_streams,
            float[][][]             _xt,
            float[][]               _wt,
            byte[][]                _pt,
            boolean[][]             _et,
            ArrayList<float[][]>    _delin_locs
    )
    {

        float[][] temp_locs = new float[2][_xt.length]; // x=0,y=1

        for (int i = 0; i < _et.length; i++) // reset examined map
            for (int j = 0; j < _et[i].length; j++)
                if (_et[i][j]) _et[i][j] = false;

        _delin_locs.clear();

        boolean all_examined;

        do {

            float max_score = Float.NEGATIVE_INFINITY;
            int max_idx = -1;

            for (int j = 0; j < _wt[_wt.length-1].length; j++) { // loop the last iteration - loop all the states/weights from the distribution
                if (!_et[_wt.length-1][j]) {
                    if (_wt[_wt.length-1][j] > max_score) {
                        max_score = _wt[_wt.length-1][j];
                        max_idx = j;
                    }
                }
            }

            all_examined = false;

            if (max_idx!=-1) { // max was found in the loop at one index

                boolean is_overlapping = false;
                boolean is_close = false;
                int sidx = max_idx;

                for (int iidx = _xt.length-1; iidx >=0 ; iidx--) {

                    if (!is_close && iidx>merging_margin_idx) {

                        if (_et[iidx][sidx]) {
                            is_overlapping = true;
                        }
                        else {
                            for (int t = 0; t < _xt[iidx].length; t++) {
                                if (_et[iidx][t] && t!=sidx && Math.pow(_xt[iidx][t][0]-_xt[iidx][sidx][0],2)+Math.pow(_xt[iidx][t][1]-_xt[iidx][sidx][1],2)<=Math.pow(merging_margin,2)) {
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
                    }

                    sidx = _pt[iidx][sidx]&0xff;

                }

                if (!is_overlapping && !is_close) { // will occur up to 4 times

                    float[][] tt = new float[temp_locs.length][temp_locs[0].length];
                    for (int i = 0; i < temp_locs.length; i++) {
                        for (int j = 0; j < temp_locs[i].length; j++) {
                            tt[i][j] = temp_locs[i][j];
                        }
                    }

                    _delin_locs.add(tt);

                }

            }
            else all_examined = true;

        }
        while (_delin_locs.size()<_nr_streams && !all_examined);

    }

    public static void extract(
            int                _nr_streams,
            float[][][]        _xt,
            float[][]          _wt,
            byte[][]           _pt,
            boolean[][]        _et,
            ArrayList<Float>   _scores,
            ArrayList<float[]> _locEnd
    )
    {

    }

    private static Overlay viz_delin(ArrayList<float[][]> _delin_locs)
    {
        Overlay ov = new Overlay();

        for (int i = 0; i < _delin_locs.size(); i++) {

            float[] plgx = _delin_locs.get(i)[0].clone(); // new float[delin.get(i).size()];
            float[] plgy = _delin_locs.get(i)[1].clone(); // new float[delin.get(i).size()];

            for (int j = 0; j < plgx.length; j++) {

//                OvalRoi or = new OvalRoi(delin.get(i).get(j)[0]-1f+.5f, delin.get(i).get(j)[1]-1f+.5f, 2f, 2f);
//                or.setFillColor(Color.RED);
//                or.setStrokeColor(Color.RED);

//                PointRoi pt = new PointRoi(delin.get(i).get(j)[0]+.5f, delin.get(i).get(j)[1]+.5f);
//                pt.setFillColor(Color.RED);
//                pt.setStrokeColor(Color.RED);

                plgx[j] += .5f;
                plgy[j] += .5f;

            }

            PolygonRoi plg = new PolygonRoi(plgx, plgy, Roi.FREELINE);
            plg.setStrokeColor(Color.RED);
            plg.setFillColor(Color.RED);
            plg.setStrokeWidth(.25f);
            ov.add(plg);

        }

        return ov;

    }

    private static Overlay viz_xt(float[][][] xt, float[][] wt, byte[][] pt, boolean show_tracks)
    {

        float rad = .5f;

        Overlay ov = new Overlay();
        for (int i = 0; i < xt.length; i++) { // iterations

            Stat.min_max_normalize(wt[i]);
            for (int j = 0; j <xt[i].length; j++) {

                OvalRoi p = new OvalRoi(xt[i][j][0]-rad+.5, xt[i][j][1]-rad+.5, 2*rad, 2*rad);

                if (i==xt.length-1) // last one in red
                    p.setFillColor(new Color(1f,0f,0f,wt[i][j]));
                else
                    p.setFillColor(new Color(1f,1f,1f,wt[i][j]));

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

    private static Overlay viz_xt

}
