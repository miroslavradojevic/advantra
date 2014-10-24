package tracing2d;

import aux.Interpolator;
import aux.ReadDET;
import aux.Stat;
import conn.Find_Connected_Regions;
import ij.*;
import ij.gui.*;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by miroslav on 6-10-14.
 * takes critical point regions and applies the bayesian tracking between the points
 * ISBI15 submission code
 * inputs are .soma and .det to define critical regions and Bayesian Tracking (copied from Bayesian Tracking class)
 * is used to trace between the points
 */
public class TreeReconstructor2D implements PlugIn {

    BayesianTracer2D    btracer;
    ReadDET             rdet;

    String _image_path;
    String _det_path;    // path to the detection image
    String _soma_path;   // path to the soma binary image file

    ImagePlus curr_img;
    int W, H;

    // tracer variables - each direction will point outwards from the CP region and will contribute with one trace
    // list will have the length equivalent to the # of traces
    int region_label = 0;   // keeps most up to date CP region label
    int soma_label = 0;     // keeps most up to date soma label
    ArrayList<float[]>      seed_xy  = new ArrayList<float[]>();            //
    ArrayList<float[]>      seed_vxy = new ArrayList<float[]>();            //
    ArrayList<float[]>      centroids_xy = new ArrayList<float[]>();        // centers of critical points
    ArrayList<Integer>      label_beg = new ArrayList<Integer>();
    ArrayList<Integer>      label_end = new ArrayList<Integer>();
    ArrayList<float[][][]>  tracer2d_Xt_xy = new ArrayList<float[][][]>();  // (trace#)[trace_len][elements][xy]
    ArrayList<float[][]>    tracer2d_wt_xy = new ArrayList<float[][]>();    // (trace#)[trace_len][elements]
    ArrayList<float[][]>    tracer2d_est_xy= new ArrayList<float[][]>();    // (trace#)[trace_len][xy]


    public void run(String s) {

        // load the image through the menu
        String in_folder = Prefs.get("id.folder", System.getProperty("user.home"));
        OpenDialog.setDefaultDirectory(in_folder);
        OpenDialog dc = new OpenDialog("Select image file");
        in_folder = dc.getDirectory();
        _image_path = dc.getPath();
        if (_image_path==null) return;
        Prefs.set("id.folder", in_folder);

        curr_img = new ImagePlus(_image_path);
        W = curr_img.getWidth();
        H = curr_img.getHeight();

        // variables to read through the menu
        if (Macro.getOptions()==null) {

            _det_path 		= 			Prefs.get("critpoint.tracing2d.det_path", "");
            _soma_path 		= 			Prefs.get("critpoint.tracing2d.soma_path", "");

            GenericDialog gd = new GenericDialog("TreeReconstructor2D");
            gd.addStringField("DetectionFilePath", 			_det_path, 50);
            gd.addStringField("SomaFilePath", 				_soma_path, 50);
            gd.showDialog();
            if (gd.wasCanceled()) return;
            _det_path       	= gd.getNextString(); 				Prefs.set("critpoint.tracing2d.det_path", _det_path);
            _soma_path       	= gd.getNextString(); 				Prefs.set("critpoint.tracing2d.soma_path", _soma_path);

        }
        else {
            _det_path   = Macro.getValue(Macro.getOptions(),  "DetectionFilePath", "");
            _soma_path  = Macro.getValue(Macro.getOptions(),  "SomaFilePath", "");

        }

        // load the critical point regions: .det
        rdet = new ReadDET(_det_path);
//        rdet.print();

        // extract seed points and the directions from .det file

        // create critical point region map, fill it with zeros
        int[][] region_map = new int[W][H];

        label_beg.clear();
        label_end.clear();
        seed_xy.clear();
        seed_vxy.clear();
        centroids_xy.clear();

//        int label = 1;

        region_label = 1;

        // loop locations from the det file and fill the labels in the region_map and generate seed vectors
        int nr_reg = rdet.x.size(); // numner of x coords

        for (int i = 0; i < nr_reg; i++) {

            int xc = Math.round(rdet.x.get(i));
            int yc = Math.round(rdet.y.get(i));
            int rc = (int) Math.ceil(rdet.r.get(i));

            int nr_seeds = rdet.v.get(i).length;

            // add seed points depending on the number of the directions
            // define the new seed points that start at the boundaries of the critical point region, defined with r
            // and point towards the direction defined previously in the critical point detection stage
            // .det will provide (x, y, vx, vy, r)  -> 1: (seedx, seedy, vx, vy)
            //                                      -> 2: (seedx, seedy, vx, vy)
            //                                      ...
            //                                      -> 4: (seedx, seedy, vx, vy)
            for (int j = 0; j < nr_seeds; j++) {

                float[] cxy = new float[]{
                        rdet.x.get(i),
                        rdet.y.get(i)
                };

                float[] sxy = new float[]{
                        rdet.x.get(i) + rdet.r.get(i) * rdet.v.get(i)[j][0],
                        rdet.y.get(i) + rdet.r.get(i) * rdet.v.get(i)[j][1]
                };

                float[] vxy = new float[]{
                        rdet.v.get(i)[j][0],
                        rdet.v.get(i)[j][1]
                };

                seed_xy.add(sxy);
                seed_vxy.add(vxy);
                centroids_xy.add(cxy);
                label_beg.add(region_label);
                label_end.add(Integer.MAX_VALUE); // default when the tag is not concrete

            }

            // region map
            // loop throughout the sphere to get the locations that are within
            for (int xreg = xc-rc; xreg <=xc+rc; xreg++) {
                for (int yreg = yc-rc; yreg <=yc+rc; yreg++) {

                    if (xreg>=0 && xreg<W && yreg>=0 && yreg < H) { // the point is within the dims

                        if ((xreg-xc)*(xreg-xc)+(yreg-yc)*(yreg-yc) <= rc*rc)  {
                            region_map[xreg][yreg] = region_label;
                        }

                    }

                }
            }

            region_label++; // for the next CP region

        }

        // done with the region labels
        System.out.println("active region label: " + region_label);



        // load soma regions: binary .tif file where somas are marked in white
        ImageProcessor somap = new ImagePlus(_soma_path).getProcessor();

        //check
        if ((somap.getWidth()!=W) || (somap.getHeight()!=H))
        {
            System.out.println("soma image error: dimensions do not match with the original");;
            return;
        }

        // convert it to a binary byte processor
        ByteProcessor score = new ByteProcessor(somap.getWidth(), somap.getHeight());
        for (int ii=0; ii<somap.getWidth()*somap.getHeight(); ii++) if (somap.getf(ii) > 0) score.set(ii, 255);
        Find_Connected_Regions conn_reg = new Find_Connected_Regions(new ImagePlus("", score), true);  // true means save locations
        conn_reg.run("");


        ArrayList<ArrayList<int[]>> soma_conn =  conn_reg.getConnectedRegions();


        // add somas top the map - each connected region will be assigned a negative label: -1, -2, -3...
        soma_label = -1;
        for (int soma_idx = 0; soma_idx < soma_conn.size(); soma_idx++) {

            for (int loc_idx = 0; loc_idx < soma_conn.get(soma_idx).size(); loc_idx++) {
                int xcoord = soma_conn.get(soma_idx).get(loc_idx)[1];
                int ycoord = soma_conn.get(soma_idx).get(loc_idx)[0];
                region_map[xcoord][ycoord] = soma_label;
            }

            soma_label = soma_label - 1;
        }

        System.out.println("active soma label: " + soma_label);

        // form trace starting points

        // tracker parameters
        String          scales_list         = "3,5,7,9"; // these parameters can be added to the dialog
        float           D                   = 5f;
        float           prior_sigma_deg     = 65f;
        int             Nt                  = 50;
        int             MAX_ITER            = 100;
        int             MIN_ITER            = 5;
        boolean         USE_ORIGINAL        = true;

        btracer = new BayesianTracer2D(
                scales_list,
                curr_img,
                D,
                prior_sigma_deg,
                Nt,
                MAX_ITER,
                MIN_ITER,
                USE_ORIGINAL
        );

        tracer2d_Xt_xy.clear();
        tracer2d_wt_xy.clear();
        tracer2d_est_xy.clear();

        int count_bdy = 0; // counts those that were reaching the body (label is -1, -2, -3...)
        int count_cpr = 0; // counts those that were reaching the other CP region (label is >0)
        int count_out = 0; // counts those that were either reaching the background or looped till the end without reaching anything concrete or they reached the border

        /*
        tracing from each seed
         */

        for (int seed_idx = 0; seed_idx < seed_xy.size(); seed_idx++) {


            //////////////////////////////////
            System.out.print("|--- trace " + seed_idx + "\t");
            int out_label = btracer.run(        // bayesian recursive tracing
                                                // to start the trace
                                                seed_xy.get(seed_idx)[0],
                                                seed_xy.get(seed_idx)[1],
                                                seed_vxy.get(seed_idx)[0],
                                                seed_vxy.get(seed_idx)[1],

                                                label_beg.get(seed_idx),
                                                region_map     // to stop the trace

            );
            label_end.set(seed_idx, out_label);
            System.out.println(" ------> " + out_label);
            //////////////////////////////////




            if (out_label<0)                            count_bdy++;
            else if (out_label==Integer.MAX_VALUE)      count_out++;
            else                                        count_cpr++;





            /*
            add the traced variables
             */




            // *** btracer.Xt_xy -> tracer2d_Xt_xy
            float[][][] _Xt_xy = new float[btracer.Xt_xy.size()][][];

            for (int i = 0; i < btracer.Xt_xy.size(); i++) {

                _Xt_xy[i] = new float[btracer.Xt_xy.get(i).length][2];

                for (int j = 0; j < btracer.Xt_xy.get(i).length; j++) {
                    _Xt_xy[i][j][0] = btracer.Xt_xy.get(i)[j][0];
                    _Xt_xy[i][j][1] = btracer.Xt_xy.get(i)[j][1];
                }

            }

            tracer2d_Xt_xy.add(_Xt_xy);





            // *** btracer.w_xy -> tracer2d_wt_xy
            float[][] _wt_xy = new float[btracer.wt_xy.size()][];

            for (int i = 0; i < btracer.wt_xy.size(); i++) {

                _wt_xy[i] = new float[btracer.wt_xy.get(i).length];

                for (int j = 0; j < btracer.wt_xy.get(i).length; j++) {
                    _wt_xy[i][j] = btracer.wt_xy.get(i)[j];
                }

            }

            tracer2d_wt_xy.add(_wt_xy);





            // *** btracer.est_xy -> tracer2d_est_xy
            float[][] _est_xy = new float[btracer.est_xy.size()][];

            for (int i = 0; i < btracer.est_xy.size(); i++) {

                _est_xy[i] = new float[btracer.est_xy.get(i).length];

                for (int j = 0; j < btracer.est_xy.get(i).length; j++) {
                    _est_xy[i][j] = btracer.est_xy.get(i)[j];
                }

            }

            tracer2d_est_xy.add(_est_xy);

        } // done tracing

        System.out.println("CPRG: " + count_cpr);
        System.out.println("SOMA: " + count_bdy);
        System.out.println("NONE: " + count_out);

//        Overlay dbgov = new Overlay();
//        ImagePlus rmap = new ImagePlus("region_map", new FloatProcessor(region_map));
//        rmap.show();

        /*
        correction for the branches that were NONE (those that ended with Integer.MAX_VALUE)
         */
        ImageStack is_traces = new ImageStack(528,255); // for plots
        for (int i = 0; i < label_end.size(); i++) {
            if (label_end.get(i)==Integer.MAX_VALUE) { // if it was finished ?

                float[] trackx = extract(tracer2d_est_xy.get(i), 0);
                float[] tracky = extract(tracer2d_est_xy.get(i), 1);

                float[] extr = extract_intensities(btracer.likelihood_xy, trackx, tracky);
                float[] cost = get_costs(extr);

                int break_idx = find_end_point_index(extr); // >0

                if (break_idx>1) {

                    // crop it
                    int length_to_keep = break_idx+1;
                    float[][] broken_tarce = new float[length_to_keep][2];
                    System.arraycopy(tracer2d_est_xy.get(i), 0, broken_tarce, 0, length_to_keep); // src, src_pos, dest, dest_pos, len
                    // add new valid CP label for the end point

                    tracer2d_est_xy.set(i, broken_tarce);
                    label_end.set(i, region_label);
                    region_label++;


                }
                else { // the length was not big enough, set it as null (overlay output will not plot it in that case)

                    tracer2d_est_xy.set(i, null);
                    label_end.set(i, Integer.MAX_VALUE); // was there before as well

                }






                float[] extrx = new float[extr.length];
                for (int j = 0; j < extr.length; j++) extrx[j] = j;

                Plot p = new Plot("", "", "", extrx, extr, Plot.BOX);
                p.setLimits(0, extrx.length, Math.min(Stat.get_min(extr), Stat.get_min(cost)), Math.max(Stat.get_max(extr), Stat.get_max(cost)) );
                p.setColor(Color.RED);
                p.addPoints(extrx, extr, Plot.LINE);
                p.draw();
                p.setColor(Color.GREEN);
                p.addPoints(extrx,cost, Plot.LINE);
                is_traces.addSlice("Break.Idx="+break_idx, p.getProcessor());

            }
        }

        new ImagePlus("", is_traces).show();


        // show the trace CP labels
        for (int i = 0; i < ; i++) {
            
        }



        // show stack with intensities along
        ImagePlus orig = curr_img.duplicate();
        orig.setTitle("orig.");

        Overlay ovl = new Overlay();

        Overlay trace_Xt_xy = getTrace_Xt_xy();
        for (int i = 0; i < trace_Xt_xy.size(); i++) ovl.add(trace_Xt_xy.get(i));

        Overlay trace_est_xy = getTrace_est_xy();
        for (int i = 0; i < trace_est_xy.size(); i++) ovl.add(trace_est_xy.get(i));

        Overlay trace_seeds = getSeeds();
        for (int i = 0; i < trace_seeds.size(); i++) ovl.add(trace_seeds.get(i));
        
        orig.setOverlay(ovl);
        orig.show();
        IJ.setTool("hand");

    }

    private static float[] extract(float[][] _trajectory_N_xy, int _idx_xt) {

        // rows are points along the trajectory
        float[] track = new float[_trajectory_N_xy.length];

        // columns are x,y, _idx_xt will choose which of those
        for (int i = 0; i < track.length; i++) track[i] = _trajectory_N_xy[i][_idx_xt];

        return track;

    }

    private float[] extract_intensities(float[][] _source_xy, float[] _path_x, float[] _path_y) {

        // allocate array of the size
        float[] out = new float[_path_x.length];
        
        // loop x,y and sample interpolated values
        for (int i = 0; i < _path_x.length; i++) {
            out[i] = Interpolator.interpolateAt(_path_x[i], _path_y[i], _source_xy);
        }

        return out;

    }

    private float[] get_costs(float[] _centerline_intensities) {

        int len = _centerline_intensities.length;

        // total sum
        float summ = 0;
        for (int i = 0; i < len; i++) {
            summ += _centerline_intensities[i];
        }

        float sum_before = _centerline_intensities[0];
        float sum_after = summ - sum_before;

        float[] score = new float[len];
        score[0] = 0;
        score[len-1] = 0;

        // loop once to get cummulative sums
        for (int i = 1; i < len-1; i++) {
            sum_before += _centerline_intensities[i];
            sum_after = summ - sum_before;
            score[i] = sum_before / (float)(i+1) - sum_after / (len-(i+1));
        }

        return score;

    }

    private int find_end_point_index(float[] _centerline_intensities){ // output is the index where the stream ends

        int len = _centerline_intensities.length;

        // total sum
        float summ = 0;
        for (int i = 0; i < len; i++) {
            summ += _centerline_intensities[i];
        }

        float sum_before = _centerline_intensities[0];
        float sum_after = summ - sum_before;

        int break_idx = 0;
        float max_score = sum_before / 1f - sum_after / (len-1f);
        
        // loop once to get cummulative sums
        for (int i = 1; i < len-1; i++) {
            sum_before += _centerline_intensities[i];
            sum_after = summ - sum_before;
            float curr_score = sum_before / (float)(i+1) - sum_after / (len-(i+1));
            if (curr_score>max_score) {
                max_score = curr_score;
                break_idx = i;
            }
        }
        
        return break_idx;

    }

    private Overlay getTrace_Xt_xy(){

        Overlay ov = new Overlay();

        // read tracer2d_Xt_xy and output the states

        for (int i = 0; i < tracer2d_Xt_xy.size(); i++) { // i - trace

            if (label_end.get(i)!=Integer.MAX_VALUE){

                for (int j = 0; j < tracer2d_Xt_xy.get(i).length; j++) {  // j - iteration along the trace


                    // weights at this iteration will define the radiuses
                    float[] radiuses = tracer2d_wt_xy.get(i)[j].clone();
                    Stat.normalize(radiuses); // max. radius becomes 1

                    for (int k = 0; k < tracer2d_Xt_xy.get(i)[j].length; k++) { // k - one particle state

                        Color colr = new Color(1f,1f,1f,radiuses[k]);
                        OvalRoi overlay_Xt = new OvalRoi(
                                tracer2d_Xt_xy.get(i)[j][k][0]+.5f - 0.5f*radiuses[k],
                                tracer2d_Xt_xy.get(i)[j][k][1]+.5f - 0.5f*radiuses[k],
                                radiuses[k], radiuses[k]);
                        overlay_Xt.setStrokeColor(colr);
                        overlay_Xt.setFillColor(colr);

                        ov.add(overlay_Xt);

                    }

                }

            }

        }

        return ov;

    }

    private Overlay getTrace_est_xy(){

        Overlay ov = new Overlay();

        for (int i = 0; i < tracer2d_est_xy.size(); i++) { // i - trace

            if (tracer2d_est_xy.get(i)!=null) { // if it is null then -> nothing



            float[] trc_x = new float[tracer2d_est_xy.get(i).length];
            float[] trc_y = new float[tracer2d_est_xy.get(i).length];

            for (int j = 0; j < tracer2d_est_xy.get(i).length; j++) { // j - iteration along the trace

                trc_x[j] = tracer2d_est_xy.get(i)[j][0]   + .5f;
                trc_y[j] = tracer2d_est_xy.get(i)[j][1]   + .5f;

            }

            if (label_end.get(i)!=Integer.MAX_VALUE){
                // paint it yellow if ok
                PolygonRoi trc = new PolygonRoi(trc_x, trc_y, PolygonRoi.FREELINE);
                trc.setStrokeWidth(1f);
                trc.setStrokeColor(Color.YELLOW);
                trc.setFillColor(Color.YELLOW);
                ov.add(trc);
            }
            else {
                // paint it blue if wrong
                PolygonRoi trc = new PolygonRoi(trc_x, trc_y, PolygonRoi.FREELINE);
                trc.setStrokeWidth(1f);
                trc.setStrokeColor(Color.BLUE);
                trc.setFillColor(Color.BLUE);
                ov.add(trc);
            }

            }

        }


        return ov;

    }

    private Overlay getSeeds(){

        float arrow_length = 5f;
        float arrow_width = 2f;

        Overlay ov = new Overlay();

        for (int i = 0; i < seed_xy.size(); i++) { // loop all the seeds

            Line seed_arrow = new Line(
                    seed_xy.get(i)[0],
                    seed_xy.get(i)[1],
                    seed_xy.get(i)[0] + arrow_length * seed_vxy.get(i)[0],
                    seed_xy.get(i)[1] + arrow_length * seed_vxy.get(i)[1]
            );


            seed_arrow.setFillColor(Color.RED);
            seed_arrow.setStrokeColor(Color.RED);
            seed_arrow.setStrokeWidth(arrow_width);
            ov.add(seed_arrow);

        }

        return ov;
    }

    private void getSwc(){
        // labels will be used to distinct the nodes
    }

}
