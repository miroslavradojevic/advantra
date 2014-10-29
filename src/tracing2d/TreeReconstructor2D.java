package tracing2d;

import aux.Interpolator;
import aux.ReadDET;
import aux.Stat;
import com.sun.org.apache.bcel.internal.generic.FLOAD;
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
import java.io.*;
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

    String 		image_dir;
    String		image_name;

    ImagePlus curr_img;
    int W, H;

    // tracer variables - each direction will point outwards from the CP region and will contribute with one trace
    // list will have the length equivalent to the # of traces
    int region_label = 0;   // keeps most up to date CP region label - used for logging regions both CP and soma

    ArrayList<float[]>      seed_xy  = new ArrayList<float[]>();            //
    ArrayList<float[]>      seed_vxy = new ArrayList<float[]>();            //

    ArrayList<float[]>      cpregion_xy = new ArrayList<float[]>();        // as many as there are critical points

    ArrayList<Integer>      label_beg = new ArrayList<Integer>();
    ArrayList<Integer>      label_end = new ArrayList<Integer>();

    int soma_label = 0;     // keeps most up to date soma label
    ArrayList<float[]>      soma_xyr = new ArrayList<float[]>();
    ArrayList<Integer>      label_soma = new ArrayList<Integer>();

    ArrayList<float[][][]>  tracer2d_Xt_xy      = new ArrayList<float[][][]>();  // (trace#)[trace_len][elements][xy]
    ArrayList<float[][]>    tracer2d_wt_xy      = new ArrayList<float[][]>();    // (trace#)[trace_len][elements]
    ArrayList<float[][]>    tracer2d_est_xy     = new ArrayList<float[][]>();    // (trace#)[trace_len][xy]

    ArrayList<float[]>      tracer2d_vals       = new ArrayList<float[]>();      // (trace#)[trace_len] values sampled from the image along the trace
    ArrayList<Float>        tracer2d_vals_med   = new ArrayList<Float>();        // (trace#) median of each of the values

    // inputs - parameters
    String _image_path;  // path to the original image (selected throughout menu)
    String _det_path;    // path to the detection file (with marked CP regions - taken from experiments from previous works)
    String _soma_path;   // path to the soma binary image file (somas are also in separate binary .tif files - extracted using soma.ijm script)

    // tracker parameters
    float           _D                   ;    // will define the step in Bayesian tracking scheme - size of the sampling sphere
    String          _scales_list         ;    // scales used when calculating neuriteness
    float           _prior_sigma_deg     ;    // gaussian sigma for the angular prediction
    int             _Nt                  ;    // number of samples to maintain throughout recursive estimation
    int             _MAX_ITER            ;    // limit number of iterations when tracing recursively
    int             _MIN_ITER            ;    // minimum amount of iterations in recursive tracing to elapse before considering joining some other region
    boolean         _USE_ORIGINAL        ;    // use original as input measurement instead of neuriteness

    // some constants (used to correct after tracking)
    int             _NR_SAMPLES_MARGIN   ;

    // interface parameters
    boolean         _show_it             ;

    public void run(String s) {

        // load the image
        String in_folder = Prefs.get("id.folder", System.getProperty("user.home"));
        OpenDialog.setDefaultDirectory(in_folder);
        OpenDialog dc = new OpenDialog("Select image file");
        in_folder = dc.getDirectory();
        _image_path = dc.getPath();
        if (_image_path==null) return;
        Prefs.set("id.folder", in_folder);

        curr_img = new ImagePlus(_image_path);

        image_dir   = curr_img.getOriginalFileInfo().directory; //  + File.separator  + image_name
        image_name  = curr_img.getShortTitle();

        W = curr_img.getWidth();
        H = curr_img.getHeight();

        if (Macro.getOptions()==null) { // variables to read through the menu

//            _det_path 		= 			Prefs.get("critpoint.tracing2d.det_path", "");
//            _soma_path 		= 			Prefs.get("critpoint.tracing2d.soma_path", "");
            _D              =   (float) Prefs.get("critpoint.tracing2d.d", 5f);
            _scales_list    =           Prefs.get("critpoint.tracing2d.scales_list", "3,5,7,9");
            _prior_sigma_deg=   (float) Prefs.get("critpoint.tracing2d.prior_sigma_deg", 65f);
            _Nt             =   (int)   Prefs.get("critpoint.tracing2d.Nt", 50);
            _MAX_ITER       =   (int)   Prefs.get("critpoint.tracing2d.max_iter", 200);
            _MIN_ITER       =   (int)   Prefs.get("critpoint.tracing2d.min_iter", 5);
            _USE_ORIGINAL   =           Prefs.get("critpoint.tracing2d.use_original", true);

            _NR_SAMPLES_MARGIN = (int)  Prefs.get("critpoint.tracing2d.margin_iter", 5);

            _show_it        =           Prefs.get("critpoint.tracing2d.show_it", true);

            GenericDialog gd = new GenericDialog("TreeReconstructor2D");

            gd.addStringField( "DetectionFilePath", 		image_dir + image_name + ".det",          60);
            gd.addStringField( "SomaFilePath", 				image_dir + image_name + ".soma.tif",          60);

            gd.addMessage("---> BAYESIAN TRACKING <---");

            gd.addNumericField("D", _D, 2, 10, "[pix.]");
            gd.addStringField("ScalesList", _scales_list, 30);
            gd.addNumericField("AngularPriorSigma",        _prior_sigma_deg,    2,  10, "[deg.]");
            gd.addNumericField("Nt",                       _Nt,                 0,  10, "count");
            gd.addNumericField("MAX_ITER",                 _MAX_ITER,           0,  15, "");
            gd.addNumericField("MIN_ITER",                 _MIN_ITER,           0,  15, "");
            gd.addCheckbox(    "USE_ORIGINAL", 		       _USE_ORIGINAL);

            gd.addMessage("\tCORRECTION");
            gd.addNumericField("MARGIN_ITER",              _NR_SAMPLES_MARGIN,  0,  15, "");

            gd.addMessage("\tINTERFACE");
            gd.addCheckbox(    "SHOW_IT", 		           _show_it);

            gd.showDialog();
            if (gd.wasCanceled()) return;

            _det_path       	= gd.getNextString(); 				//Prefs.set("critpoint.tracing2d.det_path",        _det_path);
            _soma_path       	= gd.getNextString(); 				//Prefs.set("critpoint.tracing2d.soma_path",       _soma_path);
            ////------------
            _D                  = (float) gd.getNextNumber();       Prefs.set("critpoint.tracing2d.d",               _D);
            _scales_list        =         gd.getNextString();       Prefs.set("critpoint.tracing2d.scales_list",     _scales_list);
            _prior_sigma_deg    = (float) gd.getNextNumber();       Prefs.set("critpoint.tracing2d.prior_sigma_deg", _prior_sigma_deg);
            _Nt                 = (int)   gd.getNextNumber();       Prefs.set("critpoint.tracing2d.Nt",              _Nt);
            _MAX_ITER           = (int)   gd.getNextNumber();       Prefs.set("critpoint.tracing2d.max_iter",        _MAX_ITER);
            _MIN_ITER           = (int)   gd.getNextNumber();       Prefs.set("critpoint.tracing2d.min_iter",        _MIN_ITER);
            _USE_ORIGINAL       =         gd.getNextBoolean();      Prefs.set("critpoint.tracing2d.use_original",    _USE_ORIGINAL);
            ////------------
            _NR_SAMPLES_MARGIN  = (int)   gd.getNextNumber();       Prefs.set("critpoint.tracing2d.margin_iter",     _NR_SAMPLES_MARGIN);
            ////------------
            _show_it            =         gd.getNextBoolean();      Prefs.set("critpoint.tracing2d.show_it",         _show_it);

        }
        else {
            _det_path   =                          Macro.getValue(Macro.getOptions(),  "detectionfilepath", "");
            _soma_path  =                          Macro.getValue(Macro.getOptions(),  "somafilepath", "");

            _D                  = Float.valueOf(   Macro.getValue(Macro.getOptions(),  "d",                String.valueOf(5f)));
            _scales_list        =                  Macro.getValue(Macro.getOptions(),  "scaleslist", "");
            _prior_sigma_deg    = Float.valueOf(   Macro.getValue(Macro.getOptions(),  "angularpriorsigma", String.valueOf(65f)));
            _Nt                 = Integer.valueOf( Macro.getValue(Macro.getOptions(),  "nt",                String.valueOf(50)));
            _MAX_ITER           = Integer.valueOf( Macro.getValue(Macro.getOptions(),  "max_iter",          String.valueOf(200)));
            _MIN_ITER           = Integer.valueOf( Macro.getValue(Macro.getOptions(),  "min_iter",          String.valueOf(5)));
            _USE_ORIGINAL       = Boolean.valueOf( Macro.getValue(Macro.getOptions(),  "use_original",      String.valueOf(true)));
            _NR_SAMPLES_MARGIN  = Integer.valueOf( Macro.getValue(Macro.getOptions(),  "margin_iter",       String.valueOf(5)));
            _show_it            = Boolean.valueOf( Macro.getValue(Macro.getOptions(),  "show_it",           String.valueOf(true)));

        }

//        System.out.println("D="+_D);
//        System.out.println("SCALES_LIST="+_scales_list);
//        System.out.println("AngularPriorSigma="+_prior_sigma_deg);
//        System.out.println("Nt="+_Nt);
//        System.out.println("MAX_ITER="+_MAX_ITER);
//        System.out.println("MIN_ITER="+_MIN_ITER);
//        System.out.println("USE_ORIGINAL="+_USE_ORIGINAL);
//        System.out.println("MARGIN_ITER="+_NR_SAMPLES_MARGIN);
//        System.out.println("SHOW_IT="+_show_it);
//        if (true) return;

        rdet = new ReadDET(_det_path); // load the critical point regions: .det

        /*
            definiton of the BayesianTracer - it will contain method to trace startign from the specified location
            and contain Xt_xy, wt_xy, and est_xy variables that will have states, weight and estimates obtained
            during probabilistic tracing
         */

        btracer = new BayesianTracer2D(
                _scales_list,
                curr_img,
                _D,
                _prior_sigma_deg,
                _Nt,
                _MAX_ITER,
                _MIN_ITER,
                _USE_ORIGINAL
        );

        ////////////////////////////////////////////////////////

        int[][] region_map = new int[W][H];  // create critical point region map, fill it with zeros

        label_beg.clear();
        label_end.clear();
        seed_xy.clear();
        seed_vxy.clear();
        cpregion_xy.clear();

        region_label = 1;

        // loop locations from the det file and fill the labels in the region_map and generate seed vectors: generate seeds, fill the labels...
        for (int i = 0; i < rdet.x.size(); i++) {

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

                label_beg.add(region_label);
                label_end.add(Integer.MAX_VALUE); // default when the tag is not concrete

            }

            // region map : loop throughout the sphere to get the locations that are within
            for (int xreg = xc-rc; xreg <=xc+rc; xreg++) {
                for (int yreg = yc-rc; yreg <=yc+rc; yreg++) {

                    if (xreg>=0 && xreg<W && yreg>=0 && yreg < H) { // the point is within the dims

                        if ((xreg-xc)*(xreg-xc)+(yreg-yc)*(yreg-yc) <= rc*rc)  {
                            region_map[xreg][yreg] = region_label;
                        }

                    }

                }
            }


            float[] cxy = new float[]{rdet.x.get(i), rdet.y.get(i)};

            cpregion_xy.add(cxy);

            region_label++; // for the next CP region

        } // done extracting the seeds

        ////////////////////////////////////////////////////////////////

        // load soma regions: binary .tif file where somas are marked in white
        ImageProcessor somap = new ImagePlus(_soma_path).getProcessor();

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
        soma_xyr.clear();
        label_soma.clear();

        for (int soma_idx = 0; soma_idx < soma_conn.size(); soma_idx++) {

            // find centroid and assign labels to the region map
            float xs = 0, ys = 0;

            for (int loc_idx = 0; loc_idx < soma_conn.get(soma_idx).size(); loc_idx++) {

                int xcoord = soma_conn.get(soma_idx).get(loc_idx)[1];
                int ycoord = soma_conn.get(soma_idx).get(loc_idx)[0];

                xs += xcoord;
                ys += ycoord;

                region_map[xcoord][ycoord] = soma_label;

            }

            // centroid
            xs = xs / soma_conn.get(soma_idx).size();
            ys = ys / soma_conn.get(soma_idx).size();

            // once more to estimate radius

            float rs = 0;

            for (int loc_idx = 0; loc_idx < soma_conn.get(soma_idx).size(); loc_idx++) {

                int xcoord = soma_conn.get(soma_idx).get(loc_idx)[1];
                int ycoord = soma_conn.get(soma_idx).get(loc_idx)[0];

                float radius = (float) Math.sqrt( (xcoord-xs)*(xcoord-xs) + (ycoord-ys)*(ycoord-ys) );

                rs += radius;

            }

            rs = rs / soma_conn.get(soma_idx).size();

            soma_xyr.add(new float[]{xs, ys, rs});

            label_soma.add(soma_label);

            soma_label = soma_label - 1;

        }

//        System.out.println("----------SOMA------------:");
//        for (int i = 0; i < label_soma.size(); i++)
//            System.out.print(label_soma.get(i) + "\t" + Arrays.toString(soma_xyr.get(i)));
//        System.out.println("\n---------------------------");

        ///////////////////////////////////////////////////////////////////

        tracer2d_Xt_xy.clear();
        tracer2d_wt_xy.clear();
        tracer2d_est_xy.clear();

        /*
        tracing from each seed
         */
        System.out.print("TRACING:\n[");
        for (int seed_idx = 0; seed_idx < seed_xy.size(); seed_idx++) {


            //////////////////////////////////
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
            System.out.print(".");
            //////////////////////////////////

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
        System.out.println("]");

        // those  traces that failed ar the ones with label_end.get(i) == Integer.MAX_VALUE (next step will do some corrections)

        ///////////////////////////////////////////////////////////////////

        /*
        tracer2d_vals
        tracer2d_vals_med
         */
        tracer2d_vals.clear();
        tracer2d_vals_med.clear();
        float MIN_MED = Float.POSITIVE_INFINITY;
        for (int i = 0; i < label_end.size(); i++) {
            if (label_end.get(i)!=Integer.MAX_VALUE) { // it was ok
                float[] trackx = extract(tracer2d_est_xy.get(i), 0);
                float[] tracky = extract(tracer2d_est_xy.get(i), 1);
                float[] extr = extract_intensities(btracer.likelihood_xy, trackx, tracky);
                tracer2d_vals.add(extr);
                float curr_med = Stat.median(extr);
                tracer2d_vals_med.add(curr_med);
                if(curr_med<MIN_MED) MIN_MED = curr_med;
            }
            else {
                tracer2d_vals.add(null);
                tracer2d_vals_med.add(Float.NaN);
            }
        }

        /*
        correction for the branches that were NONE (those that ended with Integer.MAX_VALUE)
         */
            //ImageStack is_traces = new ImageStack(528,255); // for plots
        for (int i = 0; i < label_end.size(); i++) {
            if (label_end.get(i)==Integer.MAX_VALUE) { // if it was finished ?

                float[] trackx = extract(tracer2d_est_xy.get(i), 0);
                float[] tracky = extract(tracer2d_est_xy.get(i), 1);
                float[] extr = extract_intensities(btracer.likelihood_xy, trackx, tracky);
                float[] cost = get_costs(extr);
                int break_idx = find_end_point_index(extr);

                if (break_idx>_NR_SAMPLES_MARGIN && break_idx<extr.length-_NR_SAMPLES_MARGIN) {

                    int length_to_keep = break_idx;
                    float[] broken_vals = new float[length_to_keep];
                    System.arraycopy(extr, 0, broken_vals, 0, length_to_keep);
                    float broken_vals_med = Stat.median(broken_vals);

                    if (broken_vals_med>MIN_MED) {

                        // add it - all the conditions were there
                        float[][] broken_tarce = new float[length_to_keep][2];
                        System.arraycopy(tracer2d_est_xy.get(i), 0, broken_tarce, 0, length_to_keep); // src, src_pos, dest, dest_pos, len


                        float[] cxy = new float[]{tracer2d_est_xy.get(i)[break_idx][0], tracer2d_est_xy.get(i)[break_idx][1]}; // last one before it is replaced

                        // add new valid CP label for the end point
                        tracer2d_est_xy.set(i, broken_tarce);
                        tracer2d_vals.set(i, broken_vals);
                        tracer2d_vals_med.set(i, broken_vals_med);
                        label_end.set(i, region_label); // update labels

                        cpregion_xy.add(cxy);

                        region_label++;

                    }
                    else { // median was below the loewst, set it as null
                        //tracer2d_vals, tracer2d_vals_med will stay null, NaN
                        tracer2d_est_xy.set(i, null);
                        label_end.set(i, Integer.MAX_VALUE);
                    }

                }
                else { // the length was not big enough, set it as null (overlay output will not plot it in that case)
                    //tracer2d_vals, tracer2d_vals_med will stay null, NaN
                    tracer2d_est_xy.set(i, null);
                    label_end.set(i, Integer.MAX_VALUE); // was there before as well
                }

//                    float[] extrx = new float[extr.length];
//                    for (int j = 0; j < extr.length; j++) extrx[j] = j;
//                    Plot p = new Plot("", "", "", extrx, extr, Plot.BOX);
//                    p.setLimits(0, extrx.length, Math.min(Stat.get_min(extr), Stat.get_min(cost)), Math.max(Stat.get_max(extr), Stat.get_max(cost)) );
//                    p.setColor(Color.RED);
//                    p.addPoints(extrx, extr, Plot.LINE);
//                    p.draw();
//                    p.setColor(Color.GREEN);
//                    p.addPoints(extrx,cost, Plot.LINE);
//                    is_traces.addSlice("Break.Idx="+break_idx, p.getProcessor());

            }
        }

        // show the outcome of the correction
        int total = 0;
        for (int i = 0; i < label_end.size(); i++) {
            if (label_end.get(i)!=Integer.MAX_VALUE) {
                total++;
            }
        }
        System.out.println("" + total + " / " + label_end.size() + " after correction");

//        new ImagePlus("", is_traces).show();

//        System.out.println("----------TRACES------------:");
//        for (int i = 0; i < label_beg.size(); i++)
//            System.out.println(label_beg.get(i) + "\t" + ((label_end.get(i)==Integer.MAX_VALUE)?"NA":label_end.get(i)) + "\t\t" +
//                    Arrays.toString(cpregion_xy.get(label_beg.get(i)-1)) + "\t" + ((tracer2d_est_xy.get(i)!=null)?tracer2d_est_xy.get(i).length:"NULL") + " pts \t" ); // + .length
//        System.out.println("\n---------------------------");

        ///////////////////////////////////////////////////////////////////
        // export swc: maintain the list with active labels
        generate_swc();
        ///////////////////////////////////////////////////////////////////

        if (_show_it) {

            // show stack with intensities along
            ImagePlus orig = curr_img.duplicate();
            orig.setTitle("TreeReconstructor2D");

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

    private void generate_swc()
    {

        int nr_somas = -soma_label -1; // decrease one because of the counter
        int nr_regions = region_label -1; // counter had one more count

        String output_dir_name = image_dir+String.format(
            "REC.D.scalesList.priorSigmaDeg.Nt.maxIter.minIter.useOrig.nrSampMargin_"+"%.2f_%s_%.2f_%d_%d_%d_%b_%d",
                _D,
                _scales_list,
                _prior_sigma_deg,
                _Nt,
                _MAX_ITER,
                _MIN_ITER,
                _USE_ORIGINAL,
                _NR_SAMPLES_MARGIN,
                _show_it
        );

        // create output directory
        File f = new File(output_dir_name);
        if (!f.exists()) f.mkdirs();

        String det_path = output_dir_name + File.separator + image_name+".swc";

        PrintWriter logWriter = null;
        try {
            logWriter = new PrintWriter(det_path);
            logWriter.print("");
            logWriter.close();
        } catch (FileNotFoundException ex) {}

        try {
            logWriter = new PrintWriter(new BufferedWriter(new FileWriter(det_path, true)));
            logWriter.println("# REC2D...\n");
        } catch (IOException e) {}



        ArrayList<Boolean>      tracer2d_addedswc = new ArrayList<Boolean>(); // tag saying if the trace was added
        tracer2d_addedswc.clear();
        for (int i = 0; i < label_beg.size(); i++) tracer2d_addedswc.add(false);

        /*
        swc format:
        ID  TYPE    X   Y   Z   R   MOTHER_ID
         */

        ArrayList<Integer> frontier = new ArrayList<Integer>();
        ArrayList<Integer> examined = new ArrayList<Integer>();

        int[] somaLab2swcId = new int[nr_somas];
        Arrays.fill(somaLab2swcId, Integer.MAX_VALUE);

        int[] regLab2swcId = new int[nr_regions];
        Arrays.fill(regLab2swcId, Integer.MAX_VALUE);

        int swc_id = 1;

        for (int soma_idx = 0; soma_idx < soma_xyr.size(); soma_idx++) {


            somaLab2swcId[-label_soma.get(soma_idx)-1] = swc_id;
            logWriter.println(String.format("%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f -1", swc_id, 1, soma_xyr.get(soma_idx)[0],soma_xyr.get(soma_idx)[1], 0f, soma_xyr.get(soma_idx)[2]));
            swc_id++;

            for (int i = 0; i < label_end.size(); i++) { // loops all TRACES here!

                if (label_end.get(i)==label_soma.get(soma_idx)) {

                    int mother_idx = somaLab2swcId[-label_soma.get(soma_idx)-1];

                    for (int j = tracer2d_est_xy.get(i).length-1; j >=0; j--) { // loop the trace end-beg

                        logWriter.println(String.format("%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f %-4d", swc_id, 2, tracer2d_est_xy.get(i)[j][0],tracer2d_est_xy.get(i)[j][1], 0f,1f, mother_idx));
                        mother_idx = swc_id;
                        swc_id++;

                    }

                    tracer2d_addedswc.set(i, true);

                    logWriter.println(String.format("%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f %-4d", swc_id, 2, cpregion_xy.get(label_beg.get(i)-1)[0],cpregion_xy.get(label_beg.get(i)-1)[1], 0f,1f, mother_idx));



                    // add the last point and add it to the swc table & frontier
                    regLab2swcId[label_beg.get(i)-1] = swc_id;
                    frontier.add(label_beg.get(i));




                    swc_id++;




                }
            }

            examined.add(label_soma.get(soma_idx));

            ///// structure thats built upon soma

            while (frontier.size()>0) {

                // pick an element from it
                int frontier_label = frontier.get(0);
                frontier.remove(0);

                for (int i = 0; i < label_beg.size(); i++) { // loop the rest of the TRACES to see if there are any that begin with this one

                    if (
                            label_beg.get(i)==frontier_label &&
                            label_end.get(i)!=Integer.MAX_VALUE &&
                            label_end.get(i) != label_beg.get(i) &&
                            !tracer2d_addedswc.get(i))
                    { // begin is in the frontier list

                        int mother_idx = regLab2swcId[label_beg.get(i)-1];

                        for (int j = 0; j <tracer2d_est_xy.get(i).length; j++) { // now loop the trace beg-end

                            logWriter.println(String.format("%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f %-4d", swc_id, 2, tracer2d_est_xy.get(i)[j][0],tracer2d_est_xy.get(i)[j][1], 0f,1f, mother_idx));
                            mother_idx = swc_id;
                            swc_id++;


                        }

                        tracer2d_addedswc.set(i, true);

                        logWriter.println(String.format("%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f %-4d", swc_id, 2, cpregion_xy.get(label_end.get(i)-1)[0],cpregion_xy.get(label_end.get(i)-1)[1], 0f,1f, mother_idx));

                        int endlab = label_end.get(i);
                        boolean isnew = true;

                        for (int j = 0; j < frontier.size(); j++) {
                            if (frontier.get(j)==endlab) {
                                isnew = false;
                                break;
                            }
                        }

                        for (int j = 0; j < examined.size(); j++) {
                            if (examined.get(j)==endlab) {
                                isnew = false;
                                break;
                            }
                        }

                        if (isnew) {
                            regLab2swcId[label_end.get(i)-1] = swc_id;
                            frontier.add(label_end.get(i));
                        }

                        swc_id++;

                    }
                    else if (
                            label_end.get(i)==frontier_label &&
                            label_end.get(i) != label_beg.get(i) &&
                            !tracer2d_addedswc.get(i))
                    { // end is in the frontier list, trace backwards


                        int mother_idx = regLab2swcId[label_end.get(i)-1];
                        for (int j = tracer2d_est_xy.get(i).length-1; j >=0; j--) { // loop the trace end-beg


                            logWriter.println(String.format("%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f %-4d", swc_id, 2, tracer2d_est_xy.get(i)[j][0],tracer2d_est_xy.get(i)[j][1], 0f,1f, mother_idx));
                            mother_idx = swc_id;
                            swc_id++;

                        }

                        tracer2d_addedswc.set(i, true);

                        logWriter.println(String.format("%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f %-4d", swc_id, 2, cpregion_xy.get(label_beg.get(i)-1)[0],cpregion_xy.get(label_beg.get(i)-1)[1], 0f,1f, mother_idx));

                        int endlab = label_beg.get(i);
                        boolean isnew = true;

                        for (int j = 0; j < frontier.size(); j++) {
                            if (frontier.get(j)==endlab) {
                                isnew = false;
                                break;
                            }
                        }

                        for (int j = 0; j < examined.size(); j++) {
                            if (examined.get(j)==endlab) {
                                isnew = false;
                                break;
                            }
                        }

                        if (isnew) {

                            regLab2swcId[label_beg.get(i)-1] = swc_id;
                            frontier.add(label_beg.get(i));

                        }

                        swc_id++;

                    }

                } // done checking all traces

                // append it to examined
                examined.add(frontier_label); // place it to examined once it is executed

            } // done with particular soma


        } // done with looping somas

        /*
            count legal traces that were not added to swc
         */
        int get_remaining=Integer.MAX_VALUE;
        for (int i = 0; i < label_end.size(); i++)
            if (label_end.get(i) != Integer.MAX_VALUE && !tracer2d_addedswc.get(i)) {
                get_remaining=label_beg.get(i);
                break;
            }

        while (get_remaining != Integer.MAX_VALUE) {

            // take the first one and expand the tree in the same way as it was expanded form soma

            frontier.add(get_remaining);

            while (frontier.size()>0) {

                // pick an element from it
                int frontier_label = frontier.get(0);
                frontier.remove(0);

                for (int i = 0; i < label_beg.size(); i++) { // loop the rest of the TRACES to see if there are any that begin with this one

                    if (
                            label_beg.get(i)==frontier_label &&
                                    label_end.get(i)!=Integer.MAX_VALUE &&
                                    label_end.get(i) != label_beg.get(i) &&
                                    !tracer2d_addedswc.get(i))
                    { // begin is in the frontier list

                        int mother_idx = regLab2swcId[label_beg.get(i)-1];

                        for (int j = 0; j <tracer2d_est_xy.get(i).length; j++) { // now loop the trace beg-end

                            logWriter.println(String.format("%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f %-4d", swc_id, 2, tracer2d_est_xy.get(i)[j][0],tracer2d_est_xy.get(i)[j][1], 0f,1f, mother_idx));
                            mother_idx = swc_id;
                            swc_id++;


                        }

                        tracer2d_addedswc.set(i, true);

                        logWriter.println(String.format("%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f %-4d", swc_id, 2, cpregion_xy.get(label_end.get(i)-1)[0],cpregion_xy.get(label_end.get(i)-1)[1], 0f,1f, mother_idx));

                        // see what was the end in this case - if it was examined before or not
                        // check if it belongs fo frontier or finished lists and if so it is eitherdouble-tracing or tracing some branch that was skipped but should have been done

                        int endlab = label_end.get(i);
                        boolean isnew = true;

                        for (int j = 0; j < frontier.size(); j++) {
                            if (frontier.get(j)==endlab) {
                                isnew = false;
                                break;
                            }
                        }

                        for (int j = 0; j < examined.size(); j++) {
                            if (examined.get(j)==endlab) {
                                isnew = false;
                                break;
                            }
                        }

                        if (isnew) {
                            regLab2swcId[label_end.get(i)-1] = swc_id;
                            frontier.add(label_end.get(i));
                        }

                        swc_id++;

                    }
                    else if (
                            label_end.get(i)==frontier_label &&
                                    label_end.get(i) != label_beg.get(i) &&
                                    !tracer2d_addedswc.get(i))
                    { // end is in the frontier list, trace backwards


                        int mother_idx = regLab2swcId[label_end.get(i)-1];
                        for (int j = tracer2d_est_xy.get(i).length-1; j >=0; j--) { // loop the trace end-beg


                            logWriter.println(String.format("%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f %-4d", swc_id, 2, tracer2d_est_xy.get(i)[j][0],tracer2d_est_xy.get(i)[j][1], 0f,1f, mother_idx));
                            mother_idx = swc_id;
                            swc_id++;

                        }

                        tracer2d_addedswc.set(i, true);

                        logWriter.println(String.format("%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f %-4d", swc_id, 2, cpregion_xy.get(label_beg.get(i)-1)[0],cpregion_xy.get(label_beg.get(i)-1)[1], 0f,1f, mother_idx));

                        int endlab = label_beg.get(i);
                        boolean isnew = true;

                        for (int j = 0; j < frontier.size(); j++) {
                            if (frontier.get(j)==endlab) {
                                isnew = false;
                                break;
                            }
                        }

                        for (int j = 0; j < examined.size(); j++) {
                            if (examined.get(j)==endlab) {
                                isnew = false;
                                break;
                            }
                        }

                        if (isnew) {

                            regLab2swcId[label_beg.get(i)-1] = swc_id;
                            frontier.add(label_beg.get(i));

                        }

                        swc_id++;

                    }

                } // done checking all traces

                examined.add(frontier_label);

            }

            // count remaining once again
            get_remaining=Integer.MAX_VALUE;
            for (int i = 0; i < label_end.size(); i++)
                if (label_end.get(i) != Integer.MAX_VALUE && !tracer2d_addedswc.get(i)) {
                    get_remaining=label_beg.get(i);
                    break;
                }


        }

        logWriter.close();
        System.out.println("done.\nExported .SWC: " + det_path);

    }

}
