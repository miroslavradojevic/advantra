package tracing2d;

import aux.Interpolator;
import aux.Stat;
import detection2d.Masker2D;
import detection2d.ProfileSpread2D;
import detection2d.Profiler2D;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.*;
import ij.plugin.PlugIn;
import ij.util.Tools;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by miroslav on 28-8-14.
 *
 */
public class BayesianTracking implements PlugIn, MouseListener {

    ImagePlus       curr_img;
    ImageCanvas     curr_can;

//    float[][] 	    inimg_xy;               			// store input image as an array
    float[][]       likelihood_xy;



    String          midresults_dir = "";
    String 		    image_dir;
    String		    image_name;

    Sphere2D        sph2d;
    Overlay         curr_ov;
    int             count_click;
    int           curr_center_x;
    int           curr_center_y;

    int             CPU_NR;

    // parameters
    float           D                   = 5f;
    float           prior_sigma_deg     = 35f;
    int             Nt = 100;
    boolean         save_midresults     = true;
    int             MAX_ITER = 500;
    boolean         ADD_PARTICLES = false;
    int             MARGIN = 20;

    // bayesian tracking iteration components
    ArrayList<float[][]>  Xt_xy = new ArrayList<float[][]>(); // this one is recursively updated, starts with 1 but keeps N elements distribution
    ArrayList<float[]>    wt_xy    = new ArrayList<float[]>(); // weights of the states that describe tube configuration
    ArrayList<float[]>  est_xy = new ArrayList<float[]>();

    private static void init(int _start_x, int _start_y,
                             ArrayList<float[][]> _xt_xy,
                             ArrayList<float[]> _wt_xy,
                             ArrayList<float[]>  _est_xy
    )
    {
        _xt_xy.clear();
        float[][] start_xy = new float[1][2];
        start_xy[0][0] = _start_x;
        start_xy[0][1] = _start_y;
        _xt_xy.add(start_xy);

        _wt_xy.clear();
        float[] start_w = new float[1];
        start_w[0] = 1f;
        _wt_xy.add(start_w);

        _est_xy.clear();
        float[] start_est = new float[2];
        start_est[0] = _start_x;
        start_est[1] = _start_y;
        _est_xy.add(start_est);

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


    public void run(String s) {
        System.out.println("started...");
        sph2d = new Sphere2D(D, 1.0f);
        IJ.open();
        curr_img = IJ.getImage();

        String scales_list = "3,5,7,9";

        String[] scales1 = scales_list.split(","); if (scales1.length==0) return;
        float[] scales = new float[scales1.length];
        for (int i=0; i<scales1.length; i++) scales[i] = Float.valueOf(scales1[i]);

        ImagePlus imout = Calc.neuriteness(curr_img, scales);
//        imout.show();

        likelihood_xy = new float[imout.getWidth()][imout.getHeight()]; 	// x~column, y~row

        float[] read = (float[]) imout.getProcessor().getPixels();
        for (int idx=0; idx<read.length; idx++) {
            likelihood_xy[idx%imout.getWidth()][idx/imout.getWidth()] = read[idx];
        }
        System.out.println("ready!");

        count_click = 0; // reset the counts

        curr_center_x = -1;//Float.NaN;
        curr_center_y = -1;//Float.NaN;

        curr_ov = new Overlay();
        curr_can = curr_img.getCanvas();
        curr_can.addMouseListener(this);
        curr_img.setOverlay(curr_ov);

        IJ.setTool("hand");

//
//        /********************************************************************/
//        CPU_NR = Runtime.getRuntime().availableProcessors() + 1;
//        image_dir = curr_img.getOriginalFileInfo().directory;
//        image_name = curr_img.getShortTitle();
//        midresults_dir = image_dir+image_name + "_midresults" + File.separator;
//        saveMidresults(true);
//        /********************************************************************/
//        System.out.print("loading image...");
//        inimg_xy = new float[curr_img.getWidth()][curr_img.getHeight()]; 	// x~column, y~row
//        if (curr_img.getType()== ImagePlus.GRAY8) {
//            byte[] read = (byte[]) curr_img.getProcessor().getPixels();
//            for (int idx=0; idx<read.length; idx++) {
//                inimg_xy[idx%curr_img.getWidth()][idx/curr_img.getWidth()] = (float) (read[idx] & 0xff);
//            }
//        }
//        else if (curr_img.getType()==ImagePlus.GRAY32) {
//            float[] read = (float[]) curr_img.getProcessor().getPixels();
//            for (int idx=0; idx<read.length; idx++) {
//                inimg_xy[idx%curr_img.getWidth()][idx/curr_img.getWidth()] = read[idx];
//            }
//        }
//        else {
//            IJ.log("image type not recognized");
//            return;
//        }
//        System.out.println("done.");
//        /********************************************************************/
//        System.out.print("Masker2D...");
//        float new_masker_radius = 1.5f*sph2d.getOuterRadius();   	// important that it is outer radius of the sphere
//        float new_masker_percentile = 20;                   		// used to have these two as argument but not necessary
//        Masker2D.loadTemplate(
//                inimg_xy,
//                (int) Math.ceil(new_masker_radius),
//                new_masker_radius,
//                new_masker_percentile); //image, margin, check, percentile
//        int totalLocs = inimg_xy.length * inimg_xy[0].length;
//        Masker2D ms_jobs[] = new Masker2D[CPU_NR];
//        for (int i = 0; i < ms_jobs.length; i++) {
//            ms_jobs[i] = new Masker2D(i*totalLocs/CPU_NR,  (i+1)*totalLocs/CPU_NR);
//            ms_jobs[i].start();
//        }
//        for (int i = 0; i < ms_jobs.length; i++) {
//            try {
//                ms_jobs[i].join();
//            } catch (InterruptedException e) {
//                e.printStackTrace();
//            }
//        }
//        Masker2D.defineThreshold();
//        Masker2D.formRemainingOutputs();
//        if (save_midresults) {
//            ImagePlus mask = new ImagePlus("mask", Masker2D.getMask());
//            IJ.saveAs(mask, "Tiff", midresults_dir+"mask_"+D+".tif");
//        }
//        /********************************************************************/
//        System.out.print("Profiler2D...");
//        Profiler2D.loadTemplate(
//                sph2d,
//                Masker2D.i2xy,
//                Masker2D.xy2i,
//                inimg_xy);
//        int totalProfileComponents = sph2d.getProfileLength();
//        Profiler2D pf_jobs[] = new Profiler2D[CPU_NR];
//        for (int i = 0; i < pf_jobs.length; i++) {
//            pf_jobs[i] = new Profiler2D(i*totalProfileComponents/CPU_NR,  (i+1)*totalProfileComponents/CPU_NR);
//            pf_jobs[i].start();
//        }
//        for (int i = 0; i < pf_jobs.length; i++) {
//            try {
//                pf_jobs[i].join();
//            } catch (InterruptedException e) {
//                e.printStackTrace();
//            }
//        }
//        /********************************************************************/
//        System.out.print("ProfileSpread2D...");
//        ProfileSpread2D.loadTemplate(
//                Masker2D.i2xy,
//                Masker2D.xy2i,
//                Profiler2D.prof2,
//                inimg_xy.length,
//                inimg_xy[0].length);
//        int totalProfileLocations = Profiler2D.prof2.length;
//        ProfileSpread2D pv_jobs[] = new ProfileSpread2D[CPU_NR];
//        for (int i = 0; i < pv_jobs.length; i++) {
//            pv_jobs[i] = new ProfileSpread2D(i*totalProfileLocations/CPU_NR,  (i+1)*totalProfileLocations/CPU_NR);
//            pv_jobs[i].start();
//        }
//        for (int i = 0; i < pv_jobs.length; i++) {
//            try {
//                pv_jobs[i].join();
//            } catch (InterruptedException e) {
//                e.printStackTrace();
//            }
//        }
//        ProfileSpread2D.threshold();
//        if (save_midresults) {
//            ImagePlus testip = new ImagePlus("", ProfileSpread2D.getMask());
//            IJ.saveAs(testip, "Tiff", midresults_dir+"mask_profile_"+D+".tif");
//        }
//        System.out.print(" " + IJ.d2s((ProfileSpread2D.getNrCritpointCandidates() * 100f) / (inimg_xy.length * inimg_xy[0].length), 2) + " % candidates...");
//        /********************************************************************/
//
//        System.out.println("\ndone.");

    }

//    public void saveMidresults(boolean flag)
//    {
//        // will put the flag on and create the output folder for midresults
//        save_midresults = flag;
//
//        if (save_midresults) {  // if set to true create the folder as well
//            File f1 = new File(midresults_dir);
//            if (!f1.exists()) {
//                f1.mkdirs();
//            }
//            else { // clean it, reset with each exec
//                File[] files = f1.listFiles();
//                for (File file : files)
//                {
//                    if (!file.delete())
//                        System.out.println("Failed to delete " + file);
//                }
//            }
//        }
//
//    }

    public void mouseClicked(MouseEvent e)
    {

        int x = e.getX();
        int y = e.getY();
        int offscreenX = curr_can.offScreenX(x);
        int offscreenY = curr_can.offScreenY(y);

//        int idx = Masker2D.xy2i[offscreenX][offscreenY];
//        if (idx == -1 ) {
//            System.out.println("\nit's backgound! click again...");
//            return;
//        }

        count_click++;

        // add center marker
        float r = 3.0f;
        OvalRoi central = new OvalRoi(offscreenX-r+.5, offscreenY-r+.5f, 2*r, 2*r);
        central.setStrokeColor(Color.WHITE);
        central.setFillColor(Color.WHITE);
        curr_ov.add(central);

        if (count_click==1) {

            curr_center_x = offscreenX; // only one starting state (expect more in sequential tracing)
            curr_center_y = offscreenY;

//            /*
//                plot prediction ring
//              */
//            float[] transition_x = new float[sph2d.N];
//            float[] transition_y = new float[sph2d.N];
//
//            for (int i = 0; i < sph2d.N; i++) {
//                transition_x[i] = offscreenX + sph2d.locsXY.get(i)[0]+.5f;
//                transition_y[i] = offscreenY + sph2d.locsXY.get(i)[1]+.5f;
//            }
//
//            PolygonRoi transition_ring = new PolygonRoi(transition_x, transition_y, PolygonRoi.NORMAL);
//            transition_ring.setStrokeWidth(.5);
//            curr_ov.add(transition_ring);

        }
        else if (count_click==2) { /// bayesian tracking iteration


            float[] vxy = new float[]{offscreenX-curr_center_x, offscreenY-curr_center_y};

            init(curr_center_x, curr_center_y, Xt_xy, wt_xy, est_xy);

            int iter = 1;
            while (iter<=MAX_ITER) {

                Overlay track_iter = new Overlay();

                if (iter==1) {

                    track_iter = bayesian_iteration(Nt, sph2d, prior_sigma_deg, likelihood_xy,
                            Xt_xy, wt_xy, est_xy,
                            vxy); // priors with manually set direction (used at the start)

                }
                else {

                    track_iter = bayesian_iteration(Nt, sph2d, prior_sigma_deg, likelihood_xy,
                            Xt_xy, wt_xy, est_xy
                            ); // priors with manually set direction (used at the start)

                }

                // check if it is in/out by some margin
                float last_x = est_xy.get(est_xy.size()-1)[0];
                float last_y = est_xy.get(est_xy.size()-1)[1];

                float x1 = last_x - 0;
                float x2 = curr_img.getWidth() - last_x;

                float y1 = last_y - 0;
                float y2 = curr_img.getHeight() - last_y;

                if (x1<MARGIN || x2<MARGIN || y1<MARGIN || y2<MARGIN) {
                    System.out.println("stop, reached the margin.");
                    break;
                }

                if (ADD_PARTICLES) {
                    for (int i = 0; i < track_iter.size(); i++) curr_ov.add(track_iter.get(i));
                }

                System.out.println("iter " + iter + " ->   " + est_xy.get(est_xy.size()-1)[0] + " , " + est_xy.get(est_xy.size()-1)[1]);

                iter++;

            }


            // input: curr_center_x,y,  offscreenX,Y -> vxy,  sph2d
            // output:  transition_ring_x,y, priors[], likelihoods[], posteriors[]
//            // priors
//            float[] priors = new float[sph2d.theta.size()]; // get priors based on the selected direction
//            sph2d.getPriors(vxy, prior_sigma_deg, priors);
//            double[] min_max_priors = Tools.getMinMax(priors);
//            // likelihoods
//            float[] likelihoods = new float[sph2d.N];//[Profiler2D.prof2[0].length];
//            idx = Masker2D.xy2i[curr_center_x][curr_center_y];
//            float likelihoods_min = Float.POSITIVE_INFINITY;
//            float likelihoods_max = Float.NEGATIVE_INFINITY;
//            for (int i=0; i<likelihoods.length; i++) {
//                likelihoods[i] = 1; //((Profiler2D.prof2[idx][i] & 0xffff) / 65535f) * 255f; // retrieve the profile
//                if (likelihoods[i]<likelihoods_min)
//                    likelihoods_min = likelihoods[i];
//                if (likelihoods[i]>likelihoods_max)
//                    likelihoods_max = likelihoods[i];
//            }
//            // posteriors
//            float[] posteriors = new float[priors.length];
//            float posteriors_sum = 0;
//            for (int i = 0; i < posteriors.length; i++) {
//                posteriors[i] = priors[i] * likelihoods[i] * 1;
//                posteriors_sum += posteriors[i];
//            }
//            for (int i = 0; i < posteriors.length; i++) {
//                posteriors[i] /= posteriors_sum;
//            }
//            double[] min_max_posteriors = Tools.getMinMax(posteriors);


            /*
                draw direction line
             */
            Line ln = new Line(offscreenX+.5f, offscreenY+.5f, curr_center_x+.5f, curr_center_y+.5f);
            ln.setStrokeColor(Color.YELLOW);
            ln.setFillColor(Color.YELLOW);
            curr_ov.add(ln);

            /*
                draw estimated trace
             */
            float[] est_x = new float[est_xy.size()];
            float[] est_y = new float[est_xy.size()];

            for (int i = 0; i < est_xy.size(); i++) {

                est_x[i] = est_xy.get(i)[0];
                est_y[i] = est_xy.get(i)[1];

//                System.out.println("drawing: " + est_x[i] +" , "+ est_y[i]);

            }
            PolygonRoi pr = new PolygonRoi(est_x, est_y, PolygonRoi.FREELINE);
            pr.setFillColor(Color.YELLOW);
            pr.setStrokeWidth(2);
            curr_ov.add(pr);

//            /*
//                draw prior vales ring (blue)
//              */
//            float[] priors_x = new float[sph2d.theta.size()];
//            float[] priors_y = new float[sph2d.theta.size()];
//
//            for (int i = 0; i < sph2d.locsXY.size(); i++) {
//                float sc = (float) (priors[i] * (D/min_max_priors[1]));
//                priors_x[i] = curr_center_x + sph2d.locsXY.get(i)[0] + sc * sph2d.vxy.get(i)[0];
//                priors_y[i] = curr_center_y + sph2d.locsXY.get(i)[1] + sc * sph2d.vxy.get(i)[1];
//            }
//
//            PolygonRoi prior_ring = new PolygonRoi(priors_x, priors_y, PolygonRoi.NORMAL);
//            prior_ring.setStrokeColor(Color.BLUE);
//            prior_ring.setStrokeWidth(.5);
//            curr_ov.add(prior_ring);
//
//            /*
//                add likelihood plot
//             */
//            float[] likelihoods_x = new float[sph2d.theta.size()];
//            float[] likelihoods_y = new float[sph2d.theta.size()];
//
//            for (int i = 0; i < sph2d.theta.size(); i++) {
//                float sc = (float) (likelihoods[i] * D);
//                likelihoods_x[i] = curr_center_x + sph2d.locsXY.get(i)[0] + sc * sph2d.vxy.get(i)[0];
//                likelihoods_y[i] = curr_center_y + sph2d.locsXY.get(i)[1] + sc * sph2d.vxy.get(i)[1];
//            }
//
//            PolygonRoi likelihood_ring = new PolygonRoi(likelihoods_x, likelihoods_y, PolygonRoi.NORMAL);
//            likelihood_ring.setStrokeColor(Color.GREEN);
//            likelihood_ring.setStrokeWidth(.5);
//
//            curr_ov.add(likelihood_ring);
//
//
//
//            /*
//                add posterior values ring (red) - bayesian - only first iteration
//              */
//            float[] posteriors_x = new float[sph2d.theta.size()];
//            float[] posteriors_y = new float[sph2d.theta.size()];
//
//            for (int i = 0; i < sph2d.theta.size(); i++) {
//                float sc = (float) (posteriors[i] * (D/min_max_posteriors[1]));
//                posteriors_x[i] = curr_center_x + sph2d.locsXY.get(i)[0] + sc * sph2d.vxy.get(i)[0];
//                posteriors_y[i] = curr_center_y + sph2d.locsXY.get(i)[1] + sc * sph2d.vxy.get(i)[1];
//            }
//
//            PolygonRoi posterior_ring = new PolygonRoi(posteriors_x, posteriors_y, PolygonRoi.NORMAL);
//            posterior_ring.setStrokeColor(Color.RED);
//            posterior_ring.setStrokeWidth(.5);
//
//            curr_ov.add(posterior_ring);



        }
        else if (count_click==3) {

            count_click=0;
            curr_ov.clear();

        }

        curr_img.draw();


    }

    public void mousePressed(MouseEvent e) {

    }

    public void mouseReleased(MouseEvent e) {

    }

    public void mouseEntered(MouseEvent e) {

    }

    public void mouseExited(MouseEvent e) {

    }
}
