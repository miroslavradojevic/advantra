package detection2d;

import aux.Interpolator;
import aux.Tools;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.*;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;
import tracing2d.BayesianTracerMulti;
import tracing2d.Extractor;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.util.ArrayList;

/**
 * Created by miroslav on 11-11-14.
 */
public class NS_Delineation implements PlugIn, MouseListener, MouseMotionListener {

    // input
    ImageCanvas cnv;
    String      image_path;
    float[][]   likelihood_xy;

    // bayesian filtering
    int R               = 4;     // will cover radiuses R * BayesianTracerMulti.sstep
    int Ns              = 50;    // nr samples
    float sigma_deg     = 100;   // degrees standard deviation
    int Nstreams        = 4;     //
    BayesianTracerMulti.Expansion expan = BayesianTracerMulti.Expansion.INNER;

    // foreground extraction
    int percentile = 90;    // how many to keep as the foreground/background, neighbourhood will be defined with radius and
    int nbhood = (int) Math.ceil(2*R*BayesianTracerMulti.sstep[BayesianTracerMulti.sstep.length-1]);
    int CPU_NR = Runtime.getRuntime().availableProcessors() + 1;

    // output
//    ArrayList<ArrayList<float[]>> delin = new ArrayList<ArrayList<float[]>>(); // will store the spatial delineation

//    float[][][] xc = new float[Ni][Ns][2];    // mean-shift convergence of the bayesian filtered states (x, y)
//    float[][][] xcc = new float[Ni][Nc][2];    // after clustering the converged values - take up to Nc clusters
//    float[][]   dsts = new float[Ns][Ns];       // inter distances - used as a clustering criteria
//    int[]       cllab = new int[Ns];                // clustering labels at each iteration
//    byte[][]     tag = new byte[Ni][Nc];       // fmm tags (0, 1, 2)
//    byte[][]     par = new byte[Ni][Nc];       // fmm pointers (-2, -1, 0, 1, 2, 3)
//    float[][]    cst = new float[Ni][Nc];      // fmm scores
    // aux
//    float[]     xn = new float[2];            // new one for the mean-shift (auxilliary variable)
//    ArrayList<int[]>    heap_vals = new ArrayList<int[]>();      // rank, Ni idx, N curr idx, N prev idx
//    ArrayList<Float>    heap_scrs = new ArrayList<Float>();          // scores

    Overlay ov = new Overlay();

    public void run(String s)
    {

        // parent map is in bytes - can cover only 256 values (0-255) - that is the limit in number of the states
        if (Ns >= 255) { IJ.log("Parent map (BYTE vals) supports up to 255 states. Current # " + Ns + ". Exiting..."); return; }

         // later will be used to trace, here is initialised to set the scale and the templates to match

        // load the image through the menu
        String in_folder = Prefs.get("id.folder", System.getProperty("user.home"));
        OpenDialog.setDefaultDirectory(in_folder);
        OpenDialog dc = new OpenDialog("Select file");
        in_folder = dc.getDirectory();
        image_path = dc.getPath();
        if (image_path==null) return;
        Prefs.set("id.folder", in_folder);

        ImagePlus img = new ImagePlus(image_path);
        if(img==null) { IJ.log("Input image was null."); return; }

//        stream_map = new byte[img.getWidth()][img.getHeight()];

        likelihood_xy = new float[img.getWidth()][img.getHeight()]; 	// x~column, y~row
        if (img.getType()== ImagePlus.GRAY8) {

            System.out.println("loading image...");
            float mn = Float.POSITIVE_INFINITY;
            float mx = Float.NEGATIVE_INFINITY;

            byte[] read = (byte[]) img.getProcessor().getPixels();
            for (int idx=0; idx<read.length; idx++) {
                float val = (float) (read[idx] & 0xff) / 255f; // normalize 0-1
                likelihood_xy[idx%img.getWidth()][idx/img.getWidth()] = val;
                if (val<mn) mn = val;
                if (val>mx) mx = val;
            }

            System.out.println("done. range=" + mn + " -- " + mx);

        }
//        else if (img.getType()==ImagePlus.GRAY32) {
//            float[] read = (float[]) img.getProcessor().getPixels();
//            for (int idx=0; idx<read.length; idx++) {
//                likelihood_xy[idx%img.getWidth()][idx/img.getWidth()] = read[idx];
//            }
//        }
        else {
            IJ.log("image type not recognized");
            return;
        }






        /*
            extract the foreground -> Masker2D.xy2i, Masker2D.i2xy
         */
        System.out.print("extracting foreground...");
        Masker2D.loadTemplate(
                likelihood_xy, nbhood, nbhood, percentile); //image, margin, check, percentile
        int totalLocs = likelihood_xy.length * likelihood_xy[0].length;
        Masker2D ms_jobs[] = new Masker2D[CPU_NR];
        for (int i = 0; i < ms_jobs.length; i++) {
            ms_jobs[i] = new Masker2D(i*totalLocs/CPU_NR,  (i+1)*totalLocs/CPU_NR);
            ms_jobs[i].start();
        }
        for (int i = 0; i < ms_jobs.length; i++) {
            try {
                ms_jobs[i].join();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }
        Masker2D.defineThreshold();
        Masker2D.formRemainingOutputs();
        System.out.println(" done. ");

//        ImagePlus mask = new ImagePlus("mask", Masker2D.getMask());
//        mask.show();//      IJ.saveAs(mask, "Tiff", midresults_dir+"mask_"+D[didx]+".tif");










        // show it and get the canvas
        img.show();
        cnv = img.getCanvas();
        cnv.addMouseListener(this);
        cnv.addMouseMotionListener(this);
        IJ.setTool("hand");

        System.out.println("ready...");

//        vizTag.show();
//        vizTag.getCanvas().zoom100Percent();

    }

    public static Overlay viz_xc(float[][][] xc)
    {

//        float rad = .5f;

        Overlay ov = new Overlay();
        for (int i = 0; i < xc.length; i++) { // iterations

//            Stat.min_max_normalize(wt[i]);

            for (int j = 0; j <xc[i].length; j++) {

                PointRoi p = new PointRoi(xc[i][j][0]+.5, xc[i][j][1]+.5);
//                p.setFillColor(new Color(1f,1f,1f,wt[i][j]));
//                Line l = new Line(xt[i][j][0]+.5, xt[i][j][1]+.5, xt[i][j][0]+xt[i][j][2]+.5, xt[i][j][1]+xt[i][j][3]+.5);

                ov.add(p);
//                ov.add(l);

            }
        }

        return ov;
    }

    public static Overlay viz_xcc(float[][][] xcc) // Ni*4*2
    {

//        float rad = .5f;

        Overlay ov = new Overlay();
        for (int i = 0; i < xcc.length; i++) { // iterations

//            Stat.min_max_normalize(wt[i]);

            for (int j = 0; j <xcc[i].length; j++) {

                if (!Float.isNaN(xcc[i][j][0])) {
                    PointRoi p = new PointRoi(xcc[i][j][0]+.5, xcc[i][j][1]+.5);
                    p.setFillColor(Color.GREEN);
                    p.setStrokeColor(Color.GREEN);
                    p.setStrokeWidth(2);
                    ov.add(p);
                }

//                p.setFillColor(new Color(1f,1f,1f,wt[i][j]));
//                Line l = new Line(xt[i][j][0]+.5, xt[i][j][1]+.5, xt[i][j][0]+xt[i][j][2]+.5, xt[i][j][1]+xt[i][j][3]+.5);
//                ov.add(l);

            }
        }

        return ov;
    }

    public void mouseClicked(MouseEvent e)
    {
        int clickX = cnv.offScreenX(e.getX());
        int clickY = cnv.offScreenY(e.getY());

        System.out.print("(x,y) -> "+clickX+" , " +clickY+" -> ");

        long t1 = System.currentTimeMillis();

        Extractor.loadTemplate(2, R, Ns, sigma_deg, Nstreams);
        boolean show_tracks = false;
        Overlay oo = Extractor.extractAt(clickX, clickY, likelihood_xy, show_tracks, expan);

        // bayesian filtering - take cate that it is initilized before

//        extract(xt, wt, pt, Nstreams, stream_map, delin); // states, weights, parent_index, nr_streams, overlap_margin, output_delineation

            // mean-shift for estimation (local maxima detection)
//            for (int i = 0; i < xc.length; i++) { // initialize xc with corresponding elements from xt
//                for (int s = 0; s < xc[0].length; s++) {
//                    xc[i][s][0] = xt[i + 1][s][0];
//                    xc[i][s][1] = xt[i + 1][s][1];
//                }
//            }
//            for (int i = 0; i < Ni; i++) { // do convergence at each iteration
//                meanshift_xy_convg(xt[i + 1], wt[i + 1], xc[i], xn, max_iter, epsilon, r2); // xc will hold the converged values
//                // clustering xc->xcc
//                dists2(xc[i], dsts); // Ns*Ns distances at each iteration
//                clustering(dsts, r2, cllab);
//                extracting(cllab, xc[i], 1, Nc, xcc[i]);
//            }
//            fmm(likelihood_xy, new float[]{clickX, clickY}, xcc, tag, par, cst, heap_vals, heap_scrs);
//            extract_streamlines(clickX, clickY, tag, par, cst, xcc, 4, delin); // xcc will be modified

        long t2 = System.currentTimeMillis();
        float tt = (t2-t1)/1000f;
        System.out.println("     "+ tt + "sec. whole image would take " + (tt*Masker2D.i2xy.length));

//        vizTag.setProcessor(getTag());
//        vizTag.updateAndDraw();
//        IJ.run(vizTag, "Enhance Contrast", "saturated=0.35");
//        new ImagePlus("", getTag()).show();

//        Overlay ov_xc   = viz_xc(xc);
//        Overlay ov_xcc  = viz_xcc(xcc);

        // comment here if you need different visualisations
//        ov.clear();
//        for (int i = 0; i <ov_xc.size(); i++) ov.add(ov_xc.get(i));
//        for (int i = 0; i <ov_xcc.size(); i++) ov.add(ov_xcc.get(i));

        cnv.setOverlay(oo);

    }

    public void mouseMoved(MouseEvent e)
    {
        mouseClicked(e);
    }

    private static void fmm(
                            float[][] inimg_xy,         // input data
                            float[] start_xy,           // input data
                            float[][][] xcc,            // input to fast march through Ni*4*2
                            byte[][] tag,               // output Ni*4 byte (0,1,2)
                            byte[][] par,               // output Ni*4 byte (-2,-1,0,1,2,3)
                            float[][] cst,              // output Ni*4 float - cost
                            ArrayList<int[]> heap_vals, // auxilliary heap
                            ArrayList<Float> heap_scrs  // auxilliary heap
    )
    {

        // some paramters used for calculating the costs
        float LAMBDA = 1;
        float I = 0.4f;
        float ALFA = (float) (Math.PI/2);
        float P = 3; // multiple of radius should be parametrized

        // reset the matrices with outputs
        for (int i = 0; i < tag.length; i++)
            for (int j = 0; j < tag[i].length; j++) {
                tag[i][j] = (byte) 0;  // 0 (none), 1 (narrow band), 2 (frozen)
                par[i][j] = (byte) -2; // -2 (none), -1 (root xy), 0,1,2,3 (index of the element in previous iteration)
                cst[i][j] = (byte) Float.POSITIVE_INFINITY;
            }
        // heap architecture: heap_vals:<rank, Ni idx, N curr idx, N parent idx> heap_scrs:<score>
        heap_vals.clear();
        heap_scrs.clear();

        // xcc contains the centroids Ni*4*2d centroids (set to NaN,NaN when there is no centroid)
        // centroids are iterating outwards with the nighest index being the most distant from the center
        // fmm finds the optimal path through the set of centroids
        // each next iteration is considered when searching for the neighbours
        // initialization part of the algorithm
        int startII = 0; // start iteration index fmm will start backwards from the first iteration

        // tag, par, and cst will be filled up while marching till the end
        for (int i = 0; i < xcc[startII].length; i++) { // check 4 of neighbours at iteration=0 (the first one)

            int par_idx     = -1;
            int ni_idx      = startII;
            int pos_idx     = i;

            if (!Float.isNaN(xcc[startII][i][0])) { // if there is a centroid at this iteration

                float curr_x = xcc[0][i][0];
                float curr_y = xcc[0][i][1];

                // cost calculation
                float I1 = Interpolator.interpolateAt(curr_x, curr_y, inimg_xy);
                float I0 = Interpolator.interpolateAt(start_xy[0], start_xy[1], inimg_xy);

                float CI = (float) Math.exp(LAMBDA * Math.pow( (I1-I0)/I,2));
                float CP = (float) Math.exp(LAMBDA * ((Math.pow(curr_x-start_xy[0],2)+Math.pow(curr_y-start_xy[1],2))/Math.pow(P,2)));
                float CALPHA = 1;

                float C = CI * CP * CALPHA;

                tag[0][i] = (byte) 1;           // tag it as narrow band
                par[0][i] = (byte) par_idx;     // parent is root xy
                cst[0][i] = C;                  // cost should accumulate the previous cost which is 0

                if (heap_vals.size()==0) { // just insert into the heap
                    heap_vals.add(new int[]{1, ni_idx, pos_idx, par_idx}); // rank, Ni idx, N curr idx, N parent idx
                    heap_scrs.add(C);
                }
                else{ // append it and update the rank depending on where the cost was
                    int cnt = 0;
                    for (int j = 0; j < heap_scrs.size(); j++) if (heap_scrs.get(j) > C) {cnt++; heap_vals.get(j)[0]++;} // count them & decrease ranking

                    int new_rank = heap_scrs.size() - cnt + 1; // ranking starts from 1

                    // add the new one
                    heap_scrs.add(C);
                    heap_vals.add(new int[]{new_rank, ni_idx, pos_idx, par_idx});
                }

            }

        } // end initialization

        while (heap_vals.size()>0) {

//            for (int i = 0; i <heap_vals.size(); i++) System.out.print(heap_vals.get(i)[0] + " , ");

            // get the top score (this can be optimized so that there is no loop)
            int get_lowest_score = -1;
            for (int i = 0; i < heap_vals.size(); i++) // loop through vals list
                if (heap_vals.get(i)[0] == 1)
                   get_lowest_score = i;
                else
                   heap_vals.get(i)[0]--; // decrease ranking

            if (get_lowest_score==-1) {
                System.out.println("wrong !!!! cannot find rank 1");
                System.exit(0);
            }

            // read heap_vals and heap_scrs first...
            int at_i = heap_vals.get(get_lowest_score)[1];  // iteration index
            int at_p = heap_vals.get(get_lowest_score)[2];  // curr idx within iteration
            int at_pp = heap_vals.get(get_lowest_score)[3]; // parent idx
            float score = heap_scrs.get(get_lowest_score);  // score

            // remove from heap list
            heap_scrs.remove(get_lowest_score);
            heap_vals.remove(get_lowest_score);

            // freeze it using the taken list values
            tag[at_i][at_p] = (byte) 2; // tag it as freezed
            par[at_i][at_p] = (byte) at_pp;
            cst[at_i][at_p] = score;

            // loop neighbours vn (at_i+1,:) of the taken one v: (at_i,at_p) (the following iteration) if there is such possibility
            if (at_i+1<xcc.length) { // if the neighbours exist

                for (int i = 0; i < xcc[at_i+1].length; i++){ // looping the neighbours

                    if (!Float.isNaN(xcc[at_i+1][i][0])) { // there was a cluster at the neighbour position

                        if (tag[at_i+1][i]!=(byte)2){ // it is not frozen

                            // (x1,y1), (x0,y0), (xp,yp)
                            float x1 = xcc[at_i+1][i][0]; // neighbour of the element taken from the heap
                            float y1 = xcc[at_i+1][i][1];

                            float x0 = xcc[at_i][at_p][0]; // just taken from the heap
                            float y0 = xcc[at_i][at_p][1];

                            // see what is xp, yp
                            float xp= Float.NaN, yp=Float.NaN;

//                            System.out.println("test " + at_i + "    --- " + at_p);

                            if (par[at_i][at_p]==(byte)-1) {      // if it was root point
                                xp = start_xy[0];
                                yp = start_xy[1];
                            }
                            else if (       par[at_i][at_p]>=(byte)0 &&
                                            par[at_i][at_p]<(byte)xcc[at_i].length
                                    )
                            {
                                xp = xcc[at_i-1][par[at_i][at_p]&0xff][0];
                                yp = xcc[at_i-1][par[at_i][at_p]&0xff][1];

                                // check
                                if (Float.isNaN(xp) || Float.isNaN(yp)) {
                                    System.out.println("parent was NaN!!!");
                                    System.exit(0);
                                }

                            }
                            else { // it was not a legal parent index
                                System.out.println("wrong parent index when calculating the score " + par[at_i][at_p] + "  --- " + (int)(par[at_i][at_p]&0xff));
                                System.exit(0);
                            }

                            // (v1x,y, v2x,y)
                            float v1x = x1-x0;
                            float v1y = y1-y0;
                            float v1 = (float) Math.sqrt(v1x*v1x+v1y*v1y);
                            v1x = v1x/v1;
                            v1y = v1y/v1;

                            float v2x = x0-xp;
                            float v2y = y0-yp;
                            float v2 = (float) Math.sqrt(v2x*v2x+v2y*v2y);
                            v2x = v2x/v2;
                            v2y = v2y/v2;

                            // ang
                            float ang = v1x*v2x+v1y*v2y;//(float) Math.acos(v1x*v2x+v1y*v2y);
                            ang = (ang>1)? 1 : (ang<-1)? -1 : ang; // to avoid NaNs
                            ang = (float) Math.acos(ang); // so that the angle is not NaN - takes the range [0 - pi]

                            // cost calculation
                            float I1 = Interpolator.interpolateAt(x1, y1, inimg_xy);
                            float I0 = Interpolator.interpolateAt(x0, y0, inimg_xy);

                            float CI = (float) Math.exp(LAMBDA * Math.pow( (I1-I0)/I,2));
                            float CP = (float) Math.exp(LAMBDA * ((Math.pow(x1-x0,2)+Math.pow(y1-y0,2))/Math.pow(P,2)));
                            float CALPHA = (float) Math.exp(LAMBDA * Math.pow(ang/ALFA,2)); //2 - cosang; //(float) (Math.pow(curr_x-prev_x,2)+Math.pow(curr_y-prev_y,2)); // directional cost of the neighbour

                            float C = CI * CP * CALPHA + cst[at_i][at_p];

//                            float eI = 1;// (float) Math.exp(5*Math.pow(1-Icurr/1f,2));
//                            float eIprev = (float) Math.exp(10*Math.pow(1-Iprev/1f,2));
//                            float e = eD * eI * eD1 + cst[at_i][at_p];// new cost is cummulative + ( (eIcurr+eIcurr)/2f );

                            if (tag[at_i+1][i]!=(byte)1) {// neighbour (at_i+1, i) is NOT in narrow band

                                tag[at_i+1][i] = (byte) 1; // tag it as narrow band
                                par[at_i+1][i] = (byte) at_p;
                                cst[at_i+1][i] = C;

                                // calculate the rank depending on where the cost was in the heap score list
                                int cnt = 0;
                                for (int j = 0; j < heap_scrs.size(); j++) if (heap_scrs.get(j) > C) {cnt++; heap_vals.get(j)[0]++;}
                                int new_rank = heap_scrs.size() - cnt + 1; // ranking starts from 1

                                // append it to the heap
                                heap_scrs.add(C);
                                heap_vals.add(new int[]{new_rank, at_i+1, i, at_p});

                            }
                            else {

//                                System.out.println("shoud be 1 " + (tag[at_i+1][i]&0xff));
                                // neighbour IS in narrow band already - means that it is already added to the heap
                                // find the corresponding heap location and update it

                                // [look for at_i+1, i]
                                boolean found = false;
                                int idx_found = -1;
                                for (int j = 0; j < heap_vals.size(); j++) {
                                    if (heap_vals.get(j)[1]==at_i+1 && heap_vals.get(j)[2]==i) { // element in the heap equals (at_i+1,i)
                                        found = true;
                                        idx_found = j;
                                        break;
                                    }
                                }

                                if (!found) {System.out.println("#### problem finding in the HEAP ####");System.exit(0);}

                                float found_score = heap_scrs.get(idx_found);

                                if (C<found_score) {

                                    /*
                                        heap updates
                                     */

                                    // update rankings
                                    float old_val = found_score;
                                    float new_val = C; // new_val is smaller!!

                                    int new_rank = 1;
                                    for (int j = 0; j < heap_scrs.size(); j++) {

                                        if (heap_scrs.get(j) >= new_val && heap_scrs.get(j)<old_val)
                                            heap_vals.get(j)[0]++;

                                        if (heap_scrs.get(j)<new_val) new_rank++;

                                    }
                                    heap_vals.get(idx_found)[0] = new_rank;
                                    heap_vals.get(idx_found)[3] = at_p;
                                    // replace the score with the cheaper one once the rankings are corrected
                                    heap_scrs.set(idx_found, new_val);

                                    /*
                                        table updates
                                     */
                                    // update parent table - means updating the score and the parent, labels do not need correction - value stays the same
                                    par[at_i+1][i] = (byte) at_p;
                                    cst[at_i+1][i] = C;

                                }

                            }

                        }

                    }

                }

            }

        }

    }

    private static void extract_streamlines(
            float                           x_root,
            float                           y_root,
            byte[][]                        tag, // input
            byte[][]                        par, // input
            float[][]                       cst, // input
            float[][][]                     xcc, // side output, extracted path will be excluded
            int                             max_streamlines,
            ArrayList<ArrayList<float[]>>   streamlines
    )
    {

        int Ni_margin = 2; // todo parameter

//        int nr_streams = 4; // to capture enough to cover all types of critpoints
        byte[][] trace_map = new byte[tag.length][tag[0].length]; // size of the tag todo: outside as auxilliary not to be allocated each time

        // will be filld with ids of the traces as they are checked (0, 1, 2, 3)
        for (int i = 0; i < trace_map.length; i++)
            for (int j = 0; j < trace_map[i].length; j++)
                trace_map[i][j] = 0; // 0 means the centroid at this location does not belong to any trace

        streamlines.clear(); // take the frozen ones from the last iteration and check their backtraces

        boolean[] examined = new boolean[tag[tag.length-1].length]; // all are false at the beginning

        int cnt_stream = 0;

        for (int i = 0; i < max_streamlines; i++) { // pick the one with minimum cost that does not overlap with any earlier taken one max_streamlines times

            // extract streamlines will contain prunning so that the non-overlapping ones with lowest cost are extracted gradually
            float min_score = Float.POSITIVE_INFINITY;
            int min_idx = -1;

            for (int j = 0; j < tag[tag.length-1].length; j++) {

                if (!examined[j]) {
                    if (tag[tag.length-1][j]==(byte)2) { // if it is frozen (fmm reached it)

                        // consider the corresponding cost as min
                        if (cst[tag.length-1][j]<min_score) {
                            min_score = cst[tag.length-1][j];
                            min_idx = i;
                        }

                    }
                    else {
                        // make it examined it can be only 0 (missing)
                        if (tag[tag.length-1][j]==(byte)1) {
                            System.out.println("sth was wrong");
                            System.exit(0);
                        }

                        examined[j] = true; // if it was not frozen consider it was examined, can be only NaN, so non-existing centroid, empty space

                    }
                }

            }

            // all are examined
            if (min_idx!=-1) { // min was found in the loop at one index

                cnt_stream++;

                // min was found - no need to examine it next iteration whether it was added in the end or not
                examined[min_idx] = true;

                // back trace to see whether it should be added
                // no need to search for the minimum now - take all that there is, the fact that it is frozen means that the cost is not infinity
                ArrayList<float[]> out_streamline = new ArrayList<float[]>(); // build it up gradually

                boolean is_overlapping = false;

                int found_j = min_idx;

                for (int loop_streamiline = tag.length-1; loop_streamiline >=-1 ; loop_streamiline--) { // add all the points of the stram plus the central one (loop_streamline==-1)

                    // before adding the element - check if it overlaps with some previous path
                    if (loop_streamiline==-1) {

                        out_streamline.add(new float[]{x_root, y_root});

                    }
                    else {

                        if (trace_map[loop_streamiline][found_j]!=(byte)0 && loop_streamiline>=Ni_margin) { // check if it overlaps with the earlier trace
                            is_overlapping = true;
                            break; // reached another stream
                        }

                        out_streamline.add(xcc[loop_streamiline][found_j].clone());
                        trace_map[loop_streamiline][found_j] = (byte) cnt_stream;
                        found_j = par[loop_streamiline][found_j]&0xff; // last one will give 255 (-1 in byte) but it won't be referenced

                    }

                }

                if (!is_overlapping) { // it it turned out to overlap
                    streamlines.add(out_streamline);
                }
            }
            else break;

        } // end looping the streams

//        for (int i = 0; i < tag[tag.length-1].length; i++) { // take everything there is that reached the last iteration - lowest cost first
//            if (tag[tag.length-1][i]==(byte) 2) {}
//        }


/*
        // use tag and cst to find initial thread point - most distant frozen that had the least cost
        boolean found = false;
        int found_i = -1;
        int found_j = -1;
        float min_cost;

        for (int i = tag.length-1; i >= 0; i--) { //  tag.length-1-margin_end

            min_cost = Float.POSITIVE_INFINITY; // searching the centroids at concentric circle - the one with min cost

            for (int j = 0; j < tag[i].length; j++) {
                if (tag[i][j]==(byte)2) { // if it was frozen
                    if (cst[i][j]<min_cost) { min_cost = cst[i][j]; found_i = i; found_j = j; }
                }
            }

            if (!Float.isInfinite(min_cost)) {
                found = true;
                break; // found get out as soon as it is found at one iteration
            }

        } // take the latest with the lowest cost

        if (found) { // found_i, found_j

            ArrayList<float[]> out = new ArrayList<float[]>();

            for (int loop_streamiline = found_i; loop_streamiline >=0 ; loop_streamiline--) {

                out.add(xcc[loop_streamiline][found_j].clone()); // start from found_i,j

                // cancel the one that was added - set it as NaN for future fast marching iteration
//                if (loop_streamiline>=2) {
//                    xcc[loop_streamiline][found_j][0] = Float.NaN;
//                    xcc[loop_streamiline][found_j][1] = Float.NaN;
//                }

                found_j = par[loop_streamiline][found_j]&0xff; // recursively go backwards

            }

            return out;

        }
        else {

            return null;

        }

*/

    }

    private static void	meanshift_xy_convg(
                                            float[][]   xy_inpt,        // one row, particular iteration
                                            float[]     w_inpt,
                                            float[][]   xy_conv,
                                            float[]     xy_next,        // auxiliary variable for iteration (to avoid allocation inside the loop)
                                            int         max_iter,
                                            float       epsilon,
                                            float       rad2
    )
    {

        for (int i = 0; i < xy_inpt.length; i++) { // converge each distribution element separately

            int iter = 0;
            float d;

            do {

                runone(xy_conv[i], xy_next, rad2, xy_inpt, w_inpt);

                d = (float) (Math.pow(xy_next[0]-xy_conv[i][0], 2) + Math.pow(xy_next[1]-xy_conv[i][1], 2));

                xy_conv[i][0] = xy_next[0];
                xy_conv[i][1] = xy_next[1];

                iter++;

            }
            while (iter < max_iter && d > epsilon*epsilon);

        }

    }

    private static void runone(float[] curr_xy, float[] next_xy, float rad2, float[][] template_xy, float[] template_w)
    {


//        System.out.println(Arrays.toString(curr_xy)+" :: " + rad2 + " :: ");

        float sum 	= 0;
        next_xy[0] 	= 0;
        next_xy[1] 	= 0;
//        int kk = 0;

        for (int l = 0; l < template_xy.length; l++) { // loop all the rest from the distribution
            if (Math.pow(template_xy[l][0]-curr_xy[0],2)+Math.pow(template_xy[l][1]-curr_xy[1],2)<=rad2) {
                sum += template_w[l];
                next_xy[0] += template_w[l] * template_xy[l][0]; //
                next_xy[1] += template_w[l] * template_xy[l][1]; //
//                kk++;
            }
        }

//        System.out.println("cases = " + kk + " / " + template_xy.length);

        if (sum>=Float.MIN_VALUE) {
            next_xy[0] /= sum;
            next_xy[1] /= sum;
        }
        else {
            System.out.println("this was not supposed to happen!");
            next_xy[0] = curr_xy[0];
            next_xy[1] = curr_xy[1];
        }

    }

    private static void dists2(float[][] vec_xy, float[][] cross_dists2)
    {

        for (int i = 0; i < vec_xy.length; i++) {
            for (int j = i; j < vec_xy.length; j++) {
                if (i==j) {
                    cross_dists2[i][j] = 0;
                }
                else {
                    cross_dists2[i][j] = (float) (Math.pow(vec_xy[i][0]-vec_xy[j][0],2)+Math.pow(vec_xy[i][1]-vec_xy[j][1],2));
                    cross_dists2[j][i] = cross_dists2[i][j];
                }
            }
        }

    }

    private static void clustering(float[][] dists2, float threshold_dists2, int[] labels)
    {

//        int[] labels = new int[idxs.length];
        for (int i = 0; i < dists2.length; i++) labels[i] = i;

        for (int i = 0; i < dists2.length; i++) {

            // one versus the rest
            for (int j = 0; j < dists2.length; j++) {

                //
                if (i != j) {

                    if (dists2[i][j]<=threshold_dists2) {

                        if (labels[j] != labels[i]) {

                            int currLabel = labels[j];
                            int newLabel  = labels[i];

                            labels[j] = newLabel;

                            //set all that also were currLabel to newLabel
                            for (int k = 0; k < labels.length; k++)
                                if (labels[k]==currLabel)
                                    labels[k] = newLabel;

                        }

                    }

                }

            }

        }

//        return labels;

    }

    private static void extracting(int[] labels, float[][] vec_xy, int min_count, int _Nc, float[][] out4x2)
    {

        boolean[] checked = new boolean[labels.length];

        ArrayList<float[]> out = new ArrayList<float[]>();
        ArrayList<Integer> cnt = new ArrayList<Integer>();    		// to make sure that it outputs sorted list

        for (int i = 0; i < labels.length; i++) {
            if (!checked[i]) {

                float centroid_x = vec_xy[i][0]; // idxs[i]
                float centroid_y = vec_xy[i][1]; // idxs[i]
                int count = 1;
                checked[i] = true;

                // check the rest
                for (int j = i+1; j < labels.length; j++) {
                    if (!checked[j]) {
                        if (labels[j]==labels[i]) {

                            centroid_x += vec_xy[j][0]; // idxs[j]
                            centroid_y += vec_xy[j][1]; // idxs[j]
                            count++;
                            checked[j] = true;

                        }
                    }
                }

                if (count >= min_count) { // consider as clusters all that have at least min_count elements
                    out.add(new float[]{centroid_x/count, centroid_y/count});
                    cnt.add(count);
                }

            }
        }


        // print before sorting
//		for (int ii = 0; ii < cnt.size(); ii++) {
//			IJ.log(ii + " : " + cnt.get(ii) + " points,  at " + Arrays.toString(out.get(ii)));
//		}

        // take top _Nc of all those that qualified

//        int nr_clusters_after_clustering = out.size();

        // sort by the counts of the converged elements in one cluster (take from Sorting.java) descending indexes
        int[] desc_idx = Tools.descending(cnt);      // it will change cnt list as well!!!!, indexes start with 0
//        int clusters_nr = (desc_idx.length>4)? 4 : desc_idx.length ;

//        ArrayList<float[]> out_sorted = new ArrayList<float[]>(clusters_nr); // top _Nc if there are as many
        for (int i = 0; i < _Nc; i++) { // extract top _Nc clusters

            if (i<out.size()) { // covered by ms-clustering
                out4x2[i][0] = out.get(desc_idx[i])[0];
                out4x2[i][1] = out.get(desc_idx[i])[1];
            }
            else {
                out4x2[i][0] = Float.NaN;
                out4x2[i][1] = Float.NaN;
            }

        }

//		out_sorted = new float[clusters_nr][2];
//		System.out.println(clusters_nr+" clusters allocated");
//		int[] clustered_counts = new int[clusters_nr];

//        for (int ii=0; ii<clusters_nr; ii++) {
//
//            float vx 	= out.get(desc_idx[ii])[0];
//            float vy 	= out.get(desc_idx[ii])[1];
//            float vcnt 	= cnt.get(ii);  // because cnt is already sorted
//
//            out_sorted.add(new float[]{vx, vy, vcnt}); // add top 1,2,3 or 4 directions based on the count
//        }
//
//        return out_sorted;

    }

    public void mousePressed(MouseEvent e) {}
    public void mouseReleased(MouseEvent e) {}
    public void mouseEntered(MouseEvent e) {}
    public void mouseExited(MouseEvent e) {}
    public void mouseDragged(MouseEvent e) {}

}
