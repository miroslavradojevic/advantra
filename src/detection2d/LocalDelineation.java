package detection2d;

import aux.Interpolator;
import aux.Stat;
import aux.Tools;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.*;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import tracing2d.BayesianTracer2D;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.geom.Arc2D;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by miroslav on 11-11-14.
 */
public class LocalDelineation implements PlugIn, MouseListener, MouseMotionListener {

    // input
//    ImagePlus   img;
    ImageCanvas cnv;
    String      image_path;
    float[][] likelihood_xy;

    // classes used for processing
    SemiCircle scirc;

    // params
    // bayesian filtering
    float radius = 1.2f;
    int Ni = 10;
    int Ns = 50;
    float sigma_deg = 40;
    // mean shift convergence
    float r = 1*radius;
    float r2 = r*r;
    int max_iter = 15;
    float epsilon = 0.001f;


    // output
    float[][][] xt = new float[Ni+1][Ns][4];    // bayesian filtering states (x,y, vx, vy)
    float[][]   wt = new float[Ni+1][Ns];       // bayesian filtering states

    float[][][] xc = new float[Ni][Ns][2];    // mean-shift convergence of the bayesian filtered states (x, y)
    float[]     xn = new float[2];            // new one for the mean-shift

    float[][][] xcc = new float[Ni][4][2];    // after clustering the converged values - take up to 4 clusters
    float[][]   dsts = new float[Ns][Ns];       // inter distances - used as a clustering criteria
    int[]       cllab = new int[Ns];                // clustering labels at each iteration

    byte[][]     tag = new byte[Ni][4];       // fmm tags (0, 1, 2)
    byte[][]     par = new byte[Ni][4];       // fmm pointers (-2, -1, 0, 1, 2, 3)
    float[][]    cst = new float[Ni][4];      // fmm scores

    ArrayList<ArrayList<float[]>> delin = new ArrayList<ArrayList<float[]>>();

    // aux
    ArrayList<int[]>    heap_vals = new ArrayList<int[]>();      // rank, Ni idx, N curr idx, N prev idx
    ArrayList<Float>    heap_scrs = new ArrayList<Float>();          // scores

    // interface - visualisations
    ImagePlus vizTag = new ImagePlus("TAG", new ByteProcessor(4,Ni,new byte[4*Ni]));
    ImagePlus vizPar = new ImagePlus("PARENT_IDX", new ByteProcessor(4,Ni,new byte[4*Ni]));
    ImagePlus vizScr = new ImagePlus("SCORE", new FloatProcessor(4,Ni,new float[4*Ni]));

    Overlay ov = new Overlay();

    public void run(String s) {

        // load the image through the menu
        String in_folder = Prefs.get("id.folder", System.getProperty("user.home"));
        OpenDialog.setDefaultDirectory(in_folder);
        OpenDialog dc = new OpenDialog("Select file");
        in_folder = dc.getDirectory();
        image_path = dc.getPath();
        if (image_path==null) return;
        Prefs.set("id.folder", in_folder);

        ImagePlus img = new ImagePlus(image_path);
        if(img==null) return;

        // set image as float[][]
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

        scirc = new SemiCircle(radius);

        // show it and get the canvas
        img.show();
        cnv = img.getCanvas();
        cnv.addMouseListener(this);
        cnv.addMouseMotionListener(this);
        IJ.setTool("hand");

//        vizTag.show();
//        vizTag.getCanvas().zoom100Percent();
//        vizTag.getCanvas().zoomIn(0,0);
//        vizTag.getCanvas().zoomIn(0,0);
//        vizTag.getCanvas().zoomIn(0,0);
//        vizTag.getCanvas().zoomIn(0,0);
//        vizTag.getCanvas().zoomIn(0,0);
//        vizTag.getCanvas().zoomIn(0,0);
//        vizTag.getCanvas().zoomIn(0,0);
//        vizTag.getCanvas().zoomIn(0,0);
//        vizTag.getCanvas().zoomIn(0,0);


    }

    public void mouseClicked(MouseEvent e)
    {
        int clickX = cnv.offScreenX(e.getX());
        int clickY = cnv.offScreenY(e.getY());

        System.out.print("(x,y) -> "+clickX+" , " +clickY+" -> ");

        long t1 = System.currentTimeMillis();

        // bayesian filtering
        BayesianTracer2D.run(clickX, clickY,
                likelihood_xy, scirc, sigma_deg, Ni, Ns, xt, wt);

        // mean-shift for estimation (local maxima detection)
        for (int i = 0; i < xc.length; i++) { // initialize xc with corresponding elements from xt
            for (int s = 0; s < xc[0].length; s++) {
                xc[i][s][0] = xt[i+1][s][0];
                xc[i][s][1] = xt[i+1][s][1];
            }
        }

        for (int i = 0; i < Ni; i++) {
            meanshift_xy_convg(xt[i+1], wt[i+1], xc[i], xn, max_iter, epsilon, r2); // xc will hold the converged values
            // clustering xc->xcc
            dists2(xc[i], dsts); // Ns*Ns distances at each iteration
            clustering(dsts, r2, cllab);
            extracting(cllab, xc[i], 1, xcc[i]);
        }

        // fmm - will need to march 4 times
//        for (int i = 0; i < Ni; i++) for (int j = 0; j < 4; j++) tag[i][j] = (byte) ( 0); // initialize tags
//        for (int i = 0; i < Ni; i++) for (int j = 0; j < 4; j++) par[i][j] = (byte) (-2); // initialize pointers
//        for (int i = 0; i < Ni; i++) for (int j = 0; j < 4; j++) cst[i][j] = Float.POSITIVE_INFINITY; // initialize scores to the highest
//        heap_scrs.clear();
//        heap_vals.clear();
        delin.clear();

        for (int i = 0; i < 4; i++) {

            fmm(likelihood_xy, new float[]{clickX, clickY}, xcc, tag, par, cst, heap_vals, heap_scrs);

            ArrayList<float[]> strm = extract_streamline(3, tag, par, cst, xcc); // xcc will be modified

            if (strm==null) {
                break;
            }
            else {
                delin.add(strm);
            }

        }

        long t2 = System.currentTimeMillis();

        System.out.println("     "+ ((t2-t1)/1000f) + "sec.");

//        vizTag.setProcessor(getTag());
//        vizTag.updateAndDraw();
//        IJ.run(vizTag, "Enhance Contrast", "saturated=0.35");
//        new ImagePlus("", getTag()).show();
////        System.out.println("TAG:");
//        for (int i = 0; i < tag.length; i++) {
//            for (int j = 0; j < tag[i].length; j++) {
//                System.out.print(tag[i][j]+"\t"+par[i][j]+"\t"+cst[i][j]+" | ");
//            }
//            System.out.println();
//        }
        Overlay ov_xt   = viz_xt(xt, wt);
        Overlay ov_xc   = viz_xc(xc);
        Overlay ov_xcc  = viz_xcc(xcc);
        Overlay ov_delin = viz_delin(delin);

        ov.clear();
//        for (int i = 0; i <ov_xt.size(); i++) ov.add(ov_xt.get(i));
        for (int i = 0; i <ov_xc.size(); i++) ov.add(ov_xc.get(i));
        for (int i = 0; i <ov_xcc.size(); i++) ov.add(ov_xcc.get(i));
        for (int i = 0; i <ov_delin.size(); i++) ov.add(ov_delin.get(i));

        cnv.setOverlay(ov);

    }

    private ByteProcessor getTag()
    {
        byte[] out = new byte[4*Ni];
        int cnt = 0;
        for (int i = 0; i <Ni; i++) for (int j = 0; j < 4; j++) out[cnt++] = tag[i][j];
        return new ByteProcessor(4, Ni, out);
    }

    public static Overlay viz_xt(float[][][] xt, float[][] wt)
    {

        float rad = .5f;

        Overlay ov = new Overlay();
        for (int i = 0; i < xt.length; i++) { // iterations

            Stat.min_max_normalize(wt[i]);

            for (int j = 0; j <xt[i].length; j++) {

                OvalRoi p = new OvalRoi(xt[i][j][0]-rad+.5, xt[i][j][1]-rad+.5, 2*rad, 2*rad);
                p.setFillColor(new Color(1f,1f,1f,wt[i][j]));
//                Line l = new Line(xt[i][j][0]+.5, xt[i][j][1]+.5, xt[i][j][0]+xt[i][j][2]+.5, xt[i][j][1]+xt[i][j][3]+.5);

                PointRoi pt = new PointRoi(xt[i][j][0]+.5, xt[i][j][1]+.5);
                pt.setStrokeColor(new Color(1f,1f,1f,wt[i][j]));
                pt.setFillColor(new Color(1f,1f,1f,wt[i][j]));
                ov.add(p);
                ov.add(pt);
//                ov.add(l);

            }
        }

        return ov;
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

        float rad = .5f;

        Overlay ov = new Overlay();
        for (int i = 0; i < xcc.length; i++) { // iterations

//            Stat.min_max_normalize(wt[i]);

            for (int j = 0; j <xcc[i].length; j++) {

                if (!Float.isNaN(xcc[i][j][0])) {
                    PointRoi p = new PointRoi(xcc[i][j][0]+.5, xcc[i][j][1]+.5);
                    p.setFillColor(Color.GREEN);
                    p.setStrokeColor(Color.GREEN);
                    ov.add(p);
                }

//                p.setFillColor(new Color(1f,1f,1f,wt[i][j]));
//                Line l = new Line(xt[i][j][0]+.5, xt[i][j][1]+.5, xt[i][j][0]+xt[i][j][2]+.5, xt[i][j][1]+xt[i][j][3]+.5);
//                ov.add(l);

            }
        }

        return ov;
    }

    public static Overlay viz_delin(ArrayList<ArrayList<float[]>> delin)
    {
        Overlay ov = new Overlay();

        for (int i = 0; i < delin.size(); i++) {

            // add the polygonroi that corresponds to the extracted streamline
            float[] plgx = new float[delin.get(i).size()];
            float[] plgy = new float[delin.get(i).size()];

            for (int j = 0; j < delin.get(i).size(); j++) {

//                OvalRoi or = new OvalRoi(delin.get(i).get(j)[0]-1f+.5f, delin.get(i).get(j)[1]-1f+.5f, 2f, 2f);
//                or.setFillColor(Color.RED);
//                or.setStrokeColor(Color.RED);

//                PointRoi pt = new PointRoi(delin.get(i).get(j)[0]+.5f, delin.get(i).get(j)[1]+.5f);
//                pt.setFillColor(Color.RED);
//                pt.setStrokeColor(Color.RED);

                plgx[j] = delin.get(i).get(j)[0]+.5f;
                plgy[j] = delin.get(i).get(j)[1]+.5f;

            }

            PolygonRoi plg = new PolygonRoi(plgx, plgy, Roi.FREELINE);
            if (i==0)       plg.setStrokeColor(Color.RED);
            else if (i==1)  plg.setStrokeColor(Color.ORANGE);
            else if (i==2)  plg.setStrokeColor(Color.BLUE);
            else if (i==3)  plg.setStrokeColor(Color.YELLOW);
            plg.setStrokeWidth(.25f);
            ov.add(plg);

        }

        return ov;

    }

    public void mousePressed(MouseEvent e) {}
    public void mouseReleased(MouseEvent e) {}
    public void mouseEntered(MouseEvent e) {}
    public void mouseExited(MouseEvent e) {}
    public void mouseDragged(MouseEvent e) {}

    public void mouseMoved(MouseEvent e)
    {
        mouseClicked(e);
    }

    // fast marching method algorithm: xcc -> par,tag,scr
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
        float I = 0.5f;
        float ALFA = (float) (Math.PI/2);
        float P = 10; // multiple of D should be parametrized

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

                float I1 = Interpolator.interpolateAt(curr_x, curr_y, inimg_xy);
                float I0 = Interpolator.interpolateAt(start_xy[0], start_xy[1], inimg_xy);

                float CI = (float) Math.exp(LAMBDA * Math.pow( (I1-I0)/I,2));
                float CP = (float) Math.exp(LAMBDA * ((Math.pow(curr_x-start_xy[0],2)+Math.pow(curr_y-start_xy[1],2))/Math.pow(P,2)));
                float CALPHA = 1;

                float C = CI * CP * CALPHA; // ( (eIcurr+eIcurr)/2 );

                tag[0][i] = (byte) 1;           // tag it as narrow band
                par[0][i] = (byte) par_idx;     // parent is root xy
                cst[0][i] = C;                  // cost should accumulate the previous cost which is 0

                if (heap_vals.size()==0) { // just insert into the heap
                    heap_vals.add(new int[]{1, ni_idx, pos_idx, par_idx}); // rank, Ni idx, N curr idx, N parent idx
                    heap_scrs.add(C);
                }
                else{ // append it and update the rank depending on where the cost was
                    int cnt = 0;
                    for (int j = 0; j < heap_scrs.size(); j++) if (heap_scrs.get(j) > e) {cnt++; heap_vals.get(j)[0]++;} // count them & decrease ranking

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

                            float x1 = xcc[at_i+1][i][0]; // neighbour of the element taken from the heap
                            float y1 = xcc[at_i+1][i][1];

                            float x0 = xcc[at_i][at_p][0]; // just taken from the heap
                            float y0 = xcc[at_i][at_p][1];

                            // see what is xp, yp
                            float xp= Float.NaN, yp=Float.NaN;
                            if (par[at_i][at_p]==(byte)-1) {      // if it was root point
                                xp = start_xy[0];
                                yp = start_xy[1];
                            }
                            else if (                       // if it was legal
                                    par[at_i][at_p]==(byte)0 ||
                                    par[at_i][at_p]==(byte)1 ||
                                    par[at_i][at_p]==(byte)2 ||
                                    par[at_i][at_p]==(byte)3
                                    )
                            {
                                xp = xcc[at_i-1][par[at_i][at_p]&0xff][0];
                                yp = xcc[at_i-1][par[at_i][at_p]&0xff][1];
                            }
                            else { // it was not a legal parent index
                                System.out.println("wrong parent index when calculating the score");
                                System.exit(0);
                            }

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



                            float cosang = v1x*v2x+v1y*v2y;// / (*Math.sqrt(v2x*v2x+v2y*v2y)));
                            float eD = 2 - cosang; //(float) (Math.pow(curr_x-prev_x,2)+Math.pow(curr_y-prev_y,2)); // directional cost of the neighbour
                            float eD1 = v1;

                            float Icurr = Interpolator.interpolateAt(x1, y1, inimg_xy);
//                            float Iprev = Interpolator.interpolateAt(prev_x, prev_y, inimg_xy);

                            float eI = 1;// (float) Math.exp(5*Math.pow(1-Icurr/1f,2)); // todo use intensity closeness can be a measure
//                            float eIprev = (float) Math.exp(10*Math.pow(1-Iprev/1f,2));

                            float e = eD * eI * eD1 + cst[at_i][at_p];// new cost is cummulative + ( (eIcurr+eIcurr)/2f );

                            if (tag[at_i+1][i]!=(byte)1) {// neighbour (at_i+1, i) is NOT in narrow band

                                tag[at_i+1][i] = (byte) 1; // tag it as narrow band
                                par[at_i+1][i] = (byte) at_p;
                                cst[at_i+1][i] = e;

                                // calculate the rank depending on where the cost was
                                int cnt = 0;
                                for (int j = 0; j < heap_scrs.size(); j++) if (heap_scrs.get(j) > e) {cnt++; heap_vals.get(j)[0]++;}
                                int new_rank = heap_scrs.size() - cnt + 1; // ranking starts from 1

                                // append it to the heap
                                heap_scrs.add(e);
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

                                if (e<found_score) {

                                    /*
                                        heap updates
                                     */

                                    // update rankings
                                    float old_val = found_score;
                                    float new_val = e; // new_val is smaller!!

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
                                    cst[at_i+1][i] = e;

//                                    System.out.println("neighbour WAS in narrow band and replaced : "+new_val);
//                                    System.out.print("ranks:: ");
//                                    for (int ii = 0; ii <heap_vals.size(); ii++) {
//                                        System.out.print(heap_vals.get(ii)[0] + " , ");
//                                    }
//                                    System.out.println("-----");

                                }

                            }

                        }

                    }

                }

            }
//            else {}// reacherd the end, no neighbours

        }

//        System.out.println("HEAP:");
//        for (int i = 0; i < heap_scrs.size(); i++) {
//            System.out.println(heap_scrs.get(i) + " | " + Arrays.toString(heap_vals.get(i))+ "      " + Arrays.toString(tag[0]));
//        }

    }

    // tag,par,scr -> xcc
    private static ArrayList<float[]> extract_streamline(
            int margin_end,                                 // will limit the point up to where we search when tracing back
//            int margin_begin,                               // will limit the
            byte[][] tag, // input
            byte[][] par, // input
            float[][] cst, // input
            float[][][] xcc // side output, extracted path will be excluded
    )
    {

        // use tag and cst to find initial thread point - most distant frozen that had the least cost
        boolean found = false;
        int found_i = -1;
        int found_j = -1;
        float min_cost;

        for (int i = tag.length-1; i > tag.length-1-margin_end; i--) {

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

        }

        if (found) { // found_i, found_j

            ArrayList<float[]> out = new ArrayList<float[]>();

            for (int loop_streamiline = found_i; loop_streamiline >=0 ; loop_streamiline--) {

                out.add(xcc[loop_streamiline][found_j].clone()); // start from found_i,j

                // cancel the one that was added - set it as NaN for future fast marching iteration
                //if (loop_streamiline>=margin_begin) {
                    xcc[loop_streamiline][found_j][0] = Float.NaN;
                    xcc[loop_streamiline][found_j][1] = Float.NaN;
                //}

                found_j = par[loop_streamiline][found_j]&0xff; // recursively go backwards

            }

            return out;

        }
        else {

            return null;

        }

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

//            System.out.print("\nelement " + i + " -> " + Arrays.toString(xy_conv[i]));

            int iter = 0;
            float d;

            do {

//                System.out.println("iter " + iter + " :: " + Arrays.toString(xy_conv[i]));

                runone(xy_conv[i], xy_next, rad2, xy_inpt, w_inpt);

                d = (float) (Math.pow(xy_next[0]-xy_conv[i][0], 2) + Math.pow(xy_next[1]-xy_conv[i][1], 2));

                xy_conv[i][0] = xy_next[0];
                xy_conv[i][1] = xy_next[1];

                iter++;

            }
            while (iter < max_iter && d > epsilon*epsilon); //

//            System.out.println("-> " + Arrays.toString(xy_conv[i]) + " in " + iter + " iterations");

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

    private static void extracting(int[] labels, float[][] vec_xy, int min_count, float[][] out4x2)
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

                if (count >= min_count) {
                    out.add(new float[]{centroid_x/count, centroid_y/count});
                    cnt.add(count);
                }

            }
        }


        // print before sorting
//		for (int ii = 0; ii < cnt.size(); ii++) {
//			IJ.log(ii + " : " + cnt.get(ii) + " points,  at " + Arrays.toString(out.get(ii)));
//		}

        // sort by the counts (take from Sorting.java) descending indexes
        int[] desc_idx = Tools.descending(cnt);      // it will change cnt list as well!!!!
        int clusters_nr = (desc_idx.length>4)? 4 : desc_idx.length ;

//        ArrayList<float[]> out_sorted = new ArrayList<float[]>(clusters_nr); // top 4  if there are as many
        for (int i = 0; i < 4; i++) {

            if (i<desc_idx.length) {
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

//    private static void take_streamline(
//            byte[][]  tag_map,
//            float[][] cst_map,
//            byte[][]  par_map,
//            float[][][] xcc,
//            ArrayList<float[]> out_stream // output
//    )
//    {
//    }

}
