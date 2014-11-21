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
import tracing2d.BayesianTracer2D;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
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
    float radius = 1.5f;
    int Ni = 10;
    int Ns = 50;
    float sigma_deg = 40;
    // mean shift convergence
    float r = 1*radius;
    float r2 = r*r;
    int max_iter = 15;
    float epsilon = 0.001f;


    // output
    float[][][] xt = new float[Ni][Ns][4];    // bayesian filtering states (x,y, vx, vy)
    float[][]   wt = new float[Ni][Ns];       // bayesian filtering states

    float[][][] xc = new float[Ni][Ns][2];    // mean-shift convergence of the bayesian filtered states (x, y)
    float[]     xn = new float[2];            // new one for the mean-shift

    float[][][] xcc = new float[Ni][4][2];    // after clustering the converged values - take up to 4 clusters
    float[][]   dsts = new float[Ns][Ns];       // inter distances - used as a clustering criteria
    int[]       cllab = new int[Ns];                // clustering labels at each iteration

    byte[][]     tag = new byte[Ni][4];       // fmm tags (0, 1, 2)
    byte[][]     par = new byte[Ni][4];       // fmm pointers (-2, -1, 0, 1, 2, 3)
    float[][]    cst = new float[Ni][4];      // fmm scores

    ArrayList<int[]>    heap_vals = new ArrayList<int[]>();      // rank, Ni idx, N curr idx, N prev idx
    ArrayList<Float>    heap_scrs = new ArrayList<Float>();          // scores

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


    }

    public void mouseClicked(MouseEvent e) {
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
                xc[i][s][0] = xt[i][s][0];
                xc[i][s][1] = xt[i][s][1];
            }
        }

        for (int i = 0; i < Ni; i++) {
            meanshift_xy_convg(xt[i], wt[i], xc[i], xn, max_iter, epsilon, r2); // xc will hold the converged values
            // clustering xc->xcc
            dists2(xc[i], dsts); // Ns*Ns distances at each iteration
            clustering(dsts, r2, cllab);
            extracting(cllab, xc[i], 1, xcc[i]);
        }

        // fmm - will need 4 times
        // frozen node will have score, tag and the pointer calculated
        for (int i = 0; i < Ni; i++) for (int j = 0; j < 4; j++) tag[i][j] = (byte) ( 0); // initialize tags
        for (int i = 0; i < Ni; i++) for (int j = 0; j < 4; j++) par[i][j] = (byte) (-2); // initialize pointers
        for (int i = 0; i < Ni; i++) for (int j = 0; j < 4; j++) cst[i][j] = Float.POSITIVE_INFINITY; // initialize scores to the highest
        heap_scrs.clear();
        heap_vals.clear();

        fmm(
                likelihood_xy,
                new float[]{clickX, clickY},
                xcc,  // Ni*4*2
                tag,
                par,
                cst,
                heap_vals,
                heap_scrs
        );

        long t2 = System.currentTimeMillis();

        System.out.println("     "+ ((t2-t1)/1000f) + "sec.");

//        System.out.println("TAG:");
        for (int i = 0; i < tag.length; i++) {
            for (int j = 0; j < tag[i].length; j++) {
                System.out.print(tag[i][j]+"\t"+par[i][j]+"\t"+cst[i][j]+" | ");
            }
            System.out.println();
        }

        Overlay ov = new Overlay();
        Overlay ov_xt   = viz_xt(xt, wt);
        Overlay ov_xc   = viz_xc(xc);
        Overlay ov_xcc  = viz_xcc(xcc);

        for (int i = 0; i <ov_xt.size(); i++) ov.add(ov_xt.get(i));
        for (int i = 0; i <ov_xc.size(); i++) ov.add(ov_xc.get(i));
        for (int i = 0; i <ov_xcc.size(); i++) ov.add(ov_xcc.get(i));

        cnv.setOverlay(ov);

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

                ov.add(p);
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

    public void mousePressed(MouseEvent e) {}
    public void mouseReleased(MouseEvent e) {}
    public void mouseEntered(MouseEvent e) {}
    public void mouseExited(MouseEvent e) {}
    public void mouseDragged(MouseEvent e) {}

    public void mouseMoved(MouseEvent e)
    {
        mouseClicked(e);
    }

    private static void fmm(
                            float[][] inimg_xy,     // input data
                            float[] start_xy,       //
                            float[][][] xcc,        // Ni*4*2
                            byte[][] tag,           //
                            byte[][] par,
                            float[][] cst,
                            ArrayList<int[]> heap_vals,
                            ArrayList<Float> heap_scrs
    )
    {

        int startII = xcc.length-1; // start iteration index fmm will start backwards from the last iteration
        // xcc contains the

        // fast marching method algorithm
        for (int i = 0; i < xcc[startII].length; i++) { // check 4 of neighbours at iteration=Ni (last one)

            if (!Float.isNaN(xcc[startII][i][0])) { // if there is a centroid at this iteration
// todo...
                float curr_x = xcc[0][i][0];
                float curr_y = xcc[0][i][1];

                float eD = 0;//(float) (Math.pow(curr_x-start_xy[0],2)+Math.pow(curr_y-start_xy[1],2));

                float Icurr = Interpolator.interpolateAt(curr_x, curr_y, inimg_xy);
                float Iprev = Interpolator.interpolateAt(start_xy[0], start_xy[1], inimg_xy);

                float eIcurr = (float) Math.exp(10*Math.pow(1-Icurr/1f,2));
                float eIprev = (float) Math.exp(10*Math.pow(1-Iprev/1f,2));

                float e = eD * ( (eIcurr+eIprev)/2 );

                tag[0][i] = (byte) 1; // tag it as narrow band
                par[0][i] = (byte) (-1);
                cst[0][i] = e;

                if (heap_vals.size()==0) {
                    // just insert into the heap
                    heap_vals.add(new int[]{1, 0, i, -1}); // rank, Ni idx, N curr idx, N parent idx
                    heap_scrs.add(e);
                }
                else{

                    // append it and update the rank depending on where the cost was
                    int cnt = 0;
                    for (int j = 0; j < heap_scrs.size(); j++) if (heap_scrs.get(j) > e) {cnt++; heap_vals.get(j)[0]++;}

                    int new_rank = heap_scrs.size() - cnt + 1; // ranking starts from 1

                    // add the new one
                    heap_scrs.add(e);
                    heap_vals.add(new int[]{new_rank, 0, i, -1});

                }

            }

        }

//        System.out.println("looping FMM...");

        boolean do_it = true;

        while (heap_vals.size()>0 && do_it) {

//            System.out.println("while again");
//            for (int i = 0; i <heap_vals.size(); i++) {
//                System.out.print(heap_vals.get(i)[0] + " , ");
//            }
//            System.out.println("-----");

            // get the top score
            int get_lowest_score = -1;
            for (int i = 0; i < heap_vals.size(); i++)
                if (heap_vals.get(i)[0] == 1)
                   get_lowest_score = i;
                else
                   heap_vals.get(i)[0]--;

            if (get_lowest_score==-1) {
                System.out.println("wrong !!!! cannot find rank 1");
//                for (int i = 0; i <heap_vals.size(); i++) {
//                    System.out.println(heap_vals.get(i)[0]);
//                }
                System.exit(0);
            }

            int at_i = heap_vals.get(get_lowest_score)[1]; // iteration index
            int at_p = heap_vals.get(get_lowest_score)[2]; // curr idx within iteration
            int at_pp = heap_vals.get(get_lowest_score)[3]; // parent idx

            // freeze it
            tag[at_i][at_p] = (byte) 2; // tag it as freezed
            par[at_i][at_p] = (byte) (at_pp);
            cst[at_i][at_p] = heap_scrs.get(get_lowest_score);

            // remove from heap list
            heap_scrs.remove(get_lowest_score);
            heap_vals.remove(get_lowest_score);

            // loop neighbours of at_i, at_p (some from the following iteration) if there is such
            if (at_i<xcc.length-1) { // if it is within the boundaries

//                int at_i_nbr = (at_i)? : ;

                for (int i = 0; i < xcc[at_i+1].length; i++){ // looping the neighbours

                    if (!Float.isNaN(xcc[at_i+1][i][0])) { // there was a cluster

//                        System.out.println("" + (at_i+1) + " :: " + i);

                        if (tag[at_i+1][i]!=2){ // it is not frozen

                            // compute dist

                            float curr_x = xcc[at_i+1][i][0];
                            float curr_y = xcc[at_i+1][i][1];

                            float prev_x = xcc[at_i][at_p][0];
                            float prev_y = xcc[at_i][at_p][1];

                            float eD = (float) (Math.pow(curr_x-prev_x,2)+Math.pow(curr_y-prev_y,2));

                            float Icurr = Interpolator.interpolateAt(curr_x, curr_y, inimg_xy);
                            float Iprev = Interpolator.interpolateAt(prev_x, prev_y, inimg_xy);

                            float eIcurr = (float) Math.exp(10*Math.pow(1-Icurr/1f,2));
                            float eIprev = (float) Math.exp(10*Math.pow(1-Iprev/1f,2));

                            float e = eD * ( (eIcurr+eIprev)/2 );

                            if (tag[at_i+1][i]!=1) {// neighbour (at_i+1, i) is NOT in narrow band

                                tag[at_i+1][i] = 1; // tag it as narrow band
                                par[at_i+1][i] = (byte) at_p;
                                cst[at_i+1][i] = e;

                                // append it to the heap
                                // append it and update the rank depending on where the cost was
                                int cnt = 0;
                                for (int j = 0; j < heap_scrs.size(); j++) if (heap_scrs.get(j) > e) {cnt++; heap_vals.get(j)[0]++;}

                                int new_rank = heap_scrs.size() - cnt + 1; // ranking starts from 1

                                // add the new one
                                heap_scrs.add(e);
                                heap_vals.add(new int[]{new_rank, at_i+1, i, at_p});

//                                System.out.println("neighbour was not in narrow band, added rank " + new_rank);

                            }
                            else {
                                // neighbour IS in narrow band already - means that it is already added to the heap - so update -
                                // find the corresponding heap location and update it

                                boolean found = false;
                                int idx_found = -1;
                                for (int j = 0; j < heap_vals.size(); j++) {
                                    if (heap_vals.get(j)[1]==at_i+1 && heap_vals.get(j)[2]==i) {
                                        found = true;
                                        idx_found = j;
                                        break;
                                    }
                                }

                                if (!found) {System.out.println("#### problem finding in the HEAP ####");System.exit(0);}

                                float found_score = heap_scrs.get(idx_found);

                                if (e<found_score) {

                                    // update parent
                                    heap_vals.get(idx_found)[3] = at_p;

                                    // update rankings
                                    float old_val = found_score;
                                    float new_val = e; // new_val is smaller!!

                                    int new_rank = 1;
                                    for (int j = 0; j < heap_scrs.size(); j++) {

                                        if (heap_scrs.get(j)<old_val && heap_scrs.get(j) >= new_val)
                                            heap_vals.get(j)[0]++;

                                        if (heap_scrs.get(j)<new_val) new_rank++;

                                    }
                                    heap_vals.get(idx_found)[0] = new_rank;

                                    // replace the score with the cheaper one once the rankings are corrected
                                    heap_scrs.set(idx_found, new_val);

//                                    System.out.println("neighbour WAS in narrow band and replaced : "+new_val);
//                                    System.out.print("ranks:: ");
//                                    for (int ii = 0; ii <heap_vals.size(); ii++) {
//                                        System.out.print(heap_vals.get(ii)[0] + " , ");
//                                    }
//                                    System.out.println("-----");


                                }
                                else{
//                                    System.out.println("neighbour WAS in narrow band and NOT replaced");
                                }// nothing the cheaper one is already there

                            }

                        }
                        else{
                            // it is frozen - do nothing, that would make closed loop
                        }


                    }

                }

            }
            else {}// reacherd the end, no neighbours

        }


//        System.out.println("HEAP:");
//        for (int i = 0; i < heap_scrs.size(); i++) {
//            System.out.println(heap_scrs.get(i) + " | " + Arrays.toString(heap_vals.get(i))+ "      " + Arrays.toString(tag[0]));
//        }

//        System.out.println("TAG");

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

//
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

    private static void take_streamline(
            byte[][]  tag_map,
            float[][] cst_map,
            byte[][]  par_map,
            float[][][] xcc,
            ArrayList<float[]> out_stream // output
    )
    {
        //


    }

}
