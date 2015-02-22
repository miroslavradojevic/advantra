package reconstruction2d;

import aux.Tools;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Line;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.process.FloatProcessor;

import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by miroslav on 15-2-15.
 */
public class Profiler2D extends Thread {

    private int begN, endN;

    public static float[][]             inimg_xy;                   // input image
    public static Zncc2D                zncc2D_module;              // module class for computing zero-mean-normalized-cross-corr.
    public static int[][]               i2xy;                       // locations mapping
    public static int[][]               xy2i;                       // coordinate mapping
    public static ArrayList<float[]>    soma_list;                  // list of soma center of masses


    // consider only those above this correlation level when creating the queue
    // all those below this relative level will not be added to queue
    // those that are within the queue are sorted by the distance towards soma so that those further away are taken first

    // output
    public static float[]   i2zncc;                   // correlations per indexed foreground location
    public static float[]   i2sigma;                  // cross-section gaussian sigma per foreground location (the one with largest zncc score)
    public static float[][] i2vxy;                  // local orientation per foreground location (the one with largest zncc score)
    public static int[]     queue;                   // all the points will be assigned with weight - those with highest weight will start first
    // weight = zncc (0-1) + min/max normalized distance towards soma center (0-1)

    public Profiler2D(int n0, int n1)
    {
        this.begN = n0;
        this.endN = n1;
    }

    public static void loadTemplate(
            Zncc2D _zncc2D,
            int[][] _i2xy,
            int[][] _xy2i,
            ArrayList<float[]> _soma_centers,
            float[][] _inimg_xy
    )
    {

        zncc2D_module = _zncc2D;
        i2xy = _i2xy;
        xy2i = _xy2i;
        soma_list = _soma_centers;
        inimg_xy = _inimg_xy;

        // allocate outputs
        i2zncc = new float[i2xy.length];
        i2sigma = new float[i2xy.length];
        i2vxy = new float[i2xy.length][2];
        // int[] queue will be allocated in createQueue() method
    }

    public void run() { // will be threaded

        for (int loc_index = begN; loc_index < endN; loc_index++)
            zncc2D_module.extract(loc_index, i2xy, inimg_xy, i2zncc, i2sigma, i2vxy); // i2zncc[loc_index]

    }

    public static void createQueue(float boundary){

        // take weights and the corresponding spatial index
        ArrayList<Integer>  sel_idxs    = new ArrayList<Integer>();

        // make the list of those that will be considered
        for (int i = 0; i < i2zncc.length; i++)
            if (i2zncc[i]>boundary) {
                sel_idxs.add(i);
//                sel_scrs.add(i2zncc[i]);
            }

//        queue = new int[sel_idxs.size()]; // allocate

        ArrayList<Float>    sel_scrs    = new ArrayList<Float>(sel_idxs.size());

        float min_dist = Float.POSITIVE_INFINITY;
        float max_dist = Float.NEGATIVE_INFINITY;

        // add the distance to the scores of selected locations
        for (int i = 0; i < sel_idxs.size(); i++) {

            float currx = i2xy[sel_idxs.get(i)][0];
            float curry = i2xy[sel_idxs.get(i)][1];

            // distance towards the first soma (it has to exist)
            float dist = (float) (
                    Math.pow(currx-soma_list.get(0)[0],2) +
                    Math.pow(curry-soma_list.get(0)[1],2));

            // remaining somas
            for (int j = 1; j < soma_list.size(); j++) {

                float dist1 = (float) (
                                Math.pow(currx-soma_list.get(j)[0],2) +
                                Math.pow(curry-soma_list.get(j)[1],2));
                if (dist1<dist) dist = dist1;
            }

            sel_scrs.add(i, dist);

            if (dist<min_dist) min_dist = dist;
            if (dist>max_dist) max_dist = dist;

        }
        
        // loop through all the points once they're calculated
        // append the zncc score to the min/max normalized distances to soma by now
        for (int i = 0; i < sel_scrs.size(); i++) {
            float normalized_dist = (sel_scrs.get(i)-min_dist)/(max_dist-min_dist);
            sel_scrs.set(i, normalized_dist * i2zncc[sel_idxs.get(i)]);
        }

        queue = Tools.descendingFloat(sel_scrs);
        for (int i = 0; i < queue.length; i++) {
            queue[i] = sel_idxs.get(queue[i]);
        }

    }

    public static ImagePlus getQueueMap(){

        int img_w = inimg_xy.length;
        int img_h = inimg_xy[0].length;

        Overlay ov = new Overlay();

        for (int i = 0; i < 1; i++) { // queue.length

            int x = i2xy[queue[i]][0];
            int y = i2xy[queue[i]][1];
            int idx = y * img_w + x;
            float rr = 1.6f;//i2sigma[queue[i]];

            // oval with scale and location
            OvalRoi ovr = new OvalRoi(x-rr+.5, y-rr+.5, 2*rr, 2*rr);
            float colr = 1 - i/(1.5f*(queue.length-1));
            Color col = new Color(colr,colr,0,colr);
            ovr.setStrokeColor(col);
            ovr.setFillColor(col);

            // directions
            Line l = new Line(x+.0, y+.0, x+.0 + 2*rr*i2vxy[queue[i]][0], y+.0 + 2*rr*i2vxy[queue[i]][1]);
            l.setStrokeColor(col);
            l.setFillColor(col);

            ov.add(ovr);
            ov.add(l);

        }

        ImagePlus img = new ImagePlus("queue", new FloatProcessor(inimg_xy));
        img.setOverlay(ov);
        return img;

    }

    public static ImagePlus getScores(){

        int img_w = inimg_xy.length;
        int img_h = inimg_xy[0].length;

        // will give stack with results of this stage
//        ImageStack is_out = new ImageStack(img_w, img_h);
        float[] corr_scrs_array = new float[img_w*img_h];
        for (int i = 0; i < i2zncc.length; i++) {

            int x = i2xy[i][0];
            int y = i2xy[i][1];
            int idx = y * img_w + x;
            corr_scrs_array[idx] = i2zncc[i];

        }

//        FloatProcessor img_orig = ;
//        ImagePlus impl_orig = new ImagePlus("orig", img_orig);
//        impl_orig.setOverlay(ov);
//        is_out.addSlice("vectors", new FloatProcessor(inimg_xy));
//        is_out.addSlice("correlations", );

        ImagePlus out = new ImagePlus("zncc", new FloatProcessor(img_w, img_h, corr_scrs_array));
//        out.setOverlay(ov);

        return out;
    }

}
