package reconstruction;

import aux.Tools;
import ij.ImagePlus;
import ij.gui.Line;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.process.FloatProcessor;

import java.awt.*;
import java.io.*;
import java.util.ArrayList;

/**
 * Created by miroslav on 15-2-15.
 */
public class Profiler2D extends Thread {

    private int begN, endN;

    public static float[][]             inimg_xy;                   // input image
    public static Zncc2D                zncc2D_module;              // module class for computing zero-mean-normalized-cross-corr.
    public static int[][]               i2xy;                       // locations mapping
    public static ArrayList<double[]>   soma_list;                  // list of soma center of masses

    // output
    public static float[]   i2zncc;                     // correlations per indexed foreground location
    public static float[]   i2sigma;                    // cross-section gaussian sigma per foreground location (the one with largest zncc score)
    public static float[][] i2vxy;                      // local orientation per foreground location (the one with largest zncc score)

    public Profiler2D(int n0, int n1)
    {
        this.begN = n0;
        this.endN = n1;
    }

    public static void loadTemplate(
            Zncc2D _zncc2D,
            int[][] _i2xy,
            ArrayList<double[]> _soma_centers,
            float[][] _inimg_xy
    )
    {

        zncc2D_module = _zncc2D;
        i2xy = _i2xy;
        soma_list = _soma_centers;
        inimg_xy = _inimg_xy;

        // allocate outputs
        i2zncc = new float[i2xy.length];
        i2sigma = new float[i2xy.length];
        i2vxy = new float[i2xy.length][2];
        // int[] queue will be allocated in createQueue() method
    }

    public void run() { // will be threaded

        float[] vals = new float[zncc2D_module.limR*zncc2D_module.limR];

        for (int loc_index = begN; loc_index < endN; loc_index+=2)
            zncc2D_module.extract(loc_index, i2xy, inimg_xy, i2zncc, i2sigma, i2vxy, vals);

    }

    public static int[] createQueue(float boundary, int weighting_type){

        ArrayList<Integer>  sel_idxs    = new ArrayList<Integer>();

        for (int i = 0; i < i2zncc.length; i++)
            if (i2zncc[i]>boundary) {
                sel_idxs.add(i);
            }

        ArrayList<Float>    sel_scrs    = new ArrayList<Float>(sel_idxs.size());

        if (weighting_type==2){
            // don't consider distances to soma
            for (int i = 0; i < sel_idxs.size(); i++) {
                sel_scrs.add(i, i2zncc[sel_idxs.get(i)]);
            }
        }
        else if (weighting_type==0 || weighting_type==1) {
            if (soma_list.size() >= 1) {
                // there is 1+ soma's

                float min_dist = Float.POSITIVE_INFINITY;
                float max_dist = Float.NEGATIVE_INFINITY;

                // add the distance to the scores of selected locations
                for (int i = 0; i < sel_idxs.size(); i++) {

                    float currx = i2xy[sel_idxs.get(i)][0];
                    float curry = i2xy[sel_idxs.get(i)][1];
//                float currz = i2xy[sel_idxs.get(i)][2];

                    // squared distance towards the first soma (it has to exist)
                    float dist = (float) (
                            Math.pow(currx - soma_list.get(0)[0], 2) +
                                    Math.pow(curry - soma_list.get(0)[1], 2) //+
//                                Math.pow(currz-soma_list.get(0)[2],2)
                    );

                    // remaining somas as much as soma_list allows
                    for (int j = 1; j < soma_list.size(); j++) {

                        float dist1 = (float) (
                                Math.pow(currx - soma_list.get(j)[0], 2) +
                                        Math.pow(curry - soma_list.get(j)[1], 2) //+
//                                    Math.pow(currz-soma_list.get(j)[2],2)
                        );
                        if (dist1 < dist) dist = dist1;
                    }

                    sel_scrs.add(i, dist);

                    if (dist < min_dist) min_dist = dist;
                    if (dist > max_dist) max_dist = dist;

                }

                // loop through all the points once they're calculated append the zncc score to the min/max normalized distances to soma by now
                for (int i = 0; i < sel_scrs.size(); i++) {

                    float normalized_dist = -1;

                    if (weighting_type == 0) { // those closer to the soma are weighted higher
                        normalized_dist = 1 - (sel_scrs.get(i) - min_dist) / (max_dist - min_dist);
                    } else if (weighting_type == 1) { // those that are the most distant from soma are weighted higher
                        normalized_dist = (sel_scrs.get(i) - min_dist) / (max_dist - min_dist);
                    }

                    sel_scrs.set(i, normalized_dist * i2zncc[sel_idxs.get(i)]);

                }

            }
            else {
                // there is no soma, make list without taking into account the distances towards soma
                for (int i = 0; i < sel_idxs.size(); i++) {
                    sel_scrs.add(i, i2zncc[sel_idxs.get(i)]);
                }
            }
        }


        int[] queue = Tools.descendingFloat(sel_scrs);
        for (int i = 0; i < queue.length; i++) {
            queue[i] = sel_idxs.get(queue[i]);
        }

        return queue;

    }

    public static boolean[][] getQueueLocationMap(int[] queue){

        // queue[] contains indexes in the foreground list i2xy[][]
        // to obtain the locations in image it is necessary to map, by x = i2xy[queue[]][0], y = i2xy[queue[]][1]
        // will mark the locations in the 2D map with the dimensions of the input image
        // T at locations (x,y) where queue element was found
        // F at locations (x,y) where there was no queue element
        // used to guide the bayesian tracing

        // original image dimensions
        int img_w = inimg_xy.length;
        int img_h = inimg_xy[0].length;

        boolean[][] queue_map_xy = new boolean[img_w][img_h]; // allocate output, initialize with falses

        for (int i = 0; i < queue.length; i++) {

            int x = i2xy[queue[i]][0];
            int y = i2xy[queue[i]][1];

            queue_map_xy[x][y] = true;

        }

        return queue_map_xy;

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

        ImagePlus out = new ImagePlus("zncc", new FloatProcessor(img_w, img_h, corr_scrs_array));

        return out;
    }

    public static void getQueueLocationSwc(int[] queue, String swc_path) {

        // queue[] contains indexes in the foreground list i2xyz[][]
        // to obtain the locations in image it is necessary to map, by x = i2xyz[queue[]][0], y = i2xyz[queue[]][1], z = i2xyz[queue[]][2]
        // will mark the locations in the 3d map with the dimensions of the input image

        // it is not just the locations but each has radius, and direction and they are used to initialize the bayesian tracing
        // also to guide it and to confirm it

        // original image dimensions
//        int img_w = inimg_xyz.length;
//        int img_h = inimg_xyz[0].length;
//        int img_l = inimg_xyz[0][0].length;

        PrintWriter logWriter = null;
        try {
            logWriter = new PrintWriter(swc_path);
            logWriter.print("");
            logWriter.close();
        } catch (FileNotFoundException ex) {}

        try {
            logWriter = new PrintWriter(new BufferedWriter(new FileWriter(swc_path, true)));
            logWriter.println("# ");
        } catch (IOException e) {}

        int cnt = 1;

        for (int i = 0; i < queue.length; i++) {

            float x = i2xy[queue[i]][0];
            float y = i2xy[queue[i]][1];

            float vx = i2vxy[queue[i]][0];
            float vy = i2vxy[queue[i]][1];

            float rr = i2sigma[queue[i]];

            logWriter.println(String.format(
                    "%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f %-4d",
                    cnt, 2, x, y, 0f, rr, -1));

            // directions
            logWriter.println(String.format(
                    "%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f %-4d",
                    cnt+1, 2, x, y, 0f, 0.1f, -1));

            logWriter.println(String.format(
                    "%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f %-4d",
                    cnt+2, 2, x+2*rr*vx, y+2*rr*vy, 0f, 0.1f, cnt+1));

            logWriter.println(String.format(
                    "%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f %-4d",
                    cnt+3, 2, x-2*rr*vx, y-2*rr*vy, 0f, 0.1f, cnt+1));

            cnt+=4;

        }

        logWriter.close();

    }


    public static Overlay getQueueLocationOverlay(int[] queue) {

        // queue[] contains indexes in the foreground list i2xy[][]
        // to obtain the locations in image it is necessary to map, by x = i2xy[queue[]][0], y = i2xy[queue[]][1]
        // will mark the locations in the 2d map with the dimensions of the input image

        // it is not just the locations but each has radius, and direction and they are used to initialize the bayesian tracing
        // also to guide it and to confirm it

        // original image dimensions
        int img_w = inimg_xy.length;
        int img_h = inimg_xy[0].length;

        Overlay ov = new Overlay();

        for (int i = 0; i < queue.length; i++) {

            int x = i2xy[queue[i]][0];
            int y = i2xy[queue[i]][1];
            float rr = 2f;//i2sigma[queue[i]];

            // oval with scale and location
            OvalRoi ovr = new OvalRoi(x-rr+.5, y-rr+.5, 2*rr, 2*rr);
            float colr = 1 - i/(1.5f*(queue.length-1)); // fade out to some low boundary as it reaches the lower weights
            Color col = new Color(colr,colr,0,colr);
            ovr.setStrokeColor(col);
            ovr.setFillColor(col);
            ov.add(ovr);

            // directions
            Line l = new Line(x+.0, y+.0, x+.0 + 2*rr*i2vxy[queue[i]][0], y+.0 + 2*rr*i2vxy[queue[i]][1]);
            l.setStrokeWidth(1.5);
            l.setStrokeColor(col);
            l.setFillColor(col);
            ov.add(l);

        }

        return ov;

    }

    public static void clean(){ // only clean those allocated that take space
        i2zncc = null;
        i2sigma = null;
        i2vxy = null;
    }
}
