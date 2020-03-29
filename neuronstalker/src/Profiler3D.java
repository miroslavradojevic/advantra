import ij.ImagePlus;
import ij.ImageStack;
import ij.process.FloatProcessor;

import java.io.*;
import java.util.ArrayList;

/**
 * Created by miroslav on 7-3-15.
 */
public class Profiler3D extends Thread {

    private int begN, endN;

    public static float[][][]           inimg_xyz;                   // input image
    public static Zncc3D                zncc3D_module;               // module class for computing zero-mean-normalized-cross-corr.
    public static int[][]               i2xyz;                       // locations mapping
    public static ArrayList<double[]>   soma_list;                   // list of soma center of masses

    // output
    public static float[]   i2zncc;                     // correlations per indexed foreground location
    public static float[]   i2sigma;                    // cross-section gaussian sigma per foreground location (the one with largest zncc score)
    public static float[][] i2vxyz;                      // local orientation per foreground location (the one with largest zncc score)

    public Profiler3D(int n0, int n1)
    {
        this.begN = n0;
        this.endN = n1;
    }

    public static void loadTemplate(
            Zncc3D _zncc3D,
            int[][] _i2xyz,
            ArrayList<double[]> _soma_centers,  // only for the queue extraction criteria
            float[][][] _inimg_xyz              // the image itself
    )
    {
        zncc3D_module = _zncc3D;
        i2xyz = _i2xyz;
        soma_list = _soma_centers;
        inimg_xyz = _inimg_xyz;

        // allocate outputs
        i2zncc = new float[i2xyz.length];
        i2sigma = new float[i2xyz.length];
        i2vxyz = new float[i2xyz.length][3];
    }

    public void run() { // will be threaded

        float[] vals = new float[zncc3D_module.limR*zncc3D_module.limR*zncc3D_module.limR];

        for (int loc_index = begN; loc_index < endN; loc_index+=4)
            zncc3D_module.extract(loc_index, i2xyz, inimg_xyz, i2zncc, i2sigma, i2vxyz, vals);

    }

    public static int[] createQueue(float boundary, int weighting_type){

        ArrayList<Integer>  sel_idxs    = new ArrayList<Integer>();

        for (int i = 0; i < i2zncc.length; i++)
            if (i2zncc[i]>boundary) {
                sel_idxs.add(i);
            }

        ArrayList<Float>    sel_scrs    = new ArrayList<Float>(sel_idxs.size());

        if (soma_list.size()>=1) {
            // there is 1+ soma's

            float min_dist = Float.POSITIVE_INFINITY;
            float max_dist = Float.NEGATIVE_INFINITY;

            // add the distance to the scores of selected locations
            for (int i = 0; i < sel_idxs.size(); i++) {

                float currx = i2xyz[sel_idxs.get(i)][0];
                float curry = i2xyz[sel_idxs.get(i)][1];
                float currz = i2xyz[sel_idxs.get(i)][2];

                // squared distance towards the first soma (it has to exist)
                float dist = (float) (
                        Math.pow(currx-soma_list.get(0)[0],2) +
                                Math.pow(curry-soma_list.get(0)[1],2) +
                                    Math.pow(currz-soma_list.get(0)[2],2)
                );

                // remaining somas as much as soma_list allows
                for (int j = 1; j < soma_list.size(); j++) {

                    float dist1 = (float) (
                            Math.pow(currx-soma_list.get(j)[0],2) +
                                Math.pow(curry-soma_list.get(j)[1],2) +
                                    Math.pow(currz-soma_list.get(j)[2],2)
                    );
                    if (dist1<dist) dist = dist1;
                }

                sel_scrs.add(i, dist);

                if (dist<min_dist) min_dist = dist;
                if (dist>max_dist) max_dist = dist;

            }

            // loop through all the points once they're calculated append the zncc score to the min/max normalized distances to soma by now
            for (int i = 0; i < sel_scrs.size(); i++) {

                float normalized_dist = -1;

                if (weighting_type==0) { // those closer to the soma are weighted higher
                    normalized_dist = 1 - (sel_scrs.get(i)-min_dist)/(max_dist-min_dist);
                }
                else if (weighting_type==1) { // those that are the most distant from soma are weighted higher
                    normalized_dist = (sel_scrs.get(i)-min_dist)/(max_dist-min_dist);
                }
                else {
                    System.out.println("wrong guidepoint weighting type");
                    System.exit(-1);
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

        int[] queue = Toolbox.descendingFloat(sel_scrs);
        for (int i = 0; i < queue.length; i++) {
            queue[i] = sel_idxs.get(queue[i]);
        }

        return queue;

    }

    public static ArrayList<boolean[][]> getQueueLocationMap(int[] queue){

        // queue[] contains indexes in the foreground list i2xyz[][]
        // to obtain the locations in image it is necessary to map, by x=i2xyz[queue[]][0], y=i2xyz[queue[]][1], z=i2xyz[queue[]][2]
        // will mark the locations in the 3D map with the dimensions of the input image
        // T at locations (x,y,z) where queue element was found, F otherwise
        // used to guide the bayesian tracing, we need to know if the trace crossed anything that was a guidepoint
        // so that the trace is nto thrown away

        // original image dimensions
        int img_w = inimg_xyz.length;
        int img_h = inimg_xyz[0].length;
        int img_l = inimg_xyz[0][0].length;

        ArrayList<boolean[][]> queue_map = new ArrayList<boolean[][]>(img_l);
        for (int i = 0; i < img_l; i++) queue_map.add(i, new boolean[img_w][img_h]);

        for (int i = 0; i < queue.length; i++) {

            int x = i2xyz[queue[i]][0];
            int y = i2xyz[queue[i]][1];
            int z = i2xyz[queue[i]][2];

            queue_map.get(z)[x][y] = true;

        }

        return queue_map;

    }

    public static ImagePlus getScores(){

        int img_w = inimg_xyz.length;
        int img_h = inimg_xyz[0].length;
        int img_l = inimg_xyz[0][0].length;

        float[][][] corr_scrs_array = new float[img_l][img_w][img_h];

        for (int i = 0; i < i2zncc.length; i++) {

            int x = i2xyz[i][0];
            int y = i2xyz[i][1];
            int z = i2xyz[i][2];

//            int idx = y * img_w + x;
            corr_scrs_array[z][x][y] = i2zncc[i];

        }

        ImageStack is_out = new ImageStack(img_w, img_h);
        for (int i = 0; i < corr_scrs_array.length; i++) {
            is_out.addSlice(new FloatProcessor(corr_scrs_array[i]));
        }

        ImagePlus out = new ImagePlus("zncc", is_out);
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

            float x = i2xyz[queue[i]][0];
            float y = i2xyz[queue[i]][1];
            float z = i2xyz[queue[i]][2];

            float vx = i2vxyz[queue[i]][0];
            float vy = i2vxyz[queue[i]][1];
            float vz = i2vxyz[queue[i]][2];

            float rr = i2sigma[queue[i]];

            logWriter.println(String.format(
                    "%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f %-4d",
                    cnt, 2, x, y, z, rr, -1));

            // directions
            logWriter.println(String.format(
                    "%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f %-4d",
                    cnt+1, 2, x, y, z, 0.1f, -1));

            logWriter.println(String.format(
                    "%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f %-4d",
                    cnt+2, 2, x+2*rr*vx, y+2*rr*vy, z+2*rr*vz, 0.1f, cnt+1));

            logWriter.println(String.format(
                    "%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f %-4d",
                    cnt+3, 2, x-2*rr*vx, y-2*rr*vy, z-2*rr*vz, 0.1f, cnt+1));

            cnt+=4;

        }

        logWriter.close();

    }

    public static void clean(){
        i2zncc = null;
        i2sigma = null;
        i2vxyz = null;
    }
}