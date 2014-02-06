package detection2d;

import conn.Find_Connected_Regions;
import ij.ImagePlus;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.process.ByteProcessor;

/**
 * Created by miroslav on 1/9/14.
 * takes the robustly estimated peaks
 * loops the locations' peak skeletons extracted using PeakAnalyzer2D (the usage of that class draws the usage of the other Thread classes as well)
 * detector checks if the robust, consistent peak recursion evolved till the proposed end, along at least three branches
 * no special unsupervised/supervised detection, just structure fit used to make the decision
 * is supposed to give the rough sketch of the areas that can be bifurcation spots, the rest is to extract the features and do thorough decision making
 * important that simple detector does not miss some bifurcations
 */

public class SimpleDetector2D extends Thread {

    private int begN, endN;

    // VARIABLES
    public static int[][] 	    i2xy;                   // selected locations
    public static int           W, H;                   // range for x and y

    // INPUT: list of extracted peaks' indexes
    public static int[][][]     delin2;             	// N x 4(max. threads) x M   every FG location with 4 selected peaks in XY format

    // PARAMETERS
    // no parameters here - just check the structure

    // OUTPUT
    public static byte[]      score2;

    public SimpleDetector2D(int n0, int n1)
    {
        this.begN = n0;
        this.endN = n1;
    }

    public static void loadTemplate(int _W, int _H, int[][] _i2xy, int[][][] _delin2) {

        W       = _W;
        H       = _H;
        i2xy    = _i2xy;
        delin2  = _delin2;

        /*
            allocate output
         */
        score2 = new byte[W*H];

    }

    public void run() {

        // loop through all the locations and decide (if at least three complete branches exist)
        for (int locationIdx = begN; locationIdx < endN; locationIdx++) {

            int atX = i2xy[locationIdx][0];
            int atY = i2xy[locationIdx][1];

            if (simpleDecision(locationIdx)) {
                score2[atY*W+atX] = (byte) 255; // row by row align
            }

        }

    }

    private boolean simpleDecision(int idx){

        int cnt = 0;

        for (int i=0; i<delin2[idx].length; i++) {

            int last_idx = delin2[idx][i].length - 1;

            if (delin2[idx][i][last_idx] != -1) {
                // means that it reached the end, converged, branch is complete
                cnt++;
            }

        }

        if (cnt >= 3) {
            return true;
        }
        else {
            return false;
        }

    }

    public static void drawDetections(){

        ByteProcessor bp = new ByteProcessor(W, H, score2);
        ImagePlus ip = new ImagePlus("DET", bp);
        ip.show();

        Overlay ov = new Overlay();

        // take detections (binary image), find connected regions, and extract out the overlay with the detections

        //ByteProcessor score = new ByteProcessor(W, H);
        //for (int ii=0; ii<W*H; ii++) if (score2[ii]) score.set(ii, 255);
        Find_Connected_Regions conn_reg = new Find_Connected_Regions(new ImagePlus("", bp), true);
        conn_reg.run("");
//        detectionOverlay = formPointOverlay(conn_reg.getConnectedRegions(), MIN_SIZE);
//        for (int ii=0; ii<imp.getWidth()*imp.getHeight(); ii++) {
//            if (score.get(ii)>=255) detectionOverlay.add(new PointRoi(ii%imp.getWidth()+.5, ii/imp.getWidth()+.5));
//        }

//        return ov;

    }

}