package detection2d;

import conn.Find_Connected_Regions;
import ij.ImagePlus;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;

import java.awt.*;
import java.util.ArrayList;

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
    public static int[][]       xy2i;                   // legend, map
    public static int           W, H;                   // range for x and y

    // INPUT: list of extracted peaks' indexes
    public static int[][][]     delin2;             	// N x 4(max. threads) x M, every FG location with 4 selected peaks in XY format
    public static float[][]     lhood2;                 // N x 5, PeakAnalyzer2D fuzzy logic output, 5 outputs with likelihood

    static float[]      Gx;                             // derivative filter weights (Sobel)
    static float[]      Gy;                             // derivative filter weights (Sobel)
    static int[][]      Gx_dXdY;                        // filter offsets
    static int[][]      Gy_dXdY;                        // filter offsets

    // PARAMETERS
    static int min_size = 1;                            // when extracting connected components
    static float min_lhood = .6f;                       // when discarding those that are not certain enough
    static int nbhood_size = 0;

    // OUTPUT
    //public static byte[]      score2;
    //public static float[]     lhoods;

    public static float[]   jun_lhood; // W x H matrix with junction (bif. and cross.) fuzzy likelihoods
    public static float[]   end_lhood; // W x H matrix with end-point fuzzy likelihoods per pixel
    public static float[]   model_grad;// intensity of the model change (expressed in terms of change of # patches)
    // fuzzy scores turn out to be instable - therefore the regions with consistent detections need to be extracted
    // consistent regions meaning those where the model stays the same and the score is high

    public SimpleDetector2D(int n0, int n1)
    {
        this.begN = n0;
        this.endN = n1;
    }

    public static void loadTemplate(int _W, int _H, int[][] _i2xy, int[][] _xy2i, int[][][] _delin2, float[][] _lhood2) {

        W       = _W;
        H       = _H;
        i2xy    = _i2xy;
        xy2i    = _xy2i;
        delin2  = _delin2;
        lhood2  = _lhood2; // 0-NON, 1-END, 2-BDY, 3,4-JUN

        // filter
        Gx = new float[6];
        Gy = new float[6];
        Gx_dXdY = new int[6][2]; // X (column), Y (row)
        Gy_dXdY = new int[6][2]; // X (column), Y (row)

        Gy[0] = 1;     Gy_dXdY[0][0] = -1; Gy_dXdY[0][1] = -1;
        Gy[1] = 2;     Gy_dXdY[1][0] =  0; Gy_dXdY[1][1] = -1;
        Gy[2] = 1;     Gy_dXdY[2][0] =  1; Gy_dXdY[2][1] = -1;
        Gy[3] = -1;    Gy_dXdY[3][0] = -1; Gy_dXdY[3][1] =  1;
        Gy[4] = -2;    Gy_dXdY[4][0] =  0; Gy_dXdY[4][1] =  1;
        Gy[5] = -1;    Gy_dXdY[5][0] =  1; Gy_dXdY[5][1] =  1;

        Gx[0] = 1;     Gx_dXdY[0][0] =  1; Gx_dXdY[0][1] = -1;
        Gx[1] = 2;     Gx_dXdY[1][0] =  1; Gx_dXdY[1][1] =  0;
        Gx[2] = 1;     Gx_dXdY[2][0] =  1; Gx_dXdY[2][1] =  1;
        Gx[3] = -1;    Gx_dXdY[3][0] = -1; Gx_dXdY[3][1] = -1;
        Gx[4] = -2;    Gx_dXdY[4][0] = -1; Gx_dXdY[4][1] =  0;
        Gx[5] = -1;    Gx_dXdY[5][0] = -1; Gx_dXdY[5][1] =  1;

        /*
            allocate output
         */
        jun_lhood   = new float[W*H];// initialize with zero l'hoods
        end_lhood   = new float[W*H];
        model_grad  = new float[W*H];
    }

    public void run() {

        // loop through all the locations and generate images:
        // likelihood for being bifurcation (fuzzy logic output)
        // likelihood for being endpoint    (fuzzy logic too)
        // delineation model change intensity (gradient of the model description - # patches)
        for (int locationIdx = begN; locationIdx < endN; locationIdx++) {

            int atX = i2xy[locationIdx][0];
            int atY = i2xy[locationIdx][1];

            if (atX<1 || atX>W-2 || atY<1 || atY>H-2) continue;

            // junction l'hood - read it directly from lhood2
            jun_lhood[atY*W+atX] = Math.max(lhood2[locationIdx][3], lhood2[locationIdx][4]);
            end_lhood[atY*W+atX] = lhood2[locationIdx][1];

            float GradX = 0;
            float GradY = 0;
            // model gradient calculate
            for (int kk = 0; kk<6; kk++) {
                GradX += Gx[kk] * modelValue(atX+Gx_dXdY[kk][0], atY+Gx_dXdY[kk][1]);
                GradY += Gy[kk] * modelValue(atX+Gy_dXdY[kk][0], atY+Gy_dXdY[kk][1]);
            }

            model_grad[atY*W+atX] = (float) Math.sqrt(GradX*GradX+GradY*GradY);

            // check the neighbours to make majority voting decision
            // take the best decision from the local neighbourhood
//            int     curr_label = -1; // 0-4
//            float   curr_lhood = -1;

//            for (int xx=atX-nbhood_size; xx<=atX+nbhood_size; xx++) {
//                for (int yy=atY-nbhood_size; yy<=atY+nbhood_size; yy++) {
//                    if (xy2i[xx][yy]!=-1) { // point is in the foreground
//
//                        float[] nbr_lhood_list = lhood2[xy2i[xx][yy]];
//
//                        for (int kk=0; kk<nbr_lhood_list.length; kk++) {
//
//                            if (nbr_lhood_list[kk]>curr_lhood) {
//                                curr_lhood = nbr_lhood_list[kk];
//                                curr_label = kk;
//                            }
//
//                        }
//                    }
//                }
//            }


//            if (simpleDecision(locationIdx)) {
//                score2[atY*W+atX] = (byte) 255; // row by row align
//            }

        }

    }

    private float modelValue(int atX, int atY)
    {

        int idx = xy2i[atX][atY];

        if (idx!=-1) {

            float score = 0;
            // read the model from delin2
            // count the # patches here
            for (int ii=0; ii<delin2[idx].length; ii++) {
                for (int jj=0; jj<delin2[idx][ii].length; jj++) {
                    if (delin2[idx][ii][jj]!=-1) {
                        score += .5f;
                    }
                }
            }

            return score;

        }
        else {
            // it is in the background
            return 0;
        }

    }

//    private boolean simpleDecision(int idx){
//
//        int cnt = 0;
//
//        for (int i=0; i<delin2[idx].length; i++) {
//
//            int last_idx = delin2[idx][i].length - 1;
//
//            if (delin2[idx][i][last_idx] != -1) {
//                // means that it reached the end, converged, branch is complete
//                cnt++;
//            }
//
//        }
//
//        if (cnt >= 3) {
//            return true;
//        }
//        else {
//            return false;
//        }
//
//    }

    public static Overlay drawDetections(){

        // eliminate those lower than min_lhood
        byte[] out_labels = new byte[W*H];

//        for (int ii=0; ii<W*H; ii++) {
//
//            if (lhoods[ii]>min_lhood) {
//                if (score2[ii]==(byte)1) out_labels[ii] = (byte)127;
//                if (score2[ii]==(byte)3 || score2[ii]==(byte)4) out_labels[ii] = (byte)255;
//            }
//
//        }

        ByteProcessor bp = new ByteProcessor(W, H, out_labels);
        ImagePlus ip = new ImagePlus("DET", bp);
        ip.show();

        // take detections (binary image), find connected regions, and extract out the overlay with the detections
        //ByteProcessor score = new ByteProcessor(W, H);
        //for (int ii=0; ii<W*H; ii++) if (score2[ii]) score.set(ii, 255);
        Find_Connected_Regions conn_reg = new Find_Connected_Regions(new ImagePlus("", bp), true);
        conn_reg.run("");
        conn_reg.showLabels().show();

        Overlay ov = formPointOverlay(conn_reg.getConnectedRegions(), min_size);

        ip.setOverlay(ov);

        return ov;

    }

    private static Overlay formPointOverlay(ArrayList<ArrayList<int[]>> regs, int minSize)
    {

        Overlay detections = new Overlay();

        for (int i=0; i<regs.size(); i++) {
            if (regs.get(i).size()>minSize) {

                float Cx=0, Cy=0, R= (float) Math.sqrt((float)regs.get(i).size()/Math.PI);
                R = (R<1)? 1 : R ;

                for (int aa=0; aa<regs.get(i).size(); aa++) {
                    Cx += regs.get(i).get(aa)[1];
                    Cy += regs.get(i).get(aa)[0];
                }
                Cx /= regs.get(i).size();
                Cy /= regs.get(i).size();

                OvalRoi ovroi = new OvalRoi(Cx-R+.5, Cy-R+.5, 2*R, 2*R);
                ovroi.setStrokeWidth(2);


                int firstY = regs.get(i).get(0)[0];
                int firstX = regs.get(i).get(0)[1];
//                if (score2[firstY*W+firstX]==(byte)1) ovroi.setStrokeColor(Color.YELLOW);
//                else ovroi.setStrokeColor(Color.RED);
                detections.add(ovroi);

            }
        }

        return detections;

    }

    public static ImagePlus showJunctionLhood()
    {
        FloatProcessor fp = new FloatProcessor(W, H, jun_lhood);
        ImagePlus im_out = new ImagePlus("JunctionLikelihood", fp);
        return im_out;
    }

    public static ImagePlus showEndLhood()
    {
        FloatProcessor fp = new FloatProcessor(W, H, end_lhood);
        ImagePlus im_out = new ImagePlus("EndLikelihood", fp);
        return im_out;
    }

    public static ImagePlus showModelGradient()
    {
        FloatProcessor fp = new FloatProcessor(W, H, model_grad);
        ImagePlus im_out = new ImagePlus("ModelGradient", fp);
        return im_out;
    }

}