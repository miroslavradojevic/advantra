package detection2d;

import ij.IJ;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.gui.PointRoi;

import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by miroslav on 1/6/14.
 * Will associate the peaks of the profiles, parallel threaded implementation
 * expands each branch recursively, M steps, outputs the indexes of the skeleton locations
 * has to satisfy geometry (expansion angle and expansion length)
 */
public class PeakAnalyzer2D extends Thread {

    private int begN, endN;

    // VARIABLES (mainly used as a table)
    public static int[][] 	    i2xy;                   // selected locations
    public static int[][]     	xy2i;                   // need for recursion

    // INPUT: list of extracted peaks
    public static int[][][]     peaks_xy;             	// N x 4(max. threads) x 2   every FG location with 4 selected peaks in XY format

    // PARAMETERS
    public static int M = 2;                            // how much it expands recursively from the center
    public static float minCos = 0.6f;                  // allowed derail

    // OUTPUT: associate the peaks and link follow-up points
    public static int[][][] delin2;                     // N(foreground locs.) x 4(max. threads) x M(follow-up locs) contains index for each location

    public PeakAnalyzer2D(int n0, int n1)
    {
        this.begN = n0;
        this.endN = n1;
    }

    public static void loadTemplate(int[][] _i2xy, int[][] _xy2i, int[][][] _peaks_xy, int _M, float _minCos) {

        i2xy = _i2xy;
        xy2i = _xy2i;
        peaks_xy = _peaks_xy;
        M = _M;
        minCos = _minCos;

        // allocate output -> set to -1
        delin2 = new int[i2xy.length][4][M];
        for (int ii=0; ii<delin2.length; ii++) {
            for (int jj=0; jj<delin2[0].length; jj++) {
                for (int kk=0; kk<delin2[0][0].length; kk++) {
                    delin2[ii][jj][kk] = -1;
                }
            }
        }

    }

    public void run()
    {
        for (int locationIdx = begN; locationIdx < endN; locationIdx++) {

            int atX = i2xy[locationIdx][0];
            int atY = i2xy[locationIdx][1];

            int[][] peaks_at_loc = peaks_xy[locationIdx];

            // check neighbouring peaks
            // access individual peaks at this point
            for (int pp = 0; pp<peaks_at_loc.length; pp++) {  // loop 4 allocated branches

                if (peaks_at_loc[pp][0] != -1) { // if the peak exists (X coord. is != -1)

                    // store the index value of this peak

                    int pkX = peaks_at_loc[pp][0];
                    int pkY = peaks_at_loc[pp][1];

                    int indexValue = xy2i[pkX][pkY];

                    delin2[locationIdx][pp][0] = indexValue; // m=0

                    int curr_index, prev_index, next_index;

                    curr_index = indexValue;
                    prev_index = locationIdx;

                    for (int m=1; m<M; m++) { // follow the rest of the indexes

                        // recursion : prev+curr->next index

                        next_index = getNext(prev_index, curr_index); // next follow-up will be calculated and sorted

                        if (next_index!=-1) { // -1 will be if the next one is not found

//                            next_index = locIndexZXY[nextXYZ[0]][nextXYZ[1]][nextXYZ[2]];
                            // store it in output matrix
                            delin2[locationIdx][pp][m] = next_index;

                        }
                        else { // follow-up does not exist, break looping m (extending further) but continue looping branches
                            break;
                        }

                    }

                }
                else { // it was -1 and the rest are not found
                    break; // stop for() loop - there are no more branches (already aligned)
                }

            }

        }

    }

    private static int getNext(int prev_index, int curr_index){

        // these are stacked as XY - take care on that
        int prevX = i2xy[prev_index][0];    // X
        int prevY = i2xy[prev_index][1];    // Y

        int currX = i2xy[curr_index][0];    // X
        int currY = i2xy[curr_index][1];    // Y

        // check peaks at curr
        int[][] pks4xXY = peaks_xy[curr_index];
        for (int pkIdx = 0; pkIdx<pks4xXY.length; pkIdx++) { // loops them by rank - those with highest weight first, to pick the first one that points outwards

            if (pks4xXY[pkIdx][0] != -1) {
                // there is a peak to check - needs to be pointing outwards

                int nextX = pks4xXY[pkIdx][0];
                int nextY = pks4xXY[pkIdx][1];

                double cosAng =
                        (
                        	(currX-prevX)*(nextX-currX) + (currY-prevY)*(nextY-currY)                               // + (currZ-prevZ)*(nextZ-currZ)
                        )
                        /
                        (
                        	Math.sqrt( Math.pow(currX-prevX, 2) + Math.pow(currY-prevY, 2) ) *             //  + Math.pow(currZ-prevZ, 2)
                            Math.sqrt( Math.pow(nextX-currX, 2) + Math.pow(nextY-currY, 2) )        //  + Math.pow(nextZ-currZ, 2)
                        );

                if (cosAng>minCos) {
                    return xy2i[nextX][nextY]; // it is aligned - add it, find its index and return as output
                }
                else {
                    // if not pointing outwards, continue further down the rank till it reaches -1 or checks all
                }
            }
            else {
                return -1; // no more peaks to search
            }

        }

        return -1; // checked all

    }

    /*
        outputs
     */
    public static Overlay getDelineation(int atX, int atY)
    {
        Overlay ov = new Overlay();
        float R = 2;
        OvalRoi ovalroi = new OvalRoi(atX-(R/2)+.5f, atY-(R/2)+.5f, R, R);
        ov.add(ovalroi);

        // read extracted peaks at this location
        int idx = Masker2D.xy2i[atX][atY];

        /*
            show selected points
         */

        if (idx!=-1) {

            int[][] delin_at_loc = PeakAnalyzer2D.delin2[idx];

            // show locs  (debug)
            IJ.log("\n_____");
            for (int a=0; a<delin_at_loc.length; a++) IJ.log(Arrays.toString(delin_at_loc[a]));
            IJ.log("_____\n");

            for (int b = 0; b<delin_at_loc.length; b++) {           // loop 4 branches

                boolean complete = true;

                for (int m=0; m<M; m++) {

                    if (delin_at_loc[b][m] != -1) {

                        // there is a point to add

                        int pt_idx = delin_at_loc[b][m];

                        int pt_x = i2xy[pt_idx][0];
                        int pt_y = i2xy[pt_idx][1];

                        ovalroi = new OvalRoi(pt_x-(R/2)+.5f, pt_y-(R/2)+.5f, R, R); // add the point to the overlay
                        ov.add(ovalroi);

                    }
                    else {
                        complete = false;
                        break; // stop with the branch
                    }

                }

                // finished along the branch
                if (complete) {

                    IJ.log("branch"+b+". -> complete");
                    // TODO: put the lines in between for the fully delineated branches that reached the end


                }

            }

        }

        /*
            show cloud of peak points
         */

        if (idx != -1) {

            // establish the overlay with the cloud of recursively traced profile peaks
            // trace works M times recursively (that is the full choice)

            /*
                initialize
             */
            int m = 1;
            ArrayList<int[]> next_pts_xy = new ArrayList<int[]>();
            next_pts_xy.clear();
            next_pts_xy.add(i2xy[idx]);     // xy

            while (m<=4) {

                // add those from the next_pts_xy list to the overlay
                for (int kk=0; kk<next_pts_xy.size(); kk++) {
                    int curr_x = next_pts_xy.get(kk)[0];
                    int curr_y = next_pts_xy.get(kk)[1];
                    PointRoi proi = new PointRoi(curr_x+.5f, curr_y+.5f);
                    ov.add(proi);
                }

                // redefine next_pts_xy as follow-up points
                ArrayList<int[]> temp = (ArrayList<int[]>) next_pts_xy.clone();
                next_pts_xy.clear();
                for (int ii=0; ii<temp.size(); ii++) {

                    int read_x = temp.get(ii)[0];
                    int read_y = temp.get(ii)[1];
                    int read_i = xy2i[read_x][read_y];

                    int[][] pks_at_loc = PeakExtractor2D.peaks_xy[read_i];

                    for (int jj = 0; jj<pks_at_loc.length; jj++) {

                        int new_x = pks_at_loc[jj][0];
                        int new_y = pks_at_loc[jj][1];

                        if (new_x != -1) {
                            next_pts_xy.add(new int[]{new_x, new_y});
                        }
                        else {
                            break;
                        }

                    }
                }

                m++;

            } // loop

        }

        return ov;

    }

}
