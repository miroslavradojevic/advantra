package detection2d;

import aux.Stat;
import com.sun.org.apache.bcel.internal.generic.FLOAD;
import detection.Interpolator;
import ij.IJ;
import ij.gui.Line;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.gui.PointRoi;

import java.awt.*;
import java.io.*;
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

    // INPUT:
    public static int[][][]     peaks_xy;             	// list of extracted peaks: N x 4(max. threads) x 2   every FG location with 4 selected peaks in XY format
	public static float[][]		inimg_xy;				// input image (necessary for feature extraction)
	public static byte[][]		backg_xy;				// background estimation (for feature extraction)

    // PARAMETERS
    public static int M = 2;                            // how much it expands recursively from the center
    public static float minCos = 0.6f;                  // allowed derail
    public static float scatterDist = 5;                // allowed scatter dist, upper limit, half of the neighbouring peaks should be within


    // OUTPUT: associate the peaks and link follow-up points
    public static int[][][] delin2;                     // N(foreground locs.) x 4(max. threads) x M(follow-up locs) contains index for each location
	public static float[][]	feat2;						// N x 7 (6+1) features

    public PeakAnalyzer2D(int n0, int n1)
    {
        this.begN = n0;
        this.endN = n1;
    }

    public static void loadTemplate(int[][] _i2xy, int[][] _xy2i, int[][][] _peaks_xy, float[][] _inimg_xy, byte[][] _backgr_xy, int _M, float _minCos, float _scatterDist) {

        i2xy = _i2xy;
        xy2i = _xy2i;
        peaks_xy = _peaks_xy;
		inimg_xy = _inimg_xy;
		backg_xy = _backgr_xy;
        M = _M;
        minCos = _minCos;
        scatterDist = _scatterDist;

        // allocate output -> set to -1
        delin2 = new int[i2xy.length][4][M];
        for (int ii=0; ii<delin2.length; ii++) {
            for (int jj=0; jj<delin2[0].length; jj++) {
                for (int kk=0; kk<delin2[0][0].length; kk++) {
                    delin2[ii][jj][kk] = -1;
                }
            }
        }

		feat2 = new float[i2xy.length][7]; // there will be 2*3+1 feature
		// set to zero - all

    }

    public void run()
    {
        for (int locationIdx = begN; locationIdx < endN; locationIdx++) {

            //int atX = i2xy[locationIdx][0];
            //int atY = i2xy[locationIdx][1];

            int[][] peaks_at_loc = peaks_xy[locationIdx];

            // check neighbouring peaks
            // access individual peaks at this point
            for (int pp = 0; pp<peaks_at_loc.length; pp++) {  // loop 4 allocated branches

                if (peaks_at_loc[pp][0] != -1) { // if the peak exists (X coord. is != -1)

                    // store the index value of this peak

                    int pkX = peaks_at_loc[pp][0];
                    int pkY = peaks_at_loc[pp][1];

                    int indexValue = xy2i[pkX][pkY];

                    if (isRobust(indexValue, locationIdx, scatterDist)) {    // locationIdx is mother index in this case
                        delin2[locationIdx][pp][0] = indexValue; // STORE IT! m=0
                    }

                    int curr_index, prev_index, next_index;

                    curr_index = indexValue;
                    prev_index = locationIdx;

                    for (int m=1; m<M; m++) { // follow the rest of the indexes

                        // recursion : prev+curr->next index

                        next_index = getNext(prev_index, curr_index); // next follow-up will be calculated and sorted

                        if (next_index!=-1) { // -1 will be if the next one is not found

                            if (isRobust(next_index, curr_index, scatterDist)) {
                                delin2[locationIdx][pp][m] = next_index;     // store it in output matrix
                            }
                            else {
                                break; // stop going further if it is not robust
                            }

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

			// delin2[locationIdx] is fully formed

            /*
            calc_feats
             */
			// extract features at this spot theta = median(line_beg_end) - backgr(line_end)
			// take top 3 branches (highest min.)  out of calc_feats
			float[][] calc_feats = new float[delin2[locationIdx].length][2]; // take just min/max, 4x2, 4x(min,max)

			for (int ii=0; ii<calc_feats.length; ii++) {

				float min_val = Float.MAX_VALUE;
				float max_val = Float.MIN_VALUE;
				float curr_val;

				for (int jj=0; jj<delin2[locationIdx][ii].length; jj++) { // along each branch peak

					// calculate curr_val
					if (delin2[locationIdx][ii][jj] != -1) {

						int prev_i, prev_x, prev_y;
						int curr_i, curr_x, curr_y;

						if (jj==0) { // jj marks current index
							prev_i = locationIdx;
						}
						else {
							prev_i = jj-1;
						}

						prev_x = i2xy[prev_i][0];
						prev_y = i2xy[prev_i][1];
						curr_i = delin2[locationIdx][ii][jj];
						curr_x = i2xy[curr_i][0];
						curr_y = i2xy[curr_i][1];

						float med_along_line = medianAlongLine(prev_x, prev_y, curr_x, curr_y, inimg_xy);
						int back_at_loc = backg_xy[curr_x][curr_y] & 0xff;
						curr_val = med_along_line - back_at_loc;
						curr_val = (curr_val>0)? curr_val : 0 ;

					}
					else {
						curr_val = 0;
					}

					if (curr_val>max_val) {
						max_val = curr_val;
					}
					if (curr_val<min_val) {
						min_val = curr_val;
					}

				}

				// store min/max val for this branch
				calc_feats[ii][0] = min_val;
				calc_feats[ii][1] = max_val;

			}

            /*
            calc_feat, central
             */
            int center_x = i2xy[locationIdx][0];
            int center_y = i2xy[locationIdx][1];
            float calc_feat = medianAtPoint(center_x, center_y, inimg_xy) - (backg_xy[center_x][center_y] & 0xff);
            calc_feat = (calc_feat>0)? calc_feat : 0 ;

            feat2[locationIdx][6] = calc_feat;  // central one at index 6

			// store top 3 rows of calc_feats in feat2[locationIdx]  (criteria: min_value, first index)
			// loop 3 times and each time take the currently highest min
            boolean[] checked = new boolean[calc_feats.length];

            for (int bb=0; bb<3; bb++) { // loop three times and pick the one with highest available min

                float max_to_beat = Float.NEGATIVE_INFINITY;
                int max_idx = -1;

                for (int aa=0; aa<calc_feats.length; aa++) {

                    if (calc_feats[aa][0]>max_to_beat && !checked[aa]) {
                        max_to_beat = calc_feats[aa][0];
                        max_idx = aa;
                    }

                }

                checked[max_idx] = true;
                feat2[locationIdx][2*bb]    = calc_feats[max_idx][0];   // bb = 0, 1, 2
                feat2[locationIdx][2*bb+1]  = calc_feats[max_idx][1];

            }

        }

    }

 /*
 * the "spatial consistency" is checked to account for the peak robustness (it's neighbours have to point to the spatially close location)
 * to avoid having some outlier peaks included but have groups of peaks pointing together, agreeing on the same information
 * isRobust() gives out the filtered version of delin2 where only the stable ones are taken
 */

    private static boolean isRobust(int test_idx, int mother_idx, float scatter_th) {

        // check for robustness
        // check if the mother-peak's neighbours agree with the follow-up
        // take the 4 or 8 neighbourhood of the mother peak and take
        // 4 or 8 euclidean-wise closest peaks to the one being checked
        // if the median of them is close enough (if they're not too scattered)
        // then consider the peak follow-up spatially robust

        int test_x = i2xy[test_idx][0];
        int test_y = i2xy[test_idx][1];

        // check how many agree-points (scatter) there can be

        // 4 neighbours
        int[][] dx_dy = new int[][]{
                {-1, 0},
                { 0,-1},
                { 1, 0},
                { 0, 1}
        };
        float[] scatter_dists = new float[dx_dy.length];

        int mother_x = i2xy[mother_idx][0];
        int mother_y = i2xy[mother_idx][1];

        for (int i=0; i<dx_dy.length; i++) {

            int neigbr_x = mother_x + dx_dy[i][0];
            int neigbr_y = mother_y + dx_dy[i][1];
            int neigbr_i = xy2i[neigbr_x][neigbr_y];

            scatter_dists[i] = xy2i.length; // max dist

            if (neigbr_i != -1) {

                int[][] get_peaks_nbr = peaks_xy[neigbr_i]; // peak signature of the neighbour

                // pick the closest one
                for (int k=0; k<get_peaks_nbr.length; k++) {

                    if (get_peaks_nbr[k][0] != -1) {

                        int scatter_x = get_peaks_nbr[k][0];
                        int scatter_y = get_peaks_nbr[k][1];

                        float d = dist(test_x, test_y, scatter_x, scatter_y);

                        if (d<scatter_dists[i]) {
                            scatter_dists[i] = d;
                        }

                    }
                    else {
                        break;
                    }

                }

            }

        }

        if (Stat.median(scatter_dists)<=scatter_th) {
            return true;
        }
        else {
            return false;
        }

    }

    private static float dist(int a_x, int a_y, int b_x, int b_y){
        return (float) Math.sqrt( Math.pow(a_x-b_x, 2) + Math.pow(a_y-b_y, 2) );
    }

    private static int getNext(int prev_index, int curr_index) {

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

    private static float medianAtPoint(int x, int y, float[][] _inimg_xy) {

        float[] nhood = new float[9];

        if (x>0 && x<_inimg_xy.length-1 && y>0 && y<_inimg_xy[0].length-1) {

            int cnt = 0;
            for (int dx=-1; dx<=1; dx++) {
                for (int dy=-1; dy<=1; dy++) {

                    nhood[cnt] = _inimg_xy[x+dx][y+dy];
                    cnt++;

                }
            }

        }

        return Stat.median(nhood);

    }

	private static float medianAlongLine(int x1, int y1, int x2, int y2, float[][] inimg_xy) {

		float increment_length = .7f;

		float dist = (float) Math.sqrt(Math.pow(x2-x1, 2) + Math.pow(y2-y1, 2)); //  + Math.pow(z2lay-z1lay, 2)

		int elementsInLine = (int) (dist / increment_length);  // how many increment can safely fit between
		float[] valuesAlongLine = new float[elementsInLine];

		float dx = (x2 - x1) / dist;
		float dy = (y2 - y1) / dist;
		// [dx, dy] is unit vector

		dx *= increment_length;
		dy *= increment_length;

		for (int cc = 0; cc<elementsInLine; cc++) {

			float atX = x1      + cc * dx;
			float atY = y1      + cc * dy;
//			float atZ = z1lay   + cc * dz;

			if (atX<0 || atX>inimg_xy.length-1 || atY<0 || atY>inimg_xy[0].length-1) {
				System.out.println(String.format("\n%5d.(%4d) element (X, Y) was set wrong: (%5.2f, %5.2f) when sampling values between (%5.2f, %5.2f)->(%5.2f, %5.2f)\n", cc, elementsInLine, atX, atY, x1, y1, x2, y2));
			}

			valuesAlongLine[cc] = Interpolator.interpolateAt(atX, atY, inimg_xy);

		}

		return Stat.median(valuesAlongLine);

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
            //IJ.log("\n_____");
            //for (int a=0; a<delin_at_loc.length; a++) IJ.log(Arrays.toString(delin_at_loc[a]));
            //IJ.log("_____\n");

            for (int b = 0; b<delin_at_loc.length; b++) {           // loop 4 branches

                boolean complete = true;

                for (int m=0; m<M; m++) {

                    if (delin_at_loc[b][m] != -1) {

                        // there is a point to add

                        int pt_idx = delin_at_loc[b][m];

                        int pt_x = i2xy[pt_idx][0];
                        int pt_y = i2xy[pt_idx][1];

                        ovalroi = new OvalRoi(pt_x-(R/2)+.5f, pt_y-(R/2)+.5f, R, R); // add the point to the overlay
                        ovalroi.setStrokeColor(java.awt.Color.RED);
                        ovalroi.setStrokeWidth(2);
                        ov.add(ovalroi);

                    }
                    else {
                        complete = false;
                        break; // stop with the branch
                    }

                }

                // finished along the branch
                if (complete) {    // add lines along complete lines

                    for (int m=0; m<M; m++) {
                        int curr_i = delin_at_loc[b][m];
                        int curr_x = i2xy[curr_i][0];
                        int curr_y = i2xy[curr_i][1];

                        int prev_i, prev_x, prev_y;

                        if (m==0) {
                            prev_x = atX;
                            prev_y = atY;
                        }
                        else{
                            prev_i = delin_at_loc[b][m-1];
                            prev_x = i2xy[prev_i][0];
                            prev_y = i2xy[prev_i][1];
                        }

                        Line l = new Line(prev_x, prev_y, curr_x, curr_y);
                        l.setStrokeColor(Color.RED);
                        l.setStrokeWidth(2);
                        ov.add(l);

                    }

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

    public static void exportFeatsCsv(String file_path) {

        IJ.log("exporting extracted features...");

        PrintWriter logWriter = null; //initialize writer

        try {
            logWriter = new PrintWriter(file_path);
            logWriter.print("");
            logWriter.close();
        } catch (FileNotFoundException ex) {}   // empty the file before logging...

        try {                                   // initialize detection log file
            logWriter = new PrintWriter(new BufferedWriter(new FileWriter(file_path, true)));
        } catch (IOException e) {}

        //main loop
        for (int ii=0; ii<feat2.length; ii++) {
            for (int jj=0; jj<feat2[ii].length; jj++) {
                logWriter.print(String.format("%6.2f", feat2[ii][jj]));
                if (jj<feat2[ii].length-1) {
                    logWriter.print(",\t");
                }
                else {
                    logWriter.print("\n");
                }
            }
        }

        logWriter.close(); // close log
        IJ.log("Saved in "+file_path);

    }

    public static void exportFeatsLegend(String file_path) {

        IJ.log("exporting feature legend...");

        PrintWriter logWriter = null; //initialize writer

        try {
            logWriter = new PrintWriter(file_path);
            logWriter.print("");
            logWriter.close();
        } catch (FileNotFoundException ex) {}   // empty the file before logging...

        try {                                   // initialize detection log file
            logWriter = new PrintWriter(new BufferedWriter(new FileWriter(file_path, true)));
        } catch (IOException e) {}


        logWriter.println("strength: branch 0 > branch 1 > branch 2");
        logWriter.println("score_at_point = median_along_line_ending at_point(or in the point neighbourhood) - background_estimate_at_point");

        logWriter.println("feature 0: \tbranch 0 min score");
        logWriter.println("feature 1: \tbranch 0 max score");

        logWriter.println("feature 2: \tbranch 1 min score");
        logWriter.println("feature 3: \tbranch 1 max score");

        logWriter.println("feature 4: \tbranch 2 min score");
        logWriter.println("feature 5: \tbranch 2 max score");

        logWriter.println("feature 6: \tcenter       score");

        logWriter.close(); // close log
        IJ.log("Saved in "+file_path);

    }

}
