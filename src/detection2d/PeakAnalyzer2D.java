package detection2d;

import aux.Stat;
import detection.Interpolator;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.*;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.awt.*;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by miroslav on 1/6/14.
 *
 * TASK 1:
 * Will associate the peaks of the profiles, parallel threaded implementation
 * expands each branch recursively, M steps, outputs the indexes of the skeleton locations
 * has to satisfy geometry (expansion angle and expansion length) (delin2)
 * TASK 2:
 * Calculates and exports the features at every location after frame delineation (feat2)
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

    private static float samplingStep = 0.9f;           // when sampling image values to extract features

    // OUTPUT
    // associate the peaks and link follow-up points
    public static int[][][] delin2;                     // N(foreground locs.) x 4(max. threads) x M(follow-up locs) contains index for each location
	// extract features
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

            int[][] peaks_at_loc = peaks_xy[locationIdx];

            // access individual peaks at this point
            for (int pp = 0; pp<peaks_at_loc.length; pp++) {  // loop 4 allocated branches

                if (peaks_at_loc[pp][0] != -1) { // if the peak exists (X coord. is != -1)

                    int pkX = peaks_at_loc[pp][0]; // store the index value of this peak
                    int pkY = peaks_at_loc[pp][1];
                    int indexValue = xy2i[pkX][pkY];
                    delin2[locationIdx][pp][0] = indexValue; // STORE IT! m=0, don't check "robustness" here, as later on

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
                                break; // stop going further with delineating if it is not robust
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

            } // delin2[locationIdx] is fully formed

            /*	calc_feats	*/
			// extract features at this spot, feat. value = median(line_beg_end) - backgr(line_end)
			// take just min/max along one thread, 4xM -> 4x2(min,max)
			float[][] calc_feats = takeMinMax(delin2[locationIdx], locationIdx, false);

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

			/*	calc_feat, central	*/
			int center_x = i2xy[locationIdx][0];
			int center_y = i2xy[locationIdx][1];
			float calc_feat = medianAtPoint(center_x, center_y, inimg_xy) - (backg_xy[center_x][center_y] & 0xff);
			calc_feat = (calc_feat>0)? calc_feat : 0 ;

			feat2[locationIdx][6] = calc_feat;  // central one at index 6

		}

    }

	private static float[][] takeMinMax(int[][] delineation_idxs, int root_idx, boolean verbose) {

        //if (verbose) { IJ.log(""); } // just in case

		float[][] calc_feats = new float[delineation_idxs.length][2];

		for (int ii=0; ii<calc_feats.length; ii++) {

			float min_val = Float.POSITIVE_INFINITY;
			float max_val = Float.NEGATIVE_INFINITY;
			float curr_val;

			for (int jj=0; jj<delineation_idxs[ii].length; jj++) { // along each branch peak

				if (delineation_idxs[ii][jj] != -1) { // calculate curr_val

					int prev_i, prev_x, prev_y;
					int curr_i, curr_x, curr_y;

					prev_i = (jj==0)? root_idx : delineation_idxs[ii][jj-1];
                    if (prev_i == -1) {
                        System.out.println("wrong:");
                        for (int kk=0; kk<delineation_idxs.length; kk++) {
                            System.out.println(Arrays.toString(delineation_idxs[kk])+"|"+Arrays.toString(delin2[root_idx][kk])+"|" +Arrays.toString(peaks_xy[root_idx][kk]));
                        }
                    }
					prev_x = i2xy[prev_i][0];
					prev_y = i2xy[prev_i][1];
					curr_i = delineation_idxs[ii][jj];
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

			calc_feats[ii][0] = min_val; // store min/max val for this branch
			calc_feats[ii][1] = max_val;

		}

		return calc_feats;

	}

	// debug
	public static void print(int atX, int atY) {

		int atLoc = xy2i[atX][atY];

        IJ.log(String.format("/**** LOC (%5d, %5d) [%10d] ****/", atX, atY, atLoc));

		if (atLoc != -1) {

			IJ.log("DELINEATION INDEXES:");
            for (int ii=0; ii<delin2[atLoc].length; ii++) {
				IJ.log("-> "+Arrays.toString(delin2[atLoc][ii]));
			}

			float[][] test = takeMinMax(delin2[atLoc], atLoc, true);
			IJ.log("CALCULATED FEATURES:");
			for (int ii=0; ii<test.length; ii++) {
				IJ.log("MIN/MAX-> "+Arrays.toString(test[ii]));
			}

            IJ.log("FEATURE VECTOR(TH11, TH12, TH21, TH22, TH31, TH32, C):");
            IJ.log(Arrays.toString(feat2[atLoc]));

        }
		else {
			IJ.log("background point, no data!");
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

	// it is public only because Sphere2D uses it for weighting peaks while expanding
//	public static float medianAlongLine(int x1, int y1, int x2, int y2, float[][] inimg_xy) {
//
//		float increment_length = .7f;
//
//		float dist = (float) Math.sqrt(Math.pow(x2-x1, 2) + Math.pow(y2-y1, 2)); //  + Math.pow(z2lay-z1lay, 2)
//
//		int elementsInLine = (int) (dist / increment_length);  // how many increment can safely fit between
//		float[] valuesAlongLine = new float[elementsInLine];
//
//		float dx = (x2 - x1) / dist;
//		float dy = (y2 - y1) / dist;
//		// [dx, dy] is unit vector
//
//		dx *= increment_length;
//		dy *= increment_length;
//
//		for (int cc = 0; cc<elementsInLine; cc++) {
//
//			float atX = x1      + cc * dx;
//			float atY = y1      + cc * dy;
////			float atZ = z1lay   + cc * dz;
//
//			if (atX<0 || atX>inimg_xy.length-1 || atY<0 || atY>inimg_xy[0].length-1) {
//				System.out.println(String.format("\n%5d.(%4d) element (X, Y) was set wrong: (%5.2f, %5.2f) when sampling values between (%5.2f, %5.2f)->(%5.2f, %5.2f)\n", cc, elementsInLine, atX, atY, x1, y1, x2, y2));
//			}
//
//			valuesAlongLine[cc] = Interpolator.interpolateAt(atX, atY, inimg_xy);
//
//		}
//
//		return Stat.median(valuesAlongLine);
//
//	}

    /*
        outputs
     */
    public static Overlay getDelineationOverlay(int atX, int atY)
    {

        // return the delineated local structure Overlay for visualization
        Overlay ov = new Overlay();
        float R = 1;
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

            for (int b = 0; b<delin_at_loc.length; b++) {           // loop 4 branches, b index defines the strength

                boolean complete = true;

                for (int m=0; m<M; m++) {

                    if (delin_at_loc[b][m] != -1) {

                        // there is a point to add

                        int pt_idx = delin_at_loc[b][m];

                        int pt_x = i2xy[pt_idx][0];
                        int pt_y = i2xy[pt_idx][1];

                        ovalroi = new OvalRoi(pt_x-(R/2)+.5f, pt_y-(R/2)+.5f, R, R); // add the point to the overlay
                        ovalroi.setStrokeColor(java.awt.Color.RED);
                        ovalroi.setStrokeWidth(1);
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

                        // add the points for the line connecting prev_xy and curr_xy
                        ArrayList<OvalRoi> line_pts = localLineLocs(prev_x, prev_y, curr_x, curr_y);
                        for (int aa = 0; aa < line_pts.size(); aa++) {
                            ov.add(line_pts.get(aa));
                        }
//                      Line l = new Line(prev_x+.5f, prev_y+.5f, curr_x+.5f, curr_y+.5f);

                        ArrayList<PointRoi> pts = localPatchLocs(prev_x, prev_y, curr_x, curr_y);
						for (int aa=0; aa<pts.size(); aa++) {
							ov.add(pts.get(aa));
						}

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

    public static ImageStack getDelineationPatches(int atX, int atY)
	{

        int patch_size = 256;

        // create new ImageStack with every layer corresponding to one patch
        // involved in modelling the local structure
        ImageStack isOut = new ImageStack(patch_size,patch_size);

        int idx = Masker2D.xy2i[atX][atY]; // read extracted peaks at this location

        if (idx!=-1) {

            int[][] delin_at_loc = PeakAnalyzer2D.delin2[idx];

            for (int b = 0; b<delin_at_loc.length; b++) {           // loop 4 branches, b index defines the strength

                if (delin_at_loc[b][M-1] != -1) { // if the last one is there - it is complete

                    for (int m = 0; m<M; m++) {

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

                        float[] vals = localPatchVals(prev_x, prev_y, curr_x, curr_y, inimg_xy);
                        int N = (int) Math.sqrt(vals.length);
                        FloatProcessor patch = new FloatProcessor(N, N, vals);
                        // necessary to add the scaled copy to be consistent with stack size
                        isOut.addSlice("("+b+","+m+")", patch.resize(patch_size, patch_size));

                    }

                }

            }

        }

        if (isOut.getSize()==0) {
            isOut.addSlice(new FloatProcessor(patch_size,patch_size));
        }

        return isOut;

    }

	public static ImageStack plotDelineationLines(int atX, int atY) {

        // plot of foreground versus local background estimates at each patch line
        ImageStack isOut = new ImageStack(528, 255);

        int idx = Masker2D.xy2i[atX][atY]; // read extracted peaks at this location

        if (idx!=-1) {

            int[][] delin_at_loc = PeakAnalyzer2D.delin2[idx];

            for (int b = 0; b<delin_at_loc.length; b++) {           // loop 4 branches, b index defines the strength

                if (delin_at_loc[b][M-1] != -1) { // if the last one is there - it is complete

                    for (int m = 0; m<M; m++) {

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

                        // get values sampled from the local patch
                        float[] vals_patch  = localPatchVals(prev_x, prev_y, curr_x, curr_y, inimg_xy);
                        // get values sampled from the local line
                        float[] vals_line   = localLineVals(prev_x, prev_y, curr_x, curr_y, inimg_xy);

                        int L = vals_line.length;

                        float medn_patch    = Stat.median(vals_patch);
                        float medn_line     = Stat.median(vals_line);

                        float vals_line_max = Float.NEGATIVE_INFINITY;
                        float vals_line_min = Float.POSITIVE_INFINITY;
                        for (int aa=0; aa<L; aa++) {
                            if (vals_line[aa]>vals_line_max) {
                                vals_line_max = vals_line[aa];
                            }
                            if (vals_line[aa]<vals_line_min) {
                                vals_line_min = vals_line[aa];
                            }
                        }

                        float[] medn_patch_plot = new float[L];
                        Arrays.fill(medn_patch_plot, medn_patch);
                        float[] medn_line_plot  = new float[L];
                        Arrays.fill(medn_line_plot, medn_line);

                        float[] vals_line_x = new float[vals_line.length];
                        for (int aa=0; aa<vals_line.length; aa++) {
                            vals_line_x[aa] = aa;
                        }

                        Plot plt = new Plot("BCK="+medn_patch+",LINE="+medn_line, "", "");
                        plt.setLimits(0, L-1, Math.min(medn_patch, vals_line_min), Math.max(medn_patch, vals_line_max));
                        plt.addPoints(vals_line_x, vals_line, Plot.CIRCLE);

                        plt.draw();
                        plt.setColor(Color.RED);
                        plt.setLineWidth(2);
                        plt.addPoints(vals_line_x, medn_line_plot, Plot.LINE);
                        plt.draw();
                        plt.setColor(Color.BLUE);
                        plt.setLineWidth(2);
                        plt.addPoints(vals_line_x, medn_patch_plot, Plot.LINE);

                        // necessary to add the scaled copy to be consistent with stack size
                        isOut.addSlice("("+b+","+m+"),ON="+medn_line+",OFF="+medn_patch, plt.getProcessor());

                    }

                }

            }

        }

        if (isOut.getSize()==0) {
            isOut.addSlice(new ByteProcessor(528,255));
        }

        return isOut;

	}

	public static float medianAlongLine(float x1, float y1, float x2, float y2, float[][] inimg_xy) {

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

//			if (atX<0 || atX>inimg_xy.length-1 || atY<0 || atY>inimg_xy[0].length-1) {
//				System.out.println(String.format("\n%5d.(%4d) element (X, Y) was set wrong: (%5.2f, %5.2f) when sampling values between (%5.2f, %5.2f)->(%5.2f, %5.2f)\n", cc, elementsInLine, atX, atY, x1, y1, x2, y2));
//			}

			valuesAlongLine[cc] = Interpolator.interpolateAt(atX, atY, inimg_xy);

		}

		return Stat.median(valuesAlongLine);

	}

	public static void exportFeatsCsv(String file_path) {

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

    }

    public static void exportFrames(String file_path) {

        PrintWriter logWriter = null; //initialize writer

        try {
            logWriter = new PrintWriter(file_path);
            logWriter.print("");
            logWriter.close();
        } catch (FileNotFoundException ex) {}   // empty the file before logging...

        try {                                   // initialize
            logWriter = new PrintWriter(new BufferedWriter(new FileWriter(file_path, true)));
        } catch (IOException e) {}

        //main loop
        for (int ii=0; ii<delin2.length; ii++) {
            for (int jj=0; jj<delin2[ii].length; jj++) { // loop streamlines
                for (int kk=0; kk<delin2[ii][jj].length; kk++) { // print M elements in the same streamline
                    logWriter.print(String.format("%6d", delin2[ii][jj][kk]));
                    if (kk<delin2[ii][jj].length-1) {
                        logWriter.print(",\t");
                    }
                    else {
                        logWriter.print("\n");
                    }
                }
            }
        }

        logWriter.close(); // close log

    }

    public static void exportFeatsLegend(String file_path) {

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

    }

    /*
        methods that deal with local image patch - (rectangle defined with prev_xy and curr_xy)
        extract set of PointRoi-s used for sampling,    localPatchLocs()   (vizualization only)
        extract array with sampled image values,        localPatchVals()   (vizualization only)
        or final median estimate of the patch values    localPatchMedn()   (real thing)
        (median serves as the estimate of the background level)
     */
    private static ArrayList<PointRoi> localPatchLocs(float x1, float y1, float x2, float y2) { // , float[][] inimg_xy

        float l = (float) Math.sqrt(Math.pow(x2-x1, 2)+Math.pow(y2-y1,2));   // l will vary depending on the chosen pixel (it's not always perfectly symmetric)
        int N = (int) Math.ceil( l / (samplingStep*2) );
        float vx, vy, wx, wy; // vectors that cover the square

        vx = (x2-x1)/l;
        vy = (y2-y1)/l;
        wx = vy;
        wy = -vx;

        // allocate outputs
        ArrayList<PointRoi> pts = new ArrayList<PointRoi>((2*N+1)*(2*N+1));

        for (int ii=0; ii<=2*N; ii++) { // loops vector v

            for (int jj=-N; jj<=N; jj++) { // loops vector w

                float curr_x = x1 + ii * samplingStep * vx + jj * samplingStep * vy;
                float curr_y = y1 + ii * samplingStep * wx + jj * samplingStep * wy;

                PointRoi pt = new PointRoi(curr_x+.5f, curr_y+.5f);
                pt.setFillColor(Color.DARK_GRAY);
                pts.add(pt);

            }

        }

        return pts;

    }

    private static float[] localPatchVals(float x1, float y1, float x2, float y2, float[][] inimg_xy) {
        // estimate the background as median or 1st quartile of image values from the square limited with (x1, y1) and (x2, y2)
        float l = (float) Math.sqrt(Math.pow(x2-x1, 2)+Math.pow(y2-y1,2));   // l will vary depending on the chosen pixel (it's not always perfectly symmetric)
        int N = (int) Math.ceil( l / (samplingStep*2) );
        float vx, vy, wx, wy; // vectors that cover the square

        vx = (x2-x1)/l;
        vy = (y2-y1)/l;
        wx = vy;
        wy = -vx;

        // allocate outputs
        float[] vals            = new float[(2*N+1)*(2*N+1)];

        int cnt = 0;

        for (int ii=0; ii<=2*N; ii++) { // loops vector v

            for (int jj=-N; jj<=N; jj++) { // loops vector w

                float curr_x = x1 + ii * samplingStep * vx + jj * samplingStep * vy;
                float curr_y = y1 + ii * samplingStep * wx + jj * samplingStep * wy;

                vals[cnt] = Interpolator.interpolateAt(curr_x, curr_y, inimg_xy);
                cnt++;

            }

        }

        return vals;

    }

    private static float localPatchMedn(float x1, float y1, float x2, float y2, float[][] inimg_xy) {
        // estimate the background as median or 1st quartile of image values from the square limited with (x1, y1) and (x2, y2)
        float l = (float) Math.sqrt(Math.pow(x2-x1, 2)+Math.pow(y2-y1,2));   // l will vary depending on the chosen pixel (it's not always perfectly symmetric)
        int N = (int) Math.ceil( l / (samplingStep*2) );
        float vx, vy, wx, wy; // vectors that cover the square

        vx = (x2-x1)/l;
        vy = (y2-y1)/l;
        wx = vy;
        wy = -vx;

        // allocate outputs
        float[] vals            = new float[(2*N+1)*(2*N+1)];

        int cnt = 0;

        for (int ii=0; ii<=2*N; ii++) { // loops vector v

            for (int jj=-N; jj<=N; jj++) { // loops vector w

                float curr_x = x1 + ii * samplingStep * vx + jj * samplingStep * vy;
                float curr_y = y1 + ii * samplingStep * wx + jj * samplingStep * wy;

                vals[cnt] = Interpolator.interpolateAt(curr_x, curr_y, inimg_xy);
                cnt++;

            }

        }

        return Stat.median(vals);

    }

	/*
		methods that deal with local line (defined with prev_xy and curr_xy)
	 */
    private static ArrayList<OvalRoi> localLineLocs(float x1, float y1, float x2, float y2) {

        float R = 0.5f; // radius of the oval circles

        float dist = (float) Math.sqrt(Math.pow(x2-x1, 2) + Math.pow(y2-y1, 2)); //  + Math.pow(z2lay-z1lay, 2)

        int elementsInLine = (int) (dist / samplingStep);  // how many increment can safely fit between
//        float[] valuesAlongLine = new float[elementsInLine];

        float dx = (x2 - x1) / dist;
        float dy = (y2 - y1) / dist;
        // [dx, dy] is unit vector

        dx *= samplingStep;
        dy *= samplingStep;

        ArrayList<OvalRoi> pts = new ArrayList<OvalRoi>(elementsInLine);

        for (int cc = 0; cc<elementsInLine; cc++) {

            float atX = x1      + cc * dx;
            float atY = y1      + cc * dy;
//			float atZ = z1lay   + cc * dz;

//            valuesAlongLine[cc] = Interpolator.interpolateAt(atX, atY, inimg_xy);

            OvalRoi ovroi = new OvalRoi(atX-(R/2)+.5f, atY-(R/2)+.5f, R, R);
            ovroi.setStrokeWidth(R/2);
            ovroi.setStrokeColor(Color.BLUE);
            pts.add(ovroi);

        }

        return pts;
    }

    private static float[] localLineVals(float x1, float y1, float x2, float y2, float[][] inimg_xy) {

        float dist = (float) Math.sqrt(Math.pow(x2-x1, 2) + Math.pow(y2-y1, 2)); //  + Math.pow(z2lay-z1lay, 2)

        int elementsInLine = (int) (dist / samplingStep);  // how many increment can safely fit between
        float[] valuesAlongLine = new float[elementsInLine];

        float dx = (x2 - x1) / dist;
        float dy = (y2 - y1) / dist;
        // [dx, dy] is unit vector

        dx *= samplingStep;
        dy *= samplingStep;

        for (int cc = 0; cc<elementsInLine; cc++) {

            float atX = x1      + cc * dx;
            float atY = y1      + cc * dy;
//			float atZ = z1lay   + cc * dz;

            valuesAlongLine[cc] = Interpolator.interpolateAt(atX, atY, inimg_xy);

        }

        return valuesAlongLine;

    }

    private static float localLineMedn(float x1, float y1, float x2, float y2, float[][] inimg_xy) {

        float dist = (float) Math.sqrt(Math.pow(x2-x1, 2) + Math.pow(y2-y1, 2)); //  + Math.pow(z2lay-z1lay, 2)

        int elementsInLine = (int) (dist / samplingStep);  // how many increment can safely fit between
        float[] valuesAlongLine = new float[elementsInLine];

        float dx = (x2 - x1) / dist;
        float dy = (y2 - y1) / dist;
        // [dx, dy] is unit vector

        dx *= samplingStep;
        dy *= samplingStep;

        for (int cc = 0; cc<elementsInLine; cc++) {

            float atX = x1      + cc * dx;
            float atY = y1      + cc * dy;
//			float atZ = z1lay   + cc * dz;

            valuesAlongLine[cc] = Interpolator.interpolateAt(atX, atY, inimg_xy);

        }

        return Stat.median(valuesAlongLine);

    }

}