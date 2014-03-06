package detection2d;

import aux.Stat;
import detection.Interpolator;
import fit.Fitter1D;
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
 *
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
    public static float D           = 4f;
    public static int M             = 2;                    // how much it expands recursively from the center
    public static float minCos      = 0.6f;                 // allowed derail
    public static float scatterDist = 5;                    // allowed scatter dist, upper limit, half of the neighbouring peaks should be within


    private static float    samplingStep = 0.6f;            // when sampling image values to extract features
	private static int		L = 4;                          // will define how many are taken along the diameter, in radial direction
	private static int      dim;                            // cross profile length
	private static int      dim_half;                       // cross profile half-length
    private static int      radial_step = 3;

    private static boolean useVarScore = false;

	public static float threshold = 10f;				// to normalize the scores in feature calculation (fuzzy inputs 0-1)

    // OUTPUT
    // associate the peaks and link follow-up points
    public static int[][][] delin2;                     // N(foreground locs.) x 4(max. threads) x M(follow-up locs) contains index for each location
	// extract features
    public static float[][]	feat2;						// N(foreground locs.) x 7 (6+1) features

    // PROCESSING UNITS
    private static Fitter1D fitter;

    public PeakAnalyzer2D(int n0, int n1)
    {
        this.begN = n0;
        this.endN = n1;
    }

    public static void loadTemplate(int[][] _i2xy, int[][] _xy2i, int[][][] _peaks_xy, float[][] _inimg_xy, byte[][] _backgr_xy,
									int _M, float _minCos, float _scatterDist, float _threshold, float _D) {

        D = _D;
        dim = (int) Math.ceil( D / (samplingStep*2) );
        dim_half = dim;
        dim = 2*dim + 1;

        fitter = new Fitter1D(dim, false); // dim = profile width with current samplingStep, verbose = false
        fitter.showTemplates();

        i2xy = _i2xy;
        xy2i = _xy2i;
        peaks_xy = _peaks_xy;
		inimg_xy = _inimg_xy;
		backg_xy = _backgr_xy;
        M = _M;
        minCos = _minCos;
        scatterDist = _scatterDist;
		threshold = _threshold;

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
			// extract features at this spot, feat. value = median(line_beg_end) - median(patch_values) = foreground - background
			// take just min/max along one streamline thread, 4xM -> 4x2(min,max)
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



			float fg_estimate 	= medianAtPoint(center_x, center_y, inimg_xy);
			float bg_estimate 	= backg_xy[center_x][center_y] & 0xff;
			float theta 		= fg_estimate - bg_estimate;

			float curr_val = score_bias(theta, threshold); // score calculation (normalized, limited...)
			feat2[locationIdx][6] = curr_val;  // central one at index 6

		}

    }

	private static float[][] takeMinMax(int[][] delineation_idxs, int root_idx, boolean verbose) {

        //if (verbose) { IJ.log(""); } // just in case

		float[][] calc_feats = new float[delineation_idxs.length][2];   // output to fill in
		int M = delineation_idxs[0].length;

		for (int ii=0; ii<calc_feats.length; ii++) { // loop branches

			if (delineation_idxs[ii][M-1] == -1) { // not complete line, last is not there
				calc_feats[ii][0] = 0;
				calc_feats[ii][1] = 0;
			}
			else {  // it is complete, store min/max val for this branch streamline

				float min_val = Float.POSITIVE_INFINITY;
				float max_val = Float.NEGATIVE_INFINITY;
				float curr_val;

				for (int jj=0; jj<delineation_idxs[ii].length; jj++) { // along each branch streamline

					if (delineation_idxs[ii][jj] != -1) { // calculate curr_val

						int prev_i, prev_x, prev_y;
						int curr_i, curr_x, curr_y;

						prev_i = (jj==0)? root_idx : delineation_idxs[ii][jj-1];
						prev_x = i2xy[prev_i][0];
						prev_y = i2xy[prev_i][1];
						curr_i = delineation_idxs[ii][jj];
						curr_x = i2xy[curr_i][0];
						curr_y = i2xy[curr_i][1];

						/*
							feature calculation  - on one segment of the streamline
						*/

						float[] fg_vals = localLineVals(prev_x, prev_y, curr_x, curr_y, inimg_xy);
						float[] bg_vals = localPatchVals(prev_x, prev_y, curr_x, curr_y);

						float fg_estimate 	= Stat.median(fg_vals);
						float bg_estimate 	= Stat.median(bg_vals);
						float sigma 		= Stat.std(fg_vals, Stat.average(fg_vals));
						float theta 		= fg_estimate - bg_estimate;

						curr_val = score_bias(theta, threshold); // score calculation (normalized, limited...)
						if (useVarScore) {
							curr_val = curr_val * score_var(theta, sigma);
						}

					}
					else {
						curr_val = 0;    // this is now redundant because it won't happen
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

		}

		return calc_feats;

	}

	// debug
//	public static void print(int atX, int atY) {
//
//		int atLoc = xy2i[atX][atY];
//
//        IJ.log(String.format("/**** LOC (%5d, %5d) [%10d] ****/", atX, atY, atLoc));
//
//		if (atLoc != -1) {
//
//			IJ.log("DELINEATION INDEXES:");
//            for (int ii=0; ii<delin2[atLoc].length; ii++) {
//				IJ.log("-> "+Arrays.toString(delin2[atLoc][ii]));
//			}
//
//			float[][] test = takeMinMax(delin2[atLoc], atLoc, true);
//			IJ.log("CALCULATED FEATURES:");
//			for (int ii=0; ii<test.length; ii++) {
//				IJ.log("MIN/MAX-> "+Arrays.toString(test[ii]));
//			}
//
//            IJ.log("FEATURE VECTOR(TH11, TH12, TH21, TH22, TH31, TH32, C):");
//			String plot_feats = "[";
//			for (int kk = 0; kk<feat2[atLoc].length; kk++) {
//				plot_feats += IJ.d2s(feat2[atLoc][kk], 2);
//				if(kk<feat2[atLoc].length-1) {
//					plot_feats += ", ";
//				}
//				else {
//					plot_feats += "]";
//				}
//			}
//			IJ.log(plot_feats); // Arrays.toString(feat2[atLoc])
//
//        }
//		else {
//			IJ.log("background point, no data!");
//		}
//
//	}

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

    /*
        outputs (to visualize delineation and extracted features)
     */
    public static Overlay getDelineationOverlay(int atX, int atY)
    {

        // return the delineated local structure Overlay for visualization
        Overlay ov = new Overlay();
        float R = 1; // radius of the circles written
        OvalRoi ovalroi = new OvalRoi(atX-(R/2)+.5f, atY-(R/2)+.5f, R, R);
        ov.add(ovalroi);

        // read extracted peaks at this location
        int idx = Masker2D.xy2i[atX][atY];

        /*
            show selected points
         */

        if (idx!=-1) {

            int[][] delin_at_loc = PeakAnalyzer2D.delin2[idx];

            for (int b = 0; b<delin_at_loc.length; b++) {           // loop 4 branches, b index defines the branch "strength"

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
                        ArrayList<OvalRoi> line_pts = localPatchCrossProfilesLocs(prev_x, prev_y, curr_x, curr_y);
                        for (int aa = 0; aa < line_pts.size(); aa++) ov.add(line_pts.get(aa));

                        // add the patch locations
                        ArrayList<PointRoi> pts = localPatchValsLocs(prev_x, prev_y, curr_x, curr_y);
                        for (int aa=0; aa<pts.size(); aa++) ov.add(pts.get(aa));

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


        // create new ImageStack with every layer corresponding to one patch
        // involved in modelling the local structure
        ImageStack isOut = new ImageStack(dim, M*dim);//(patch_size,patch_size);

        int idx = Masker2D.xy2i[atX][atY]; // read extracted peaks at this location

        if (idx!=-1) {

            int[][] delin_at_loc = PeakAnalyzer2D.delin2[idx];      //

            for (int b = 0; b<delin_at_loc.length; b++) {           // loop 4 branches, b index defines the strength

                if (delin_at_loc[b][M-1] != -1) {                   // if the last one is there - it is complete

					float[] concat = new float[M*dim*dim];
					int cnt = 0;

                    for (int m = M-1; m>=0; m--) { // align them backwards

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

                        float[] vals = localPatchVals(prev_x, prev_y, curr_x, curr_y);

						for (int aa=0; aa<vals.length; aa++) { // append it
							concat[cnt] = vals[aa];
                            cnt++;
						}

                    }

					// add slice after looping
					FloatProcessor patch = new FloatProcessor(dim, M*dim, concat);
					isOut.addSlice("bch "+b+")", patch); // .resize(patch_size, patch_size)

				}

            }

        }

        if (isOut.getSize()==0) {
            isOut.addSlice(new FloatProcessor(dim, dim));//(patch_size,patch_size));
        }

        return isOut;

    }

	public static ImageStack plotDelineationProfiles(int atX, int atY)
    {

        ImageStack isOut = new ImageStack(528, 255);

        int idx = Masker2D.xy2i[atX][atY]; // read extracted peaks at this location

        if (idx!=-1) {

            int[][] delin_at_loc = PeakAnalyzer2D.delin2[idx];

            for (int b = 0; b<delin_at_loc.length; b++) {           // loop 4 branches, b index defines the strength

                if (delin_at_loc[b][M-1] != -1) { // if the last one is there - it is complete

                    ArrayList<float[]> profiles_along = new ArrayList<float[]>(); // list of profiles along the branch

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

                        // get cross-profile values sampled from the local patch (aligned with the patch)
                        ArrayList<float[]> vals = localPatchCrossProfiles(prev_x, prev_y, curr_x, curr_y);

                        // append "vals" to "profiles_along" list
                        for (int aa = 0; aa < vals.size(); aa++) profiles_along.add(vals.get(aa));

                    }

                    // profiles_along with M*L profiles
                    float[] xx = new float[dim];  // xaxis plot
                    for (int aa=0; aa<dim; aa++) xx[aa] = aa;

                    Plot plt = new Plot("", "", "");
                    plt.setLimits(0, dim-1, 0, 1);
                    for (int aaa=0; aaa<profiles_along.size(); aaa++)
                        plt.addPoints(xx, profiles_along.get(aaa), Plot.LINE);
                    plt.draw();

                    // fitting gaussians
                    plt.setColor(Color.RED);
                    plt.setLineWidth(2);
                    System.out.println("bch " + b + " :");
                    for (int aaa=0; aaa<profiles_along.size(); aaa++) { // add the fittings to the plot
                        float[] out_idx_scr = fitter.fit(profiles_along.get(aaa), "NSSD");
                        float[] curr_fit = fitter.getTemplate((int)out_idx_scr[0]);
                        plt.addPoints(xx, curr_fit, Plot.LINE);
                    }
                    plt.draw();
                    isOut.addSlice("bch " + b + "", plt.getProcessor());

                }

            }

        }

        if (isOut.getSize()==0) {
            isOut.addSlice(new ByteProcessor(528,255));
        }

        return isOut;

	}

    public static ImageProcessor plotDelineationFeatures(int atX, int atY)
    {

        Plot feature_plot = new Plot("", "", "");
        feature_plot.setLimits(0, 5.5, 0, 1);
        // 4x(MxL) features, xx is x axis plot, yy is y axis plot
        float[] xx = new float[4*M*L];
        float[] yy = new float[4*M*L];
        for (int i=0; i<(4*M*L); i++) yy[i]=1;
        for (int i=0; i<4; i++) {
            for (int j=0; j<(M*L); j++) {
                xx[i*(M*L)+j] = i*1.5f + j * (1 / (float)(M*L-1));
            }
        }

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

                        // get cross-profile values sampled from the local patch (aligned with the patch)
                        float[] vals = localPatchCrossProfileFitScores(prev_x, prev_y, curr_x, curr_y);

                        // append L "vals" to the feature vector 4*(M*L)
                        for (int l = 0; l < vals.length; l++) yy[b*(M*L)+m*L+l] = vals[l];

                    }

                }

            }

        }

        IJ.log("------");
        String[] leg = new String[4];
        for (int ii=0; ii<4; ii++) {
            leg[ii] = "";
            for (int jj=0; jj<(L*M); jj++) {
                leg[ii] += IJ.d2s(yy[ii*(L*M)+jj], 2)+"\t";
            }
            IJ.log(leg[ii]);
        }
        IJ.log("------");

        feature_plot.addPoints(xx, yy, Plot.BOX);
        return feature_plot.getProcessor();

    }

    public static float[] getDelineationFeatures(int atX, int atY)
    {
        // will calculate the features (fitting scores of the gaussian profiles along the delineated branch, M*L cross-sections)
        // & store them in (M*L) dimensional vector, with the first one being the closest to the central root location
        // 4*(M*L) fit scores
        // 4*(M*L) overlap scores (one for each fit)
        // every cross section line has one score for fit and one for overlap with the highest one from the other branches
		// there are max 4 branches in 2D, in case they are missing - the rest of the features are filled with modelled values - modeling bad scores
        float[][] scores = new float[2][4*M*L]; // 4 would cover both cross sections and junctions

		// first row - normalized fitting scores
        for (int i=0; i<scores[0].length; i++) {
			scores[0][i] = 1f; // missing values modelled as bad fit
		}

		// second row - overlap scores
		for (int i=0; i<scores[1].length; i++) {
			scores[1][i] = 0f; // missing values modelled as no overlap
		}

        int idx = Masker2D.xy2i[atX][atY]; // read extracted peaks at this location

        if (idx!=-1) {

            int[][] delin_at_loc = PeakAnalyzer2D.delin2[idx];

            for (int b = 0; b<delin_at_loc.length; b++) {           // loop 4 branches, b index defines the strength

                if (delin_at_loc[b][M-1] != -1) { // if the last one is there - it is complete

                    for (int m = 0; m<M; m++) {      				// loop patches

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

                        // get cross-profile values sampled from the local patch (aligned with the patch)
                        float[] vals = localPatchCrossProfileFitScores(prev_x, prev_y, curr_x, curr_y);  // L scores per patch

                        // append L "vals" from the patch to the feature vector 4*(M*L)
                        for (int l = 0; l < vals.length; l++) scores[0][b*(M*L)+m*L+l] = vals[l];

                    }

                }

            }

        }

        return feats;

    }

	/*
		score calculation
	 */
	private static float score_bias(float theta, float threshold) { // theta is foreground-background

		if (theta<=0) {
			return 0;
		}
		else if (theta>=threshold) {
			return 1;
		}
		else {
			return (float) (1 / (1 + Math.exp( -(theta-threshold/2) / (threshold/10) )));
		}

	}

	private static float score_var(float theta, float sigma) {

		float theta1 = 1, sigma1 = 1;

		if (theta>1) {
			theta1 = theta;
		}

		if (sigma>1) {
			sigma1 = sigma;
		}

		float sigma_over_theta = sigma1 / theta1;

		float out = (float) (1f - 1f/(1f+Math.exp(-(sigma_over_theta - 1f/2) / (1f/10)) ));

//		IJ.log("s/d="+sigma_over_theta+", out="+out);

		return out;

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
         localPatchAvgProfile() (average of the profiles along)
        (median serves as the estimate of the background level)
     */

    private static ArrayList<PointRoi> localPatchValsLocs(float x1, float y1, float x2, float y2)
    {

        float l = (float) Math.sqrt(Math.pow(x2-x1, 2)+Math.pow(y2-y1,2));
        float vx = (x2-x1)/l;
		float vy = (y2-y1)/l;
		float wx = vy;
		float wy = -vx;

        ArrayList<PointRoi> pts = new ArrayList<PointRoi>(dim*dim);

        for (int ii=0; ii<=2*dim_half; ii++) { // loops vector v
            for (int jj=-dim_half; jj<=dim_half; jj++) { // loops vector w
                float curr_x = x2 - ii * samplingStep * vx + jj * samplingStep * vy;
                float curr_y = y2 - ii * samplingStep * wx + jj * samplingStep * wy;
                PointRoi pt = new PointRoi(curr_x+.5f, curr_y+.5f);
                pt.setFillColor(Color.DARK_GRAY);
                pts.add(pt);
            }
        }

        return pts;

    }

    private static float[] localPatchVals(float x1, float y1, float x2, float y2)
    {

        float l = (float) Math.sqrt(Math.pow(x2-x1, 2)+Math.pow(y2-y1, 2));
        float vx = (x2-x1)/l;
        float vy = (y2-y1)/l;
        float wx = vy;
        float wy = -vx;

        float[] vals = new float[dim*dim];
        int cnt = 0;

        for (int ii=0; ii<=2*dim_half; ii++) { // loops vector v
            for (int jj=-dim_half; jj<=dim_half; jj++) { // loops vector w
                float curr_x = x2 - ii * samplingStep * vx + jj * samplingStep * vy;
                float curr_y = y2 - ii * samplingStep * wx + jj * samplingStep * wy;
                vals[cnt] = Interpolator.interpolateAt(curr_x, curr_y, inimg_xy);
                cnt++;
            }
        }

        return vals;

    }

    private static ArrayList<float[]> localPatchCrossProfiles(float x1, float y1, float x2, float y2)
    {

        float l = (float) Math.sqrt(Math.pow(x2-x1, 2)+Math.pow(y2-y1, 2));
        float vx = (x2-x1)/l;
        float vy = (y2-y1)/l;
        float wx = vy;
        float wy = -vx;

        ArrayList<float[]> vals            = new ArrayList<float[]>();

		float samplingStepRadial = D / (float)(L-1);

        for (int ii=0; ii<L; ii++) { // loops L of them in radial direction

            float[] val = new float[dim];
            int cnt = 0;
			float val_min = Float.POSITIVE_INFINITY;
			float val_max = Float.NEGATIVE_INFINITY;

            for (int jj=-dim_half; jj<=dim_half; jj++) { // loops vector w
                float curr_x = x2 - ii * samplingStepRadial * vx + jj * samplingStep * vy;
                float curr_y = y2 - ii * samplingStepRadial * wx + jj * samplingStep * wy;

                val[cnt] = Interpolator.interpolateAt(curr_x, curr_y, inimg_xy);

				if (val[cnt]>val_max) val_max = val[cnt];
				if (val[cnt]<val_min) val_min = val[cnt];

                cnt++;
            }

			// normalize min-max so that they're from 0 to 1
			for (int iii = 0; iii<val.length; iii++){
				val[iii] = (val[iii]-val_min)/(val_max-val_min);
			}

            vals.add(val);

        }

        return vals;

    }

	private static ArrayList<OvalRoi> localPatchCrossProfilesLocs(float x1, float y1, float x2, float y2)
    {

		float R = 1f;

		float l = (float) Math.sqrt(Math.pow(x2-x1, 2)+Math.pow(y2-y1, 2));
		float vx = (x2-x1)/l;
		float vy = (y2-y1)/l;
		float wx = vy;
		float wy = -vx;

		ArrayList<OvalRoi> pts = new ArrayList<OvalRoi>(dim*L);

		float samplingStepRadial = D / (float)(L-1);

		for (int ii=0; ii<L; ii++) { // loops L of them in radial direction

//			float[] val = new float[dim];
//			int cnt = 0;
//			float val_min = Float.POSITIVE_INFINITY;
//			float val_max = Float.NEGATIVE_INFINITY;

			for (int jj=-dim_half; jj<=dim_half; jj++) { // loops vector w
				float curr_x = x2 - ii * samplingStepRadial * vx + jj * samplingStep * vy;
				float curr_y = y2 - ii * samplingStepRadial * wx + jj * samplingStep * wy;

//				val[cnt] = Interpolator.interpolateAt(curr_x, curr_y, inimg_xy);

//				if (val[cnt]>val_max) val_max = val[cnt];
//				if (val[cnt]<val_min) val_min = val[cnt];

//				cnt++;

				OvalRoi pt = new OvalRoi(curr_x+.5f-(R/2), curr_y+.5f-(R/2), R, R);
				pt.setFillColor(Color.BLUE);
				pts.add(pt);

			}

			// normalize min-max so that they're from 0 to 1
//			for (int iii = 0; iii<val.length; iii++){
//				val[iii] = (val[iii]-val_min)/(val_max-val_min);
//			}

//			pts.add(val);

		}

		return pts;


	}

    private static float[] localPatchCrossProfileFitScores(float x1, float y1, float x2, float y2)
    {

        float l = (float) Math.sqrt(Math.pow(x2-x1, 2)+Math.pow(y2-y1, 2));
        float vx = (x2-x1)/l;
        float vy = (y2-y1)/l;
        float wx = vy;
        float wy = -vx;

        float[] fit_scores = new float[L];
        float[] dummy;

        float samplingStepRadial = D / (float)(L-1);

        for (int ii=0; ii<L; ii++) { // loops L of them in radial direction

            float[] val = new float[dim];
            int cnt = 0;
            float val_min = Float.POSITIVE_INFINITY;
            float val_max = Float.NEGATIVE_INFINITY;

            for (int jj=-dim_half; jj<=dim_half; jj++) { // loops vector w
                float curr_x = x2 - ii * samplingStepRadial * vx + jj * samplingStep * vy;
                float curr_y = y2 - ii * samplingStepRadial * wx + jj * samplingStep * wy;

                val[cnt] = Interpolator.interpolateAt(curr_x, curr_y, inimg_xy);

                if (val[cnt]>val_max) val_max = val[cnt];
                if (val[cnt]<val_min) val_min = val[cnt];

                cnt++;
            }

            // normalize min-max so that they're from 0 to 1
            for (int iii = 0; iii<val.length; iii++){
                val[iii] = (val[iii]-val_min)/(val_max-val_min);
            }

            // fit the normalized profile
            dummy = fitter.fit(val, "NSSD");
            fit_scores[ii] = dummy[1];

        }

        return fit_scores;

    }
	/*
		methods that deal with local line (defined with prev_xy and curr_xy)
	 */
//    private static ArrayList<OvalRoi> localLineLocs(float x1, float y1, float x2, float y2) {
//
//        float R = 0.5f; // radius of the oval circles
//
//        float dist = (float) Math.sqrt(Math.pow(x2-x1, 2) + Math.pow(y2-y1, 2)); //  + Math.pow(z2lay-z1lay, 2)
//
//        int elementsInLine = (int) (dist / samplingStep);  // how many increment can safely fit between
////        float[] valuesAlongLine = new float[elementsInLine];
//
//        float dx = (x2 - x1) / dist;
//        float dy = (y2 - y1) / dist;
//        // [dx, dy] is unit vector
//
//        dx *= samplingStep;
//        dy *= samplingStep;
//
//        ArrayList<OvalRoi> pts = new ArrayList<OvalRoi>(elementsInLine);
//
//        for (int cc = 0; cc<elementsInLine; cc++) {
//
//            float atX = x1      + cc * dx;
//            float atY = y1      + cc * dy;
////			float atZ = z1lay   + cc * dz;
//
////            valuesAlongLine[cc] = Interpolator.interpolateAt(atX, atY, inimg_xy);
//
//            OvalRoi ovroi = new OvalRoi(atX-(R/2)+.5f, atY-(R/2)+.5f, R, R);
//            ovroi.setStrokeWidth(R/2);
//            ovroi.setStrokeColor(Color.BLUE);
//            pts.add(ovroi);
//
//        }
//
//        return pts;
//    }

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

}