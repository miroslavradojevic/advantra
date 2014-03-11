package detection2d;

import aux.Stat;
import detection.Interpolator;
import fit.Fitter1D;
import ij.IJ;
import ij.ImageStack;
import ij.gui.*;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;

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
 * Calculates and exports the fit scores (NCC) features at every location after frame delineation (feat2)
 * feat2:
 *      there are L such NCC values per patch, and M patches per branch (up to 4 branches): every location 4*M*L elements
 * ratio2:
 *      scores for each branch 4*1 (M*L elements sum up to one number telling the ratio of those above the threshold)
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
    public static float     D           = 4f;
    public static int       M           = 2;              // how much it expands recursively from the center
    public static float     minCos      = 0.6f;             // allowed derail
    public static float     scatterDist = 5;                // allowed scatter dist, upper limit, half of the neighbouring peaks should be within


    private static float    samplingStep = 0.6f;            // when sampling image values to extract features
    private static int		L;                              // will define how many are taken along the diameter, in radial direction
    private static float    samplingStepLongitudinal;       // sampling along the streamline of patches
    private static int      dim;                            // cross profile length
	private static int      dim_half;                       // cross profile half-length

	public static float     thresholdNCC = .85f;			// to threshold the NCC scores in two sets (ON, OFF)

    // OUTPUT
    // associate the peaks and link follow-up points
    public static int[][][]     delin2;                     // N(foreground locs.) x 4(max. threads) x M(follow-up locs) contains index for each location
	// features
    public static float[][][]   feat2;						// N(foreground locs.) x 4 x (M x L) fit scores
    public static float[][]     ratio2;                     // N(foreground locs.) x 4 + 1 ratio of those above the threshold and score for the point
    public static float[][]     lhood2;

    // PROCESSING UNITS
    private static Fitter1D fitter;                         // class used to calculate fitting scores
    private static Fuzzy2D  fuzzy_logic_system;             // class used to calculate fuzzy score

    public PeakAnalyzer2D(int n0, int n1)
    {
        this.begN = n0;
        this.endN = n1;
    }

    public static void loadTemplate(int[][] _i2xy, int[][] _xy2i, int[][][] _peaks_xy, float[][] _inimg_xy, byte[][] _backgr_xy,
									int _M, float _minCos, float _scatterDist, float _thresholdNCC, float _D)
    {

        D = _D;
        dim = (int) Math.ceil( D / (samplingStep*2) );
        dim_half = dim;
        dim = 2*dim + 1;

        L = Math.round(1.5f*D);
        samplingStepLongitudinal = D / (float)(L-1);

        fitter = new Fitter1D(dim, false); // dim = profile width with current samplingStep, verbose = false
        fitter.showTemplates();

        fuzzy_logic_system = new Fuzzy2D(4*5+1, 0.25f);
        fuzzy_logic_system.showFuzzification();
        fuzzy_logic_system.showDefuzzification();

        i2xy = _i2xy;
        xy2i = _xy2i;
        peaks_xy = _peaks_xy;
		inimg_xy = _inimg_xy;
		backg_xy = _backgr_xy;
        M = _M;
        minCos = _minCos;
        scatterDist = _scatterDist;
		thresholdNCC = _thresholdNCC;

        // allocate output -> set to -1
        delin2 = new int[i2xy.length][4][M];
        for (int ii=0; ii<delin2.length; ii++)
            for (int jj = 0; jj < delin2[0].length; jj++)
                for (int kk = 0; kk < delin2[0][0].length; kk++)
                    delin2[ii][jj][kk] = -1;

		feat2   = new float[i2xy.length][4][M*L]; // fitting scores
        ratio2  = new float[i2xy.length][4+1];      // ratio of ON NCC scores (higher than threshold)
        lhood2  = new float[i2xy.length][fuzzy_logic_system.L];

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

                            if (true) { // isRobust(next_index, curr_index, scatterDist)
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

            /*	calculate features	*/
            int atX = i2xy[locationIdx][0];
            int atY = i2xy[locationIdx][1];
            getDelineationFeatures(atX, atY, feat2[locationIdx], ratio2[locationIdx], lhood2[locationIdx]);

		}

    }

	public static void print(int atX, int atY)
	{

		int atLoc = xy2i[atX][atY];

        IJ.log(String.format("/**** LOC (%5d, %5d) [%10d] ****/", atX, atY, atLoc));

		if (atLoc != -1) {
            String printout = "";

            printout += "\nDELINEATION:\n";
            for (int ii=0; ii<delin2[atLoc].length; ii++) {
                printout += ii+"\t->\t";
                for (int jj=0; jj<delin2[atLoc][ii].length; jj++) {
                    printout += (delin2[atLoc][ii][jj]!=-1)? IJ.d2s(delin2[atLoc][ii][jj], 2) : "NONE";
                    if (jj==delin2[atLoc][ii].length-1) printout += "\n";
                    else printout += ", ";
                }
            }

            printout += "\nNCC:\n";
            for (int ii=0; ii<feat2[atLoc].length; ii++) {
                printout += ii+"\t->\t";
                for (int jj=0; jj<feat2[atLoc][ii].length; jj++) {
                    printout += IJ.d2s(feat2[atLoc][ii][jj], 2);
                    if (jj==feat2[atLoc][ii].length-1) printout += "\n";
                    else printout += ", ";
                }
            }

            printout += "\nRATIO (NCC ABOVE "+IJ.d2s(thresholdNCC, 2)+"):\n";
            for (int ii=0; ii<ratio2[atLoc].length; ii++) {
                printout += (ii+1) +"\t->\t" + IJ.d2s(ratio2[atLoc][ii], 2) + "\n";
            }

            printout += "\nFUZZY LIKELIHOODS:\n";
            printout += "NON -> " + IJ.d2s(lhood2[atLoc][0], 2) + "\n";
            printout += "END -> " + IJ.d2s(lhood2[atLoc][1], 2) + "\n";
            printout += "BDY -> " + IJ.d2s(lhood2[atLoc][2], 2) + "\n";
            printout += "BIF -> " + IJ.d2s(lhood2[atLoc][3], 2) + "\n";
            printout += "CRS -> " + IJ.d2s(lhood2[atLoc][4], 2) + "\n";

            IJ.log(printout);

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
        for (int pkIdx = 0; pkIdx<pks4xXY.length; pkIdx++) { // loops them by rank - those with highest weight first, to pick the first one with good angle

            if (pks4xXY[pkIdx][0] != -1) {
                // there is a peak to check - needs to be pointing outwards more than defined

                int nextX = pks4xXY[pkIdx][0];
                int nextY = pks4xXY[pkIdx][1];

                double cosAng =
                        (
                        	(currX-prevX)*(nextX-currX) + (currY-prevY)*(nextY-currY)           // + (currZ-prevZ)*(nextZ-currZ)
                        )
                        /
                        (
                        	Math.sqrt( Math.pow(currX-prevX, 2) + Math.pow(currY-prevY, 2) ) *  //  + Math.pow(currZ-prevZ, 2)
                            Math.sqrt( Math.pow(nextX-currX, 2) + Math.pow(nextY-currY, 2) )    //  + Math.pow(nextZ-currZ, 2)
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
        outputs (to visualize delineation and extracted features)
     */
    public static Overlay getDelineationOverlay(int atX, int atY)
    {

        // return the delineated local structure Overlay for visualization
        Overlay ov = new Overlay();
        float R = 0.5f; // radius of the circles written
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
                    for (int aaa=0; aaa<profiles_along.size(); aaa++) { // add the fittings to the plot
                        float[] out_idx_scr = fitter.fit(profiles_along.get(aaa), "NCC");
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

    public static ImageStack plotDelineationFeatures(int atX, int atY)
    {

		// just read from the feat2 and ratio2
        ImageStack is_out = new ImageStack(528, 255);
        int loc_idx = xy2i[atX][atY];

		int nrFeats = 5;
        int perBin = 40;
        float[] xaxis2 = new float[nrFeats*perBin];
        float[] yaxis2 = new float[nrFeats*perBin];
        for (int i=0; i<xaxis2.length; i++) {
            xaxis2[i] = i * (nrFeats / (float)(nrFeats*perBin));
            if (Math.abs(xaxis2[i]-Math.round(xaxis2[i]))<0.1) {
                yaxis2[i] = 0;
            }
            else {
                yaxis2[i] = ratio2[loc_idx][i/perBin];
            }
        }

        Plot feature_plot1 = new Plot("", "", "ratio");
        feature_plot1.setLimits(0, nrFeats, 0, 1);
        feature_plot1.addPoints(xaxis2, yaxis2, Plot.LINE);
        feature_plot1.setColor(Color.BLACK);
        feature_plot1.draw();
        feature_plot1.setLineWidth(3);
        feature_plot1.addPoints(xaxis2, yaxis2, Plot.LINE);
        feature_plot1.draw();
        is_out.addSlice("ratioON", feature_plot1.getProcessor());

        // plot fit scores for each branch + one plot with ratios
        float[] xaxis1 = new float[M*L];
        for (int i=0; i<M*L; i++) xaxis1[i] = i;
        float[] thNCC = new float[M*L];
        for (int i=0; i<M*L; i++) thNCC[i] = thresholdNCC;

        for (int b= 0; b < 4; b++) {
            Plot feature_plot = new Plot("", "#", "NCC");
            feature_plot.setLimits(0, M*L-1, 0, 1);
            feature_plot.addPoints(xaxis1, feat2[loc_idx][b], Plot.LINE);
            feature_plot.draw();
            feature_plot.setColor(Color.RED);
            feature_plot.setLineWidth(3);
            feature_plot.addPoints(xaxis1, thNCC, Plot.LINE);
            is_out.addSlice("bch="+b, feature_plot.getProcessor());
        }

        return is_out;

    }

    public static void getDelineationFeatures(int atX, int atY, float[][] fitNCC, float[] ratioON, float[] lhood) // feat2[locIdx] ratio2[locIdx]
    {
        // will calculate the features (fitting scores of the gaussian profiles along the delineated branch, M*L cross-section fits)
        // & store them in (M*L) dimensional vector, with the first one being the closest to the central root location
        // 4*(M*L) fit scores
        // 4 ratio scores + 1 value saying how high above the background
        // every cross section line has one score for fit and one for overlap with the highest one from the other branches
		// there are max 4 branches in 2D, in case they are missing - the rest of the features are filled with modelled values - modeling bad scores

        int idx = Masker2D.xy2i[atX][atY]; // read extracted peaks at this location

        if (idx!=-1) {       // location is in foreground, there is a delineation there, change fitNCC and ratioON

            int[][] delin_at_loc = PeakAnalyzer2D.delin2[idx];

            // loop branches twice:
            // 1 - to extract cross-profiles' locations and to extract their fit scores
            // 2 - to extract ratio of those higher than threshold

            // 1
            for (int b = 0; b<delin_at_loc.length; b++) {           // loop 4 branches, b index defines the strength

                if (delin_at_loc[b][M-1] != -1) {                   // if the last one is there - it is complete

                    for (int m = 0; m<M; m++) {      				// loop patches outwards

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

                        // get cross-profile values sampled from the local patch (aligned with the patch) store them in fitNCC
                        localPatchCrossProfileFitScores(prev_x, prev_y, curr_x, curr_y, b, m, fitNCC);

                    }

                }

            }

            // 2 calculate ratios in every branch
            for (int b=0; b<4; b++) {
                int cntON = 0;
                for (int l=0; l<fitNCC[b].length; l++) {
                    if (fitNCC[b][l]>=thresholdNCC) cntON++;
                }
                ratioON[b] = cntON / (float) fitNCC[0].length;
            }

			// last feature
			ratioON[4] = Masker2D.fg_score[atX][atY];

            fuzzy_logic_system.critpointScores(ratioON[0], ratioON[1], ratioON[2], ratioON[3], ratioON[4], lhood);


        }

    }

	/*
		score calculation
	 */

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

	public static void exportNcc(String file_path)
	{
        PrintWriter logWriter = null;
        try {
            logWriter = new PrintWriter(file_path);
            logWriter.print("");
            logWriter.close();
        } catch (FileNotFoundException ex) {}

        try {
            logWriter = new PrintWriter(new BufferedWriter(new FileWriter(file_path, true)));
        } catch (IOException e) {}

        for (int ii=0; ii<feat2.length; ii++) {
            for (int jj=0; jj<feat2[ii].length; jj++) {
                for (int kk=0; kk<feat2[ii][jj].length; kk++) {
                    logWriter.print(String.format("%1.2f", feat2[ii][jj][kk]));
                    if (jj==feat2[ii].length-1 && kk==feat2[ii][jj].length-1) logWriter.print("\n");
                    else logWriter.print(",\t");
                }
            }
        }

        logWriter.close();
	}

	public static void exportRatios(String file_path)
	{
        PrintWriter logWriter = null;
        try {
            logWriter = new PrintWriter(file_path);
            logWriter.print("");
            logWriter.close();
        } catch (FileNotFoundException ex) {}

        try {
            logWriter = new PrintWriter(new BufferedWriter(new FileWriter(file_path, true)));
        } catch (IOException e) {}

        for (int ii=0; ii<ratio2.length; ii++) {
            for (int jj=0; jj<ratio2[ii].length; jj++) {
                    logWriter.print(String.format("%1.2f", ratio2[ii][jj]));
                    if (jj==ratio2[ii].length-1) logWriter.print("\n");
                    else logWriter.print(",\t");
            }
        }

        logWriter.close();
	}

    public static void exportFrames(String file_path)
	{

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

    public static void exportFeatsLegend(String file_path)
	{

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
        methods that deal with local image patch - (rectangle defined with prev_x,y and curr_x,y)
     */

    // (VIZ) extract set of PointRoi-s used for local patch sampling
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

    // (VIZ) extract array with local patch sampled values
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

    // (VIZ) list of arrays containing the local patch cross profiles that were fitted
    private static ArrayList<float[]> localPatchCrossProfiles(float x1, float y1, float x2, float y2)
    {

        float l = (float) Math.sqrt(Math.pow(x2-x1, 2)+Math.pow(y2-y1, 2));
        float vx = (x2-x1)/l;
        float vy = (y2-y1)/l;
        float wx = vy;
        float wy = -vx;

        ArrayList<float[]> vals            = new ArrayList<float[]>();

        for (int ii=0; ii<L; ii++) { // loops L of them in radial direction

            float[] val = new float[dim];
            int cnt = 0;
			float val_min = Float.POSITIVE_INFINITY;
			float val_max = Float.NEGATIVE_INFINITY;

            for (int jj=-dim_half; jj<=dim_half; jj++) { // loops vector w
                float curr_x = x2 - ii * samplingStepLongitudinal * vx + jj * samplingStep * vy;
                float curr_y = y2 - ii * samplingStepLongitudinal * wx + jj * samplingStep * wy;

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

    // (VIZ) list of of local patch cross profile line segments described with dots
	private static ArrayList<OvalRoi> localPatchCrossProfilesLocs(float x1, float y1, float x2, float y2)
    {

		float l = (float) Math.sqrt(Math.pow(x2-x1, 2)+Math.pow(y2-y1, 2));
		float vx = (x2-x1)/l;
		float vy = (y2-y1)/l;
		float wx = vy;
		float wy = -vx;

        float R = .35f;

		ArrayList<OvalRoi> pts = new ArrayList<OvalRoi>(dim*L);

		for (int ii=0; ii<L; ii++) { // loops L of them in radial direction

			for (int jj=-dim_half; jj<=dim_half; jj++) { // loops vector w
				float curr_x = x2 - ii * samplingStepLongitudinal * vx + jj * samplingStep * vy;
				float curr_y = y2 - ii * samplingStepLongitudinal * wx + jj * samplingStep * wy;

				OvalRoi pt = new OvalRoi(curr_x+.5f-(R/2), curr_y+.5f-(R/2), R, R);
				pt.setFillColor(Color.BLUE);
				pts.add(pt);

			}

		}

		return pts;

	}

    // (CALC) calculation of the fitting scores for every local patch cross profile
    private static void localPatchCrossProfileFitScores(float x1, float y1, float x2, float y2, int branch_idx, int patch_idx, float[][] fit_scores_out)
    {

        float l = (float) Math.sqrt(Math.pow(x2-x1, 2)+Math.pow(y2-y1, 2));
        float vx = (x2-x1)/l;
        float vy = (y2-y1)/l;
        float wx = vy;
        float wy = -vx;

        float[] dummy;
        float[] val = new float[dim];

        for (int ii=0; ii<L; ii++) { // loops L of them in radial direction

            int cnt = 0;
            float val_min = Float.POSITIVE_INFINITY;
            float val_max = Float.NEGATIVE_INFINITY;

            for (int jj=-dim_half; jj<=dim_half; jj++) { // loops vector w
                float curr_x = x2 - ii * samplingStepLongitudinal * vx + jj * samplingStep * vy;
                float curr_y = y2 - ii * samplingStepLongitudinal * wx + jj * samplingStep * wy;

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
            dummy = fitter.fit(val, "NCC");
            fit_scores_out[branch_idx][patch_idx*L + ii] = dummy[1];

        }

    }

    // (CALC) calculation of the line segment geometry for every local patch cross profile
    // TODO gave this up
//    private static float[][] localPatchCrossProfileGeometry(float x1, float y1, float x2, float y2)
//    {
//        // extract [px py rx ry] describing cross profile at this patch
//        float l = (float) Math.sqrt(Math.pow(x2-x1, 2)+Math.pow(y2-y1, 2));
//        float vx = (x2-x1)/l;
//        float vy = (y2-y1)/l;
//        float wx = vy;
//        float wy = -vx;
//
//        float[][] geom = new float[L][4];
//
//        float samplingStepRadial = D / (float)(L-1);
//
//        for (int ii=0; ii<L; ii++) { // loops L of them in radial direction
//
//            float[] val = new float[dim];
//            int cnt = 0;
//            //float val_min = Float.POSITIVE_INFINITY;
//            //float val_max = Float.NEGATIVE_INFINITY;
//
//            float px = x2 - ii * samplingStepRadial * vx + (-dim_half) * samplingStep * vy;
//            float py = 99;
//            float rx = 99;
//            float ry = 99;
//
//            for (int jj=-dim_half; jj<=dim_half; jj++) { // loops vector w
//
//                float curr_x = x2 - ii * samplingStepRadial * vx + jj * samplingStep * vy;
//                float curr_y = y2 - ii * samplingStepRadial * wx + jj * samplingStep * wy;
//
////                if (val[cnt]>val_max) val_max = val[cnt];
////                if (val[cnt]<val_min) val_min = val[cnt];
//
//            }
//
//            // normalize min-max so that they're from 0 to 1
////            for (int iii = 0; iii<val.length; iii++){
////                val[iii] = (val[iii]-val_min)/(val_max-val_min);
////            }
//
//            // fit the normalized profile
////            dummy = fitter.fit(val, "NSSD");
////            fit_scores[ii] = dummy[1];
//
//        }
//
//
//
//        return geom;
//
//
//
//    }

    // method that calculates overlap of two line segments

	/*
		methods that deal with local line (defined with prev_xy and curr_xy)    // TODO put it to misc
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

//    private static float[] localLineVals(float x1, float y1, float x2, float y2, float[][] inimg_xy) {
//
//        float dist = (float) Math.sqrt(Math.pow(x2-x1, 2) + Math.pow(y2-y1, 2)); //  + Math.pow(z2lay-z1lay, 2)
//
//        int elementsInLine = (int) (dist / samplingStep);  // how many increment can safely fit between
//        float[] valuesAlongLine = new float[elementsInLine];
//
//        float dx = (x2 - x1) / dist;
//        float dy = (y2 - y1) / dist;
//        // [dx, dy] is unit vector
//
//        dx *= samplingStep;
//        dy *= samplingStep;
//
//        for (int cc = 0; cc<elementsInLine; cc++) {
//
//            float atX = x1      + cc * dx;
//            float atY = y1      + cc * dy;
////			float atZ = z1lay   + cc * dz;
//
//            valuesAlongLine[cc] = Interpolator.interpolateAt(atX, atY, inimg_xy);
//
//        }
//
//        return valuesAlongLine;
//
//    }

}