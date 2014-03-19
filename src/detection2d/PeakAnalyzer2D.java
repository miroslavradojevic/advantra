package detection2d;

import aux.Stat;
import aux.Interpolator;
import fit.Fitter1D;
import ij.IJ;
import ij.ImageStack;
import ij.gui.*;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;

import java.awt.*;
import java.io.*;
import java.util.ArrayList;

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
    public static int[][]       peaks_i;             	// list of extracted peaks: N x 4 every FG location with 4 extracted peaks in indexed format
	public static int[][]		peaks_w;                // weight assigned to each peak (for expansion)
	public static float[][]		inimg_xy;				// input image (necessary for feature extraction)
	public static byte[][]		backg_xy;				// background estimation (for feature extraction)

    // PARAMETERS
    public static float     D;
    public static int       M;              // how much it expands recursively from the center
    public static float     minCos;             // allowed derail
//    public static float     scatterDist = 5;                // allowed scatter dist, upper limit, half of the neighbouring peaks should be within
    private static float    samplingStep = 0.6f;            // when sampling image values to extract features
    private static int		L;                              // will define how many are taken along the diameter, in radial direction
    private static float    samplingStepLongitudinal;       // sampling along the streamline of patches
    private static int      dim;                            // cross profile length
	private static int      dim_half;                       // cross profile half-length

//	public static float     thresholdNCC = .85f;			// to threshold the NCC scores in two sets (ON, OFF)
    private static String mode = "MSE";
    private static int nr_points = 4*5+1;    // todo agruments of loadTemplate() because they are likelihood extraction parameteres
    private static float mu_ON = 0;
    private static float sig_ON = 0.25f;
    private static float mu_OFF = 0.5f;
    private static float sig_OFF = 0.25f;


    // OUTPUT
    // associate the peaks and link follow-up points
    public static int[][][]     delin2;                     // N(foreground locs.) x 4(max. threads) x M(follow-up locs) contains index for each location
	// FEATURES, F. DESCRIPTOR, OUT LIKELIHOODS
    public static float[][][]   feat2;						// N(foreground locs.) x 4 x ((1..M) x L)       fit scores
    public static float[][][]   desc2;                      // N(foreground locs.) x 4 x 2 (mean,variance)  description
    public static float[][]     lhood2;                     // N(foreground locs.) x 5 (NON..CRS) fuzzy logic output is stored here

    // PROCESSING UNITS
    private static Fitter1D fitter;                         // class used to calculate fitting scores
    private static Fuzzy2D  fuzzy_logic_system;             // class used to calculate fuzzy score

    public PeakAnalyzer2D(int n0, int n1)
    {
        this.begN = n0;
        this.endN = n1;
    }

    public static void loadTemplate(int[][] _i2xy, int[][] _xy2i, int[][] _peaks_i, int[][] _peaks_w, float[][] _inimg_xy, byte[][] _backgr_xy,
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

        fuzzy_logic_system = new Fuzzy2D(nr_points, mu_ON, sig_ON, mu_OFF, sig_OFF);
        fuzzy_logic_system.showFuzzification();
        fuzzy_logic_system.showDefuzzification();

        i2xy = _i2xy;
        xy2i = _xy2i;
        peaks_i = _peaks_i;
		peaks_w = _peaks_w;
		inimg_xy = _inimg_xy;
		backg_xy = _backgr_xy;
        M = _M;
        minCos = _minCos;
//        scatterDist = _scatterDist;

        // allocate output -> set to -1
        delin2 = new int[i2xy.length][4][M];
        for (int ii=0; ii<delin2.length; ii++)
            for (int jj = 0; jj < delin2[0].length; jj++)
                for (int kk = 0; kk < delin2[0][0].length; kk++)
                    delin2[ii][jj][kk] = -1;

		feat2   = new float[i2xy.length][4][];                      // fitting scores
        desc2  = new float[i2xy.length][4][];                       // description of the fit scores
        lhood2  = new float[i2xy.length][fuzzy_logic_system.L];     // fuzzy likelihood output

    }

    public void run()
    {
        for (int locationIdx = begN; locationIdx < endN; locationIdx++) {

            int[] peaks_at_loc = peaks_i[locationIdx];

            // access individual peaks at this point
            for (int pp = 0; pp<peaks_at_loc.length; pp++) {  // loop 4 allocated branches

                if (peaks_at_loc[pp] != -1) { // if the peak exists

                    int indexValue = peaks_at_loc[pp];//xy2i[pkX][pkY];

                    delin2[locationIdx][pp][0] = indexValue;

                    int curr_index, prev_index, next_index;

                    curr_index = indexValue;
                    prev_index = locationIdx;

                    for (int m=1; m<M; m++) { // follow the rest of the indexes

                        // recursion : prev+curr->next index
                        next_index = getNext(prev_index, curr_index); // next follow-up will be calculated and sorted

                        if (next_index!=-1) { // -1 will be if the next one is not found

                            delin2[locationIdx][pp][m] = next_index;     // store it in output matrix

                        }
                        else { // follow-up does not exist, break looping m (extending further) but continue looping branches
                            break;
                        }

                    }

                }
                else { // it was -1 and the rest are not found
                    break; // stop for() loop - there are no more branches
                }

            } // delin2[locationIdx] is fully formed

            /*	calculate features, descriptors and likelihoods	*/
            getDelineationFeatures(locationIdx, mode);

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

            printout += "\nFEATURES:\n";
            for (int b=0; b<feat2[atLoc].length; b++) {
                printout += b+"\t->\t";

                if (feat2[atLoc][b] != null) {

                    for (int l=0; l<feat2[atLoc][b].length; l++) {
                        printout += IJ.d2s(feat2[atLoc][b][l], 2);
                        if (l==feat2[atLoc][b].length-1) printout += "\n";
                        else printout += ", ";
                    }

                }
                else {
                    printout += "NONE\n";
                }

            }

            printout += "\nDESCRIPTORS:\n";
            for (int b=0; b<desc2[atLoc].length; b++) {
                printout += (b+1) +"\t->\t"; //+ IJ.d2s(ratio2[atLoc][ii], 2) + "\n"

                if (desc2[atLoc][b] != null) {

                    for (int l=0; l<desc2[atLoc][b].length; l++) {
                        printout += IJ.d2s(desc2[atLoc][b][l], 2);
                        if (l==desc2[atLoc][b].length-1) printout += "\n";
                        else printout += ",\t";
                    }

                }
                else {
                    printout += "NONE\n";
                }

            }

            printout += "\nFUZZY LIKELIHOODS:\n";
            printout += "NON -> " + IJ.d2s(lhood2[atLoc][0], 2) + " "; // "\n";
            printout += "END -> " + IJ.d2s(lhood2[atLoc][1], 2) + " "; //"\n";
            printout += "BDY -> " + IJ.d2s(lhood2[atLoc][2], 2) + " ";
            printout += "BIF -> " + IJ.d2s(lhood2[atLoc][3], 2) + " ";
            printout += "CRS -> " + IJ.d2s(lhood2[atLoc][4], 2) + " ";

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

    /*
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

                int[][] get_peaks_nbr = new int[4][2];// peaks_xy[neigbr_i]; // peak signature of the neighbour
                for (int ii=0; ii<4; ii++) {
                    get_peaks_nbr[ii][0] = i2xy[peaks_i[neigbr_i][ii]][0];
                    get_peaks_nbr[ii][1] = i2xy[peaks_i[neigbr_i][ii]][1];
                }

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
    */

//    private static float dist(int a_x, int a_y, int b_x, int b_y){
//        return (float) Math.sqrt( Math.pow(a_x-b_x, 2) + Math.pow(a_y-b_y, 2) );
//    }

    private static int getNext(int prev_index, int curr_index) {

        // these are stacked as XY - take care on that
        int prevX = i2xy[prev_index][0];    // X
        int prevY = i2xy[prev_index][1];    // Y

        int currX = i2xy[curr_index][0];    // X
        int currY = i2xy[curr_index][1];    // Y

        // check peaks at curr
        int[]   pks4xI  = peaks_i[curr_index]; // list all indexes
		int[] 	pks4xW	= peaks_w[curr_index];

		// take the one with highest weight that is above minCos
		int next_index = -1;
		int next_weight = -1;

		for (int p = 0; p<pks4xI.length; p++) {
			if (pks4xI[p]!=-1) {
				// peak in foreground
				int check_I = pks4xI[p];
				int next_X 	= i2xy[check_I][0];
				int next_Y	= i2xy[check_I][1];

				double cosAng =
						(
								(currX-prevX)*(next_X-currX) + (currY-prevY)*(next_Y-currY)           // + (currZ-prevZ)*(nextZ-currZ)
						)
								/
								(
										Math.sqrt( Math.pow(currX-prevX, 2) + Math.pow(currY-prevY, 2) ) *  //  + Math.pow(currZ-prevZ, 2)
												Math.sqrt( Math.pow(next_X-currX, 2) + Math.pow(next_Y-currY, 2) )    //  + Math.pow(nextZ-currZ, 2)
								);

				if (cosAng>minCos) {
					// peak is outwards pointing
					if (pks4xW[p]>next_weight) {
						next_weight = pks4xW[p];
						next_index 	= pks4xI[p];
					}
				}

			}
		}

		return next_index;

	}

    /*
        outputs (to visualize delineation and extracted features)
     */
    public static Overlay getDelineationOverlay(int atX, int atY)
    {

        // return the delineated local structure Overlay for visualization
        Overlay ov = new Overlay();

        PeakExtractor2D.getPeaks(atX, atY, 3, ov);

		float Rd = 1.5f; // radius of the circles written for delineation
        Color cd = Color.RED;
        float wd = 0.25f;

		// central location
        OvalRoi ovalroi = new OvalRoi(atX-(Rd/2)+.5f, atY-(Rd/2)+.5f, Rd, Rd);
		ovalroi.setFillColor(cd);
		ovalroi.setStrokeWidth(wd);
        ov.add(ovalroi);

		// read extracted peaks at this location
        int idx = xy2i[atX][atY];

        // show delineation
        if (idx!=-1) {

            int[][] delin_at_loc = delin2[idx];

            for (int b = 0; b<delin_at_loc.length; b++) {           // loop 4 branches, b index defines the branch "strength"

                for (int m=0; m<M; m++) {

                    if (delin_at_loc[b][m] != -1) {

                        int pt_idx = delin_at_loc[b][m]; // there is a point to add
                        int pt_x = i2xy[pt_idx][0];
                        int pt_y = i2xy[pt_idx][1];

                        ovalroi = new OvalRoi(pt_x-(Rd/2)+.5f, pt_y-(Rd/2)+.5f, Rd, Rd); // add the point to the overlay
                        ovalroi.setStrokeColor(cd);
                        ovalroi.setStrokeWidth(wd);
                        ov.add(ovalroi);

						// find previous indexes
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

						// add local cross profile plots and patch locations
						ArrayList<OvalRoi> line_pts = localPatchCrossProfilesLocs(prev_x, prev_y, pt_x, pt_y);
						for (int aa = 0; aa < line_pts.size(); aa++) ov.add(line_pts.get(aa));
						//ArrayList<PointRoi> pts = localPatchValsLocs(prev_x, prev_y, pt_x, pt_y);
						//for (int aa=0; aa<pts.size(); aa++) ov.add(pts.get(aa));

                    }

                }

            }

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

        int idx = xy2i[atX][atY]; // read extracted peaks at this location

        if (idx!=-1) {

            int[][] delin_at_loc = delin2[idx];
            ArrayList<float[]> profiles_along = new ArrayList<float[]>(); // list of profiles for all branches

            for (int b = 0; b<delin_at_loc.length; b++) {           // loop 4 branches, b index defines the strength

                for (int m = 0; m<M; m++) {

                    if (delin_at_loc[b][m] != -1) {

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

                } // l loop

            } // b loop


            // plot
            float[] xx = new float[dim];  // xaxis plot
            for (int aa=0; aa<dim; aa++) xx[aa] = aa;

            for (int aaa=0; aaa<profiles_along.size(); aaa++) {
                Plot plt = new Plot("", "", "");
                plt.setLimits(0, dim-1, 0, 1);
                plt.addPoints(xx, profiles_along.get(aaa), Plot.LINE);
                plt.draw();
                // fitting templates
                plt.setColor(Color.RED);
                plt.setLineWidth(2);
                float[] out_idx_scr = fitter.fit(profiles_along.get(aaa), mode);
                float[] curr_fit = fitter.getTemplate((int)out_idx_scr[0]);
                plt.addPoints(xx, curr_fit, Plot.LINE);
                plt.draw();
                isOut.addSlice("feat="+IJ.d2s(out_idx_scr[1], 2), plt.getProcessor());
            }

        }

        if (isOut.getSize()==0) {
            isOut.addSlice(new ByteProcessor(528,255));
        }

        return isOut;

	}

    public static ImageStack plotDelineationFeatures(int atX, int atY)
    {

		// just read from the feat2 and desc2
        ImageStack is_out = new ImageStack(528, 255);
        int loc_idx = xy2i[atX][atY];

        // calculate total nr. feats
        if (loc_idx != -1) {

            int totalFeats = 0;
            for (int b=0; b<feat2[loc_idx].length; b++) {

                if (feat2[loc_idx][b] != null) {
                    totalFeats += feat2[loc_idx][b].length;
                }

            }

            float[] xaxis   = new float[totalFeats];
            float[] yaxis1  = new float[totalFeats];
            float[] yaxis2  = new float[totalFeats];

            int cnt = 0;
            for (int b=0; b<feat2[loc_idx].length; b++) {
                if (feat2[loc_idx][b] != null) {
                    for (int l=0; l<feat2[loc_idx][b].length; l++) {
                        xaxis[cnt]  = cnt;
                        yaxis1[cnt] = feat2[loc_idx][b][l];
                        yaxis2[cnt] = desc2[loc_idx][b][0];
                        cnt++;
                    }
                }
            }

            Plot p = new Plot("", "", "");
            p.setLimits(0, cnt-1, mu_ON, mu_OFF);
            p.addPoints(xaxis, yaxis1, Plot.X);
            p.draw();
            p.setColor(Color.RED);
            p.setLineWidth(3);
            p.addPoints(xaxis, yaxis2, Plot.LINE);
            is_out.addSlice("features, descriptors", p.getProcessor());

        }
        else {
            is_out.addSlice(new ByteProcessor(528,255));
        }
        return is_out;
    }

//        // plot fit scores for each branch + one plot with ratios
//        float[] xaxis1 = new float[M*L];
//        for (int i=0; i<M*L; i++) xaxis1[i] = i;
//        float[] thNCC = new float[M*L];
//        for (int i=0; i<M*L; i++) thNCC[i] = thresholdNCC;
//
//        for (int b= 0; b < 4; b++) {
//            Plot feature_plot = new Plot("", "#", "NCC");
//            feature_plot.setLimits(0, M*L-1, 0, 1);
//            feature_plot.addPoints(xaxis1, feat2[loc_idx][b], Plot.LINE);
//            feature_plot.draw();
//            feature_plot.setColor(Color.RED);
//            feature_plot.setLineWidth(3);
//            feature_plot.addPoints(xaxis1, thNCC, Plot.LINE);
//            is_out.addSlice("bch="+b, feature_plot.getProcessor());
//        }
//
//        return is_out;


    public static void getDelineationFeatures(int loc_idx, String mode) // feat2[locIdx] ratio2[locIdx]
    {
        if (loc_idx!=-1) {       // location is in foreground, there is a delineation there, fill 'features' and 'descriptors' up

            int[][] delin_at_loc = delin2[loc_idx];

            // loop branches twice:
            // 1 - to extract cross-profiles' geometry and to extract their fit scores (features)
            // 2 - to extract descriptors

            // 1
            for (int b = 0; b<delin_at_loc.length; b++) {           // loop 4 branches, b index defines the strength

                if (delin_at_loc[b][0]==-1) {
                    // whole streamline is missing - first one was -1, recursion was stopped
                    feat2[loc_idx][b] = null;  // float[] is null
                }
                else {

                    // there is at least one patch in the streamline
                    // loop the rest to count how many patches there are to allocate the array
                    int count_patches = 1;
                    for (int m=1; m<M; m++) {
                        if (delin_at_loc[b][m] != -1) count_patches++;
                        else break;
                    }

                    // allocate
                    feat2[loc_idx][b] = new float[count_patches*L];

                    for (int m = 0; m<M; m++) {      				// loop patches outwards, from the beginning

                        if (delin_at_loc[b][m]!=-1) {

                            // there is a patch, add the features to the matrix
                            int curr_i = delin_at_loc[b][m];
                            int curr_x = i2xy[curr_i][0];
                            int curr_y = i2xy[curr_i][1];

                            int prev_i, prev_x, prev_y;

                            if (m==0) {
                                prev_x = i2xy[loc_idx][0];  // cnetral locaiton
                                prev_y = i2xy[loc_idx][1];
                            }
                            else{
                                prev_i = delin_at_loc[b][m-1];
                                prev_x = i2xy[prev_i][0];
                                prev_y = i2xy[prev_i][1];
                            }

                            // get cross-profile values sampled from the local patch (aligned with the patch) store them in fitNCC
                            int init = m * L;
                            localPatchFeatures(prev_x, prev_y, curr_x, curr_y, mode, feat2[loc_idx][b], init); // will fill L values starting from init index

                        }
                        else break;

                    }

                }

            }

            // 2 calculate feature descriptors in every branch at loc_idx location
            ArrayList<Integer> branches_found = new ArrayList<Integer>(4);
            for (int b=0; b<delin_at_loc.length; b++) {

                if (feat2[loc_idx][b] != null) {

                    desc2[loc_idx][b]       = new float[2];
                    desc2[loc_idx][b][0]    = Stat.average(feat2[loc_idx][b]);
                    desc2[loc_idx][b][1]    = Stat.var(feat2[loc_idx][b], desc2[loc_idx][b][0]);

                    branches_found.add(b);

                }
                else { // it is null there is nothing

                    desc2[loc_idx][b] = null;

                }

            }

            // calculate fuzzy likelihoods
            if (branches_found.size()==4) {
                // desc2[loc][branch][mean value descriptor index]
                fuzzy_logic_system.critpointScores(
                        desc2[loc_idx][0][0],
                        desc2[loc_idx][1][0],
                        desc2[loc_idx][2][0],
                        desc2[loc_idx][3][0],
                        lhood2[loc_idx]);
            }
            else if (branches_found.size()==3) {
                fuzzy_logic_system.critpointScores(
                        desc2[loc_idx][branches_found.get(0)][0],
                        desc2[loc_idx][branches_found.get(1)][0],
                        desc2[loc_idx][branches_found.get(2)][0],
                        lhood2[loc_idx]);
            }
            else if (branches_found.size()==2) {
                fuzzy_logic_system.critpointScores(
                        desc2[loc_idx][branches_found.get(0)][0],
                        desc2[loc_idx][branches_found.get(1)][0],
                        lhood2[loc_idx]);
            }
            else if (branches_found.size()==1) {
                fuzzy_logic_system.critpointScores(
                        desc2[loc_idx][branches_found.get(0)][0],
                        lhood2[loc_idx]);
            }

        }

    }

	/*
		score calculation
	 */


	public static void exportFeats(String file_path)
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

                if (feat2[ii][jj] != null) {

                    for (int kk=0; kk<feat2[ii][jj].length; kk++) {
                        logWriter.print(String.format("%1.2f", feat2[ii][jj][kk]));
                        if (jj==feat2[ii].length-1 && kk==feat2[ii][jj].length-1) logWriter.print("\n");
                        else logWriter.print(", ");
                    }

                }

            }
        }

        logWriter.close();
	}

	public static void exportDescripts(String file_path)
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

        for (int ii=0; ii<desc2.length; ii++) {
            for (int jj=0; jj<desc2[ii].length; jj++) {

                if (desc2[ii][jj] != null) {

                    for (int kk=0; kk<desc2[ii][jj].length; kk++) {
                        logWriter.print(String.format("%1.2f, ", desc2[ii][jj][kk]));
                        //if (jj==desc2[ii].length-1) logWriter.print("\n");
                        //else logWriter.print(", ");
                    }

                }


            }
            logWriter.print("\n");
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

//        logWriter.println("strength: branch 0 > branch 1 > branch 2");
//        logWriter.println("score_at_point = median_along_line_ending at_point(or in the point neighbourhood) - background_estimate_at_point");
//
//        logWriter.println("feature 0: \tbranch 0 min score");
//        logWriter.println("feature 1: \tbranch 0 max score");
//
//        logWriter.println("feature 2: \tbranch 1 min score");
//        logWriter.println("feature 3: \tbranch 1 max score");
//
//        logWriter.println("feature 4: \tbranch 2 min score");
//        logWriter.println("feature 5: \tbranch 2 max score");
//
//        logWriter.println("feature 6: \tcenter       score");

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

        float R = .5f;

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

    // (CALC) NaN fitting scores for local patch cross profile
    private static void localPatchNaNScores(int branch_idx, int patch_idx, float[][] fit_scores_out)
    {
        for (int ii=0; ii<L; ii++) {
            fit_scores_out[branch_idx][patch_idx*L + ii] = Float.NaN;
        }
    }

    // (CALC) calculation of the fitting scores for every local patch cross profile
    private static void localPatchFeatures(float x1, float y1, float x2, float y2, String mode, float[] feats_out, int init_index) // L values per patch in feat_out
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

            // WARNING - tricky normalization if they are all completely the same
            if (Math.abs(val_max-val_min)<1) {

                // in case they are as close as 1 measure of 8bit level, then give the predefined profile
                // filled with zeros
                for (int iii = 0; iii<val.length; iii++){
                    val[iii] = 0f;
                }

            }
            else {
                // normalize min-max so that they're from 0 to 1
                for (int iii = 0; iii<val.length; iii++){
                    val[iii] = (val[iii]-val_min)/(val_max-val_min);
                }
            }

            // fit the normalized profile
            dummy = fitter.fit(val, mode);
            feats_out[init_index + ii] = dummy[1];

        }

    }

    // (CALC) calculation of the line segment geometry for every local patch cross profile
    // TODO gave this up, add it to misc
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
//        return geom;
//    }

    // method that calculates overlap of two line segments

}