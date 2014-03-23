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
public class Delineator2D extends Thread {

    private int begN, endN;

    // VARIABLES (mainly used as a table)
    public static int[][] 	    i2xy;                       // selected locations
    public static int[][]     	xy2i;                       // need for recursion

    // INPUT:
    public static int[][]       peaks_i;             	    // list of extracted peaks: N x 4 every FG location with 4 extracted peaks in indexed format
	public static int[][]		peaks_w;                    // weight assigned to each peak (for expansion)
	public static float[][]		inimg_xy;				    // input image (necessary for feature extraction)

    // PARAMETERS
    public static float     D;
    public static int       M;                              // how much it expands recursively from the center
    public static float     minCos;                         // allowed derail
    private static float    samplingStep;                   // when sampling image values to extract features
    private static int		L;                              // will define how many are taken along the diameter, in radial direction
    private static float    samplingStepLongitudinal;       // sampling along the streamline of patches
    private static int      dim;                            // cross profile length
	private static int      dim_half;                       // cross profile half-length
    private static String   mode;                           // fitting score
    private static int      nr_points = 4*1+1;              // number of the points given in fuzzy logic system
    private static float    mu_ON;
    private static float    sig_ON;
    private static float    mu_OFF;
    private static float    sig_OFF;

    // OUTPUT
    // associate the peaks and link follow-up points
    public static int[][][]     delin2;                     // N(foreground locs.) x 4(max. threads) x (1..M) (follow-up locs.) contains index for each location
    public static float[][][][] delin_refined_locs2;        // N(foreground locs.) x 4(max. threads) x 2 x ((1..M) x L) (2 dimensional)
    public static float[][][]   delin_refined_dirs2;        // N(foreground locs.) x 4(max. threads) x (2x2) [x0, y0, vx, vy]



//    public static float[][][][] delin_refined_vecs2;        // N(foreground locs.) x 4(max. threads) x 2 x ((1..M) x L) (2 dimensional)

	// FEATURES, F. DESCRIPTOR, OUT LIKELIHOODS
    public static float[][][]   offset2;						// N(foreground locs.) x 4 x ((1..M) x L) (offset)
//    // temporary
//    public static float[][][]   fitsco2;						// N(foreground locs.) x 4 x ((1..M) x L) (fit score)
//    public static float[][][]   varian2;						// N(foreground locs.) x 4 x ((1..M) x L) (variance)
//    public static float[][][]   desc2;                          // N(foreground locs.) x 4 x 3   (averages)
//    public static float[][]     lhood2;                         // N(foreground locs.) x 5 (NON..CRS) fuzzy logic output is stored here

    public Delineator2D(int n0, int n1)
    {
        this.begN = n0;
        this.endN = n1;
    }

    public static void loadTemplate(int[][] _i2xy, int[][] _xy2i, int[][] _peaks_i, int[][] _peaks_w, float[][] _inimg_xy,// byte[][] _backgr_xy,
									int     _M,
                                    float   _minCos,
                                    float   _D,
                                    float   _mu_ON,
                                    float   _mu_OFF,
                                    float   _sig_ON,
                                    float   _sig_OFF,
                                    int     _L,                         // will define sampling longitudinal
                                    float   _sampling_crosswise,
                                    String  _mode)
    {

        M               = _M;
        minCos          = _minCos;
        D               = _D;
        L               = _L;
        mu_ON           = _mu_ON;
        mu_OFF          = _mu_OFF;
        sig_ON          = _sig_ON;
        sig_OFF         = _sig_OFF;
        samplingStep    = _sampling_crosswise;
        mode            = _mode;

        dim_half = (int) Math.ceil( D / (samplingStep*2) );
        dim = 2*dim_half + 1;      // number of samples in a cross-profile

        samplingStepLongitudinal = D / (float)(L-1);

        i2xy = _i2xy;
        xy2i = _xy2i;
        peaks_i = _peaks_i;
		peaks_w = _peaks_w;
		inimg_xy = _inimg_xy;

        // allocate output
        delin2 				= new int[i2xy.length][4][M];
        // initialize with -2 (same meaning as with peaks_i)
        for (int i=0; i<delin2.length; i++)
            for (int j=0; j<delin2[i].length; j++)
                for (int k=0; k<delin2[i][j].length; k++)
                    delin2[i][j][k] = -2;

		delin_refined_locs2 = new float[i2xy.length][4][2][];       // (x,y)
        delin_refined_dirs2 = new float[i2xy.length][4][];          // (x,y)
//		delin_refined_vecs2 = new float[i2xy.length][4][2][];       // (x,y)

        offset2   = new float[i2xy.length][4][];                    // offset
//        // temporary
//        fitsco2   = new float[i2xy.length][4][];                    // fit score
//        varian2   = new float[i2xy.length][4][];                    // variance
//
//        desc2  	= new float[i2xy.length][4][];                      // description of the fit scores
//        lhood2  = new float[i2xy.length][];     					// fuzzy likelihood output

    }

    public void run()
    {
        for (int locationIdx = begN; locationIdx < endN; locationIdx++) {

			/*
				processing classes for this run
				they will interfere with outputs
				need to be allocated independently
				for this run interval
			*/

			Fitter1D fitter  			= new Fitter1D(dim, false); // dim = profile width with current samplingStep, verbose = false
//			fitter.showTemplates();
            /*
            central FLS
             */
			Fuzzy2D  fls 				= new Fuzzy2D(nr_points, mu_ON, sig_ON, mu_OFF, sig_OFF); // class used to calculate fuzzy score
            /*
            stream FLS
             */

			// todo add another Fuzzy class for the cascade
//            SegmentFLS
//			fuzzy_logic_system.showFuzzification();
//			fuzzy_logic_system.showDefuzzification();

            int[] peaks_at_loc = peaks_i[locationIdx];

            // access individual peaks at this point (-2:not exist, -1:background, >=0:foreground)

            /*
                delin2[locationIdx] calculation
             */
            for (int pp = 0; pp<peaks_i[locationIdx].length; pp++) {  // loop 4 allocated branches (1st generation)

                if (peaks_i[locationIdx][pp] >= 0) {

                    //if the peak exists in the foreground
                    delin2[locationIdx][pp][0] = peaks_i[locationIdx][pp];

                    int curr_index, prev_index, next_index;
                    curr_index = peaks_i[locationIdx][pp];
                    prev_index = locationIdx;

                    for (int m=1; m<M; m++) { // follow the recursion for the rest of the indexes

                        next_index = getNext(prev_index, curr_index); // recursion : prev+curr->next index
                        // >=0 if there was at least one in foreground
                        // will give -1 if the next one was in the background and there was no other foreground peak
                        // -2 if it there was no other foreground peak and background either

                        if (next_index>=0) {
                            delin2[locationIdx][pp][m] = next_index;     // store it in output matrix
                        }
                        else if (next_index==-1){ // not found but found in the background
                            delin2[locationIdx][pp][m] = -1;
                            // fill the rest with -1 once one fell out
                            for (int m_aux=m+1; m_aux<M; m_aux++) delin2[locationIdx][pp][m_aux] = -1;
                            break; // stop looping further along the streamline here
                        }
                        else if (next_index==-2) { // not found at all
                            delin2[locationIdx][pp][m] = -2;
                            // fill the rest with -2 once one had a dead end
                            for (int m_aux=m+1; m_aux<M; m_aux++) delin2[locationIdx][pp][m_aux] = -2;
                            break; // stop looping further along the streamline here
                        }

                        // recursion - to expand further in the next iteration
                        prev_index = curr_index;
                        curr_index = next_index;

                    }

                }
                else if (peaks_i[locationIdx][pp] == -1) {

                    // 1st generation peak index was -1 (peak in the background according to the mask)
                    // this is a case of an incomplete delineation, having it's parts in the background
                    // it can be anything and mask should not influence the final decision making
                    // if the peak was discarded as it did not exist - wrong it can affect the decision
                    // solution: keep only those delineations unaffected by the mask, where the whole delineation is in foreground
                    // assign the whole delineation as incomplete here - giving -1 to all elements and finish
					for (int m=0; m<M; m++) delin2[locationIdx][pp][m] = -1;

                }
                else if (peaks_i[locationIdx][pp] == -2) {

                    // those are the places that were not filled up with peaks
                    // propagate it to delin2
                    for (int m=0; m<M; m++) delin2[locationIdx][pp][m] = -2;

                }

            }
            /*
             delin2[locationIdx] is fully formed
             */

            /*
            examples delin2[locationIdx][][] indexes start from the one defined with the first index
            (M=3):
            23  34 56
            678 -1 -1
            -2  -2 -2
            --------------
            23  -2 -2      -> still extract the features
            678 -1 -1      -> still extract the features
            -2  -2 -2
            --------------
            23  -2 -2
            678 -1 -1
            -1  -1 -1      -> peak in background from the beginning (don't use it, unsure, skip this case in detection)
            --------------
            23  90 101    |
            -2  -2 -2     |-> endpoint candidate
            -2  -2 -2     |
             */

            /*
                delin_refined_locs2[locationIdx], offset2[locationIdx], delin_refined_dirs2[locationIdx]
            */

            for (int b = 0; b<delin2[locationIdx].length; b++) { // loop 4 branches

                if (delin2[locationIdx][b][0]==-1) {
                    // whole streamline is missing: the first one was -1, recursion was stopped
                    delin_refined_locs2[locationIdx] = null;
                    delin_refined_dirs2[locationIdx] = null;
                    offset2[locationIdx] = null;
                    break; // stop looping branches further
                }
                else if (delin2[locationIdx][b][0]==-2) {
                    // no streamline here
                    delin_refined_locs2[locationIdx][b] = null;
                    delin_refined_dirs2[locationIdx][b] = null;
                    offset2[locationIdx][b] = null;
                }
                else if (delin2[locationIdx][b][0]>=0) {
                    // there is at least one patch in the streamline
                    // loop the rest to count how many patches there are to allocate the array
                    int count_patches = 1;
                    for (int m=1; m<M; m++) {
                        if (delin2[locationIdx][b][m] >= 0) count_patches++;
                        else break; // because the rest are filled up with -1 or -2 anyway
                    }

                    // allocate
                    delin_refined_locs2[locationIdx][b][0] = new float[count_patches*L]; // x coordinates allocate
                    delin_refined_locs2[locationIdx][b][1] = new float[count_patches*L]; // y coordinates allocate

                    delin_refined_dirs2[locationIdx][b]    = new float[4];               // x0,y0,vx,vy

                    offset2[locationIdx][b]                = new float[count_patches*L]; //

                    // fill  up the allocated arrays
                    for (int m = 0; m<M; m++) {      				// loop patches outwards, from the beginning

                        if (delin2[locationIdx][b][m]>=0) {

                            // there is a patch, add the features to the matrix
                            int curr_i = delin2[locationIdx][b][m];
                            int curr_x = i2xy[curr_i][0];
                            int curr_y = i2xy[curr_i][1];

                            int prev_i, prev_x, prev_y;

                            if (m==0) {
                                prev_x = i2xy[locationIdx][0];  // central location
                                prev_y = i2xy[locationIdx][1];
                            }
                            else{
                                prev_i = delin2[locationIdx][b][m-1];
                                prev_x = i2xy[prev_i][0];
                                prev_y = i2xy[prev_i][1];
                            }

                            // get refined locations sampled from the local patch (aligned with the patch)
                            if (m==0) { // first patch from the center
                                localPatchRefined(
                                        prev_x, prev_y, curr_x, curr_y,
                                        delin_refined_locs2[locationIdx][b],
                                        offset2[locationIdx][b],
                                        m*L,
                                        delin_refined_dirs2[locationIdx][b]);
                            }
                            else {
                                localPatchRefined(
                                        prev_x, prev_y, curr_x, curr_y,
                                        delin_refined_locs2[locationIdx][b],
                                        offset2[locationIdx][b],
                                        m*L);
                            }


                        }
                        else break; // because the rest are filled up with -1 or -2 anyway

                    }



                }

            }
            /*
                delin_refined_locs2[locationIdx], offset2[locationIdx] formed
            */

            /*
                delin_refined_dirs2[locationIdx]  using delin_refined_locs2[locationIdx]
            */



            /*
                delin_refined_dirs2[locationIdx] formed
            */



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

                    if (delin2[atLoc][ii][jj]==-1) {
                        printout += "BGRD";
                    }
                    else if (delin2[atLoc][ii][jj]==-2) {
                        printout += "NONE";
                    }
                    else {
                        printout += IJ.d2s(delin2[atLoc][ii][jj], 2);
                    }

                    if (jj==delin2[atLoc][ii].length-1) printout += "\n";
                    else printout += ",  ";
                }
            }

            printout += "\nREFINED LOCS:\n";
            if (delin_refined_locs2[atLoc]!=null) {

                for (int b=0; b<delin_refined_locs2[atLoc].length; b++) {

                    printout += b+"\t->\t";

                    if (delin_refined_locs2[atLoc][b]!=null) {

                        for (int l=0; l<delin_refined_locs2[atLoc][0].length; l++) {

                            printout += "("+IJ.d2s(delin_refined_locs2[atLoc][b][0][l], 2)+", "+IJ.d2s(delin_refined_locs2[atLoc][b][1][l], 2)+")";

                            if (l==delin_refined_locs2[atLoc][b][0].length-1) printout += "\n";
                            else printout += ", ";
                        }

                    }
                    else {
                        printout += "NONE\n";
                    }

                }

            }
            else {
                printout += "SKIPPED CALCULATING HERE (THERE WAS A THREAD POINTING TO BGRD)\n";
            }


//            // skip descriptors and likelihoods for now
//            printout += "\nDESCRIPTORS:\n";
//            for (int b=0; b<desc2[atLoc].length; b++) {
//                printout += (b+1) +"\t->\t"; //+ IJ.d2s(ratio2[atLoc][ii], 2) + "\n"
//
//                if (desc2[atLoc][b] != null) {
//
//                    for (int l=0; l<desc2[atLoc][b].length; l++) {
//                        printout += IJ.d2s(desc2[atLoc][b][l], 2);
//                        if (l==desc2[atLoc][b].length-1) printout += "\n";
//                        else printout += ",\t";
//                    }
//
//                }
//                else {
//                    printout += "NONE\n";
//                }
//
//            }
//
//            printout += "\nFUZZY LIKELIHOODS:\n";
//            printout += "NON->" + IJ.d2s(lhood2[atLoc][0], 2) + "\t\t";
//            printout += "END->" + IJ.d2s(lhood2[atLoc][1], 2) + "\t\t";
//            printout += "BDY->" + IJ.d2s(lhood2[atLoc][2], 2) + "\t\t";
//            printout += "BIF->" + IJ.d2s(lhood2[atLoc][3], 2) + "\t\t";
//            printout += "CRS->" + IJ.d2s(lhood2[atLoc][4], 2) + "\t\t";

            IJ.log(printout);

        }
		else {
			IJ.log("background point, no data!");
		}

	}

    private static int getNext(int prev_index, int curr_index)
    {

        // these are stacked as XY - take care on that
        int prevX = i2xy[prev_index][0];    // X
        int prevY = i2xy[prev_index][1];    // Y

        int currX = i2xy[curr_index][0];    // X
        int currY = i2xy[curr_index][1];    // Y

        // check peaks at curr
        int[]   pks4xI  = peaks_i[curr_index]; // list all indexes
		int[] 	pks4xW	= peaks_w[curr_index];

		// take the one with highest weight that is above minCos
		int next_index 	= -2;
		int next_weight = Integer.MIN_VALUE;
        boolean found_valid_followup = false;

		for (int p = 0; p<pks4xI.length; p++) {

			if (pks4xI[p]>=0) {
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
                        found_valid_followup = true;
						next_weight = pks4xW[p];
						next_index 	= pks4xI[p];
					}
				}

			}
			else {
				// peak in the background and there is nothing acceptable in foreground so far
                if (!found_valid_followup) {
                    next_index = -1;
                    //next_weight neutral
                }
                else {
                    // there is a follow up to use, ignore background point
                }

			}
		}

		return next_index;
		// will give -1 if the next one was in the background and there was no other foreground peak
        // -2 if it there was no other foreground peak and background either
        // >=0 if there was at least one in foreground

	}

    /*
        outputs (to visualize delineation and extracted features)
     */
    public static Overlay getDelineationOverlay(int atX, int atY)
    {

        // return the delineated local structure Overlay for visualization
        Overlay ov = new Overlay();

        /*
            adds generations of peaks in green
         */
        PeakExtractor2D.getPeaks(atX, atY, 3, ov);   // peaks will be stored in ov

        float Rd = 1.5f; // radius of the circles written for delineation
        Color cd = Color.RED;
        float wd = 0.25f;

		/*
		 central location
		  */
        OvalRoi ovalroi = new OvalRoi(atX-(Rd/2)+.5f, atY-(Rd/2)+.5f, Rd, Rd);
		ovalroi.setFillColor(cd);
		ovalroi.setStrokeWidth(wd);
        ov.add(ovalroi);

        /*
            add local cross-section locations (re-calculate) and patch delineation spots (delin2)
         */
        if (xy2i[atX][atY]>=0) {

            int[][] delin_at_loc = delin2[xy2i[atX][atY]];

            for (int b = 0; b<delin_at_loc.length; b++) {           // loop 4 branches

                for (int m=0; m<M; m++) {

                    if (delin_at_loc[b][m] >= 0) {

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

						ArrayList<OvalRoi> line_pts = localPatchCrossProfilesLocs(prev_x, prev_y, pt_x, pt_y);  // will be recalculated
						for (int aa = 0; aa < line_pts.size(); aa++) ov.add(line_pts.get(aa));

                    }

                }

            }

        }

        /*
        add refined streamline locations (read from delin_refined_locs2)
         */
        if (xy2i[atX][atY]>=0) {

            float[][][] local_refined_locs = delin_refined_locs2[xy2i[atX][atY]];

            for (int b = 0; b<local_refined_locs.length; b++) {

                // check whether there is a refinement at all here, add in case there is
                if (local_refined_locs[b]!=null) {

                    // loop points to add them
                    for (int l=0; l<local_refined_locs[b][0].length; l++) {

                        float refined_x = local_refined_locs[b][0][l];
                        float refined_y = local_refined_locs[b][1][l];
                        ovalroi = new OvalRoi(refined_x-(Rd/2)+.5f, refined_y-(Rd/2)+.5f, Rd, Rd);
                        ovalroi.setFillColor(Color.YELLOW);
                        ovalroi.setStrokeWidth(Rd/2);
                        ov.add(ovalroi);

                    }

                }

            }

        }

        return ov;

    }


    /*
    visualize the extracted images of square patches in separate frames
     */
    public static ImageStack getDelineationPatches(int atX, int atY)
	{

        // create new ImageStack with every layer corresponding to one streamline
        ImageStack isOut = new ImageStack(dim, M*dim);

        int idx = xy2i[atX][atY]; // read extracted peaks at this location

        if (idx>=0) { // index is foreground

            int[][] delin_at_loc = delin2[idx];

            for (int b = 0; b<delin_at_loc.length; b++) {           // loop 4 branches

					float[] concat = new float[M*dim*dim];
					int cnt = 0;

                    for (int m = 0; m<M; m++) { // outwards loop

                        if (delin_at_loc[b][m] >= 0) { // add those that exist, the rest stays blank

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

                    }

					// add slice after looping
					FloatProcessor patch = new FloatProcessor(dim, M*dim, concat);
					isOut.addSlice("bch "+b+")", patch);

            }

        }

        if (isOut.getSize()==0) {
            isOut.addSlice(new FloatProcessor(dim, dim));//(patch_size,patch_size));
        }

        return isOut;

    }

	public static ImageStack plotDelineationProfiles(int atX, int atY)
    {


		/*
			processing unit
		 */

		Fitter1D fitter = new Fitter1D(dim, false);
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

    /*
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

            float[] xaxis       = new float[totalFeats];
            float[] yaxis1      = new float[totalFeats];
            float[] yaxis2      = new float[totalFeats];
            float[] yaxisON     = new float[totalFeats];
            float[] yaxisOFF    = new float[totalFeats];

            int cnt = 0;
            for (int b=0; b<feat2[loc_idx].length; b++) {
                if (feat2[loc_idx][b] != null) {
                    for (int l=0; l<feat2[loc_idx][b].length; l++) {

                        if (l==0) {
                            if (cnt==0) {
                                xaxis[cnt]  = 0;
                            }
                            else {
                                xaxis[cnt]  = xaxis[cnt-1] + 20;
                            }
                        }
                        else {
                            xaxis[cnt] = xaxis[cnt-1] + 1;
                        }

                        yaxis1[cnt]     = feat2[loc_idx][b][l]; // features
                        yaxis2[cnt]     = desc2[loc_idx][b][0]; // mean estimate
                        yaxisON[cnt]    = mu_ON;
                        yaxisOFF[cnt]   = mu_OFF;
                        cnt++;

                    }
                }
            }

            // find limits for the features
            float min_feat = Float.MAX_VALUE;
            float max_feat = Float.NEGATIVE_INFINITY;
            float max_axis = Float.NEGATIVE_INFINITY;

            for (int ii=0; ii<yaxis1.length; ii++) {
                if (yaxis1[ii]<min_feat) {
                    min_feat = yaxis1[ii];
                }
                if (yaxis1[ii]>max_feat) {
                    max_feat = yaxis1[ii];
                }
                if (xaxis[ii]>max_axis) {
                    max_axis = xaxis[ii];
                }
            }

            Plot p = new Plot("", "", "");
            p.setLimits(0, max_axis, Math.min(min_feat, mu_ON-.1f), Math.max(max_feat, mu_OFF+.1f));
            p.addPoints(xaxis, yaxis1, Plot.BOX);
            p.draw();

            p.setLineWidth(4);
            p.setColor(Color.GREEN);
            p.addPoints(xaxis, yaxisON, Plot.DOT);
            p.draw();

            p.setLineWidth(4);
            p.setColor(Color.BLUE);
            p.addPoints(xaxis, yaxisOFF, Plot.DOT);
            p.draw();

            p.setLineWidth(4);
            p.setColor(Color.RED);
            p.addPoints(xaxis, yaxis2, Plot.DOT);
            p.draw();

            is_out.addSlice("features, descriptors", p.getProcessor());

        }
        else {
            is_out.addSlice(new ByteProcessor(528,255));
        }
        return is_out;
    }
    */

    /*
    private static void getDelineationFeatures(int loc_idx, String mode, Fitter1D _fitter, Fuzzy2D _fls)
    {
        //if (loc_idx!=-1) {       // location is in foreground, there is a delineation there, fill 'features' and 'descriptors' up

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
                                prev_x = i2xy[loc_idx][0];  // central location
                                prev_y = i2xy[loc_idx][1];
                            }
                            else{
                                prev_i = delin_at_loc[b][m-1];
                                prev_x = i2xy[prev_i][0];
                                prev_y = i2xy[prev_i][1];
                            }

                            // get cross-profiles sampled from the local patch (aligned with the patch), and their features
                            int init = m * L;
                            localPatchFeatures(prev_x, prev_y, curr_x, curr_y, mode, feat2[loc_idx][b], init, _fitter); // will fill L values starting from init index

                        }
                        else break;

                    }

                }

            }

            // 2 calculate feature descriptors in every branch at loc_idx location
            ArrayList<Integer> branches_found = new ArrayList<Integer>();
			branches_found.clear();

            for (int b=0; b<feat2[loc_idx].length; b++) {

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

			lhood2[loc_idx] = new float[_fls.L]; // allocate in the length of fls output vector

            // calculate fuzzy likelihoods
            if (branches_found.size()==4) {
                // desc2[loc][branch][mean value descriptor index]
                _fls.critpointScores(
											desc2[loc_idx][0][0],
											desc2[loc_idx][1][0],
											desc2[loc_idx][2][0],
											desc2[loc_idx][3][0],
											lhood2[loc_idx]);
            }
            else if (branches_found.size()==3) {
                _fls.critpointScores(
                        desc2[loc_idx][branches_found.get(0)][0],
                        desc2[loc_idx][branches_found.get(1)][0],
                        desc2[loc_idx][branches_found.get(2)][0],
                        lhood2[loc_idx]);
            }
            else if (branches_found.size()==2) {
                _fls.critpointScores(
                        desc2[loc_idx][branches_found.get(0)][0],
                        desc2[loc_idx][branches_found.get(1)][0],
                        lhood2[loc_idx]);
            }
            else if (branches_found.size()==1) {
                _fls.critpointScores(
                        desc2[loc_idx][branches_found.get(0)][0],
                        lhood2[loc_idx]);
            }

        //}

    }
    */

    /*
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
                        logWriter.print(String.format("%1.6f", feat2[ii][jj][kk]));
                        if (kk<feat2[ii][jj].length-1) logWriter.print(", ");
                        else logWriter.print("\n");
                    }

                }
				else {
					logWriter.print("null\n");
				}
            }
        }

        logWriter.close();
	}
	*/


    /*
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
                        logWriter.print(String.format("%1.6f, ", desc2[ii][jj][kk]));
                        if (kk<desc2[ii][jj].length-1) logWriter.print(", ");
                        else logWriter.print("\n");
                    }

                }
				else {
					logWriter.print("null\n");
				}
            }
        }

        logWriter.close();
	}
	*/

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
                    if (kk<delin2[ii][jj].length-1) logWriter.print(", ");
                    else logWriter.print("\n");
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
	public static void exportLikelihoods(String file_path)
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
		for (int ii=0; ii<lhood2.length; ii++) {
			for (int jj=0; jj<lhood2[ii].length; jj++) { // loop streamlines
					logWriter.print(String.format("%1.8f", lhood2[ii][jj]));
					if (jj<lhood2[ii].length-1) logWriter.print(", ");
					else logWriter.print("\n");
			}
		}

		logWriter.close();

	}
	*/

    /*
    public static ImageStack exportLikelihoods(int[] choice)
    {

        int w = inimg_xy.length;
        int h = inimg_xy[0].length;

        ImageStack is_out = new ImageStack(w, h);

        for (int i=0; i<choice.length; i++) {
            if (choice[i] >=0 && choice[i]<lhood2[0].length) {

                float[][] t = new float[w][h];

                for (int ii = 0; ii<lhood2.length; ii++) {

                    int x = i2xy[ii][0];
                    int y = i2xy[ii][1];

                    boolean isMax = true;

					for (int k=0; k<lhood2[0].length; k++) {
                        if (k!=choice[i] && lhood2[ii][k]>lhood2[ii][choice[i]]) {
                            //isMax = false;
                        }
                    }

                    t[x][y] = (isMax)? lhood2[ii][choice[i]] : 0;

                }

                is_out.addSlice(""+IJ.d2s(choice[i], 0), new FloatProcessor(t));
            }
        }

        return is_out;
    }
    */

    /*
        methods that deal with local image patch - (rectangle defined with prev_x,y and curr_x,y)
     */

    // (VIZ) extract array with local patch sampled values
    private static float[] localPatchVals(float x1, float y1, float x2, float y2)
    {

        float l = (float) Math.sqrt(Math.pow(x2-x1, 2)+Math.pow(y2-y1, 2));
        float vx = (x2-x1)/l;
        float vy = (y2-y1)/l;
        float wx = vy;
        float wy = -vx;

        float x_root = x2 - vx * D;
        float y_root = y2 - wx * D;

        float[] vals = new float[dim*dim];
        int cnt = 0;

        for (int ii=0; ii<=2*dim_half; ii++) { // loops vector v
            for (int jj=-dim_half; jj<=dim_half; jj++) { // loops vector w

                float curr_x = x_root + ii * samplingStep * vx + jj * samplingStep * vy;
                float curr_y = y_root + ii * samplingStep * wx + jj * samplingStep * wy;

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

        float x_root = x2 - vx * D;
        float y_root = y2 - wx * D;

        ArrayList<float[]> vals            = new ArrayList<float[]>();

        for (int ii=0; ii<L; ii++) { // loops L of them in radial direction

            float[] val = new float[dim];
            int cnt = 0;
			float val_min = Float.POSITIVE_INFINITY;
			float val_max = Float.NEGATIVE_INFINITY;

            for (int jj=-dim_half; jj<=dim_half; jj++) { // loops vector w

                float curr_x = x_root + ii * samplingStepLongitudinal * vx + jj * samplingStep * vy;
                float curr_y = y_root + ii * samplingStepLongitudinal * wx + jj * samplingStep * wy;

                val[cnt] = Interpolator.interpolateAt(curr_x, curr_y, inimg_xy);

				if (val[cnt]>val_max) val_max = val[cnt];
				if (val[cnt]<val_min) val_min = val[cnt];

                cnt++;
            }

            if (Math.abs(val_max-val_min)<Float.MIN_VALUE) Arrays.fill(val, 0f);
            else {
                for (int iii = 0; iii<val.length; iii++) val[iii] = (val[iii] - val_min) / (val_max - val_min);
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

        float x_root = x2 - vx * D;
        float y_root = y2 - wx * D;

        float R = samplingStep/2;

		ArrayList<OvalRoi> pts = new ArrayList<OvalRoi>(dim*L);

		for (int ii=0; ii<L; ii++) { // loops L of them in radial direction

			for (int jj=-dim_half; jj<=dim_half; jj++) { // loops vector w

                float curr_x = x_root + ii * samplingStepLongitudinal * vx + jj * samplingStep * vy;
				float curr_y = y_root + ii * samplingStepLongitudinal * wx + jj * samplingStep * wy;

				OvalRoi pt = new OvalRoi(curr_x+.5f-(R/2), curr_y+.5f-(R/2), R, R);
				pt.setFillColor(Color.BLUE);
				pt.setStrokeWidth(R/2);
                pts.add(pt);

			}

		}

		return pts;

	}

    /*
        refined local directions per profile
     */
    private static void localPatchRefinedVecs(float[][] refined_centerline_locs_xy, float[][] refined_centerline_vecs_xy)
    {

        // will do local average of directions (average 2 consecutive directions)

        for (int ii=0; ii<refined_centerline_locs_xy.length; ii++) {



        }


    }

    /*
        refined locations per cross profile
     */
    private static void localPatchRefined(float x1, float y1, float x2, float y2,
                                          float[][] refined_centerline_locs_xy,
                                          float[] refined_offsets,
                                          int init_index,
                                          float[] refined_centerline_dir)
    {
        /*
            standard way to loop through a patch defined with 2 points
         */
        float l = (float) Math.sqrt(Math.pow(x2-x1, 2)+Math.pow(y2-y1, 2));
        float vx = (x2-x1)/l;
        float vy = (y2-y1)/l;
        float wx = vy;
        float wy = -vx;

        float x_root = x2 - vx * D;
        float y_root = y2 - wx * D;

        for (int ii=0; ii<L; ii++) {                        // loops L of them in radial direction

            // precompute the profile in advance to seek its local maximum
            float[] patch_profile = new float[2*dim_half+1 +2 +2];

            patch_profile[0] = Float.NEGATIVE_INFINITY;
            patch_profile[1] = Float.NEGATIVE_INFINITY;

            patch_profile[patch_profile.length-1] = Float.NEGATIVE_INFINITY;
            patch_profile[patch_profile.length-2] = Float.NEGATIVE_INFINITY;

            for (int jj=-dim_half; jj<=dim_half; jj++) {
                float curr_x = x_root + ii * samplingStepLongitudinal * vx + jj * samplingStep * vy;
                float curr_y = y_root + ii * samplingStepLongitudinal * wx + jj * samplingStep * wy;
                patch_profile[jj+dim_half+2] = Interpolator.interpolateAt(curr_x, curr_y, inimg_xy);
            }

            // find local max and take the one with min dist towards center
            // initialize
            float refined_x = x_root + ii * samplingStepLongitudinal * vx + (-dim_half-1) * samplingStep * vy;
            float refined_y = y_root + ii * samplingStepLongitudinal * wx + (-dim_half-1) * samplingStep * wy;
            float offset_min = Math.abs(1 - (patch_profile.length-1)/2) * samplingStep;

            for (int loop_prof=2; loop_prof<patch_profile.length-2; loop_prof++) {
                boolean is_local_max =
                                patch_profile[loop_prof]>=patch_profile[loop_prof-1] &&
                                patch_profile[loop_prof]>=patch_profile[loop_prof-2] &&
                                patch_profile[loop_prof]>=patch_profile[loop_prof+1] &&
                                patch_profile[loop_prof]>=patch_profile[loop_prof+2];
                if (is_local_max) {
                    // measure dist
                    float dist_from_center = Math.abs(loop_prof - (patch_profile.length - 1) / 2) * samplingStep;
                    if (dist_from_center<offset_min) {
                        offset_min = dist_from_center;

                        int index_from_center = loop_prof-dim_half-2;
                        refined_x = x_root + ii * samplingStepLongitudinal * vx + index_from_center * samplingStep * vy;
                        refined_y = y_root + ii * samplingStepLongitudinal * wx + index_from_center * samplingStep * wy;

                    }
                }
            }

            refined_centerline_locs_xy[0][init_index + ii] = refined_x;
            refined_centerline_locs_xy[1][init_index + ii] = refined_y;
            refined_offsets[init_index + ii] = offset_min;

        }

        /*
        calculate the direction based on fitting a line through calculated refined_centerline_locs_xy
         */
        float l1x, l1y, l2x, l2y;
        if (refined_centerline_dir!=null) {
            for (int jj=-dim_half; jj<=dim_half; jj++) { // loop first cross-section points
                for (int kk=-dim_half; kk<=dim_half; kk++) { // loop last cross-section points

                    //

                }


            }
        }


    }

    /*
        refined locations per cross profile with outward direction estimation, based on the first patch refined points
     */
    private static void localPatchRefined(float x1, float y1, float x2, float y2,
                                          float[][] refined_centerline_locs_xy,
                                          float[] refined_offsets,
                                          int init_index,
                                          float[] refined_dir)
    {


    }

    // (CALC) calculation of the fitting scores for every local patch cross profile
    private static void localPatchFeatures(float x1, float y1, float x2, float y2, String mode, float[] feats_out, int init_index, Fitter1D _fitter) // L values per patch in feat_out
    {

        float[] dummy;
        float[] val = new float[dim];

        /*
            standard scheme to loop through the patch
         */

        float l = (float) Math.sqrt(Math.pow(x2-x1, 2)+Math.pow(y2-y1, 2));
        float vx = (x2-x1)/l;
        float vy = (y2-y1)/l;
        float wx = vy;
        float wy = -vx;

        float x_root = x2 - vx * D;
        float y_root = y2 - wx * D;

        for (int ii=0; ii<L; ii++) { // loops L of them in radial direction

            int cnt = 0;
            float val_min = Float.POSITIVE_INFINITY;
            float val_max = Float.NEGATIVE_INFINITY;

            for (int jj=-dim_half; jj<=dim_half; jj++) { // loops vector w

                float curr_x = x_root + ii * samplingStepLongitudinal * vx + jj * samplingStep * vy;
                float curr_y = y_root + ii * samplingStepLongitudinal * wx + jj * samplingStep * wy;

                val[cnt] = Interpolator.interpolateAt(curr_x, curr_y, inimg_xy);

                if (val[cnt]>val_max) val_max = val[cnt];
                if (val[cnt]<val_min) val_min = val[cnt];

                cnt++;
            }

            if (Math.abs(val_max-val_min)<Float.MIN_VALUE) Arrays.fill(val, 0f);
            else {
                for (int iii = 0; iii<val.length; iii++) val[iii] = (val[iii] - val_min) / (val_max - val_min);
            }

            // fit the normalized profile
            dummy = _fitter.fit(val, mode);
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