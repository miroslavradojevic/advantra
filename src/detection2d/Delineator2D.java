package detection2d;

import aux.Geometry;
import aux.Stat;
import aux.Interpolator;
import fit.Fitter1D;
import ij.IJ;
import ij.ImageStack;
import ij.gui.*;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;

import java.awt.*;
import java.awt.geom.Arc2D;
import java.awt.geom.CubicCurve2D;
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
	public static float[][] 	peaks_lhood;				// likelihoods for each peaks taken from the profile
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
	public static float[][][][] delin_refined_vecs2;        // N(foreground locs.) x 4(max. threads) x 2 x ((1..M) x L) (2 dimensional)
	public static float[][][]   delin_refined_dirs2;        // N(foreground locs.) x 4(max. threads) x (2x2) [x0, y0, x1, y1]
    public static float[][][]   offset2;				    // N(foreground locs.) x 4(max. threads) x (L) (offsets from the first patch only)
    // features
    public static float[][][]   fitsco2;					// N(foreground locs.) x 4(max. threads) x ((1..M) x L) (fit scores for the refined locs)
    public static float[][][]   stdev2;						// N(foreground locs.) x 4(max. threads) x ((1..M) x L) (variances of real values)
    // descriptions (calculations on features)
    public static int           NR_DESC2 = 4;               // number of descriptors (fitsco2 mean, fitsco2 median, stdev2 mean, stdev2 median)
    public static float[][]     desc2;                      // N(foreground locs.) x 4(max. threads) x 4 (nr.desc)
	public static float[][]		peaks_lhood_reg; 			// (normalized actually now 0-1)regularized likelihoods, using mean ncc stored in desc2
	// fuzzy system outputs
    public static float[][]     streamline_score2;          // N(foreground locs.) x 4(max. threads)
    public static float[][]     endpoint_score2;            // by fuzzy logic
    public static float[][]     junction_score2;            // at least two or three

    public Delineator2D(int n0, int n1)
    {
        this.begN = n0;
        this.endN = n1;
    }

    public static void loadTemplate(int[][] _i2xy, int[][] _xy2i, int[][] _peaks_i, int[][] _peaks_w, float[][] _peaks_lhood, float[][] _inimg_xy,// byte[][] _backgr_xy,
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
		peaks_lhood = _peaks_lhood;
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
		delin_refined_vecs2 = new float[i2xy.length][4][2][];       // (x,y)
        offset2             = new float[i2xy.length][4][];          // offset from the patch center
        fitsco2   			= new float[i2xy.length][4][];          // fit score
        stdev2   			= new float[i2xy.length][4][];          // variance

        desc2  				= new float[i2xy.length][];          // description of the fit scores
		peaks_lhood_reg		= new float[i2xy.length][];
        streamline_score2  	= new float[i2xy.length][];     		// critpoint likelihood output
        endpoint_score2 = new float[i2xy.length][];
        junction_score2 = new float[i2xy.length][];


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

			// aux variables used allocation
			Fitter1D fitter  			= new Fitter1D(dim, false); // dim = profile width with current samplingStep, verbose = false
//			fitter.showTemplates();
            /*
            central FLS
             */
			FuzzyPoint2D  fls 				= new FuzzyPoint2D(nr_points, mu_ON, sig_ON, mu_OFF, sig_OFF); // class used to calculate fuzzy score
            /*
            stream FLS
             */

//            SegmentFLS
//			fls.showFuzzification();
//			fls.showDefuzzification();

            int[] peaks_at_loc = peaks_i[locationIdx];

            // access individual peaks at this point (-2:not exist, -1:background, >=0:foreground)

            /* delin2[locationIdx] calculation */
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

            /* delin_refined_locs2[locationIdx], offset2[locationIdx], delin_refined_dirs2[locationIdx] */
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
                    delin_refined_dirs2[locationIdx][b]    = new float[4];               // x0,y0,x1,y1
                    offset2[locationIdx][b]                = new float[L];               // offsets along the patch

                    // fill  up the allocated arrays for each patch
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
                                localPatchRefinedLocs(
                                        prev_x, prev_y, curr_x, curr_y,
                                        m * L,                                    // start from
                                        delin_refined_locs2[locationIdx][b],
                                        offset2[locationIdx][b],
                                        delin_refined_dirs2[locationIdx][b]);
                            }
                            else {
                                localPatchRefinedLocs(
                                        prev_x, prev_y, curr_x, curr_y,
                                        m * L,
                                        delin_refined_locs2[locationIdx][b],
                                        null,                               // no offsets here
                                        null);                              // no refined dirs
                            }

                        }
                        else break; // because the rest are filled up with -1 or -2 anyway

                    }

                }

            }

            /* calculate delin_refined_vecs2[locationIdx]  using delin_refined_locs2[locationIdx]*/

            if (delin_refined_locs2[locationIdx]!=null) {
                for (int b = 0; b<delin_refined_locs2[locationIdx].length; b++) {
                    if (delin_refined_locs2[locationIdx][b]!=null) {

                        int to_allocate = delin_refined_locs2[locationIdx][b][0].length;
                        delin_refined_vecs2[locationIdx][b][0] = new float[to_allocate]; // x
                        delin_refined_vecs2[locationIdx][b][1] = new float[to_allocate]; // y

                        // calculate by averaging neighbouring directions

						int windowSize = 5;

                        for (int l=0; l<to_allocate; l++) {

                            float avg_v_x=0, avg_v_y=0;

                            for (int l_nbr=l-windowSize/2; l_nbr<=l+windowSize/2; l_nbr++) {

                            	if (l_nbr>=0 && l_nbr+1<to_allocate) {

                                    float x1 = delin_refined_locs2[locationIdx][b][0][l_nbr];
                                    float y1 = delin_refined_locs2[locationIdx][b][1][l_nbr];
                                    float x2 = delin_refined_locs2[locationIdx][b][0][l_nbr+1];
                                    float y2 = delin_refined_locs2[locationIdx][b][1][l_nbr+1];

                                    float norm = (float) Math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
                                    float vx = (x2-x1)/norm;
                                    float vy = (y2-y1)/norm;

                                    avg_v_x += vy;
                                	avg_v_y += -vx;

                                }

                            }

                            float norm = (float) Math.sqrt(avg_v_x*avg_v_x + avg_v_y*avg_v_y);
                            if (norm>Float.MIN_VALUE) {
                                delin_refined_vecs2[locationIdx][b][0][l] = avg_v_x/norm;
                                delin_refined_vecs2[locationIdx][b][1][l] = avg_v_y/norm;
                            }
                            else {
                                delin_refined_vecs2[locationIdx][b][0][l] = 0;
                                delin_refined_vecs2[locationIdx][b][1][l] = 0;
                            }

                        }

                    }
                    else {
                        delin_refined_vecs2[locationIdx][b] = null;
                    }
                }
            }
            else {
                delin_refined_vecs2[locationIdx] = null;
            }

			/* fitsco2[locationIdx] stdev2[locationIdx] */
			//
			float[] cross_profile = new float[dim];   // aux variable, can be allocated outside the loop (this way it's allocated for every location)
			float[] fit_idx_score = new float[2];     // aux variable
			//
			if (delin_refined_locs2[locationIdx]!=null) {

				for (int b = 0; b<delin_refined_locs2[locationIdx].length; b++) {

					if (delin_refined_locs2[locationIdx][b]!=null) {

						int to_allocate = delin_refined_locs2[locationIdx][b][0].length;  // length of x locs
						fitsco2[locationIdx][b] = new float[to_allocate];
						stdev2[locationIdx][b] = new float[to_allocate];

						// calculate 1-ncc fit scores at every cross-section defined with the point and the direction
						for (int l=0; l<to_allocate; l++) {  // loop points

							int cnt = 0;
							for (int shift=-dim_half; shift<=dim_half; shift++) {
								// sample the values using refined points and orthogonal direction
								float at_x = delin_refined_locs2[locationIdx][b][0][l] + shift * samplingStep * delin_refined_vecs2[locationIdx][b][0][l];
								float at_y = delin_refined_locs2[locationIdx][b][1][l] + shift * samplingStep * delin_refined_vecs2[locationIdx][b][1][l];
								cross_profile[cnt] = Interpolator.interpolateAt(at_x, at_y, inimg_xy);
								cnt++;
							}
							// profile is formed

							// extract std of the cross-profile before normalizing it
							stdev2[locationIdx][b][l] = Stat.std(cross_profile);

							// calculate normalized cross_profile[]
							Stat.normalize(cross_profile); // stays in cross_profile variable

							// calculate the fit
							fit_idx_score = fitter.fit(cross_profile, mode); // returns [0] - profile index, [1] - fit score
							fitsco2[locationIdx][b][l] = fit_idx_score[1];

						}

					}
                    else {
                        fitsco2[locationIdx][b] = null;
                    }
				}
			}
            else {
                fitsco2[locationIdx] = null;
            }

            /*  desc2[locationIdx] using fitsco2[locationIdx] and stdev2[locationIdx] */
            if (fitsco2[locationIdx]!=null) {

                desc2[locationIdx] = new float[4];

                for (int b = 0; b<fitsco2[locationIdx].length; b++) {

                    if (fitsco2[locationIdx][b]!=null) {
                        // calculate means, medians and other descriptors
                        /*
                        HERE IS HOW CALCULATION
                         */
                        desc2[locationIdx][b] = Stat.average(fitsco2[locationIdx][b]);//Stat.median(fitsco2[locationIdx][b]) Stat.average(stdev2[locationIdx][b])
//                        desc2[locationIdx][b][3] = Stat.median(stdev2[locationIdx][b]);

                    }
                    else {
                        desc2[locationIdx][b] = Float.NaN;
                    }

                }

            }
            else {
                desc2[locationIdx]=null;
            }

			/* peaks_lhood_reg[locationIdx] using desc2 and peaks_lhood (from PeakAnalyzer2D) */

			// if fixed to 4
			if (desc2[locationIdx]!=null) {
                peaks_lhood_reg[locationIdx] = new float[4];
				//float sum = 0;
//                float mn = Float.POSITIVE_INFINITY;
//                float mx = Float.NEGATIVE_INFINITY;
                for (int b=0; b<desc2[locationIdx].length; b++) {
					if (!Float.isNaN(desc2[locationIdx][b])) {
						//branches_found.add(b); // add the index of the branch found
						peaks_lhood_reg[locationIdx][b] = peaks_lhood[locationIdx][b];//* (1 - desc2[locationIdx][b][0]);
						//sum += peaks_lhood_reg[locationIdx][b];
//                        if (peaks_lhood_reg[locationIdx][b]>mx) mx = peaks_lhood_reg[locationIdx][b];
//                        if (peaks_lhood_reg[locationIdx][b]<mn) mn = peaks_lhood_reg[locationIdx][b];
					}
					else {
						peaks_lhood_reg[locationIdx][b] = Float.NaN;
					}
				}

			}
			else {
				peaks_lhood_reg[locationIdx] = null;
			}


//			// if keep total # varying
//            // check how many branches there are in desc2[locationIdx]
//            ArrayList<Integer> branches_found = new ArrayList<Integer>();
//            branches_found.clear();
//            if (desc2[locationIdx]!=null) {
//                for (int b=0; b<desc2[locationIdx].length; b++)
//                    if (desc2[locationIdx][b] != null) branches_found.add(b); // add the index of the branch found
//            }
//
//			if (desc2[locationIdx]!=null){
//				peaks_lhood_reg[locationIdx] = new float[branches_found.size()]; // allocate in the length of branches found
//				float sum = 0;
//				for(int bf=0; bf<branches_found.size(); bf++) {
//					int bch_index = branches_found.get(bf);
//					peaks_lhood_reg[locationIdx][bch_index] = peaks_lhood[locationIdx][bch_index] * (1 - desc2[locationIdx][bch_index][0]);
//					sum += peaks_lhood_reg[locationIdx][bch_index];
//				}
//				for(int bf=0; bf<branches_found.size(); bf++) {
//					int bch_index = branches_found.get(bf);
//					peaks_lhood_reg[locationIdx][bch_index] /= sum;
//				}
//			}
//			else {
//				peaks_lhood_reg[locationIdx] = null;
//			}

			/* lhood2[locationIdx] using peaks_lhood_reg[locationIdx] */
            if (peaks_lhood_reg[locationIdx]!=null) {

                float sum = 0;
                int cnt = 0;

                for (int b=0; b<desc2[locationIdx].length; b++) {
                    if (Float.isNaN(desc2[locationIdx][b])) {
                        //branches_found.add(b); // add the index of the branch found
                        sum += (1 - desc2[locationIdx][b]);//peaks_lhood_reg[locationIdx][b];
                        cnt++;
                    }
                }
                endpoint_score2[locationIdx] = new float[]{0};
                //lhood2[locationIdx] = new float[]{1};

//                lhood2[locationIdx] = new float[4];
//                float f1 = peaks_lhood_reg[locationIdx][0];
//                float f2 = peaks_lhood_reg[locationIdx][1];
//                float f3 = peaks_lhood_reg[locationIdx][2];
//                float f4 = peaks_lhood_reg[locationIdx][3];
//                extract_cemd(f1, f2, f3, f4, lhood2[locationIdx]);  // end, bdy, bif, crs, how close it is to each of these distributions

            }
            else {
                //lhood2[locationIdx] = new float[]{0f};
                endpoint_score2[locationIdx] = new float[]{0};
            }

		}

    }

    private static void extract_cemd(float f1, float f2, float f3, float f4, float[] cemd_end_bdy_bif_crs) {

        // store 4 distances
        cemd_end_bdy_bif_crs[0] = Stat.cemd(f1, f2, f3, f4, 1f, 0, 0, 0);
        cemd_end_bdy_bif_crs[1] = Math.min(Stat.cemd(f1, f2, f3, f4, 0.5f, 0.5f, 0, 0), Stat.cemd(f1, f2, f3, f4, 0.5f, 0, 0.5f, 0));
        cemd_end_bdy_bif_crs[2] = Stat.cemd(f1, f2, f3, f4, 1f/3, 1f/3, 1f/3, 0);
        cemd_end_bdy_bif_crs[3] = Stat.cemd(f1, f2, f3, f4, 1f/4, 1f/4, 1f/4, 1f/4);
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

                        for (int l=0; l<delin_refined_locs2[atLoc][b][0].length; l++) {

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

            printout += "\nREFINED VECS:\n";
            if (delin_refined_vecs2[atLoc]!=null) {

                for (int b=0; b<delin_refined_vecs2[atLoc].length; b++) {

                    printout += b+"\t->\t";

                    if (delin_refined_vecs2[atLoc][b]!=null) {

                        for (int l=0; l<delin_refined_vecs2[atLoc][b][0].length; l++) {

                            printout += "("+IJ.d2s(delin_refined_vecs2[atLoc][b][0][l], 2)+", "+IJ.d2s(delin_refined_vecs2[atLoc][b][1][l], 2)+")";

                            if (l==delin_refined_vecs2[atLoc][b][0].length-1) printout += "\n";
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


            // skip descriptors and likelihoods for now
            printout += "\nDESCRIPTORS:\n";
            if (desc2[atLoc]!=null) {
                for (int b=0; b<desc2[atLoc].length; b++) {
                    printout += (b+1) +"\t->\t"; //+ IJ.d2s(ratio2[atLoc][ii], 2) + "\n"

                    if (Float.isNaN(desc2[atLoc][b])) {

                        //for (int l=0; l<desc2[atLoc][b].length; l++) {
							float aa = desc2[atLoc][b];
							float aaa = 1-aa;//(float) Math.exp(-(aa*aa)/(2*0.1*0.1));
                            printout += IJ.d2s(aa, 2)+"(ncc="+IJ.d2s(aaa,2)+")";
                            //if (l==desc2[atLoc][b].length-1) printout += "\t(avg(fit), median(fit), avg(std), med(std))\n";
                            //else printout += ",\t\t";
                        //}

                    }
                    else {
                        printout += "NONE\n";
                    }

                }
            }
            else {
                printout += "NONE\n";
            }

//			printout += "\nLIKELIHOODS:\n";
//            printout += "END->" + IJ.d2s(lhood2[atLoc][0], 2) + "\t\t";
//            printout += "BDY->" + IJ.d2s(lhood2[atLoc][1], 2) + "\t\t";
//            printout += "BIF->" + IJ.d2s(lhood2[atLoc][2], 2) + "\t\t";
//            printout += "CRS->" + IJ.d2s(lhood2[atLoc][3], 2) + "\t\t";
//            printout += "CRS->" + IJ.d2s(lhood2[atLoc][4], 2) + "\t\t";

            IJ.log(printout);

			IJ.log("\npeaks_lhood_reg !!!\n");
			IJ.log(Arrays.toString(peaks_lhood_reg[atLoc]));
			IJ.log(Arrays.toString(peaks_lhood[atLoc]));
            IJ.log("\ndesc2 !!!\n");
            IJ.log(Arrays.toString(desc2[atLoc]));
//			float[] get_pties = peaks_lhood_reg[atLoc];
//			if (get_pties.length>0) {
//				for (int ii=0; ii<get_pties.length; ii++) {
//					IJ.log(IJ.d2s(get_pties[ii], 3));
//					//IJ.log(IJ.d2s(get_pties[ii]*(1-desc2[atLoc][ii][0]),4));	// 0 corrs. to mean
//				}
//			}
//			else {
//				IJ.log("there was no peaks");
//			}



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
        int locationIdx = xy2i[atX][atY];

        if (locationIdx!=-1) {

			if (delin_refined_locs2[locationIdx]!=null) {
				for (int b = 0; b<delin_refined_locs2[locationIdx].length; b++) {

					// check whether there is a refinement at all here, add in case there is
					if (delin_refined_locs2[locationIdx][b]!=null) {

						// loop points to add them
						for (int l=0; l<delin_refined_locs2[locationIdx][b][0].length; l++) {

							float refined_x = delin_refined_locs2[locationIdx][b][0][l];
							float refined_y = delin_refined_locs2[locationIdx][b][1][l];

							ovalroi = new OvalRoi(refined_x-(samplingStep/2)+.5f, refined_y-(samplingStep/2)+.5f, samplingStep, samplingStep);
							ovalroi.setFillColor(Color.YELLOW);
							ovalroi.setStrokeWidth(samplingStep/4);
							ov.add(ovalroi);

						}

					}

				}
			}

        }

        /*
        add lines marking the refined cross sections
         */
        if (locationIdx!=-1) {

			if (delin_refined_vecs2[locationIdx]!=null) {
				for (int b=0; b<delin_refined_vecs2[locationIdx].length; b++) {

					if (delin_refined_vecs2[locationIdx][b]!=null) {

						for (int l=0; l<delin_refined_vecs2[locationIdx][b][0].length; l++) {

							float end_1_x = delin_refined_locs2[locationIdx][b][0][l] - dim_half * samplingStep * delin_refined_vecs2[locationIdx][b][0][l];
							float end_1_y = delin_refined_locs2[locationIdx][b][1][l] - dim_half * samplingStep * delin_refined_vecs2[locationIdx][b][1][l];

							float end_2_x = delin_refined_locs2[locationIdx][b][0][l] + dim_half * samplingStep * delin_refined_vecs2[locationIdx][b][0][l];
							float end_2_y = delin_refined_locs2[locationIdx][b][1][l] + dim_half * samplingStep * delin_refined_vecs2[locationIdx][b][1][l];

							Line lne = new Line(end_1_x+.5f, end_1_y+.5f, end_2_x+.5f, end_2_y+.5f);
							lne.setStrokeColor(Color.CYAN);
							ov.add(lne);

						}

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

    public static ImageStack plotDelineationFeatures(int atX, int atY)
    {

		// just read from the fitsco2 and stdev2, and desc2
        ImageStack is_out = new ImageStack(528, 255);
        int loc_idx = xy2i[atX][atY];

        // calculate total nr. feats
        if (loc_idx!=-1 && fitsco2[loc_idx] != null) {

            int totalFeats = 0;
            for (int b=0; b<fitsco2[loc_idx].length; b++)
				if (fitsco2[loc_idx][b] != null) totalFeats += fitsco2[loc_idx][b].length;

            float[] xaxis       = new float[totalFeats];
            float[] yaxis1      = new float[totalFeats];
            float[] yaxis2      = new float[totalFeats];
            float[] yaxisON     = new float[totalFeats];
            float[] yaxisOFF    = new float[totalFeats];

            int cnt = 0;
            for (int b=0; b<fitsco2[loc_idx].length; b++) {
                if (fitsco2[loc_idx][b] != null) {
                    for (int l=0; l<fitsco2[loc_idx][b].length; l++) {

                        if (l==0) 	if (cnt == 0) xaxis[cnt] = 0;
									else xaxis[cnt] = xaxis[cnt - 1] + 20;
                        else xaxis[cnt] = xaxis[cnt - 1] + 1;

                        yaxis1[cnt]     = fitsco2[loc_idx][b][l];   // features
                        yaxis2[cnt]     = desc2[loc_idx][b];     // mean fitsco2
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
                if (yaxis1[ii]<min_feat) min_feat = yaxis1[ii];
                if (yaxis1[ii]>max_feat) max_feat = yaxis1[ii];
                if (xaxis[ii]>max_axis) max_axis = xaxis[ii];
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

	public static void exportLikelihoods(String file_path)
	{

		PrintWriter logWriter = null; //initialize writer

		try {
			logWriter = new PrintWriter(file_path);	logWriter.print(""); logWriter.close();
		} catch (FileNotFoundException ex) {}   // empty the file before logging...

		try {                                   // initialize
			logWriter = new PrintWriter(new BufferedWriter(new FileWriter(file_path, true)));
		} catch (IOException e) {}

//		//main loop
//		for (int ii=0; ii<lhood2.length; ii++) {
//			for (int jj=0; jj<lhood2[ii].length; jj++) { // loop streamlines
//					logWriter.print(String.format("%1.8f", lhood2[ii][jj]));
//					if (jj<lhood2[ii].length-1) logWriter.print(", ");
//					else logWriter.print("\n");
//			}
//		}
		logWriter.close();

	}

    public static ImageStack exportLikelihoods() // {0, 1, 2, 3, 4}
    {

        int w = inimg_xy.length;
        int h = inimg_xy[0].length;

        ImageStack is_out = new ImageStack(w, h);

        //for (int i=0; i<endpoint_score2[0].length; i++) {
            //if (choice[i] >=0 && choice[i]<4) {

                float[][] t = new float[w][h];

                for (int ii = 0; ii<endpoint_score2.length; ii++) {

                    int x = i2xy[ii][0];
                    int y = i2xy[ii][1];

                    boolean isMax = true;

//					for (int k=0; k<lhood2[0].length; k++) {
//                        if (k!=choice[i] && lhood2[ii][k]>lhood2[ii][choice[i]]) {
//                            isMax = false;
//                        }
//                    }

                    if (endpoint_score2[ii]!=null) {
                        t[x][y] = (isMax)? endpoint_score2[ii][0] : 0;
                    }


                }

                is_out.addSlice("ENDPT", new FloatProcessor(t));
            //}
        //}

        return is_out;
    }

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
        refined locations per cross profile with outward direction estimation, based on the first patch refined points
     */
    private static void localPatchRefinedLocs(float x1, float y1, float x2, float y2,
                                          int init_index,
                                          float[][] refined_centerline_locs_xy,
                                          float[] refined_offsets,               // optional
                                          float[] refined_centerline_dir)        // optional
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

            if (refined_offsets!=null) refined_offsets[0 + ii] = offset_min;

        }

        /*
        calculate the direction based on fitting a line through calculated refined_centerline_locs_xy
         */
		if(refined_centerline_dir!=null) {

			float l1x, l1y, l2x, l2y, px, py;
			float min_ssd = Float.POSITIVE_INFINITY;
			refined_centerline_dir[0] = Float.NaN;
			refined_centerline_dir[1] = Float.NaN;
			refined_centerline_dir[2] = Float.NaN;
			refined_centerline_dir[3] = Float.NaN;

			for (int jj=-dim_half; jj<=dim_half; jj++) { // loop first cross-section points
				for (int kk=-dim_half; kk<=dim_half; kk++) { // loop last cross-section points

					// ii=0
					l1x = x_root + 0 * samplingStepLongitudinal * vx + jj * samplingStep * vy;
					l1y = y_root + 0 * samplingStepLongitudinal * wx + jj * samplingStep * wy;
					// ii=L-1
					l2x = x_root + (L-1) * samplingStepLongitudinal * vx + kk * samplingStep * vy;
					l2y = y_root + (L-1) * samplingStepLongitudinal * wx + kk * samplingStep * wy;

					// calculate m.sq.distance wrt to the rest
					float curr_ssd = 0;

					for (int ii=0; ii<L; ii++) {

						px = refined_centerline_locs_xy[0][init_index + ii]; // use the ones calculated in previous loop
						py = refined_centerline_locs_xy[1][init_index + ii];

						float dist = Geometry.point_to_line(px, py, l1x, l1y, l2x, l2y);
						curr_ssd += dist * dist;

					}

					if (curr_ssd<min_ssd) {
						min_ssd = curr_ssd;
						refined_centerline_dir[0] = l1x;
						refined_centerline_dir[1] = l1y;
						refined_centerline_dir[2] = l2x;
						refined_centerline_dir[3] = l2y;
					}

				}
			}

		}



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

}