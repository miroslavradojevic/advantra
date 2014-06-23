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
 * refines the locations for cross profile extraction
 * TASK 2:
 * Calculates and exports the fit scores (NCC) features at every location after frame delineation (feat2)
 * feat2:
 *      there are L such NCC values per patch, and M patches per branch (up to 4 branches): every location 4*M*L elements
 * ratio2:
 *      scores for each branch 4*1 (M*L elements sum up to one number telling the ratio of those above the threshold)
 */
public class FeatureExtractor2D extends Thread {

    private int begN, endN;

    // INPUT:
	public static int[][] 	    i2xy;                       // selected locations
	public static int[][]     	xy2i;                       // need for recursion
	public static int[][]       peaks_i;             	    // list of extracted peaks: N x 4 every FG location with 4 extracted peaks in indexed format
	public static int[][]		peaks_w;                    // weight assigned to each peak (for expansion)
	public static float[][]		inimg_xy;				    // input image (necessary for feature extraction)

    // PARAMETERS
    public static float     D;
	public static float		cross_sigma_ratio;
    public static int       M;                              // how much it expands recursively from the center
    public static float     minCos;                         // allowed derail
    private static float    samplingStep;                   // when sampling image values to extract features
    private static int		L;                              // will define how many are taken along the diameter, in radial direction
    private static float    samplingStepLongitudinal;       // sampling along the streamline of patches
    private static int      dim;                            // cross profile length
	private static int      dim_half;                       // cross profile half-length
    private static String   ncc_estimator;                  // select ncc values estimator

	// detection params...
    private static float    ncc_high_mean;
    private static float    ncc_low_mean;
	private static float    likelihood_high_mean;
	private static float    likelihood_low_mean;

	private static float    output_sigma;

    // OUTPUT
    // associate the peaks and link follow-up points
    public static int[][][]     delin2;                     // N(foreground locs.) x 4(max. threads) x (1..M) (follow-up locs.) contains index for each location
    public static float[][][][] delin_refined_locs2;        // N(foreground locs.) x 4(max. threads) x 2 x ((1..M) x L) (2 dimensional)
	public static float[][][][] delin_refined_vecs2;        // N(foreground locs.) x 4(max. threads) x 2 x ((1..M) x L) (2 dimensional)

    // features
    public static float[][][]   fitsco2;					// N(foreground locs.) x 4(max. threads) x ((1..M) x L) (fit scores ncc for refined locs)
    public static float[][][]   stdev2;						// N(foreground locs.) x 4(max. threads) x ((1..M) x L) (variances of real values)

    // descriptions (calculations on features)
    public static float[][]     ncc2;                       // N(foreground locs.) x 4(max. threads) (ncc moment along the stream)
	public static float[][]		lhoods2; 					// N(foreground locs.) x 4(max. threads) (normalized 0-1) likelihoods, read from Peak

	// fuzzy system outputs
    public static float[][][]   streamline_score2;          // N(foreground locs.) x 3 (off,none,on)x 4(max. threads)
    public static float[][]   	critpoint_score2;           // N(foreground locs.) x 3(endpoint,nonepoint,bifpoint) given by fuzzy logic

    public FeatureExtractor2D(int n0, int n1)
    {
        this.begN = n0;
        this.endN = n1;
    }

    public static void loadTemplate(
										   	int[][] _i2xy,
											int[][] _xy2i,
											int[][] _peaks_i,
											int[][] _peaks_w,
											float[][] _peaks_lhood,  // it was necessary for the detection

											float[][] _inimg_xy,
											int     _M,
                                    		float   _minCos,

											float   _D,
											float  	_cross_sigma_ratio, // the same one used for profiler, now necessary for the fitter profiles

                                            String  _ncc_estimator,

                                    		float   _ncc_high_mean,
                                    		float   _ncc_low_mean,
											float 	_likelihood_high_mean,
											float 	_likelihood_low_mean,
											float 	_output_sigma
//											int     _L,                         // will define sampling longitudinal
//                                    		float   _sampling_crosswise
    )
    {
        M               = _M;
        minCos          = _minCos;
        D               = _D;
        L               = Math.round(_D) + 1;
		samplingStep    = .5f;
		cross_sigma_ratio = _cross_sigma_ratio;

        ncc_estimator = _ncc_estimator;

        ncc_high_mean   = _ncc_high_mean;
        ncc_low_mean    = _ncc_low_mean;
		likelihood_high_mean   = _likelihood_high_mean;
        likelihood_low_mean    = _likelihood_low_mean;
		output_sigma = _output_sigma;

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

        lhoods2 			= _peaks_lhood; // these features have been extracted before

		delin_refined_locs2 = new float[i2xy.length][4][2][];       // (x,y)
		delin_refined_vecs2 = new float[i2xy.length][4][2][];       // (x,y)
        fitsco2   			= new float[i2xy.length][4][];          // fit score
        stdev2   			= new float[i2xy.length][4][];          // variance
        ncc2  			    = new float[i2xy.length][];          	// average of the fit scores along streamline
		streamline_score2  	= new float[i2xy.length][][];     		// branch score 	fg. location x 3 (off, none, on)x 4 (direcitons)
        critpoint_score2    = new float[i2xy.length][];             // critpoint score  fg. location x 3 (end, non, jun)

	}

	public static ImageStack getFuzzyAgg(int atX, int atY){

        ImageStack is_out = new ImageStack(528, 255);

		int locationIdx = xy2i[atX][atY];
		if (locationIdx!=-1 && critpoint_score2[locationIdx]!=null) {

			Fuzzy2D fls = new Fuzzy2D(
											    // ncc
											    ncc_high_mean,// ncc_high_sigma,      // high
											    ncc_low_mean, //ncc_low_sigma,   // low
											    // lhood
											    likelihood_high_mean, //likelihood_high_sigma, // high
											    likelihood_low_mean, //likelihood_low_sigma,  // low
                                                0,
                                                0,
											    output_sigma  // std output membership functions - defines separation margin
			);

			if (ncc2[locationIdx]!=null) {

				// check how many inputs there are
				ArrayList<Integer> b_sel = new ArrayList<Integer>();
				b_sel.clear();
				for (int b=0; b<ncc2[locationIdx].length; b++)
					if (!Float.isNaN(ncc2[locationIdx][b])) b_sel.add(b); // add the index of the branch found

				// switch according to the length of valid branches found
				float ncc_1, likelihood_1, ncc_2, likelihood_2, ncc_3, likelihood_3, ncc_4, likelihood_4; // aux variables
				float[] tmp = new float[3];
				boolean show_agg = true;


				switch (b_sel.size()) {
					case 1:
						ncc_1 = ncc2[locationIdx][b_sel.get(0)];
						likelihood_1 = lhoods2[locationIdx][b_sel.get(0)];
						is_out = fls.critpointScore(ncc_1, likelihood_1, tmp, show_agg);
						break;
					case 2:
						ncc_1 = ncc2[locationIdx][b_sel.get(0)];
						likelihood_1 = lhoods2[locationIdx][b_sel.get(0)];
						ncc_2 = ncc2[locationIdx][b_sel.get(1)];
						likelihood_2 = lhoods2[locationIdx][b_sel.get(1)];
						is_out = fls.critpointScore(ncc_1, likelihood_1, ncc_2, likelihood_2, tmp, show_agg);
						break;
					case 3:
						ncc_1 = ncc2[locationIdx][b_sel.get(0)];
						likelihood_1 = lhoods2[locationIdx][b_sel.get(0)];
						ncc_2 = ncc2[locationIdx][b_sel.get(1)];
						likelihood_2 = lhoods2[locationIdx][b_sel.get(1)];
						ncc_3 = ncc2[locationIdx][b_sel.get(2)];
						likelihood_3 = lhoods2[locationIdx][b_sel.get(2)];
						is_out = fls.critpointScore(ncc_1, likelihood_1, ncc_2, likelihood_2, ncc_3, likelihood_3, tmp, show_agg);
						break;
					case 4:
						ncc_1 = ncc2[locationIdx][b_sel.get(0)];
						likelihood_1 = lhoods2[locationIdx][b_sel.get(0)];
						ncc_2 = ncc2[locationIdx][b_sel.get(1)];
						likelihood_2 = lhoods2[locationIdx][b_sel.get(1)];
						ncc_3 = ncc2[locationIdx][b_sel.get(2)];
						likelihood_3 = lhoods2[locationIdx][b_sel.get(2)];
						ncc_4 = ncc2[locationIdx][b_sel.get(3)];
						likelihood_4 = lhoods2[locationIdx][b_sel.get(3)];
						is_out = fls.critpointScore(ncc_1, likelihood_1, ncc_2, likelihood_2, ncc_3, likelihood_3, ncc_4, likelihood_4, tmp, show_agg);
						break;
					default:
						//imp_out = null;
						break;
				}

			}
			else {
                Plot p = new Plot("", "", "", new float[1], new float[1]);
				is_out.addSlice(p.getProcessor());
			}



		}
		else {
            Plot p = new Plot("", "", "", new float[1], new float[1]);
            is_out.addSlice(p.getProcessor());
		}

		return is_out;

	}

    public void run()
    {

		// processor modules for this run - each run should have it's own initialized with the same parameteres
        Fitter1D fitter  			= new Fitter1D(dim, cross_sigma_ratio, false); // dim = profile width with current samplingStep, verbose = false
		//fitter.showTemplates();
        Fuzzy2D fls = new Fuzzy2D(
                // ncc
                ncc_high_mean,// ncc_high_sigma,      // high
                ncc_low_mean, //ncc_low_sigma,   // low
                // lhood
                likelihood_high_mean, //likelihood_high_sigma, // high
                likelihood_low_mean, //likelihood_low_sigma,  // low
                0,0,
                output_sigma  // std output membership functions - defines separation margin
        );
        // aux variables
		float[] cross_profile = new float[dim];
		float[] fit_idx_score;
        float[] streamline_off_none_on = new float[3];
        float[] critpoint_endpoint_nonepoint_bifpoint = new float[3];

        for (int locationIdx = begN; locationIdx < endN; locationIdx++) {

            /* delin2[locationIdx] calculation */
            for (int pp = 0; pp<peaks_i[locationIdx].length; pp++) {  // loop 4 allocated branches (1st generation)
                // access individual peaks at this point (-2:not exist, -1:background, >=0:foreground)
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
                        else if (next_index==-1){ // not found but in the background
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



            /* delin_refined_locs2[locationIdx] */
            for (int b = 0; b<delin2[locationIdx].length; b++) { // loop 4 branches

                if (delin2[locationIdx][b][0]==-1) {
                    // whole streamline is missing: the first one was -1, recursion was stopped
                    delin_refined_locs2[locationIdx] = null;
                    break; // stop looping branches further
                }
                else if (delin2[locationIdx][b][0]==-2) {
                    // no streamline here
                    delin_refined_locs2[locationIdx][b] = null;
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
                            localPatchRefinedLocs(
                                    prev_x, prev_y, curr_x, curr_y,
                                    m * L,        // start from
                                    delin_refined_locs2[locationIdx][b]);

                        }
                        else break; // because the rest are filled up with -1 or -2 anyway

                    }

                }

            }





            /*
            *	calculate delin_refined_vecs2[locationIdx]  using delin_refined_locs2[locationIdx]
            */
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

























            /*
            *  calculate fitsco2[locationIdx]  and stdev2[locationIdx] using cross profiles at refined locations
            */
			if (delin_refined_locs2[locationIdx]!=null) {

				for (int b = 0; b<delin_refined_locs2[locationIdx].length; b++) {

					if (delin_refined_locs2[locationIdx][b]!=null) {

						int to_allocate = delin_refined_locs2[locationIdx][b][0].length;  // length of x locs
						fitsco2[locationIdx][b] = new float[to_allocate];
						stdev2[locationIdx][b] = new float[to_allocate];

						// calculate ncc fit scores at every cross-section defined with the point and the direction
						for (int l=0; l<to_allocate; l++) {  // loop points

							int cnt = 0;
							for (int shift=-dim_half; shift<=dim_half; shift++) {
								// sample the values using refined locs and refined vecs
								float at_x = delin_refined_locs2[locationIdx][b][0][l] + shift * samplingStep * delin_refined_vecs2[locationIdx][b][0][l];
								float at_y = delin_refined_locs2[locationIdx][b][1][l] + shift * samplingStep * delin_refined_vecs2[locationIdx][b][1][l];
								cross_profile[cnt] = Interpolator.interpolateAt(at_x, at_y, inimg_xy);
								cnt++;
							}
							// cross-profile is formed

							// extract std of the cross-profile before normalizing it
							stdev2[locationIdx][b][l] = Stat.std(cross_profile);

							// calculate normalized cross_profile[]
							Stat.normalize(cross_profile); // stays in cross_profile variable

							// calculate the fit
							fit_idx_score = fitter.fit(cross_profile, "NCC"); // returns [0] - profile index, [1] - fit score
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


            /*
            *  calculate ncc2[locationIdx] using fitsco2[locationIdx] as an estimate of all the ncc fit scores extracted
            */
            if (fitsco2[locationIdx]!=null) {
                ncc2[locationIdx] = new float[4];
                for (int b = 0; b<fitsco2[locationIdx].length; b++) {
                    if (fitsco2[locationIdx][b]!=null) {
                        // choose estimator of ncc scores along the stream (mean, median, percentile(quantile))
                        if (ncc_estimator.equals("AVERAGE")) {
                            ncc2[locationIdx][b] = Stat.average(fitsco2[locationIdx][b]);
                        }
                        else if (ncc_estimator.equals("MEDIAN")) {
                            ncc2[locationIdx][b] = Stat.median(fitsco2[locationIdx][b]);
                        }
                        else if (ncc_estimator.equals("QUANTILE")) {
                            ncc2[locationIdx][b] = Stat.quantile(fitsco2[locationIdx][b], 6, 20); // hardcoded 6/20
                        }
                        else {
							System.out.println("wrong estimator argument...");
							ncc2[locationIdx][b] = Float.NaN;
                        }

                    }
                    else {
                        ncc2[locationIdx][b] = Float.NaN;
                    }
                }
            }
            else {
                ncc2[locationIdx]=null;
            }




            /*
            *  calculate streamline_score[locationIdx] -> {OFF, NONE, ON} response of each individual branch
            */
            if (ncc2[locationIdx]!=null) {
                streamline_score2[locationIdx] = new float[3][4]; // {OFF, NONE, ON} X NR. BRANCHES
                for (int b=0; b<ncc2[locationIdx].length; b++) {
                    if (!Float.isNaN(ncc2[locationIdx][b])) {
                        fls.branchStrength(ncc2[locationIdx][b], lhoods2[locationIdx][b], streamline_off_none_on);
                        streamline_score2[locationIdx][0][b] = streamline_off_none_on[0];
                        streamline_score2[locationIdx][1][b] = streamline_off_none_on[1];
                        streamline_score2[locationIdx][2][b] = streamline_off_none_on[2];
                    }
                    else {
                        streamline_score2[locationIdx][0][b] = Float.NaN;
                        streamline_score2[locationIdx][1][b] = Float.NaN;
                        streamline_score2[locationIdx][2][b] = Float.NaN;
                    }
                }
            }
            else {
                streamline_score2[locationIdx] = null;
            }

            /*
            *  calculate critpoints_score[locationIdx] -> {ENDPOINT, NONEPOINT, BIFPOINT}
            */

            if (ncc2[locationIdx]!=null) {

				critpoint_score2[locationIdx] =new float[3];

                // check how many inputs there are
                ArrayList<Integer> b_sel = new ArrayList<Integer>();
                b_sel.clear();
                for (int b=0; b<ncc2[locationIdx].length; b++)
                    if (!Float.isNaN(ncc2[locationIdx][b])) b_sel.add(b); // add the index of the branch found

				// switch according to the length of valid branches found
                float ncc_1, likelihood_1, ncc_2, likelihood_2, ncc_3, likelihood_3, ncc_4, likelihood_4;
                switch (b_sel.size()) {
                    case 1:
                        ncc_1 = ncc2[locationIdx][b_sel.get(0)];
                        likelihood_1 = lhoods2[locationIdx][b_sel.get(0)];
						fls.critpointScore(ncc_1, likelihood_1, critpoint_endpoint_nonepoint_bifpoint, false);
						critpoint_score2[locationIdx][0] = critpoint_endpoint_nonepoint_bifpoint[0];
						critpoint_score2[locationIdx][1] = critpoint_endpoint_nonepoint_bifpoint[1];
						critpoint_score2[locationIdx][2] = critpoint_endpoint_nonepoint_bifpoint[2];
                        break;
                    case 2:
                        ncc_1 = ncc2[locationIdx][b_sel.get(0)];
                        likelihood_1 = lhoods2[locationIdx][b_sel.get(0)];
                        ncc_2 = ncc2[locationIdx][b_sel.get(1)];
                        likelihood_2 = lhoods2[locationIdx][b_sel.get(1)];

						fls.critpointScore(ncc_1, likelihood_1, ncc_2, likelihood_2, critpoint_endpoint_nonepoint_bifpoint, false);
						critpoint_score2[locationIdx][0] = critpoint_endpoint_nonepoint_bifpoint[0];
						critpoint_score2[locationIdx][1] = critpoint_endpoint_nonepoint_bifpoint[1];
						critpoint_score2[locationIdx][2] = critpoint_endpoint_nonepoint_bifpoint[2];
                        break;
                    case 3:
                        ncc_1 = ncc2[locationIdx][b_sel.get(0)];
                        likelihood_1 = lhoods2[locationIdx][b_sel.get(0)];
                        ncc_2 = ncc2[locationIdx][b_sel.get(1)];
                        likelihood_2 = lhoods2[locationIdx][b_sel.get(1)];
                        ncc_3 = ncc2[locationIdx][b_sel.get(2)];
                        likelihood_3 = lhoods2[locationIdx][b_sel.get(2)];

						fls.critpointScore(ncc_1, likelihood_1, ncc_2, likelihood_2, ncc_3, likelihood_3, critpoint_endpoint_nonepoint_bifpoint, false);
						critpoint_score2[locationIdx][0] = critpoint_endpoint_nonepoint_bifpoint[0];
						critpoint_score2[locationIdx][1] = critpoint_endpoint_nonepoint_bifpoint[1];
						critpoint_score2[locationIdx][2] = critpoint_endpoint_nonepoint_bifpoint[2];
                        break;
                    case 4:
                        ncc_1 = ncc2[locationIdx][b_sel.get(0)];
                        likelihood_1 = lhoods2[locationIdx][b_sel.get(0)];
                        ncc_2 = ncc2[locationIdx][b_sel.get(1)];
                        likelihood_2 = lhoods2[locationIdx][b_sel.get(1)];
                        ncc_3 = ncc2[locationIdx][b_sel.get(2)];
                        likelihood_3 = lhoods2[locationIdx][b_sel.get(2)];
                        ncc_4 = ncc2[locationIdx][b_sel.get(3)];
                        likelihood_4 = lhoods2[locationIdx][b_sel.get(3)];

						fls.critpointScore(ncc_1, likelihood_1, ncc_2, likelihood_2, ncc_3, likelihood_3, ncc_4, likelihood_4, critpoint_endpoint_nonepoint_bifpoint, false);
						critpoint_score2[locationIdx][0] = critpoint_endpoint_nonepoint_bifpoint[0];
						critpoint_score2[locationIdx][1] = critpoint_endpoint_nonepoint_bifpoint[1];
						critpoint_score2[locationIdx][2] = critpoint_endpoint_nonepoint_bifpoint[2];
                        break;
                    default:
                        critpoint_score2[locationIdx][0] = 0;
                        critpoint_score2[locationIdx][1] = 0;
                        critpoint_score2[locationIdx][2] = 0;
                        break;
                }

            }
            else {
                critpoint_score2[locationIdx] = null;
            }

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

            printout += "\nFEATURES \n";
            if (ncc2[atLoc]!=null) {
                for (int b=0; b<ncc2[atLoc].length; b++) {
                    printout += (b+1) +"\t->\t"; //+ IJ.d2s(ratio2[atLoc][ii], 2) + "\n"
                    if (!Float.isNaN(ncc2[atLoc][b])) {
						printout += "NCC " + IJ.d2s(ncc2[atLoc][b], 2)+"  \t LHOOD "+IJ.d2s(lhoods2[atLoc][b], 2)+"  \t \t  OFF "
											+IJ.d2s(streamline_score2[atLoc][0][b], 2)+"  \t  ON "+IJ.d2s(streamline_score2[atLoc][2][b], 2)+"\n";
                    }
                    else {
                        printout += "NONE\n";
                    }
                }
            }
            else {
                printout += "NONE\n";
            }

//			printout += "\nBRANCH STRENGTH \n...todo...\n";

			printout += "\nCRITPOINT \n";
			if (critpoint_score2[atLoc]!=null) {
				printout += "END "+critpoint_score2[atLoc][0]+", NONE "+critpoint_score2[atLoc][1]+", BIF "+critpoint_score2[atLoc][2];
			}
			else {
				printout += "NONE\n";
			}

            IJ.log(printout);

//			Fuzzy2D fls = new Fuzzy2D(
//											 // ncc
//											 ncc_high_mean,// ncc_high_sigma,      // high
//											 ncc_low_mean, //ncc_low_sigma,   // low
//											 // lhood
//											 likelihood_high_mean, //likelihood_high_sigma, // high
//											 likelihood_low_mean, //likelihood_low_sigma,  // low
//											 output_sigma  // std output membership functions - defines separation margin
//			);


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
		Fitter1D fitter = new Fitter1D(dim, cross_sigma_ratio, false);
		ImageStack isOut = new ImageStack(528, 255);

		int idx = xy2i[atX][atY]; // read extracted peaks at this location

		if (idx!=-1) {

			int[][] delin_at_loc = delin2[idx];
			ArrayList<float[]> profiles_along = new ArrayList<float[]>(); // list of profiles for all branches

			for (int b = 0; b<delin_at_loc.length; b++) {           // loop 4 branches, b index defines the strength

				for (int m = 0; m<M; m++) {

					if (delin_at_loc[b][m] != -1 && delin_at_loc[b][m] != -2) {

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
				float[] out_idx_scr = fitter.fit(profiles_along.get(aaa), "NCC");
				float[] curr_fit = fitter.getTemplate((int)out_idx_scr[0]);
				plt.addPoints(xx, curr_fit, Plot.LINE);
				plt.draw();
				isOut.addSlice("feat="+ IJ.d2s(out_idx_scr[1], 2), plt.getProcessor());
			}

		}

		if (isOut.getSize()==0) {
			isOut.addSlice(new ByteProcessor(528,255));
		}

		return isOut;

	}

    public static ImageStack plotDelineationFeatures(int atX, int atY)
    {

		// just read from the fitsco2 and stdev2, and ncc_avg2
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
//            float[] yaxisON     = new float[totalFeats];
//            float[] yaxisOFF    = new float[totalFeats];

            int cnt = 0;
            for (int b=0; b<fitsco2[loc_idx].length; b++) {
                if (fitsco2[loc_idx][b] != null) {
                    for (int l=0; l<fitsco2[loc_idx][b].length; l++) {

                        if (l==0) 	if (cnt == 0) xaxis[cnt] = 0;
									else xaxis[cnt] = xaxis[cnt - 1] + 50;
                        else xaxis[cnt] = xaxis[cnt - 1] + 1;

                        yaxis1[cnt]     = fitsco2[loc_idx][b][l];   // features
                        yaxis2[cnt]     = ncc2[loc_idx][b];     // mean fitsco2
//                        yaxisON[cnt]    = mu_ON;
//                        yaxisOFF[cnt]   = mu_OFF;
                        cnt++;

                    }
                }
            }

            // find limits for the features
            float min_feat = Float.MAX_VALUE;
            float max_feat = Float.NEGATIVE_INFINITY;
            float max_axis = Float.NEGATIVE_INFINITY;

            for (int ii=0; ii<yaxis1.length; ii++) {
//                if (yaxis1[ii]<min_feat) min_feat = yaxis1[ii];
//                if (yaxis1[ii]>max_feat) max_feat = yaxis1[ii];
                if (xaxis[ii]>max_axis) max_axis = xaxis[ii];
            }

            Plot p = new Plot("", "", "");
            p.setLimits(0, max_axis, 0, 1.01);//Math.min(min_feat, mu_ON-.1f), Math.max(max_feat, mu_OFF+.1f));
            p.addPoints(xaxis, yaxis1, Plot.X);
            p.draw();

//            p.setLineWidth(4);
//            p.setColor(Color.GREEN);
//            p.addPoints(xaxis, yaxisON, Plot.DOT);
//            p.draw();

//            p.setLineWidth(4);
//            p.setColor(Color.BLUE);
//            p.addPoints(xaxis, yaxisOFF, Plot.DOT);
//            p.draw();

            p.setLineWidth(2);
            p.setColor(Color.RED);
            p.addPoints(xaxis, yaxis2, Plot.LINE);
            p.draw();

            is_out.addSlice("features,descriptors", p.getProcessor());

        }
        else {
            Plot p = new Plot("", "", "", new float[1], new float[1]);
            is_out.addSlice(p.getProcessor());
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

        for (int ii=0; ii<ncc_avg2.length; ii++) {
            for (int jj=0; jj<ncc_avg2[ii].length; jj++) {

                if (ncc_avg2[ii][jj] != null) {

                    for (int kk=0; kk<ncc_avg2[ii][jj].length; kk++) {
                        logWriter.print(String.format("%1.6f, ", ncc_avg2[ii][jj][kk]));
                        if (kk<ncc_avg2[ii][jj].length-1) logWriter.print(", ");
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
//        logWriter.println("feature 0: \tbranch 0 min score");
//        logWriter.println("feature 1: \tbranch 0 max score");
//        logWriter.println("feature 2: \tbranch 1 min score");
//        logWriter.println("feature 3: \tbranch 1 max score");
//        logWriter.println("feature 4: \tbranch 2 min score");
//        logWriter.println("feature 5: \tbranch 2 max score");
//        logWriter.println("feature 6: \tcenter       score");

        logWriter.close(); // close log

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


    /*
        refined locations per cross profile with outward direction estimation, based on the first patch refined points
     */
    private static void localPatchRefinedLocs(float x1, float y1, float x2, float y2,
                                          int init_index,
                                          float[][] refined_centerline_locs_xy)
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

	public static void exportCritpointScores(float[][] out_scores, int idx_scores) { // 0,1,2 idx_scores -> END,NON,BIF

		for (int ll=0; ll<critpoint_score2.length; ll++) {
			if (critpoint_score2[ll]!=null) {

				int x = i2xy[ll][0];
				int y = i2xy[ll][1];

				out_scores[x][y] = critpoint_score2[ll][idx_scores]; // 0,1,2 corresponds to END,NON,BIF

			}

		}
	}

}