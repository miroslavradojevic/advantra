package detection2d;

import aux.Hist;
import aux.Interpolator;
import ij.ImageStack;
import ij.gui.Line;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.gui.Plot;

import java.awt.*;
import java.util.ArrayList;

/**
 * Created by miroslav on 6/16/14.
 */
public class Delineator2D extends Thread {

	// associate the peaks and link follow-up points, delineate local structure

	private int begN, endN;

	// INPUT:
	public static int[][] 	    i2xy;                       // selected locations
	public static int[][]     	xy2i;                       // need for recursion
	public static int[][]       peaks_i;             	    // list of extracted peaks: N x 4 every FG location with 4 extracted peaks in indexed format
	public static int[][]		peaks_w;                    // weight assigned to each peak (for expansion)
	public static float[][]		inimg_xy;				    // input image (necessary for feature extraction)
	public static boolean[] 	critpoint_candidate;        // if profile had enough of range (same length as i2xy)

	public static float     D;								// neuron diameter parameter
	public static int       M;                              // how much it expands recursively from the center
	public static float     minCos;                         // allowed derail
	private static int		L;                              // will define how many are taken along the diameter, in radial direction

	private static float 	dx0 = Float.NaN;				// needs to be initialized
	private static float	dx = Float.NaN;
	private static int 		windowSize = 5;                 // vxy (refined vectors) calculate by averaging neighbouring directions
	private static float 	samplingStep    = .5f;          // cross-side sampling step
	private static float    samplingStepLongitudinal;       // sampling along the streamline of patches
	private static int      dim;                            // cross profile length (number of samples in a cross-profile)
	private static int      dim_half;                       // cross profile half-length

	// OUTPUT
	public static int[][][]     delin2;                     // N(foreground locs.) x 4(max. threads) x (1..M) (follow-up locs.) contains index for each location
	public static float[][][][] xy2;        				// N(foreground locs.) x 4(max. threads) x 2 x ((1..M) x L) (2 dimensional)
	public static float[][][][] vxy2;        				// N(foreground locs.) x 4(max. threads) x 2 x ((1..M) x L) (2 dimensional)
	public static float[][]		smoothness;					// N(foreground locs.) x 4(max. threads) smoothness score

	public Delineator2D(int n0, int n1)
	{
		this.begN = n0;
		this.endN = n1;
	}

	public static void loadTemplate(
										   int[][] 		_i2xy,
										   int[][] 		_xy2i,
										   boolean[] 	_critpoint_candidate,
										   int[][] 		_peaks_i,
										   int[][] 		_peaks_w,
										   float[][] 	_inimg_xy,
										   Sphere2D		_extractionSphere,
										   int     		_M,
										   float   		_minCos
										   )
	{

		i2xy = _i2xy;
		xy2i = _xy2i;
		peaks_i = _peaks_i;
		peaks_w = _peaks_w;
		inimg_xy = _inimg_xy;
		critpoint_candidate = _critpoint_candidate;

		D 				= _extractionSphere.getNeuronDiameter();
		M               = _M;
		minCos          = _minCos;
		L               = (int) (Math.ceil(_extractionSphere.getNeuronDiameter()) + 2);  // added two because of second derivative calcualtion on the refined locations

		samplingStepLongitudinal = D / (float)(L-1);
		dim_half = (int) Math.ceil( D / (samplingStep*2) );
		dim = 2*dim_half + 1;

		dx0 			= _extractionSphere.getInitialRadialSeparation();
		dx 				= samplingStepLongitudinal;

		// allocate output
		delin2 				= new int[i2xy.length][4][M];
		// initialize with -2 (same meaning as with peaks_i)
		for (int i=0; i<delin2.length; i++)
			for (int j=0; j<delin2[i].length; j++)
				for (int k=0; k<delin2[i][j].length; k++)
					delin2[i][j][k] = -2;

		// TODO: enable allocating [4][2] afterwards...
		xy2 		= new float[i2xy.length][4][2][];       // N (foreground points) x 4 (branches) x 2 (x,y) x M*1..L (along branches) (x,y)
		smoothness 	= new float[xy2.length][4];            // N(foreground locs.)   x 4 (branches)
		vxy2 		= new float[i2xy.length][4][2][];       // N (foreground points) x 4 (branches) x 2 (vx,vy) x M*1..L (along branches) (vx,vy)


	}

	public void run()
	{

		for (int locationIdx = begN; locationIdx < endN; locationIdx++) {

			// auxilliary array to store local cross section y values
			// (side output when calculating refined locs, used for smoothness, second gradient calculation)
			float[][] yptch = new float[4][];

			if (!critpoint_candidate[locationIdx]) { // ruled out because of insignificant profile variation

				delin2[locationIdx] = null;
				xy2[locationIdx] = null;
				vxy2[locationIdx] = null;

			}
			else {

				/*
				* 	delin2[locationIdx]
				*/
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
				-2  -2 -2      -> whole streamline is empty nothing was filled up there
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
				* 	xy2[locationIdx]
				* 	y2local[locationIdx]  used temporarily to accumulate variables for the smoothness calculation
				*/
				// will simultaneously fill yptch[4][] up four times or skip using it at all
				for (int b = 0; b<delin2[locationIdx].length; b++) { // loop 4 branches

					if (delin2[locationIdx][b][0]==-1) {
						// whole streamline is missing: the first one was -1, recursion was stopped
						xy2[locationIdx] = null;
						smoothness[locationIdx] = null;
						break; 	// stop looping branches further
					}
					else if (delin2[locationIdx][b][0]==-2) {
						// no streamline here
						xy2[locationIdx][b] = null;
						smoothness[locationIdx][b] = Float.NaN;
						yptch[b] = null;
					}
					else if (delin2[locationIdx][b][0]>=0) {
						// there is at least one patch in the streamline
						// loop the rest to count how many patches there are to allocate the array
						int count_patches = 1;
						for (int m=1; m<M; m++) {
							if (delin2[locationIdx][b][m] >= 0) count_patches++;
							else break; // because the rest are filled up with -1 or -2 anyway
						}

						// allocate (now we know how much to allocate)
						xy2[locationIdx][b][0] 	= new float[count_patches*L]; // x coordinates allocate
						xy2[locationIdx][b][1] 	= new float[count_patches*L]; // y coordinates allocate
						yptch[b] 				= new float[count_patches*L]; // aux variable allocate

						smoothness[locationIdx][b] = 0; 					 // initialize

						// fill  up the allocated arrays for each patch
						for (int m = 0; m<M; m++) {      					// loop patches outwards, from the beginning

							if (delin2[locationIdx][b][m]>=0) {  // will be at least one true because it was used to count the patches in previous loop

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
								float get_smoothness = localPatchRefinedLocs(
															 prev_x, prev_y, curr_x, curr_y,
															 m * L,        // start from
															 xy2[locationIdx][b],
															 yptch[b]);

								// smoothness summed only along the first patch
								if (m==0) smoothness[locationIdx][b] = get_smoothness;

							}
							else break; // because the rest are filled up with -1 or -2 anyway

						}  // xy2 is formed for this branch

					} // branch case

				} // loop branches
				// result is yptch[4][]:
				//  null in case there was nothing to take from
				// array with values in case there was a branch



				/*
				*	smoothness[locationIdx]
				*/
				// use auxilliary variable yptch to calculate smoothness by integrating second derivative
				// this is an improved more robust version that takes into account the spatioal distribution
				// of the refined points along branches
				for (int loop_bches = 0; loop_bches < xy2[locationIdx].length; loop_bches++) {

					//
					//
					if (xy2[locationIdx][loop_bches]==null) {
						smoothness[locationIdx][loop_bches] = Float.NaN; // did it once, just in case
					}
					else {
						// refined coordinates of the branch exist, calculate second derivative

						// before that calculate the distance towards the other refined branch
						// check if it the refined point was erroneously assigned to the other branch
						// criteria: there are at least two points close enough in some other branch
						// in order to conclude that it was captured by that branch - second smallest distance
//						boolean[] use_it = new boolean[L]; // use only first L (one patch)

						ArrayList<float[]> selected_local_xy = new ArrayList<float[]>(L);

						for (int i = 0; i < xy2[locationIdx][loop_bches][0].length; i++) {


							//find smallest second minimum towards the other branch
							float all_2nd_closest = Float.POSITIVE_INFINITY;

							for (int j = 0; j < xy2[locationIdx].length; j++) { // loop other bches
								if (j!=loop_bches && xy2[locationIdx][j]!=null) { // valid other bch

									// j is a valid branch index to check, loop along it, look for 2nd closest

									float bch_2nd_closest = Float.NaN;
									float bch_1st_closest = dist(
																		xy2[locationIdx][loop_bches][0][i],
																		xy2[locationIdx][loop_bches][1][i],

																		xy2[locationIdx][j][0][0],
																		xy2[locationIdx][j][1][0]
									);
									for (int k = 0; k < xy2[locationIdx][j][0].length; k++) {
										float curr_dist = dist(
																	  xy2[locationIdx][loop_bches][0][i],
																	  xy2[locationIdx][loop_bches][1][i],

																	  xy2[locationIdx][j][0][k],
																	  xy2[locationIdx][j][1][k]
																	  );
										if (curr_dist<bch_1st_closest) {
											bch_2nd_closest = bch_1st_closest;
											bch_1st_closest = curr_dist;
										}
									}

									// got 1st and 2nd for bch j

									if (bch_2nd_closest<all_2nd_closest){
										all_2nd_closest = bch_2nd_closest; // add it to overall
									}

								}
							} // loop other bches

							// got best 2nd distance for point i in branch loop_bches
							if (all_2nd_closest<=samplingStepLongitudinal) {
								selected_local_xy.add();
							}

						}

						// finished extracting local_xy for second derivative calculation



						float[] dy = new float[L];

						float[] dx = new float[L];


						// calculate derivative and smoothness






					}

				}
				
				







				/*
				*	vxy2[locationIdx]
				*/
				if (xy2[locationIdx]!=null) {
					for (int b = 0; b<xy2[locationIdx].length; b++) {
						if (xy2[locationIdx][b]!=null) {

							int to_allocate = xy2[locationIdx][b][0].length;
							vxy2[locationIdx][b][0] = new float[to_allocate]; // vx
							vxy2[locationIdx][b][1] = new float[to_allocate]; // vy

							for (int l=0; l<to_allocate; l++) {

								float avg_v_x=0, avg_v_y=0;

								for (int l_nbr=l-windowSize/2; l_nbr<=l+windowSize/2; l_nbr++) {

									if (l_nbr>=0 && l_nbr+1<to_allocate) {

										float x1 = xy2[locationIdx][b][0][l_nbr];
										float y1 = xy2[locationIdx][b][1][l_nbr];
										float x2 = xy2[locationIdx][b][0][l_nbr+1];
										float y2 = xy2[locationIdx][b][1][l_nbr+1];

										float norm = (float) Math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
										float vx = (x2-x1)/norm;
										float vy = (y2-y1)/norm;

										avg_v_x += vy;
										avg_v_y += -vx;

									}

								}

								float norm = (float) Math.sqrt(avg_v_x*avg_v_x + avg_v_y*avg_v_y);
								if (norm>Float.MIN_VALUE) {
									vxy2[locationIdx][b][0][l] = avg_v_x/norm;
									vxy2[locationIdx][b][1][l] = avg_v_y/norm;
								}
								else {
									vxy2[locationIdx][b][0][l] = 0;
									vxy2[locationIdx][b][1][l] = 0;
								}

							}

						}
						else {
							vxy2[locationIdx][b] = null;
						}
					}
				}
				else {
					vxy2[locationIdx] = null;
				}

			} // it was critpoint candidate (selected from foreground points)

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

	private static float localPatchRefinedLocs(float x1, float y1, float x2, float y2,
											  int init_index,
											  float[][] refined_centerline_locs_xy,
											  float[] 	refined_centerline_locs_y)
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

		// varibles necessary for the smoothness calculation
		float sthness = 0; 									// integral of squared second derivatives per dx

		float yptch_3 = 0;
		float yptch_2 = 0;   								// keep the values for the recursion
		float yptch_1 = 0;

		for (int ii=0; ii<L; ii++) {                        // loops L of them in radial direction with 0(root coordinate) and L included

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
			}   // patch_profile formed

			// find local max and take the one with min dist towards center

			// initialize values to optimize in this iteration
			float refined_x = x_root + ii * samplingStepLongitudinal * vx + (-dim_half-1) * samplingStep * vy;
			float refined_y = y_root + ii * samplingStepLongitudinal * wx + (-dim_half-1) * samplingStep * wy;
			float yptch = (-dim_half-1) * samplingStep; // starts from the patch border, keep patch cross-projection along the cross-section

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

						int index_from_center = loop_prof-dim_half-2; // will range from -dim_half to +dim_half  as loop_prof changes
						refined_x = x_root + ii * samplingStepLongitudinal * vx + index_from_center * samplingStep * vy;
						refined_y = y_root + ii * samplingStepLongitudinal * wx + index_from_center * samplingStep * wy;
						yptch = index_from_center * samplingStep;
					}
				}
			}

			// store refined_x, refined_y
			refined_centerline_locs_xy[0][init_index + ii] = refined_x; // ii is radial index
			refined_centerline_locs_xy[1][init_index + ii] = refined_y;
			// store cross y values
			refined_centerline_locs_y[init_index + ii] = yptch;

			// calculate smoothness iteratively, as new values arrive
			if (ii==0) {
				yptch_1 = yptch; // just for initial conditions
			}
			else if (ii==1) {  // case when there is discrepancy in x, there will be exception when calculating
				yptch_2 = yptch;
			}
			else { // ii>=2

				yptch_3 = yptch;
				float ff = ( (yptch_3-yptch_2)/dx - (yptch_2-yptch_1)/dx ) / dx;
				sthness += Math.pow(ff, 2) * dx;

				yptch_1 = yptch_2;
				yptch_2 = yptch_3;

			}

		}  // end looping radial indexes

		return sthness;

	}

	public static ImageStack getSmoothnessDistribution(int nr_bins)
	{

		ArrayList<Float> concatenate_smoothness_scores = new ArrayList<Float>(smoothness.length*4);

		System.out.println(smoothness.length + " is amount of points");

		int cnt = 0;

		for (int i = 0; i < smoothness.length; i++) {

			if (smoothness[i]!=null) {

				for (int j = 0; j < smoothness[i].length; j++) {

					if (smoothness[i][j]!=Float.NaN) {
						concatenate_smoothness_scores.add(smoothness[i][j]);

					}
					else {
						cnt++;
					}


				}

			}


		}

		System.out.println("total " + concatenate_smoothness_scores.size() + " smoothness scores");
		System.out.println("minimum = " + Hist.find_min(concatenate_smoothness_scores));
		System.out.println("maximum = " + Hist.find_max(concatenate_smoothness_scores));

		ImageStack is_out = new ImageStack(528, 255);

		float[] bins = new float[nr_bins];
		float[] distribution = new float[nr_bins];

		Hist.getDistribution(concatenate_smoothness_scores, nr_bins, bins, distribution);
		Plot p = new Plot("", "", "", bins, distribution, Plot.LINE);
		is_out.addSlice("", p.getProcessor());
		return is_out;

	}

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

			if (delin_at_loc!=null) {

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

			} // delin2 != null
			else{
				// do nothing delineatoin was null here because profile was too flat
			}

		}

        /*
        add refined streamline locations (read from xy2)
         */
		int locationIdx = xy2i[atX][atY];

		if (locationIdx!=-1) {

			if (xy2[locationIdx]!=null) {

				for (int b = 0; b<xy2[locationIdx].length; b++) {

					// check whether there is a refinement at all here, add in case there is
					if (xy2[locationIdx][b]!=null) {




						// loop points to add them
						if (xy2[locationIdx][b][0]!=null) {

							for (int l=0; l<xy2[locationIdx][b][0].length; l++) {

								float refined_x = xy2[locationIdx][b][0][l];
								float refined_y = xy2[locationIdx][b][1][l];

								ovalroi = new OvalRoi(refined_x-(samplingStep/2)+.5f, refined_y-(samplingStep/2)+.5f, samplingStep, samplingStep);
								ovalroi.setFillColor(Color.YELLOW);
								ovalroi.setStrokeWidth(samplingStep/4);
								ov.add(ovalroi);

							}


						}
						else {

							System.out.println("xy2 null");


							float old_Rd = Rd;

							Rd = 3*Rd;
							OvalRoi ovalroi_test = new OvalRoi(atX-(Rd/2)+.5f, atY-(Rd/2)+.5f, Rd, Rd);
							ovalroi_test.setFillColor(Color.MAGENTA);
							ovalroi_test.setStrokeWidth(3);
							ov.add(ovalroi_test);
							System.out.println("ADDED");


							Rd = old_Rd;

						}



					}

				}
			}
			else {
				// do nothing, do nothing, points were not extracted because profile was flat
			}

		}

        /*
        add lines marking the refined cross sections
         */
		if (locationIdx!=-1) {

			if (vxy2[locationIdx]!=null) {

				for (int b=0; b<vxy2[locationIdx].length; b++) {

					if (vxy2[locationIdx][b]!=null) {






						// loop points to add them
						if (vxy2[locationIdx][b][0]==null) {
							System.out.println("vxy2 null");

							float old_Rd = Rd;

							Rd = 3*Rd;
							OvalRoi ovalroi_test = new OvalRoi(atX-(Rd/2)+.5f, atY-(Rd/2)+.5f, Rd, Rd);
							ovalroi_test.setStrokeColor(Color.MAGENTA);
							ovalroi_test.setStrokeWidth(3);
							ov.add(ovalroi_test);

							Rd = old_Rd;

						}
						else {

							for (int l=0; l<vxy2[locationIdx][b][0].length; l++) {

								float end_1_x = xy2[locationIdx][b][0][l] - dim_half * samplingStep * vxy2[locationIdx][b][0][l];
								float end_1_y = xy2[locationIdx][b][1][l] - dim_half * samplingStep * vxy2[locationIdx][b][1][l];

								float end_2_x = xy2[locationIdx][b][0][l] + dim_half * samplingStep * vxy2[locationIdx][b][0][l];
								float end_2_y = xy2[locationIdx][b][1][l] + dim_half * samplingStep * vxy2[locationIdx][b][1][l];

								Line lne = new Line(end_1_x+.5f, end_1_y+.5f, end_2_x+.5f, end_2_y+.5f);
								lne.setStrokeColor(Color.CYAN);
								ov.add(lne);

							}


						}












					}

				}

			}
			else {
				// profile was too flat here to extract
			}

		}

		return ov;

	}

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

}