package detection2d;

import aux.Stat;
import detection.Interpolator;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;

import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 12/17/13
 * Time: 10:37 AM
 * Computation unit for 2D processing, contains precomputed values for profile extraction and peak detection
 * indexing, neighbour list (for peaks), offset list (for profiles)
 */
public class Sphere2D {

    private static float    arcRes 	        = 0.7f;
    private static float 	arcNbhood       = 2*arcRes;
    private static float    samplingStep    = 0.7f;
    public static float     R_FULL 	        = 1.00f;
    public static float     T_HALF 	        = 0.50f;
    public static int 		weightStdRatioToD = 4;                  // could be a parameter
    public static int       MAX_ITER        = 15;
    public static int       EPSILON         = 0;

    private float 	radius;
//    private float   scale;
    private float   neuronDiameter;
    private int     N;

    private int 	limR, limT;

	private static float TWO_PI = (float) (2 * Math.PI);
	private static float ONE_PI = (float) (1 * Math.PI);

    //
    // all discretized angle (theta) values will be indexed in a list
    // masks will contain indexes of local neighbours (necessary for peak extraction and clustering)
    // offsets used in oriented filtering are precomputed for each indexed angle
    // there is also a symmetric table that stores precomputed distance differences between each indexed direction (used in clustering methods after converging the indexes)
    //

    private static ArrayList<Float>         theta = new ArrayList<Float>(); 	        // list of elements (theta) covering the circle

    private static ArrayList<int[]> 		masks = new ArrayList<int[]>(); 	    // list of list indexes of the neighbours for each list element

    private static ArrayList<float[][]> 	offstXY = new ArrayList<float[][]>(); 	// list of filter offsets for each direction

    private static float[] weights;

    private static float[][] diffs;

    /*
    *********************************************************************
     */

    public Sphere2D(float neuronDiam, float scale) {

        this.radius 	= scale * neuronDiam;
//        this.scale      = scale;
        this.neuronDiameter = neuronDiam;
        this.N 	= (int) Math.ceil( ( (2 * Math.PI * radius) / arcRes) );    // N will influence theta list size, and offstXY list size
        this.limR = (int) Math.ceil(R_FULL*neuronDiameter/samplingStep);    // how many to take radially with given sampling step
        this.limT = (int) Math.ceil(T_HALF*neuronDiameter/samplingStep);    //

        theta.clear();
        for (int i=0; i<N; i++) {
            theta.add( i * ( (float)(2*Math.PI) / N ) );
        }

        masks.clear();
        for (int ii=0; ii<theta.size(); ii++) {

            float curr_theta = theta.get(ii);

            // loop the rest of the elements to see if they are in the angNeighbourhood range
            ArrayList<Integer> nbrs = new ArrayList<Integer>();
            for (int jj = 0; jj<theta.size(); jj++) {

                if (jj!=ii) {

                    // check
                    float theta_test = theta.get(jj);

                    float arc_btw = arcBetweenDirections(curr_theta, theta_test);
                    if (arc_btw<=arcNbhood) {
                        nbrs.add(jj);
                    }

                }

            }

            // convert list to regular array and add
            int[] nbrsArray = new int[nbrs.size()];
            for (int iii=0; iii<nbrs.size(); iii++) {
                nbrsArray[iii] = nbrs.get(iii);
            }

            masks.add(nbrsArray);

        }

        offstXY.clear();
        weights = new float[(2*limT+1)*(limR+1)];
        float sumWgt = 0;

        for (int ii = 0; ii<theta.size(); ii++) {

            float curr_theta = theta.get(ii);

            /*
				form sampling (offset) points
			 */

            float[][] offsetsPerDirection = new float[(2*limT+1)*(limR+1)][2];

            int cnt = 0;

            for (int k=-limR; k<=0; k++) {

                for (int i=-limT; i<=limT; i++) {

//                    for (int j = -limT; j<=limT; j++) {

                        float px = i * samplingStep;
//                        float py = j * samplingStep;
                        float py = k * samplingStep;

                        offsetsPerDirection[cnt][0] = px;
                        offsetsPerDirection[cnt][1] = py;
//                        offsetsPerDirection[cnt][2] = pz;

                        //*** HERE IT DEFINES THE FILTER PROFILE WEIGHTS (only in first iteration, the rest are the same)
                        if (ii==0) {
                                float dstAxis = point2line(0, 0,        0, 1,       px, py);
                            weights[cnt] = (float) Math.exp(-(dstAxis*dstAxis)/(2*(neuronDiam/weightStdRatioToD)*(neuronDiam/weightStdRatioToD)));
                            sumWgt += weights[cnt];
                        }

                        cnt++;

//                    }

                }

            }

            /*
				transformation for offsets before adding
			 */
            transY(radius, offsetsPerDirection);
            //rotY(-phi+HalfPI, offsetsPerDirection);
            rotZ(curr_theta, offsetsPerDirection);
            offstXY.add(offsetsPerDirection); //store

        }

        /*
				normalize weights
	    */
        for (int iii=0; iii<weights.length; iii++) {
            weights[iii] /= sumWgt;
        }

        /*
                form table with differences (used for clustering)
         */

        diffs = new float[theta.size()][theta.size()];
        for (int i = 0; i < theta.size(); i++) {
            for (int j = i; j < theta.size(); j++) {
                if (i==j) {
                    diffs[i][j] = 0;
                }
                else {
                    float theta1 = theta.get(i);
                    float theta2 = theta.get(j);
                    float dtheta = wrap_diff(theta1, theta2);
                    diffs[i][j] = radius * dtheta;
                    diffs[j][i] = radius * dtheta;
                }
            }
        }

    }

    public ImagePlus showSampling(){

        int DIM = 2 * (int) Math.ceil(Math.sqrt(Math.pow(neuronDiameter*T_HALF,2)+Math.pow(radius, 2))) + 1;
        int CX = DIM/2;
        int CY = CX;

        ImageStack isOut = new ImageStack(DIM, DIM);
        Overlay ov = new Overlay();

        for (int i=0; i<offstXY.size(); i++) {
            isOut.addSlice(new ByteProcessor(DIM, DIM));

            // center
            OvalRoi p = new OvalRoi(CX+0.5 - .5, CY+0.5 -.5, 1, 1);
            p.setPosition(i+1);
            p.setFillColor(Color.RED);
            p.setFillColor(Color.RED);
            ov.add(p);

            // sampling
            for (int i1=0; i1<offstXY.get(i).length; i1++) {

                float offX = offstXY.get(i)[i1][0];
                float offY = offstXY.get(i)[i1][1];

                PointRoi p1 = new PointRoi(CX+offX+.5, CY+offY+.5);
                p1.setPosition(i+1);
                ov.add(p1);

            }

        }

        ImagePlus outIm = new ImagePlus("offsets", isOut);
        outIm.setOverlay(ov);
        return outIm;

    }

    public ImagePlus showWeights(){

        float sum = 0;
        for (int k=0; k<weights.length; k++) sum += weights[k];
        return new ImagePlus("weights", new FloatProcessor((2*limT+1), (limR+1), weights));

    }

    public int getProfileLength()
    {
        return offstXY.size();
    }

    public int getOuterRadius(){
        return (int) Math.ceil(Math.sqrt(Math.pow(neuronDiameter*T_HALF, 2)+Math.pow(radius, 2)));
    }

    public short extractProfile(int profileIdx, float atX, float atY, float[][] _inimg_xy) {
        // one element filter output (indexed with profileIdx)

        float value = 0;

        for (int offsetIdx=0; offsetIdx<offstXY.get(profileIdx).length; offsetIdx++) {

            float x_offs_pix = offstXY.get(profileIdx)[offsetIdx][0];
            float y_offs_pix = offstXY.get(profileIdx)[offsetIdx][1];

            float imgVal = Interpolator.interpolateAt(atX + x_offs_pix, atY + y_offs_pix, _inimg_xy); // , atZ + z_offs_lay
            value += weights[offsetIdx] * imgVal;

        }

        return  (short) ((int) ((value/255f)*65535f)); // &  0xffff

    }

	public void peakCoords_4xXY(short[] _profile,                           	// main input
                                int[] start_pts, int[] end_pts,          		// aux. arrays (to avoid allocating them inside method each time) - end_pts is also an output
                                int atX, int atY,                           	// for global coordinate outputs
                                float[][] _inimg_xy, int[][] _xy2i,         	// to calculate median() along line
                                int[][] peaks_loc_xy, float[][] peaks_ang_theta){	// outputs

        // this is the block that extracts 4 strongest peaks in global 2d coordinates (that's why atX, atY at the input) and profile indexes
        // ranking can be based on median between locations of number of points that converged to each cluster
        // in order to rank using medians - image array, look-up tables are added to be able to compute medians
        // otherwise ranking can be based on number of convergence points - then all the image data is not necessary (position is enough)
        // this implementation will have the number of convergence points but use median along the line as ranking estimate finally
        runLineSearch(
                start_pts,
                _profile,
                MAX_ITER,
                EPSILON,
                end_pts);

        // cluster end_pts together
		int[] labs = clustering(end_pts, diffs, arcNbhood);

		// extract the cluster centroids out  -> <theta, weight>
		ArrayList<float[]> clss = extracting(labs, end_pts, theta);

		appendClusters(atX, atY, _xy2i, _inimg_xy, clss, peaks_loc_xy, peaks_ang_theta);

	}

    /*
    *********************************************************************
     */

	private void appendClusters(
									   int atX, int atY,
									   int[][] 		_xy2i,  						// lookup table
									   float[][] 	_inimg_xy,                     	// input image
									   ArrayList<float[]> clusters_to_append,       // <theta, weight>
									   int[][] 		destination_locs,               // 4x2
									   float[][] 	destination_angs)               // 4x1
	{
		// take top 4 and store them according to the importance (criteria: median along the connecting line, or convergence iterations)
		float[] medAlongLin = new float[destination_locs.length];
		Arrays.fill(medAlongLin, -1f);

		for (int t=0; t<clusters_to_append.size(); t++) { // check every peak theta angle

			float peak_theta = clusters_to_append.get(t)[0];   // value in [rad]
			//float peak_theta_weight = clusters_to_append.get(t)[1];

			float x_peak_pix = atX + getX(radius, peak_theta);
			float y_peak_pix = atY + getY(radius, peak_theta);

			int x_peak_pix_base = (int) Math.floor(x_peak_pix);
			int y_peak_pix_base = (int) Math.floor(y_peak_pix);

			/*
				pick the best out of 4 neighbours (compensate location inaccuracy)
			 */

			float maxWeight = Float.MIN_VALUE;
			int[] currentPeaksXY   = new int[2];

			boolean modified = false;

			for (int ii = 0; ii <= 1; ii++) { // check around peak
				for (int jj = 0; jj <= 1; jj++) {

					int check_x = x_peak_pix_base + ii;
					int check_y = y_peak_pix_base + jj;

					if (check_x>=0 && check_x<_inimg_xy.length && check_y>=0 && check_y<_inimg_xy[0].length) {

						float currMedian = medianAlongLine(	atX, atY, check_x, check_y, _inimg_xy );

						if (_xy2i[check_x][check_y]!=-1)  {  // ensures retrieval from the list

							if (currMedian>maxWeight) {
								// set this one as peak
								currentPeaksXY[0] = check_x;
								currentPeaksXY[1] = check_y;

								modified = true;

								// update maxMedian
								maxWeight = currMedian;
							}

						}
					}
				}
			}

			if (modified) {

			    // if there was a result add it to sorted list
				// insert maxMedian and currentPeaksXYZ to the list of 4 (destination_array)
				for (int k = 0; k < 4; k++) {

					if (medAlongLin[k] == -1f) {            			// store immediately

						destination_locs[k][0]    = currentPeaksXY[0];
						destination_locs[k][1]    = currentPeaksXY[1];

						medAlongLin[k]  = maxWeight;

						destination_angs[k][0]     = peak_theta;

						break;

					}
					else if (maxWeight > medAlongLin[k]) {

						// shift the rest first
						for (int kk = 4-2; kk>=k; kk--) {   // shift them from the one before last, shift back (last one dissapears)

							destination_locs[kk+1][0] = destination_locs[kk][0];
							destination_locs[kk+1][1] = destination_locs[kk][1];

							medAlongLin[kk+1] = medAlongLin[kk];

							destination_angs[kk+1][0] = destination_angs[kk][0];

						}

						// store it at k
						destination_locs[k][0] = currentPeaksXY[0];
						destination_locs[k][1] = currentPeaksXY[1];

						medAlongLin[k] = maxWeight;

						destination_angs[k][0] = peak_theta;

						break;

					}
					else {

						// if smaller, loop further

					}

				}

			}// modified

		}


	}

    private void rotZ(float ang, float[][] coords) {
        for (int i=0; i<coords.length; i++) {
            float x_temp = coords[i][0];
            float y_temp = coords[i][1];
            coords[i][0] = x_temp * (float) Math.cos(ang) - y_temp * (float) Math.sin(ang);
            coords[i][1] = x_temp * (float) Math.sin(ang) + y_temp * (float) Math.cos(ang);
        }
    }

    private void transY(float ty, float[][] coords){
        for (int i=0; i<coords.length; i++){
            coords[i][1] += ty;
        }
    }

    private float point2line(float n1x, float n1y,  // float n1z,
                             float n2x, float n2y,  // float n2z,
                             float px,  float py    //, float pz
    )
    {

        float d = 0;

        double[] p_b = new double[2];

        //double[] n21 = new double[3];
        float n21Len = (float) Math.sqrt(Math.pow(n2x-n1x,2)+Math.pow(n2y-n1y,2)); // +Math.pow(n2z-n1z,2)
        float n21x = (n2x-n1x)/n21Len;
        float n21y = (n2y-n1y)/n21Len;
//        float n21z = (n2z-n1z)/n21Len;

        float proj = (px - n1x) * n21x + (py - n1y) * n21y; // + (pz - n1z) * n21z; // dot prod

        p_b[0] = -(px - n1x) + proj * n21x;
        p_b[1] = -(py - n1y) + proj * n21y;
//        p_b[2] = -(pz - n1z) + proj * n21z;

        return (float) Math.sqrt(p_b[0]*p_b[0] + p_b[1]*p_b[1]); // + p_b[2]*p_b[2]

    }

    private float arcBetweenDirections(float theta1, float theta2){

        float x1 = getX(radius, theta1);
        float y1 = getY(radius, theta1);

        float x2 = getX(radius, theta2);
        float y2 = getY(radius, theta2);

        return radius * (float) Math.acos( (x1*x2+y1*y2)/(radius * radius) );

    }

    private float getX(float r, float theta){return (-1) * r * (float) Math.sin(theta);}

    private float getY(float r, float theta){return (+1) * r * (float) Math.cos(theta);}

    private static int runOneMax(int curr_pos, short[] _input_profile) {
        int 	    new_pos     = curr_pos;
        int 		max	 		= Integer.MIN_VALUE;

        // curr_pos will define the set of neighbouring indexes
        int[] neighbour_idxs = masks.get(curr_pos);

        for (int i=0; i<neighbour_idxs.length; i++) {
            int neighbour_idx = neighbour_idxs[i];
            int read_value = (int) (_input_profile[neighbour_idx] & 0xffff);
            if (read_value>max) {
                max = read_value;
                new_pos = neighbour_idx;
            }
        }

        return new_pos;
    }

    private static void runLineSearch(
            int[] 	    start_idxs,
            short[] 	input_profile,
            int 		max_iter,
            int 	    epsilon,
            int[] 	    end_idxs // same length as start (this would be output)
    )
    {

        // initialize output array
        for (int i = 0; i < end_idxs.length; i++) {
            end_idxs[i] = start_idxs[i];
        }

        // iterate each of the elements of the output
        for (int i = 0; i < end_idxs.length; i++) {

            int iter = 0;
            int d;

            do{

                int new_pos = runOneMax(end_idxs[i], input_profile);
                int pre_value = input_profile[end_idxs[i]] & 0xffff;
                int new_value = input_profile[new_pos]     & 0xffff;
                d = Math.abs(new_value - pre_value);
                end_idxs[i] = new_pos;
                iter++;
            }
            while(iter < max_iter && d > epsilon);

        }

    }

    private static int[] clustering(int[] idxs, float[][] dists, float threshold_dists)   // essentially clustering of the indexes
    {
        // indxs represent indexes of values that need to be clustered
        // intended to place here indexes after the convergence
        // dists are the distances
        // threshold_dists is the distance limit
        // output is list of unique labels

        int[] labels = new int[idxs.length];
        for (int i = 0; i < labels.length; i++) labels[i] = i;  // initialize the output

        //System.out.println("INIT. LABELS:");
        //System.out.println(Arrays.toString(labels));

        for (int i = 0; i < idxs.length; i++) {

            // one versus the rest
            for (int j = 0; j < idxs.length; j++) {

                // check the rest of the values
                if (i != j) {

                    int idx_i = idxs[i]; // will be used to read diff from the table
                    int idx_j = idxs[j]; //

                    if (dists[idx_i][idx_j]<=threshold_dists) {

                        if (labels[j] != labels[i]) {
                            // propagate the label
                            int currLabel = labels[j];
                            int newLabel  = labels[i];

                            labels[j] = newLabel;

                            //set all that also were currLabel to newLabel
                            for (int k = 0; k < labels.length; k++)
                                if (labels[k]==currLabel)
                                    labels[k] = newLabel;

                        }

                    }

                }

            }

        }

        //System.out.println("OUT LABELS:");
//		for (int ii = 0; ii < labels.length; ii++)
//			System.out.print(labels[ii]+" ");
        //System.out.println(Arrays.toString(labels));

        return labels;

    }

    public static ArrayList<float[]> extracting(int[] labels, int[] idxs, ArrayList<Float> vals) {

        // loop obtained labels (labels & idxs should have the same length)

        boolean[] checked = new boolean[idxs.length];      // aux
        ArrayList<float[]> out = new ArrayList<float[]>(); // allocate the output

        for (int i = 0; i < idxs.length; i++) {
            if (!checked[i]) {
                // this is the first value
                float centroid  = vals.get(idxs[i]); // vals[ idxs[i] ];
                float shifts    = 0; // naturally 0 shift for the first one
                int count = 1;
                checked[i] = true;

                // check the rest
                for (int j = i+1; j < idxs.length; j++) {
                    if (!checked[j]) {
                        if (labels[j]==labels[i]) {

                            // clustering said they're together
							// important that vals and centroid values are all wrapped in [0, 2PI) range
                            float add_value = vals.get(idxs[j]);
							float add_diff = wrap_diff(centroid, add_value);
							if (centroid<ONE_PI) {
								// there is always pi values on right
								if (add_value>=centroid & add_value<centroid+ONE_PI) {
									// sign is (+)
								}
								else {
									// sign is (-)
									add_diff = (-1) * add_diff;
								}

							}
							else {
								// >= ONE_PI
								// there is always pi values on left
								if (add_value<=centroid & add_value>centroid-ONE_PI) {
									// sign is (-)
									add_diff = (-1) * add_diff;
								}
								else {
									// sign in (+)
								}
							}

							shifts += add_diff;
                            count++;
                            checked[j] = true;

                        }
                    }
                }

                centroid += shifts/count;
//                centroid = ;
                out.add(new float[]{wrap_0_2PI(centroid), count});  // outputs centroid in [rad]

            }
        }

        return out;

    }

    private static float wrap_diff(float theta_1, float theta_2) { // wraps the angle difference theta_1, theta_2 range [0, 2PI)

        float d = theta_1 - theta_2;

        d = (d>0)? d : (-d) ;

        d = (d>Math.PI)? (float) (2*Math.PI-d) : d ;

        return d;

    }

	private static float wrap_0_2PI(float ang) {
		float out = ang;
		while (out<0) {
			out += TWO_PI;
		}
		while (out>=TWO_PI) {
			out -= TWO_PI;
		}
		return out;
	}

	private float medianAlongLine(float x1, float y1, float x2, float y2, float[][] inimg_xy) {

		// DEBUG:
		//System.out.println(String.format("median_along_line(%4.2f, %4.2f)->(%4.2f, %4.2f)", x1, y1, x2, y2));

        float increment_length = .7f;

		int elementsInLine = (int) (radius / increment_length);  // how many increment can safely fit between
		float[] valuesAlongLine = new float[elementsInLine];

		float dist = (float) Math.sqrt(Math.pow(x2-x1, 2) + Math.pow(y2-y1, 2)); //  + Math.pow(z2lay-z1lay, 2)

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

}
