package detection2d;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.PointRoi;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 12/17/13
 * Time: 3:50 PM
 * Will extract the peaks of the profiles, parallel threaded implementation
 */
public class PeakExtractor2D extends Thread {

    private int begN, endN;

    public static Sphere2D      sph2d;                      //
    public static float[][]     inimg_xy;                   //
    public static int[][] 	    i2xy;                       // selected locations
    public static int[][]     	xy2i;                       //
    public static short[][]	    prof2;                      // profiles

    public static int       MAX_ITER        = 50;           //
    public static int       EPSILON         = 0;            //

    private static float TWO_PI = (float) (2 * Math.PI);
    private static float ONE_PI = (float) (1 * Math.PI);

    // OUTPUTS
    public static float[][][]	peaks_theta;                // N x 4 x 1    4 selected peaks in abscissa coordinates X
    public static int[][]       peaks_i;             		// N x 4        4 selected peaks in indexed format (rounded locations)
    public static int[][]       peaks_w;                    // N x 4        peak weights

	// circular statistics
	public static float[]       R_;							// N            mean resultant length
	public static float[]       theta_;						// N            mean resultant length
	public static float[]       R_;							// N            mean resultant length
	public static float[]       R_2;						// N            mean resultant length
	public static float[]       entropy;                    // N            entropy of the profile at this location

	public static float[]       kurtosis;                   // N            kurtosis of the profile at this location

    public PeakExtractor2D(int n0, int n1)
    {
        this.begN = n0;
        this.endN = n1;
    }

	public static void loadTemplate(Sphere2D _sph2d, int[][] _i2xy, short[][] _prof2, float[][] _inimg_xy, int[][] _xy2i)
	{

		sph2d           = _sph2d;
		inimg_xy        = _inimg_xy;
		i2xy      		= _i2xy;
		prof2           = _prof2;
		xy2i 			= _xy2i;

		// allocate output -> set to -1
		peaks_i  		= new int[i2xy.length][4];
		peaks_theta  	= new float[i2xy.length][4][1];
		peaks_w 		= new int[i2xy.length][4];
        entropy         = new float[i2xy.length];
        kurtosis        = new float[i2xy.length];

        // initialization legend:
        // -2 is assuming that it is not assigned, empty location (all are initialized that way)
        // -1 will be added if it is in the background
        // >=0 in case it is some foreground location
		for (int ii = 0; ii<i2xy.length; ii++) {
			for (int jj = 0; jj<4; jj++) {

                peaks_i[ii][jj] = -2;

				for (int kk=0; kk<1; kk++) {
					peaks_theta[ii][jj][kk] = -2;
				}

				peaks_w[ii][jj] = -2;

			}
		}

	}

	public void run()
	{

		// auxiliary variables, just for this thread interval
        int[] start_indexes = new int[sph2d.getProfileLength()];
        for (int i = 0; i < start_indexes.length; i++) start_indexes[i] = i;    // zero indexing is used
        int[] end_indexes = new int[sph2d.getProfileLength()];                  // zeros at the beginning
		// array to hold the p-ty distribution for every direction (angle) in profiles
		// will be redefined for every location loop
		float[] p_distr = new float[prof2[0].length];
		float profile_mass = 0;

		for (int locationIdx = begN; locationIdx < endN; locationIdx++) { // all foreground locations

			int atX = i2xy[locationIdx][0];
			int atY = i2xy[locationIdx][1];

            extractPeaks(
                    prof2[locationIdx],
                    sph2d,
                    start_indexes, end_indexes,
                    atX, atY,
                    peaks_i[locationIdx],      // indexed location
                    peaks_theta[locationIdx],  // angle in radians 0 to 2pi
					peaks_w[locationIdx]
			);

            // prepare
			profile_mass = 0;  // set it to zero at every location
            for (int ii=0; ii<prof2[locationIdx].length; ii++)  {  // one loop to calculate the mass
				p_distr[ii] = ((prof2[locationIdx][ii] & 0xffff) / 65535f) * 255f; // take the float value 0-255 range
                profile_mass += p_distr[ii];
			}
			for (int ii=0; ii<prof2[locationIdx].length; ii++) {
				p_distr[ii] /= profile_mass; // so that they sum up to 1
			}

			/* entropy[locationIdx] */
			entropy[locationIdx] = 0;
            for (int ii=0; ii<prof2[locationIdx].length; ii++) {
                entropy[locationIdx] += p_distr[ii] * Math.log(p_distr[ii]);
            }
            entropy[locationIdx] = -entropy[locationIdx];

            /* kurtosis[locationIdx] */


		}

	}

    private void extractPeaks(  short[]     _profile,           // profile input
                                Sphere2D    _profile_sphere,	// sphere used for this profile (defines resolution and the metrics)
                                int[]       start_pts,          //
                                int[]       end_pts,          	// aux. arrays (to avoid allocating them inside method each time) - end_pts is also an output
                                int         atX,                //
                                int         atY,                //
                                int[]       peaks_loc_i,        // out
                                float[][]   peaks_ang_theta,    // out [radians 0 to 2pi]
								int[]		peaks_weight        // out
    )
    {

        // this is the block that extracts 4 strongest peaks in global 2d coordinates (that's why atX, atY at the input) and profile indexes
        // ranking can be based on number of convergence points - then all the image data is not necessary (position is enough)
        // this implementation will have the number of convergence points but use median along the line as ranking estimate finally
        runLineSearch(
                start_pts,
                _profile,
                _profile_sphere,
                MAX_ITER,
                EPSILON,
                end_pts);

        int[] labs = clustering(end_pts, _profile_sphere.diffs, 2*_profile_sphere.arcNbhood);// cluster the end_pts together

        // extract the cluster centroids out  -> clss will be list of <theta[radians 0 to 2pi], weight(#convergence points)>
        ArrayList<float[]> clss = extracting(labs, end_pts, _profile_sphere.theta);

        appendClusters(atX, atY, _profile_sphere, clss,
							  peaks_loc_i,		// out
							  peaks_ang_theta,  // out
							  peaks_weight      // out
		);

    }

    private static void runLineSearch(
            int[] 	    start_idxs,
            short[] 	input_profile,
            Sphere2D    input_profile_sphere,
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

                int new_pos = runOneMax(end_idxs[i], input_profile, input_profile_sphere);
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

//      System.out.println("OUT LABELS:");
//		for (int ii = 0; ii < labels.length; ii++)
//	    System.out.print(labels[ii]+" ");
//      System.out.println(Arrays.toString(labels));

        return labels;

    }

    private void appendClusters(
            int                 atX,
            int                 atY,
            Sphere2D            sphere_atXY,
            ArrayList<float[]>  clusters_to_append,         // <theta, weight>
            int[]  		        destination_locs,           // 4x1
            float[][] 	        destination_angs,            // 4x1
			int[] 				destination_weights
    )
    {

        int[] weights = new int[destination_locs.length];
        Arrays.fill(weights, -1);

        for (int t=0; t<clusters_to_append.size(); t++) { // check every peak theta angle

            float 	peak_theta    = clusters_to_append.get(t)[0];   				// value in [rad] theta 0 to 2pi
            int 	peak_weight   = Math.round(clusters_to_append.get(t)[1]);   	// convergence score

            int x_peak_pix_rounded = Math.round(atX + sphere_atXY.getX(peak_theta));
            int y_peak_pix_rounded = Math.round(atY + sphere_atXY.getY(peak_theta));

            if (xy2i[x_peak_pix_rounded][y_peak_pix_rounded]==-1) {
                /*
                 belongs to the mask defined background
                  */

//                boolean added = false;

                // add it to the first available place if there is such, give priority to those from the foreground
                for (int k = 0; k < 4; k++) {

                    if (weights[k] == -1) { // store on first available location

                        destination_locs[k]         = -1;//xy2i[x_peak_pix_rounded][y_peak_pix_rounded];
                        destination_angs[k][0]      = -1;
                        destination_weights[k]      = peak_weight;

                        weights[k] 					= peak_weight;

//                        added = true;
                        break;
                    }

                }
                // if all th locations are filled up, just ignore the background point
                // there is enough to make the decision so that the background peaks are not disturbance

            }
            else {
                /*
                 belongs to the mask defined foreground
                  */

                boolean added = false;

                for (int k = 0; k < 4; k++) { // add it to the first available place

                    if (weights[k] == -1) { // store immediately, location is available

                        destination_locs[k]         = xy2i[x_peak_pix_rounded][y_peak_pix_rounded];
                        destination_angs[k][0]      = peak_theta;     // radians 0 to 2pi
                        destination_weights[k]      = peak_weight;

                        weights[k] 					= peak_weight;

                        added = true;
                        break;
                    }

                }

                if (!added) { // no available slot

                    // loop once more and put it instead of the one with lowest weight
                    int min_weight_idx 	= -1;
                    int min_weight 		= Integer.MAX_VALUE;

                    for (int kk=0; kk<4; kk++) {
                        if (weights[kk]<min_weight) {
                            min_weight = weights[kk];
                            min_weight_idx = kk;
                        }
                    }

                    // add it instead
                    destination_locs[min_weight_idx] = xy2i[x_peak_pix_rounded][y_peak_pix_rounded];
                    destination_angs[min_weight_idx][0] = peak_theta;
                    destination_weights[min_weight_idx] = peak_weight;

                    weights[min_weight_idx] = peak_weight;

                }

            }

            // otherwise the value of the peak stays -2 in case it was not filled up (peak did not exist, value given at initialization stays)

        }

    }

    public static ArrayList<float[]> extracting(int[] labels, int[] idxs, ArrayList<Float> vals)
    {

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
                if (count>1)
                    out.add(new float[]{wrap_0_2PI(centroid), count});  // outputs centroid in [rad]

            }
        }

        return out; // list <[float, float]> : <direction angle[rad], number of samples after converging>

    }


    public static ImageStack getProfileWithPeaks(int atX, int atY)
    {

        // reads from prof2 array
        ImageStack is_out = new ImageStack(528, 255);

        int idx = xy2i[atX][atY];
        if (idx != -1) {
            int len = prof2[0].length;
            float[] f = new float[len];
            float[] fx = new float[len];

            float profile_min = Float.MAX_VALUE;
            float profile_max = Float.MIN_VALUE;

            for (int i=0; i<len; i++) {
                f[i] = ((Profiler2D.prof2[idx][i] & 0xffff) / 65535f) * 255f; // retrieve the profile from static array
                fx[i] = (i / (float) len) * 360; // abscissa in degs

                if (f[i]>profile_max) profile_max = f[i];
                if (f[i]<profile_min) profile_min = f[i];
            }

            Plot p = new Plot("", "", "entropy="+IJ.d2s(entropy[idx],2), fx, f);

            float[][] get_thetas = peaks_theta[idx];
//            IJ.log("peaks theta:");
//            for (int ii=0; ii<get_thetas.length; ii++) IJ.log(Arrays.toString(get_thetas[ii]));
//            int[] get_weights = peaks_w[idx];
//            IJ.log("peaks wieghts:");
//            for (int ii=0; ii<get_weights.length; ii++) IJ.log(IJ.d2s(get_weights[ii], 2));
//            int[] get_idxs = peaks_i[idx];
//            IJ.log("peaks idxs:");
//            for (int ii=0; ii<get_idxs.length; ii++) IJ.log(IJ.d2s(get_idxs[ii], 2));


            for (int k=0; k<get_thetas.length; k++) {
                if (get_thetas[k][0] != -1) {
                    float curr_theta_degs = rad2deg(get_thetas[k][0]);
                    float[][] pks_abscissa = new float[][]{ {curr_theta_degs, curr_theta_degs}, {profile_min, profile_max} }; // 0 -> abscissa, 1 -> limits
                    p.setColor(Color.RED);
                    p.setLineWidth(2);
                    p.addPoints(pks_abscissa[0], pks_abscissa[1], Plot.LINE);
                    p.setColor(Color.BLACK);
                }
            }

            is_out.addSlice(p.getProcessor());
        }
        else {
            float[] fx = new float[sph2d.getProfileLength()];
            for (int i=0; i<fx.length; i++) fx[i] = ((i/(float)fx.length)*360);
            Plot p = new Plot("", "", "", fx, new float[fx.length]);
            is_out.addSlice(p.getProcessor());
        }

        return is_out;

    }

    public static Overlay getPeaks(int atX, int atY)
    {

        Overlay ov = new Overlay();

        float R = 0.5f;
        OvalRoi ovalroi = new OvalRoi(atX-(R/2)+.5f, atY-(R/2)+.5f, R, R);
        ov.add(ovalroi);

        int idx = xy2i[atX][atY];

        if (idx!=-1) {
            int[] pk_locs_i = peaks_i[idx];

            for (int i = 0; i<pk_locs_i.length; i++) {

                if (pk_locs_i[i]!=-1) { // peak exists

                    int pk_x = i2xy[pk_locs_i[i]][0];
                    int pk_y = i2xy[pk_locs_i[i]][1];

                    ovalroi = new OvalRoi(pk_x-(R/2)+.5f, pk_y-(R/2)+.5f, R, R);
                    ovalroi.setFillColor(Color.GREEN);

                }

                ov.add(ovalroi);
            }

        }

        return ov;

    }

    public static void getPeaks(int atX, int atY, int N, Overlay out_ov)
    {
        // recursively plot peaks of the peaks
        // if peak was foreground then it's plotted in (in current implementation, background peaks are not saved, memorized as -1 index)
        float R = 0.5f;
        Color c = Color.GREEN;
        float w = .25f;

        if (N==0) {

            int point_idx = xy2i[atX][atY];

            if (point_idx>=0) {
                OvalRoi pt = new OvalRoi(atX-(R/2)+.5f, atY-(R/2)+.5f, R, R);
                pt.setFillColor(c);
                pt.setStrokeWidth(w);
                out_ov.add(pt);
            }
        }
        else {

            int point_idx = xy2i[atX][atY];

            // add current point
            if (point_idx>=0) {

                OvalRoi pt = new OvalRoi(atX-(R/2)+.5f, atY-(R/2)+.5f, R, R);
                pt.setFillColor(c);
                pt.setStrokeWidth(w);
                out_ov.add(pt);

                for (int p = 0; p<peaks_i[point_idx].length; p++) {

                    int peak_idx = peaks_i[point_idx][p];

                    if (peak_idx>=0) {

                        // corresponding x, y
                        int peak_x = i2xy[peak_idx][0];
                        int peak_y = i2xy[peak_idx][1];

                        getPeaks(peak_x, peak_y, N-1, out_ov);

                    }

                }
            }
        }
    }

    public static ImagePlus getEntropy(){

        int w = inimg_xy.length;
        int h = inimg_xy[0].length;

        float[][] t = new float[w][h];

        float min_entropy = Float.NEGATIVE_INFINITY;
        float max_entropy = Float.POSITIVE_INFINITY;

        for (int ii = 0; ii<entropy.length; ii++) {
            if (entropy[ii]>max_entropy) max_entropy = entropy[ii];
            if (entropy[ii]<min_entropy) min_entropy = entropy[ii];
        }

        for (int xx=0; xx<w; xx++) {
            for (int yy=0; yy<h; yy++) {
                int idx = xy2i[xx][yy];
                if (idx!=-1) {
                    t[xx][yy] = entropy[idx];
                }
                else {
                    t[xx][yy] = max_entropy;
                }
            }
        }

        return new ImagePlus("profile_entropy", new FloatProcessor(t));

    }

    private static int runOneMax(int curr_pos, short[] _input_profile, Sphere2D _input_profile_sphere) {
        int 	    new_pos     = curr_pos;
        int 		max	 		= Integer.MIN_VALUE;

        // curr_pos will define the set of neighbouring indexes
        int[] neighbour_idxs = _input_profile_sphere.masks.get(curr_pos);

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

    private static final float rad2deg(float ang_rad){
        return (ang_rad / (float) Math.PI) * 180f;
    }

    private static final float deg2rad(float ang_deg) {
        return (ang_deg / 180f) * (float) (Math.PI);
    }

}
