package detection2d;

import ij.gui.Plot;
import ij.process.ImageProcessor;

import java.awt.*;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 12/17/13
 * Time: 3:50 PM
 * Will extract the peaks of the profiles, parallel threaded implementation
 */
public class PeakExtractor2D extends Thread {

    private int begN, endN;

    public static Sphere2D      sph2d;

    public static float[][]     inimg_xy;

    public static int[][] 	    i2xy;                       // selected locations

    public static int[][]     	xy2i;

    public static short[][]	    extracted_profiles;         // profiles

    // outputs are lists of selected peaks
    public static float[][][]	peaks_theta;                // N x (4x1)   4 selected peaks in abscissa coordinates X
    public static int[][][]     peaks_xy;             		// N x (4x2)   main output  4 selected peaks in XY format (TODO: can be only index to save memory)
    // peaks_xy are the best rounded locations found close to the spots determined with peaks_theta
    // TODO MAYBE: add another output - peak strength (useful for plot and for later tracing usage)

    public PeakExtractor2D(int n0, int n1)
    {
        this.begN = n0;
        this.endN = n1;
    }

	public static void loadTemplate(Sphere2D _sph2d, int[][] _i2xy, short[][] _extracted_profiles, float[][] _inimg_xy, int[][] _xy2i) {

		sph2d           = _sph2d;
		inimg_xy        = _inimg_xy;                        // just necessary to rank peaks
		i2xy      		= _i2xy;
		extracted_profiles = _extracted_profiles;
		xy2i 			= _xy2i;

		// allocate output -> set to -1
		peaks_xy  		= new int[i2xy.length][4][2];
		peaks_theta  	= new float[i2xy.length][4][1];

		for (int ii = 0; ii<i2xy.length; ii++) {
			for (int jj = 0; jj<4; jj++) {

				for (int kk = 0; kk<2; kk++) {
					peaks_xy[ii][jj][kk] = -1;
				}

				for (int kk=0; kk<1; kk++) {
					peaks_theta[ii][jj][kk] = -1;
				}

			}
		}

	}

	public void run()
	{
        // allocate variables just for this thread range
        int[] start_indexes = new int[sph2d.getProfileLength()];
        for (int i = 0; i < start_indexes.length; i++) start_indexes[i] = i;    // zero indexing is used
        int[] end_indexes = new int[sph2d.getProfileLength()];                  // zeros at the beginning
        // end allocation

		for (int locationIdx = begN; locationIdx < endN; locationIdx++) {

			int atX = i2xy[locationIdx][0];
			int atY = i2xy[locationIdx][1];

			sph2d.peakCoords_4xXY(
                    extracted_profiles[locationIdx],
                    start_indexes, end_indexes,
                    atX, atY,
                    inimg_xy, xy2i,
                    peaks_xy[locationIdx], peaks_theta[locationIdx]); // peaks_* will fill up (hopefully)

		}

	}

    public static ImageProcessor getProfileWithPeaks(int atX, int atY, int[][] _xy2i){

        // reads from prof2 array, Profiler2D static class member
        // reads from peaks_*, PeakExtractor2D class member

        int idx = _xy2i[atX][atY];
        if (idx != -1) {
            int len = Profiler2D.prof2[0].length;
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

            Plot p = new Plot("profile at ("+atX+","+atY+")", "", "filtered", fx, f);
            p.setSize(600, 300);

            // add detected peaks on top of the profile plot using addPoints
            float[][] get_thetas = peaks_theta[idx];  // 4 x 1

            for (int k=0; k<get_thetas.length; k++) {
                if (get_thetas[k][0] != -1) {
                    float curr_theta_degs = rad2deg(get_thetas[k][0]);
                    float[][] pks_abscissa = new float[][]{ {curr_theta_degs, curr_theta_degs}, {profile_min, profile_max} }; // 0 -> abscissa, 1 -> limits

                    Color c = Color.BLACK;
                    if(k==0) {
                        c = Color.RED;
                    }
                    else if (k==1) {
                        c = Color.YELLOW;
                    }
                    else if (k==2) {
                        c = Color.GREEN;
                    }
                    else if (k==3) {
                        c = Color.BLUE;
                    }

                    p.setColor(c);
                    p.addPoints(pks_abscissa[0], pks_abscissa[1], Plot.LINE);
                    p.setColor(Color.BLACK);

                }
            }

            return p.getProcessor();
        }
        else {
            return null;
        }

    }

    private static final float rad2deg(float ang_rad){
        return (ang_rad / (float) Math.PI) * 180f;
    }

    private static final float deg2rad(float ang_deg) {
        return (ang_deg / 180f) * (float) (Math.PI);
    }

}
