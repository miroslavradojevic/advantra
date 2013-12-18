package detection2d;

import detection3d.Sphere3D;

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
    public static int[][][]     peaks1;                     // N x (4x1)   4 selected peaks in abscissa coordinates X
    public static int[][][]     peaks2;                     // N x (4x2)   main output  4 selected peaks in XY format
    // TODO add another output - peak strength (useful for plot and for later tracing usage)

    public PeakExtractor2D(int n0, int n1)
    {
        this.begN = n0;
        this.endN = n1;
    }

	public static void loadTemplate(Sphere2D _sph2d, int[][] _i2xy, short[][] _extracted_profiles, float[][] _inimg_xy, int[][] _xy2i) {

		sph2d           = _sph2d;
		inimg_xy        = _inimg_xy;                   // just necessary to rank peaks
		i2xy      		= _i2xy;
		extracted_profiles = _extracted_profiles;
		xy2i 			= _xy2i;

		// allocate output -> set to -1
		peaks1  = new int[i2xy.length][4][1];
		peaks2  = new int[i2xy.length][4][2];

		for (int ii = 0; ii<i2xy.length; ii++) {
			for (int jj = 0; jj<4; jj++) {

				for (int kk = 0; kk<2; kk++) {
					peaks2[ii][jj][kk] = -1;
				}

				for (int kk=0; kk<1; kk++) {
					peaks1[ii][jj][kk] = -1;
				}

			}
		}

	}

	public void run()
	{
        // variables just for this thread range
        int[] start_indexes = new int[sph2d.getProfileLength()];
        for (int i = 0; i < start_indexes.length; i++) start_indexes[i] = i; // zero indenxing is used

        int[] end_indexes = new int[sph2d.getProfileLength()];        // zeros at the beginning

		for (int locationIdx = begN; locationIdx < endN; locationIdx++) {

			int atX = i2xy[locationIdx][0];
			int atY = i2xy[locationIdx][1];

			sph2d.peakCoords_4xXY(
                    extracted_profiles[locationIdx],
                    start_indexes, end_indexes,
                    atX, atY,
                    inimg_xy, xy2i,
                    peaks2[locationIdx], peaks1[locationIdx]
            ); // peaks1,2 will fill up

		}

	}

}
