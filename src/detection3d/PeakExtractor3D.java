package detection3d;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 12/6/13
 * Time: 3:42 PM
 */
public class PeakExtractor3D extends Thread {

    private int begN, endN;

    public static Sphere3D 	    sph3;                       // processing unit

    public static float[][][]   img3_zxy;                   // link to input image array

    public static int[][] 	    listLocs3D;                 // N x 3 (Z,X,Y), list of locations

    public static float         zDist;                      // to properly calibrate layer

    public static int[][][]     peaks3;                     // N x (4x3) main output

    public static short[][]     extracted_profiles;         // profiles extracted with Profiler3D

    public PeakExtractor3D(int n0, int n1)
    {
        this.begN = n0; // split oriented filters
        this.endN = n1;
    }

    public static void loadTemplate(Sphere3D sph3_init, int[][] listLocs3D_zxy, short[][] profiles_input, float[][][] img3_zxy_input, float zDist_input) {

        sph3            = sph3_init;
        img3_zxy        = img3_zxy_input;                   // just necessary to rank peaks
        listLocs3D      = listLocs3D_zxy;
        zDist           = zDist_input;
        extracted_profiles = profiles_input;

        // allocate output -> set to -1
        peaks3 = new int[listLocs3D.length][4][3];
        for (int ii = 0; ii<listLocs3D.length; ii++) {
            for (int jj = 0; jj<4; jj++) {
                for (int kk = 0; kk<3; kk++) {
                    peaks3[ii][jj][kk] = -1;
                }
            }
        }

    }

    public void run()
    {

        for (int locationIdx = begN; locationIdx < endN; locationIdx++) {

            int atZ = listLocs3D[locationIdx][0];
            int atX = listLocs3D[locationIdx][1];
            int atY = listLocs3D[locationIdx][2];

            sph3.peakCoords_4xXYZ(extracted_profiles[locationIdx], atX, atY, atZ, img3_zxy, zDist, peaks3[locationIdx]);

        }

    }

    public static int[][] getPeakCoords(int atX, int atY, int atZ, int[][][] indexTableZXY) {

        int locationIndex = indexTableZXY[atZ][atX][atY];

        return  peaks3[locationIndex]; // will return only reference
    }

}