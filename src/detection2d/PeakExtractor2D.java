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

    public static int[][][]     xy2i;

    public static short[][]	    extracted_profiles;         // profiles

    // outputs are lists of selected peaks
    public static int[][][]     peaks1;                     // N x (4x1)   4 selected peaks in abscissa coordinates X
    public static int[][][]     peaks2;                     // N x (4x2)   main output  4 selected peaks in XY format (OUTPUT)

    public PeakExtractor2D(int n0, int n1)
    {
        this.begN = n0;
        this.endN = n1;
    }





}
