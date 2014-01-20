package detection2d;

import detection3d.Sphere3D;
import ij.gui.Plot;
import ij.process.ImageProcessor;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 12/17/13
 * Time: 9:44 AM
 * Loops selected list of foreground locations and extracts profiles from those locations
 */
public class Profiler2D extends Thread {

	private int begN, endN;

	public static float[][] inimg_xy;
    public static Sphere2D  sph2d;
    public static int[][] 	i2xy;                       // selected locations

    // output
    public static short[][]	prof2;

    public Profiler2D(int n0, int n1)
    {
        this.begN = n0;
        this.endN = n1;
    }

    public static void loadTemplate(Sphere2D _sph2d, int[][] _i2xy, float[][] _inimg_xy){

        sph2d = _sph2d; // just assign link, no allocation necessary since there will be one Sphere3D instance for all
        i2xy = _i2xy;
        inimg_xy = _inimg_xy;

        // allocate output
        prof2 = new short[i2xy.length][sph2d.getProfileLength()];

    }

    public void run() {

        for (int profileComponentIdx = begN; profileComponentIdx < endN; profileComponentIdx++) {

            for (int locIdx = 0; locIdx < i2xy.length; locIdx++) {

                int atX = i2xy[locIdx][0];
                int atY = i2xy[locIdx][1];

                prof2[locIdx][profileComponentIdx] = sph2d.extractProfile(profileComponentIdx, atX, atY, inimg_xy);

            }

        }
    }

    /*
        various versions of the output (plot imageprocessor, just an array with corresponding thetas in degs)
     */

    public static ImageProcessor getProfile(int atX, int atY, int[][] _xy2i){  // reads from prof2 array class member

        int idx = _xy2i[atX][atY];
        if (idx != -1) {
            int len = prof2[0].length;
            float[] f = new float[len];
            float[] fx = new float[len];

            for (int i=0; i<len; i++) {
                f[i] = ((prof2[idx][i] & 0xffff) / 65535f) * 255f; // retrieve the profile
                //fx[i] = i; // abscissa in cnt
                fx[i] = (i / (float) len) * 360; // abscissa in degs

            }

            Plot p = new Plot("profile at ("+atX+","+atY+")", "", "filtered", fx, f);
            p.setSize(600, 300);

            return p.getProcessor();
        }
        else {
            return null;
        }

    }

}