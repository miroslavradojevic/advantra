package detection3d;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 11/13/13
 * Time: 11:40 AM
 */
public class Profiler3D extends Thread {

	private int begN, endN;

    public static float[][][] img3_zxy;                 // link to input image array
	public static float zDist;                          // zDist (goes together with 3d image)

    public static Sphere3D 	sph3;                       // processing unit
	public static int[][] 	listLocs3D;                 // N x 3 (Z,X,Y)
	public static short[][]	prof3;                      // profiles could take much of memory (therefore 'short') (OUTPUT)

	public Profiler3D(int n0, int n1)
	{
		this.begN = n0; // split oriented filters
		this.endN = n1;
	}

	public static void loadTemplate(Sphere3D sph3_init, int[][] listLocs3D_zxy, float[][][] img3_zxy_input, float zDist_input){

		sph3 = sph3_init; // just assign link, no allocation necessary since there will be one Sphere3D instance for all
		listLocs3D = listLocs3D_zxy;
        img3_zxy = img3_zxy_input;
        zDist = zDist_input;

        // allocate output
		prof3 = new short[listLocs3D.length][sph3.getProfileLength()];

	}

	public void run() {

        for (int profileComponentIdx = begN; profileComponentIdx < endN; profileComponentIdx++) {

            for (int locIdx = 0; locIdx < listLocs3D.length; locIdx++) {

                int atZ = listLocs3D[locIdx][0];
                int atX = listLocs3D[locIdx][1];
                int atY = listLocs3D[locIdx][2];

                prof3[locIdx][profileComponentIdx] = sph3.extractProfile(profileComponentIdx, atX, atY, atZ, img3_zxy, zDist);

            }

        }
    }

    public static short[][] getProfiles() {
        return  prof3;
    }

}
