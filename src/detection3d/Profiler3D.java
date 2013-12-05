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
	public static short[][]	prof3;                      // profiles could take much of memory (therefore 'short')
    public static int[][][]   peaks3;                   // better solution would be to keep the peaks and do peak extraction immediately on extracted peaks
    // peaks3  ->  N x 4 x 3

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

//		/*
//			set prof3
//		 */
//		System.out.println("allocated for profiles: "+listLocs3D.length+" locations x "+sph3.getProfileLength()+" (profile length)");
		prof3 = new short[listLocs3D.length][sph3.getProfileLength()];
		System.out.println("input "+prof3.length+" locations x "+prof3[0].length+" profiles");

        peaks3 = new int[listLocs3D.length][4][3]; // 4 points if they exist
        for (int ii = 0; ii<listLocs3D.length; ii++) {
            for (int jj = 0; jj<4; jj++) {
                for (int kk = 0; kk<3; kk++) {
                    peaks3[ii][jj][kk] = -1;
                }
            }
        }

        System.out.println("initialized peaks and profiles!");

	}

	public void run()
	{


//        if (true) {

        // profile components are split in parallel
        // prof3 is filled up in parallel

            for (int profileComponentIdx = begN; profileComponentIdx < endN; profileComponentIdx++) {

                for (int locIdx = 0; locIdx < listLocs3D.length; locIdx++) {

                    int atZ = listLocs3D[locIdx][0];
                    int atX = listLocs3D[locIdx][1];
                    int atY = listLocs3D[locIdx][2];

                    prof3[locIdx][profileComponentIdx] = sph3.extractProfile(profileComponentIdx, atX, atY, atZ, img3_zxy, zDist);

                }

            }

//        }


//        if (false) {
//
//        // extract peak locations in parallel for subset of locations - locations are split in parallel
//        // peaks3 is filled up
//        int profileSize = sph3.getProfileLength();
//        float[] profile = new float[profileSize];                   // allocate array where circular neighbourhood values will be stored
//
//        for (int locIdx=begN; locIdx<endN; locIdx++) {              // work on the range of operations
//
//            int atZ = listLocs3D[locIdx][0];
//            int atX = listLocs3D[locIdx][1];
//            int atY = listLocs3D[locIdx][2];
//
//            sph3.extractProfile(atX, atY, atZ, img3_zxy, zDist, profile); // output will be stored in profile array
//
//            sph3.peakCoords_4xXYZ(
//                    profile,
//                    atX,
//                    atY,
//                    atZ,
//                    img3_zxy,
//                    zDist,
//                    peaks3[locIdx]
//            );// detect peaks and store them in corresponding location storage
//
//        }
//
//        }

    }

	public static float[] sequence(float beginValue, float endValue, int Nelements, boolean lastIncluded)
	{

		float[] seqArray = new float[Nelements];

		if (lastIncluded) {
			float step = (endValue-beginValue)/(Nelements-1);
			for (int ii=0; ii<Nelements; ii++) seqArray[ii] = beginValue + ii * step;
		}
		else {
			float step = (endValue-beginValue)/Nelements;
			for (int ii=0; ii<Nelements; ii++) seqArray[ii] = beginValue + ii * step;
		}

		return seqArray;

	}

}
