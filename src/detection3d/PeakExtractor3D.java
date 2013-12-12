package detection3d;

import ij.IJ;
import ij.ImagePlus;

import java.util.ArrayList;

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

    public static int[][][]     lookupIdxZXY;               // dimZ x dimX x dimY

    public static float         zDist;                      // to properly calibrate layer

    public static short[][]     extracted_profiles;         // profiles extracted with Profiler3D

    // outputs are lists of selected peaks
    public static int[][][]     peaks2;                     // N x (4x2)   4 selected peaks in planisphere image coordinates XY
    public static int[][][]     peaks3;                     // N x (4x3)    main output  4 selected peaks in XYZ format (OUTPUT)

    public PeakExtractor3D(int n0, int n1)
    {
        this.begN = n0;
        this.endN = n1;
    }

    public static void loadTemplate(Sphere3D sph3_init, int[][] listLocs3D_zxy, short[][] profiles_input, float[][][] img3_zxy_input, float zDist_input, int[][][] lookupIdxZXY_input) {

        sph3            = sph3_init;
        img3_zxy        = img3_zxy_input;                   // just necessary to rank peaks
        listLocs3D      = listLocs3D_zxy;
        zDist           = zDist_input;
        extracted_profiles = profiles_input;
        lookupIdxZXY = lookupIdxZXY_input;

        // allocate output -> set to -1
        peaks3  = new int[listLocs3D.length][4][3];
        peaks2  = new int[listLocs3D.length][4][2];
        for (int ii = 0; ii<listLocs3D.length; ii++) {
            for (int jj = 0; jj<4; jj++) {

                for (int kk = 0; kk<3; kk++) {
                    peaks3[ii][jj][kk] = -1;
                }

                for (int kk=0; kk<2; kk++) {
                    peaks2[ii][jj][kk] = -1;
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

            sph3.peakCoords_4xXYZ(extracted_profiles[locationIdx], atX, atY, atZ, img3_zxy, zDist, lookupIdxZXY, peaks3[locationIdx], peaks2[locationIdx]); // 4x3 and 4x2 fill up

        }

    }

    public static ArrayList<int[]> getPeaksXYZ(int atX, int atY, int atZ, int[][][] indexTableZXY) {

        int locationIndex = indexTableZXY[atZ][atX][atY];

		ArrayList<int[]> out = new ArrayList<int[]>(4);

		if (locationIndex != -1) {
		for (int cc = 0; cc < 4; cc++) {
			if (peaks3[locationIndex][cc][0]!=-1) { // value was filled
				out.add(peaks3[locationIndex][cc]);
			}

		}
		}

        return  out; // will return only reference
    }

    public static void summary(int atX, int atY, int atZ) {

        IJ.log("### Peak cloud 4x4x4 etc... ###");
        int locationIdx = lookupIdxZXY[atZ][atX][atY];

        // center
        IJ.log("0 " + 1 + " " + (atX+0.5) + " " + (atY+0.5) + " " + (atZ+0.5) + " " + 0.5f + " -1");

		new ImagePlus("profile0", sph3.drawProfile(extracted_profiles[locationIdx])).show();

        // 1st generation 4 peaks from center
        for (int loop1=0; loop1<4; loop1++) {

            int g1_x = peaks3[locationIdx][loop1][0];
            int g1_y = peaks3[locationIdx][loop1][1];
            int g1_z = peaks3[locationIdx][loop1][2];

            if (g1_x != -1) {

                int spotIndex1 = lookupIdxZXY[g1_z][g1_x][g1_y];

                if (spotIndex1 != -1 && g1_x != -1) {

                    // spot exists
                    IJ.log("1 " + 1 + " " + (g1_x+0.5) + " " + (g1_y+0.5) + " " + (g1_z+0.5) + " " + 0.3f + " -1");
					new ImagePlus("profile"+loop1, sph3.drawProfile(extracted_profiles[spotIndex1])).show();
                    // loop it's peaks if it exists

                    for (int loop2=0; loop2<4;loop2++) {

                        int g2_x = peaks3[spotIndex1][loop2][0];
                        int g2_y = peaks3[spotIndex1][loop2][1];
                        int g2_z = peaks3[spotIndex1][loop2][2];

                        if (g2_x != -1) {


                            int spotIndex2 = lookupIdxZXY[g2_z][g2_x][g2_y];

                            if (spotIndex2 != -1 && g2_x != -1) {

                                // spot exists
                                IJ.log("2 " + 1 + " " + (g2_x+0.5) + " " + (g2_y+0.5) + " " + (g2_z+0.5) + " " + 0.3f + " -1");

                            }


                        }



                    }

                }

            }

        }


    }

    }