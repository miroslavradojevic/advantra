package detection3d;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 12/7/13
 * Time: 8:25 PM
 */
public class PeakAnalyzer3D extends Thread {

    private int begN, endN;

    public static Sphere3D 	sph3;                       // processing unit

    public static float[][][]   img3_zxy;               // link to input image array

    public static int[][] 	listLocs3D;                 // N x 3 (Z,X,Y), list of locations

    public static int[][][] listPeaks3D;                // N x (4x3) main output  4 points in XYZ format

    public static int[][][] locIndexXYZ;                // correlate xyz with location index: dimX x dimY x dimZ

    public static int M = 2;                            // how much it expands recursively from the center

    // mainly to associate the peaks and link follow-up points
    public static int[][][] delin3; // N(foreground locs.) x 4(max. threads) x M(follow-up locs) contains index for each location (OUTPUT)
//    public static float[][] theta;  // N(foreground locs.) x 3 (top 3 thetas inputs to fuzzy system) (ALTERNATIVE OUTPUT)

    public PeakAnalyzer3D(int n0, int n1){
        this.begN = n0;
        this.endN = n1;
    }

    public static void loadTemplate(Sphere3D sph3_init, int[][] listLocs3D_zxy, int[][][] listPeaks3D_input, int[][][] locIndexXYZ_input, float[][][] img3_zxy_input) {

        sph3 = sph3_init;
        listLocs3D = listLocs3D_zxy;
        img3_zxy = img3_zxy_input;
        listPeaks3D = listPeaks3D_input;                // peaks at each location already sorted by strength at lowest index
        locIndexXYZ = locIndexXYZ_input;                // lookup table

        // initialize delin3  output
        delin3 = new int[listLocs3D.length][4][M];      // will contain indexes for every 3d location (xyz can be recovered with look-up table)
        for (int ii = 0; ii<listLocs3D.length; ii++) {
            for (int jj=0; jj<4; jj++) {
                for (int kk=0; kk<M; kk++) {
                    delin3[ii][jj][kk] = -1;
                }
            }
        }

    }

    public void run() {

        for (int locationIdx = begN; locationIdx < endN; locationIdx++) {

            // access individual peaks at this point
            for (int peakAtLoc = 0; peakAtLoc<4; peakAtLoc++) {  // there are 4 peaks allocated
                if (listPeaks3D[locationIdx][peakAtLoc][0] != -1) { // x value stored as peak

                    // store the index value of this peak
                    int atX = listPeaks3D[locationIdx][peakAtLoc][0];
                    int atY = listPeaks3D[locationIdx][peakAtLoc][1];
                    int atZ = listPeaks3D[locationIdx][peakAtLoc][2];

                    int indexValue = locIndexXYZ[atX][atY][atZ];
                    delin3[locationIdx][peakAtLoc][0] = indexValue; // m=0

                    int curr_index, prev_index, next_index;

                    int[] prevXYZ = new int[3]; // TODO could be outside loop
                    int[] currXYZ = new int[3];
                    int[] nextXYZ = new int[3];

                    curr_index = indexValue;
                    prev_index = locationIdx;

                    for (int m=1; m<M; m++) {

                        // recursion : prev+curr->next index

                        // these are stacked as ZXY take care on that
                        prevXYZ[0] = listLocs3D[prev_index][1];    // X
                        prevXYZ[1] = listLocs3D[prev_index][2];    // Y
                        prevXYZ[2] = listLocs3D[prev_index][0];    // Z

                        currXYZ[0] = listLocs3D[curr_index][1];    // X
                        currXYZ[1] = listLocs3D[curr_index][2];    // Y
                        currXYZ[2] = listLocs3D[curr_index][0];    // Z

                        boolean is_found = getNext(prevXYZ, currXYZ, nextXYZ); // next follow-up will be calculated and sotred in nextXYZ

                        if (is_found) {

                            next_index = locIndexXYZ[nextXYZ[0]][nextXYZ[1]][nextXYZ[2]];

                            // store it in output matrix
                            delin3[locationIdx][peakAtLoc][m] = next_index;


                        }
                        else { // follow-up does not exist, break looping m (extending further) but continue looping peaks
                            break;
                        }

                        // turn the indexes back
                        prev_index = curr_index;
                        curr_index = next_index;

                    }


                }
                else { // it was -1 and the rest are not really peaks
                    break; // stop for() loop
                }
            }


            //int atZ = listLocs3D[locationIdx][0];
            //int atX = listLocs3D[locationIdx][1];
            //int atY = listLocs3D[locationIdx][2];

            //sph3.peakCoords_4xXYZ(extracted_profiles[locationIdx], atX, atY, atZ, img3_zxy, zDist, peaks3[locationIdx]);

        }

    }

    private static boolean getNext(int[] prevXYZ, int[] currXYZ, int[] nextXYZ){
        boolean isOK = false;
        // modify nextXYZ[0][1][2]

        // loop options - read peaks from currXYZ and see how they correlate to current and previous point


        return isOK;
    }


}
