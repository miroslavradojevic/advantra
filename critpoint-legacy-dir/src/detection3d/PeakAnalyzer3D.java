package detection3d;

import ij.IJ;

import java.util.Arrays;

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

    public static int[][][] locIndexZXY;                // correlate xyz with location index: dimX x dimY x dimZ (look up table)

    public static int M = 2;                            // how much it expands recursively from the center

	public static float minCos = 0.6f;                  // allowed derail

    // mainly to associate the peaks and link follow-up points
    public static int[][][] delin3;                     // N(foreground locs.) x 4(max. threads) x M(follow-up locs) contains index for each location (OUTPUT)
//    public static float[][] theta;                    // N(foreground locs.) x 3 (top 3 thetas inputs to fuzzy system) (ALTERNATIVE OUTPUT)

    public PeakAnalyzer3D(int n0, int n1){
        this.begN = n0;
        this.endN = n1;
    }

    public static void loadTemplate(Sphere3D sph3_init, int[][] listLocs3D_zxy, int[][][] listPeaks3D_input, int[][][] locIndexZXY_input, float[][][] img3_zxy_input, float _minCosAng) {

        sph3 = sph3_init;
        listLocs3D = listLocs3D_zxy;
        img3_zxy = img3_zxy_input;
        listPeaks3D = listPeaks3D_input;                // peaks at each location already sorted by strength at lowest index
        locIndexZXY = locIndexZXY_input;                // lookup table
		minCos = _minCosAng;

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

//            System.out.println("checking point (ZXY)-> "+Arrays.toString(listLocs3D[locationIdx])+" location "+locationIdx);

            // access individual peaks at this point
            for (int peakAtLoc = 0; peakAtLoc<4; peakAtLoc++) {  // there are 4 branches allocated

                if (listPeaks3D[locationIdx][peakAtLoc][0] != -1) { // x value stored as peak

                    // store the index value of this peak
                    int atX = listPeaks3D[locationIdx][peakAtLoc][0];
                    int atY = listPeaks3D[locationIdx][peakAtLoc][1];
                    int atZ = listPeaks3D[locationIdx][peakAtLoc][2];

                    int indexValue = locIndexZXY[atZ][atX][atY];

//                    System.out.println("peak found at ZXY "+atZ+","+atX+","+atY);

//                    if (indexValue==-1) System.out.println("FOUND IT! mask says "+locIndexZXY[atZ][atX][atY]);

//                    System.out.println(locationIdx+": peak found at ZXY "+atZ+","+atX+","+atY+" is : "+indexValue+" ::: double check (ZXY): "+ Arrays.toString(listLocs3D[indexValue]));

                    delin3[locationIdx][peakAtLoc][0] = indexValue; // m=0

                    int curr_index, prev_index, next_index;

                    curr_index = indexValue;
                    prev_index = locationIdx;

                    for (int m=1; m<M; m++) {

//                        System.out.println("prev : "+prev_index);
//                        System.out.println("curr : "+curr_index);

                        // recursion : prev+curr->next index
                        next_index = getNext(prev_index, curr_index); // next follow-up will be calculated and sotred in nextXYZ

//                        System.out.println("next : "+next_index);

                        if (next_index!=-1) { // -1 will be if the next one is not found

//                            next_index = locIndexZXY[nextXYZ[0]][nextXYZ[1]][nextXYZ[2]];
                            // store it in output matrix
                            delin3[locationIdx][peakAtLoc][m] = next_index;

                        }
                        else { // follow-up does not exist, break looping m (extending further) but continue looping branches
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

        }

    }

    private static int getNext(int prev_index, int curr_index){

        // these are stacked as ZXY take care on that
        int prevZ = listLocs3D[prev_index][0];    // Z
        int prevX = listLocs3D[prev_index][1];    // X
        int prevY = listLocs3D[prev_index][2];    // Y
        int currZ = listLocs3D[curr_index][0];    // Z
        int currX = listLocs3D[curr_index][1];    // X
        int currY = listLocs3D[curr_index][2];    // Y

        // check peaks at curr
        int[][] pks4xXYZ = listPeaks3D[curr_index];
        for (int pkIdx = 0; pkIdx<4; pkIdx++) { // loops them by rank - those with highest weight first, pick the first one that points outwards

            if (pks4xXYZ[pkIdx][0] != -1) {
                // there is a peak to check - needs to be pointing outwards

                int nextX = pks4xXYZ[pkIdx][0];
                int nextY = pks4xXYZ[pkIdx][1];
                int nextZ = pks4xXYZ[pkIdx][2];

                double cosAng =
                        (
                        (currX-prevX)*(nextX-currX) + (currY-prevY)*(nextY-currY) + (currZ-prevZ)*(nextZ-currZ)
                        )
                                /
                        (
                                Math.sqrt( Math.pow(currX-prevX, 2) +  Math.pow(currY-prevY, 2) + Math.pow(currZ-prevZ, 2) ) *
                                        Math.sqrt( Math.pow(nextX-currX, 2) + Math.pow(nextY-currY, 2) + Math.pow(nextZ-currZ, 2) )
                        );

                if (cosAng>minCos) {

                    // it is aligned - add it, find its index location and output
                    return locIndexZXY[nextZ][nextX][nextY];

                }
                else {
                    // if not pointing outwards, continue further
                }
            }
            else {
                // no more peaks to search
                return -1;
            }

        }

        return -1; // checked all

    }

    public static void summary(int atX, int atY, int atZ) {

        IJ.log("### Delineation (connecting lines) ###");
        int locationIdx = locIndexZXY[atZ][atX][atY];

        int[][] skeleton = delin3[locationIdx];

        for (int i1 = 0; i1<skeleton.length; i1++) {

			System.out.println("skeleton "+i1+""+Arrays.toString(skeleton[i1]));

			for (int i2=0; i2<skeleton[0].length; i2++) {
                if (skeleton[i1][i2] != -1) {

                    int pointId = skeleton[i1][i2];

					int pointZ = listLocs3D[pointId][0];
                    int pointX = listLocs3D[pointId][1];
                    int pointY = listLocs3D[pointId][2];

					if (i2==0) {
						// line ~ median along the line
						IJ.log(locationIdx+" " + (i2+1) + " " + (atX+0.5) +    " " + (atY+0.5) +    " " + (atZ+0.5) +    " " + 0.4f + " -1");
						IJ.log(pointId+    " " + (i2+1) + " " + (pointX+0.5) + " " + (pointY+0.5) + " " + (pointZ+0.5) + " " + 0.4f + " "+locationIdx);

					}
					else {

                        int prevPtIdx = skeleton[i1][i2-1];
                        int prevPtZ   = listLocs3D[prevPtIdx][0];
                        int prevPtX   = listLocs3D[prevPtIdx][1];
                        int prevPtY   = listLocs3D[prevPtIdx][2];

						// line
						IJ.log(prevPtIdx+" " + (i2+1) + " " + (prevPtX+0.5) +    " " + (prevPtY+0.5) +    " " + (prevPtZ+0.5) +    " " + 0.4f + " -1");
						IJ.log(pointId+    " " + (i2+1) + " " + (pointX+0.5) + " " + (pointY+0.5) + " " + (pointZ+0.5) + " " + 0.4f + " "+prevPtIdx);

					}


                }
            }
        }


    }

}
