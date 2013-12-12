package detection3d;

import detection.Interpolator;

import java.util.Arrays;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 12/12/13
 * Time: 2:30 PM
 */
public class ScoreCalculator3D extends Thread {

    private int begN, endN;

    public static float[][][]   img3_zxy;                // link to input image array

    public static byte[][][]    back3_zxy;               // background values (Masker3D output)

    public static int[][][]     lookup_idx_zxy;          // correlate xyz with location index: dimX x dimY x dimZ (look up table)

    public static int[][]       loc_zxy;                 // zxy coordinates for each location index

    public static int[][][]     delin_idxs_Nx4x3;        // list of delineated peaks around central point in at least 4 directions  N x (4 x 3)

    public static float[][]     theta3;                  // OUTPUT, 4 scores to qualify the critical point

    public ScoreCalculator3D(int n0, int n1){
        this.begN = n0;
        this.endN = n1;
    }

    public static void loadTemplate(int[][][] _delin_idxs_Nx4x3, float[][][] _img3_zxy, int[][][] _lookup_idx_zxy, int[][] _loc_zxy, byte[][][] _back3_zxy) {

        delin_idxs_Nx4x3 = _delin_idxs_Nx4x3;
        img3_zxy = _img3_zxy;
        lookup_idx_zxy = _lookup_idx_zxy;
        loc_zxy = _loc_zxy;
        back3_zxy = _back3_zxy;

        // initialize output with zeros
        theta3 = new float[delin_idxs_Nx4x3.length][4];

    }

    public void run() {

        for (int locationIdx = begN; locationIdx < endN; locationIdx++) {

            int locationX = loc_zxy[locationIdx][1];
            int locationY = loc_zxy[locationIdx][2];
            int locationZ = loc_zxy[locationIdx][0];


            // calculate theta3[location] array using delin_indxs_Nx4x3
            float[] _8Nbhood = new float[3*3*3];
            int cnt = 0;
            for (int dx=-1; dx<=1; dx++) {
                for (int dy=-1; dy<=1; dy++) {
                    for (int dz=-1; dz<=1; dz++) {

                        int x = locationX+dx;
                        int y = locationY+dy;
                        int z = locationZ+dz;

                        if (x>=0 && x<img3_zxy[0].length && y>=0 && y<img3_zxy[0][0].length && z>=0 && z<img3_zxy.length) {

                            _8Nbhood[cnt] = img3_zxy[z][x][y];

                        }

                        cnt++;

                    }
                }
            }
            theta3[locationIdx][0] = Stat.median(_8Nbhood); // theta0




            int[] idx_container = new int[delin_idxs_Nx4x3[0][0].length+1];

            for (int threadIdx=0; threadIdx<delin_idxs_Nx4x3[locationIdx].length; threadIdx++) {
                // loop along the thread

                // initialize theta for this thread
                float thetaToAdd = 0;
                // reset index container
                Arrays.fill(idx_container, -1);
                boolean completed = true;

                // first component
                idx_container[0] = locationIdx; // start line location

                for (int compInd=0; compInd<delin_idxs_Nx4x3[locationIdx][threadIdx].length; compInd++) {

                    int current_point_index = delin_idxs_Nx4x3[locationIdx][threadIdx][compInd];

                    if ( current_point_index != -1 ) {

                        idx_container[compInd+1] = delin_idxs_Nx4x3[locationIdx][threadIdx][compInd];

                    }
                    else {
                        // recursion is not complete
                        completed = false;
                        break;
                    }

                }

                if (completed) {

                    thetaToAdd = medianAlongLine(idx_container, img3_zxy, loc_zxy);

                    System.out.println("2 add: "+thetaToAdd);

                    System.out.println(locationIdx+"before: "+Arrays.toString(theta3[locationIdx]));

                    // theta 1-3 - arrange it in the array - highest first (first value is filled by now)
                    for (int q=0; q<theta3[locationIdx].length; q++) {
                        if (thetaToAdd > theta3[locationIdx][q]) {

                            // shift those below and add it to the list at q
                            for (int w=theta3[locationIdx].length-1; w>=q+1; w--) {
                                theta3[locationIdx][w] = theta3[locationIdx][w-1];
                            }

                            // add it
                            theta3[locationIdx][q] = thetaToAdd;


                        }
                    }

                    System.out.println(locationIdx+"after: "+Arrays.toString(theta3[locationIdx]));

                }
                else {
                    thetaToAdd = 0;
                }


            }

            System.out.println("(" + locationX + " , " + locationY + " , " + locationZ + ") -> " + Arrays.toString(theta3[locationIdx]));


        } // location loop

    }

    private static float medianAlongLine(int[] _idx_xyz, float[][][] img3d_zxy, int[][] _idx2zxy) { // (indexes_of_points_along_the_line, input_image, lookup_table)

        int nr_sub_lines = _idx_xyz.length - 1;

        int elements_per_line = 10;

        float[] values_along_line = new float[nr_sub_lines*elements_per_line];

        for (int sub_line_idx=1; sub_line_idx<=nr_sub_lines; sub_line_idx++) {

            int prevX = _idx2zxy[ _idx_xyz[sub_line_idx-1] ][1];
            int prevY = _idx2zxy[ _idx_xyz[sub_line_idx-1] ][2];
            int prevZ = _idx2zxy[ _idx_xyz[sub_line_idx-1] ][0];

            int currX = _idx2zxy[ _idx_xyz[sub_line_idx] ][1];
            int currY = _idx2zxy[ _idx_xyz[sub_line_idx] ][2];
            int currZ = _idx2zxy[ _idx_xyz[sub_line_idx] ][0];

            float v  = (float) Math.sqrt( Math.pow(currX-prevX, 2) + Math.pow(currY-prevY, 2) + Math.pow(currZ-prevZ, 2) );
            float dx = (currX-prevX)/v; // unit length
            dx *= ((float) (currX-prevX)/elements_per_line);
            float dy = (currY-prevY)/v;
            dy *= ((float) (currY-prevY)/elements_per_line);
            float dz = (currZ-prevZ)/v;
            dz *= ((float) (currZ-prevZ)/elements_per_line);

            int reference_idx = (sub_line_idx-1) * elements_per_line;

            // loop along the line
            for (int cc=0; cc<elements_per_line; cc++) {
                float atX = prevX   + cc * dx;
                float atY = prevY   + cc * dy;
                float atZ = prevZ   + cc * dz;
                values_along_line[reference_idx + cc] = Interpolator.interpolateAt(atX, atY, atZ, img3d_zxy);
            }

        }

//        for (int lineNr = 0; lineNr < nr_lines; lineNr++) {
//            // define borders for lineNr = 0
//            int x2 = x[lineNr+1];
//            int x1 = x[lineNr];
//
//            int y2 = y[lineNr+1];
//            int y1 = y[lineNr];
//
//            int z2lay = z_lay[lineNr+1];
//            int z1lay = z_lay[lineNr];
//
//            float v = (float) Math.sqrt(Math.pow(x2-x1, 2) + Math.pow(y2-y1, 2) + Math.pow(z2lay-z1lay, 2));
//            float dx = (x2 - x1) / v;
//            float dy = (y2 - y1) / v;
//            float dz = (z2lay - z1lay) / v;
//
//            int reference_index = lineNr * elementsInLine;
//
//            for (int cc = 0; cc < elementsInLine; cc++) {
//
//                float atX = x1      + cc * dx;
//                float atY = y1      + cc * dy;
//                float atZ = z1lay   + cc * dz;
//
//                valuesAlongLine[cc+reference_index] = Interpolator.interpolateAt(atX, atY, atZ, img3d_zxy);
//
//            }
//
//        }

        return Stat.median(values_along_line);

    }



}
