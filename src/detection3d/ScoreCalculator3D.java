package detection3d;

import aux.Stat;
import detection.Interpolator;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ByteProcessor;

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

	public static Fuzzy3D		fz3d;					 // used to calculate the score

    public static float[][]     theta3;                  // OUTPUT, 4 scores to qualify the critical point

	public static float[]		score3;					 // array of scores provided by Fuzzy3D

    public ScoreCalculator3D(int n0, int n1){
        this.begN = n0;
        this.endN = n1;
    }

    public static void loadTemplate(int[][][] _delin_idxs_Nx4x3, float[][][] _img3_zxy, int[][][] _lookup_idx_zxy, int[][] _loc_zxy, byte[][][] _back3_zxy, Fuzzy3D _fz3d) {

        delin_idxs_Nx4x3 = _delin_idxs_Nx4x3;
        img3_zxy = _img3_zxy;
        lookup_idx_zxy = _lookup_idx_zxy;
        loc_zxy = _loc_zxy;
        back3_zxy = _back3_zxy;
		fz3d = _fz3d;

        // initialize output with zeros
        theta3 = new float[delin_idxs_Nx4x3.length][4];
		score3 = new float[delin_idxs_Nx4x3.length];

    }

    public void run() {

        for (int locationIdx = begN; locationIdx < endN; locationIdx++) {

            int locationX = loc_zxy[locationIdx][1];
            int locationY = loc_zxy[locationIdx][2];
            int locationZ = loc_zxy[locationIdx][0];

            // calculate theta3[location] array using delin_indxs_Nx4x3
            float[] _8Nbhood 	= new float[3*3*3];
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

            theta3[locationIdx][0] = Stat.median(_8Nbhood) - (back3_zxy[locationZ][locationX][locationY] & 0xff); // theta0

            int[] idx_container = new int[delin_idxs_Nx4x3[0][0].length+1];

            for (int threadIdx=0; threadIdx<delin_idxs_Nx4x3[locationIdx].length; threadIdx++) { // loop along the thread

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

                        idx_container[compInd+1] = current_point_index;

                    }
                    else {
                        // recursion is not complete
                        completed = false;
                        break;
                    }

                }

                if (completed) {

                    //thetaToAdd = medianAlongLine(idx_container, img3_zxy, loc_zxy); // median of the raw image values
					thetaToAdd = medianDiffAlongLine(idx_container, img3_zxy, back3_zxy, loc_zxy); // median of differences between image values and background values along the line

//                    System.out.println("2 add: "+thetaToAdd);
//                    System.out.println(locationIdx+"before: "+Arrays.toString(theta3[locationIdx]));

                    // theta 1-3 - arrange it in the array - highest first (first value is filled by now)
                    for (int q=0; q<theta3[locationIdx].length; q++) {
                        if (thetaToAdd > theta3[locationIdx][q]) {

							if (q==theta3[locationIdx].length-1) {
								// it is the last one, just replace it, no need to break here
								theta3[locationIdx][q] = thetaToAdd; // add it
							}
							else {
								// it is not the last one
								// shift those below and add it to the list at q
								for (int w=theta3[locationIdx].length-1; w>=q+1; w--) {
									theta3[locationIdx][w] = theta3[locationIdx][w-1];
								}
								theta3[locationIdx][q] = thetaToAdd; // add it
								break;
							}

                        }
                    }
                }
            }
            //System.out.println("(" + locationX + " , " + locationY + " , " + locationZ + ") -> " + Arrays.toString(theta3[locationIdx]));

			// FLS score
			score3[locationIdx] = fz3d.bifurcationess(
															 theta3[locationIdx][0],
															 theta3[locationIdx][1],
															 theta3[locationIdx][2],
															 theta3[locationIdx][3]
															 );

        } // location loop

    }

	public static ImagePlus getThresholdedScores(float th) {
		// turn thresholded scores into an output byte image
		int W = img3_zxy[0].length;
		int H = img3_zxy[0][0].length;
		int L = img3_zxy.length;

		ImageStack out_stack = new ImageStack(W, H);
		for (int l=0; l<L; l++) {
			out_stack.addSlice(new ByteProcessor(W, H));
		}

		for (int l=0; l<score3.length; l++) {

			if (score3[l]>=th) {
				int x = loc_zxy[l][1];
				int y = loc_zxy[l][2];
				int z = loc_zxy[l][0];
				out_stack.setVoxel(x,y,z, 255);
			}

		}

		ImagePlus out_im = new ImagePlus("scores_th="+th, out_stack);
		return out_im;

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

        return Stat.median(values_along_line);

    }

	private static float medianDiffAlongLine(int[] _idx_xyz, float[][][] img3d_zxy, byte[][][] back3d_zxy, int[][] _idx2zxy) { // it is not necessary to have image and background arrays as arguments

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

			// take the background as value calculated around currX, currY, currZ  for this line
			byte background_value = back3d_zxy[currZ][currX][currY];

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
				values_along_line[reference_idx + cc] = Interpolator.interpolateAt(atX, atY, atZ, img3d_zxy) - (background_value & 0xff);
			}

		}

		return Stat.median(values_along_line);

	}

}