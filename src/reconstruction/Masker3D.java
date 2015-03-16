package reconstruction;

import aux.Stat;
import ij.IJ;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by miroslav on 1-3-15.
 */
public class Masker3D extends Thread {

    private int begN, endN;

    public static float             zDist;

    public static int 				image_width;
    public static int 				image_height;
    public static int               image_length;

    public static float[][][]		inimg_xyz;

    private static float            radiusCheck;
    private static float 			globalTh;
    private static int              marginPix;
    private static int              marginLay;

    private static float            percentile;

    // outputs
    public static float[]  		            criteria;	// scores for each voxels are stored
    public static ArrayList<boolean[][]>    mask_xyz;   // output mask
    public static int[][] 			        i2xyz;      // mapping
    public static ArrayList<int[][]>	    xyz2i;      // output tag map

    public Masker3D (int n0, int n1) {
        begN = n0;
        endN = n1;
    }

    public static void loadTemplate(
            float[][][] _inimg_xyz,
            int         _margin_xy,
            int         _margin_z,
            float       _radiusNbhoodCheck,
            float       _percentile,
            float       _zDist
    )
    {

        zDist = _zDist;

        radiusCheck = _radiusNbhoodCheck;
        percentile  = (_percentile<0)? 0 : (_percentile>100)? 100 : _percentile;
        percentile = percentile / 5;

//        marginPix = (int) Math.ceil(radiusCheck);
        marginPix = _margin_xy;//(_margin>marginPix)? _margin : marginPix ;
        marginLay = _margin_z;//(int) Math.ceil(radiusCheck/zDist);

        inimg_xyz  = _inimg_xyz;
        image_height = inimg_xyz[0].length;
        image_width = inimg_xyz.length;
        image_length = inimg_xyz[0][0].length;

        // initialize outputs
        mask_xyz = new ArrayList<boolean[][]>(image_length);
        for (int i = 0; i < image_length; i++) mask_xyz.add(new boolean[image_width][image_height]);
        criteria = new float[image_width*image_height*image_length];
        i2xyz = null;
        xyz2i = new ArrayList<int[][]>(image_length);

    }

    public static void clean() {

//        image_width = -1;
//        image_height = -1;
//        image_length = -1;
//        inimg_xyz = null;
//        radiusCheck = -1;
//        globalTh = -1;
//        marginPix = -1;
//        percentile = -1;
//        zDist = -1;

        criteria = null;

        for (int i = 0; i < mask_xyz.size(); i++) mask_xyz.set(i, null);
        mask_xyz.clear();

        i2xyz = null;

        for (int i = 0; i < xyz2i.size(); i++) xyz2i.set(i, null);
        xyz2i.clear();

    }

    public void run()
    {

        int circNeighSize = sizeCircularNbhood(radiusCheck, zDist);
        float[] circNeigh = new float[circNeighSize];

        for (int locIdx=begN; locIdx<endN; locIdx++) {

            // this is where indexed location in 3d is expressed in x,y,z
            int atX = locIdx%image_width;
            int atZ = locIdx/(image_width*image_height);
            int atY = locIdx/image_width - atZ * image_height;

            boolean processIt =
                    (atX>=marginPix) && (atY>=marginPix) && (atZ>=marginLay) &&
                            (atX<image_width-marginPix-1) && (atY<image_height-marginPix-1) && (atZ<image_length-marginLay-1);

            if (processIt) {

                extractCircularNbhood(atX, atY, atZ, radiusCheck, circNeigh, inimg_xyz, zDist); // extract from the image
                float m05 	= Stat.quantile(circNeigh, 1, 20);
                float m95 	= Stat.quantile(circNeigh, 19, 20);
                criteria[locIdx] = m95 - m05; // atY*image_width+atX

            }

        }

    }

    private static int sizeCircularNbhood(float sphereRadius, float zDist)
    {

        int rPix 		= Math.round(sphereRadius);
        int rLayers 	= (int) Math.floor(sphereRadius / zDist);

        int cnt = 0;

        for (int a=-rPix; a<=rPix; a++) {
            for (int b=-rPix; b<=rPix; b++) {
                for (int cLayers=-rLayers; cLayers<=rLayers; cLayers++) {
                    // loop in layers
                    float c = cLayers * zDist;   // back to pixels
                    if ( a*a+b*b+c*c <= rPix*rPix ) cnt++;
                }
            }
        }

        return cnt;
    }

    private static void extractCircularNbhood(
            int atX,
            int atY,
            int atZ,
            float sphereRadius,
            float[] values, // allocated in advance
            float[][][] _inimg_xyz,
            float zDist
    )
    {

        int rPix = Math.round(sphereRadius);
        int rLay = (int) Math.floor(sphereRadius / zDist);

        if (atX-rPix>=0 && atY-rPix>=0 && atX+rPix<_inimg_xyz.length && atY+rPix<_inimg_xyz[0].length && atZ-rLay>=0 && atZ+rLay < _inimg_xyz[0][0].length) {

            int cnt = 0;

            for (int xLoc=atX-rPix; xLoc<=atX+rPix; xLoc++) {
                for (int yLoc=atY-rPix; yLoc<=atY+rPix; yLoc++) {
                    for (int zLoc=atZ-rLay; zLoc<=atZ+rLay; zLoc++) { // loop in layers
                    float c = (zLoc-atZ) * zDist;   // back to pixels
                        if ( (xLoc-atX)*(xLoc-atX)+(yLoc-atY)*(yLoc-atY)+c*c <= rPix*rPix ) {
                            values[cnt] = _inimg_xyz[xLoc][yLoc][zLoc];
                            cnt++;
                        }
                    }
                }
            }
        }
        else {
            Arrays.fill(values, 0f); // set all values to zero
        }

    }

    public static void defineThreshold()
    {

        // count the number of those that have discrepancy criteria above zero, they will be thresholded using percentile threshold
        int count = 0;
        for (int i = 0; i < criteria.length; i++) {
            if (criteria[i]>Float.MIN_VALUE)
                count++;
        }

        float[] criteria_reduced = new float[count]; // useful if there are many zeros in the background

        count = 0;
        for (int i = 0; i < criteria.length; i++) {
            if (criteria[i]>Float.MIN_VALUE) {
                criteria_reduced[count] = criteria[i];
                count++;
            }
        }

        // exclude those that were zero criteria (reduced) when calculating percentile threshold
        globalTh = Stat.quantile(criteria_reduced, (int) percentile, 20);
        criteria_reduced = null; // to clear the memory faster

        for (int i = 0; i < criteria.length; i++) {
            if (criteria[i]>globalTh) {
                mask_xyz.get(i/(image_width*image_height))[i%image_width][i/image_width-(i/(image_width*image_height))*image_height] = true;
            }
            else {
                mask_xyz.get(i/(image_width*image_height))[i%image_width][i/image_width-(i/(image_width*image_height))*image_height] = false;
            }
        }

//        float[][] criteria_local_max = new float[criteria.length][criteria[0].length];
//
//        int circNeighSize = sizeCircularNbhood(radiusCheck);
//        float[] circNeigh = new float[circNeighSize];
//
//        // this is effectively dilatation
//        for (int x=0; x<criteria.length; x++) {
//            for (int y=0; y<criteria[0].length; y++) {
//                extractCircularNbhood(x, y, radiusCheck, circNeigh, criteria);
//                criteria_local_max[x][y] = Stat.get_max(circNeigh);
//            }
//        }
//
//        for (int ii=0; ii<criteria.length; ii++) {
//            for (int jj=0; jj<criteria[0].length; jj++) {
//                criteria[ii][jj] = criteria_local_max[ii][jj];
//            }
//        }
//        // this is effectively dilatation - just to change criteria
//
//        float[] criteria_temp = new float[criteria.length*criteria[0].length];
//        int		cnt = 0;
//        for (int ii=0; ii<criteria.length; ii++) {
//            for (int jj=0; jj<criteria[0].length; jj++) {
//                criteria_temp[cnt] = criteria_local_max[ii][jj];
//                cnt++;
//            }
//        }
//
//        cnt = 0;
//        for (int xx=0; xx<image_width; xx++) {
//            for (int yy=0; yy<image_height; yy++) {

//                if (criteria[xx][yy] > globalTh) {
//                    mask_xy[xx][yy] = true;
//                }
//                else {
//                    mask_xy[xx][yy] = false;
//                }

//                criteria[xx][yy] = criteria_temp[cnt];

//                cnt++;

//				if (criteria[xx][yy] > globalTh) {
//					mask_xy[xx][yy] = true;
//				}
//				else {
//					mask_xy[xx][yy] = false;
//				}
//            }
//        }

    }

    public static void formRemainingOutputs() // will fill in the remainder of the mapping arrays
    {

        xyz2i = new ArrayList<int[][]>(image_length);

        int cnt = 0;

        for (int i = 0; i < image_length; i++) {

            xyz2i.add(new int[image_width][image_height]);

            for (int xx=0; xx<image_width; xx++) {
                for (int yy=0; yy<image_height; yy++) {

                    if (mask_xyz.get(i)[xx][yy]) {
                        xyz2i.get(i)[xx][yy] = cnt;
                        cnt++;
                    }
                    else {
                        xyz2i.get(i)[xx][yy] = -1;
                    }

                }
            }

        }
        
		float perc = (cnt*100f) / (image_width*image_height*image_length);
        System.out.println(String.format("%3.1f %% vol. foreground extracted", perc));

        // component: foreground point list: locations and background estimates
        i2xyz = new int[cnt][3];

        cnt =0;

        for (int i = 0; i < image_length; i++) {
            for (int xx=0; xx<image_width; xx++) {
                for (int yy=0; yy<image_height; yy++) {

                    if (mask_xyz.get(i)[xx][yy]) {
                        i2xyz[cnt][0] = xx;
                        i2xyz[cnt][1] = yy;
                        i2xyz[cnt][2] = i;
                        cnt++;
                    }

                }
            }
        }

    }

}
