import ij.IJ;
import ij.ImagePlus;
import ij.process.FloatProcessor;

import java.util.Arrays;

/**
 * Created by miroslav on 13-2-15.
 */
public class Masker2D extends Thread {

    private int begN, endN;

    public static int 				image_width;
    public static int 				image_height;

    public static float[][]			inimg_xy;

    private static float            radiusCheck;
    private static float 			globalTh;
    private static int              marginPix;
    private static float            percentile;

    // outputs
    public static float[]  		    criteria;	// where scores for each pixel are located
    public static boolean[][]		mask_xy;    // output mask
    public static int[][] 			i2xy;       // mapping
    public static int[][]			xy2i;

    public Masker2D (int n0, int n1) {
        begN = n0;
        endN = n1;
    }

    public static void loadTemplate(
            float[][]   _inimg_xy,
            int         _margin,
            float       _radiusNbhoodCheck,
            float       _percentile
    )
    {

        radiusCheck     = _radiusNbhoodCheck;
        percentile      = (_percentile<0)? 0 : (_percentile>100)? 100 : _percentile;
        percentile      = percentile / 5;

        marginPix       = (int) Math.ceil(radiusCheck);
        marginPix       = (_margin>marginPix)? _margin : marginPix ;  // narrow down selection in XY plane

        inimg_xy        = _inimg_xy;
        image_height 	= inimg_xy[0].length;
        image_width 	= inimg_xy.length;

		// initialize outputs
        mask_xy = new boolean[image_width][image_height];
        criteria = new float[image_width*image_height];
        i2xy 	= null;//new int[1][1];  // values known after run()
        xy2i 	= null;//new int[1][1];  // values known after run()

    }

    public static void clean(){
        inimg_xy = null;
        criteria = null;
        mask_xy = null;
        i2xy = null;
        xy2i = null;
    }

    public void run()
    {

        int circNeighSize = sizeCircularNbhood(radiusCheck);
        float[] circNeigh = new float[circNeighSize];

        for (int locIdx=begN; locIdx<endN; locIdx++) {

            int atX = locIdx % image_width;
            int atY = locIdx / image_width;

            boolean processIt =
                    (atX>=marginPix) && (atY>=marginPix) && //(atZ>=marginLay) &&
                            (atX<image_width-marginPix-1) && (atY<image_height-marginPix-1);// && (atZ<image_length-marginLay-1);

            if (processIt) {

                extractCircularNbhood(atX, atY, radiusCheck, circNeigh, inimg_xy); // extract from the image
                float m05 	= Toolbox.quantile(circNeigh, 1, 20);
                float m95 	= Toolbox.quantile(circNeigh, 19, 20);
                criteria[locIdx] = m95 - m05; // atY*image_width+atX

            }

        }

    }

    private static int sizeCircularNbhood(float sphereRadius)
    {

        //int rLayers 	= Math.round(sphereRadius / zDist);
        int rPix 		= Math.round(sphereRadius);

        int cnt = 0;

        for (int a=-rPix; a<=rPix; a++) {
            for (int b=-rPix; b<=rPix; b++) {
                //for (int cLayers=-rLayers; cLayers<=rLayers; cLayers++) { // loop in layers
                //float c = cLayers * zDist;   // back to pixels
                if ( a*a+b*b <= rPix*rPix ) cnt++; // +c*c
                //}
            }
        }

        return cnt;
    }

    private static void extractCircularNbhood(
            int atX,
            int atY,
            float sphereRadius,
            float[] values,
            float[][] _inimg_xy
    )
    {
        int rPix = Math.round(sphereRadius);
        //int rLay = Math.round(sphereRadius / zDist);
        if (atX-rPix>=0 && atY-rPix>=0 && atX+rPix<_inimg_xy.length && atY+rPix<_inimg_xy[0].length) {  // && atZ-rLay>=0  // && atZ+rLay < image_length
            int cnt = 0;
            for (int xLoc=atX-rPix; xLoc<=atX+rPix; xLoc++) {
                for (int yLoc=atY-rPix; yLoc<=atY+rPix; yLoc++) {
                    //for (int zLoc=atZ-rLay; zLoc<=atZ+rLay; zLoc++) { // loop in layers
                    //float c = (zLoc-atZ) * zDist;   // back to pixels
                    if ( (xLoc-atX)*(xLoc-atX)+(yLoc-atY)*(yLoc-atY) <= rPix*rPix ) { // +c*c

                        values[cnt] = _inimg_xy[xLoc][yLoc];// instack3[atZ][xLoc][yLoc];
                        cnt++;

                    }
                    //}
                }
            }
        }
        else Arrays.fill(values, 0f); // set all values to zero
    }

    private static void extractCircularNbhoodLocs(
            int atX,
            int atY,
            float sphereRadius,
            int[][] locs_xy,
            float[][] _inimg_xy
    )
    {
        int rPix = Math.round(sphereRadius);
        //int rLay = Math.round(sphereRadius / zDist);
        if (atX-rPix>=0 && atY-rPix>=0 && atX+rPix<_inimg_xy.length && atY+rPix<_inimg_xy[0].length) {  // && atZ-rLay>=0  // && atZ+rLay < image_length
            int cnt = 0;
            for (int xLoc=atX-rPix; xLoc<=atX+rPix; xLoc++) {
                for (int yLoc=atY-rPix; yLoc<=atY+rPix; yLoc++) {
                    //for (int zLoc=atZ-rLay; zLoc<=atZ+rLay; zLoc++) { // loop in layers
                    //float c = (zLoc-atZ) * zDist;   // back to pixels
                    if ( (xLoc-atX)*(xLoc-atX)+(yLoc-atY)*(yLoc-atY) <= rPix*rPix ) { // +c*c

                        locs_xy[cnt][0] = xLoc;
                        locs_xy[cnt][1] = yLoc;
                        cnt++;

                    }
                    //}
                }
            }
        }
        else locs_xy=null; //Arrays.fill(values, 0f); // set all values to zero
    }

    public static void defineThreshold()
    {

        float MIN_CRITERIA = 0.1f; // since we are working with 8 bit

        // count the number of those that have discrepancy criteria above MIN_CRITERIA
        int count = 0;
        for (int i = 0; i < criteria.length; i++) {
            if (criteria[i]>MIN_CRITERIA)
                count++;
        }


        float[] criteria_red = new float[count]; // useful if there are many zeros in the background

        count = 0;
        for (int i = 0; i < criteria.length; i++) {
            if (criteria[i]>MIN_CRITERIA) {
                criteria_red[count] = criteria[i];
                count++;
            }
        }

        // exclude those that were zero criteria (reduced) when calculating percentile threshold
        globalTh = Toolbox.quantile(criteria_red, (int) percentile, 20);

        criteria_red = null;

        // now that it was used - set it to null
//        criteria_red = null;

        // form the mask
        if (false) {
            // version 1: percentile threshold of the original score map
            for (int i = 0; i < criteria.length; i++) {
                if (criteria[i]>globalTh) {
                    mask_xy[i%image_width][i/image_width] = true;
                }
                else {
                    mask_xy[i%image_width][i/image_width] = false;
                }
            }
        }
        else {
            // version 2: percentile threshold of the dilated score map (local max expanded in circular neighbourhood)
            // threshold is applied on radius dilated score map
            // it is the same as the criteria is expanded to the whole diameter of the neighbourhood that was used for checking
            // if the score calculated at the center is foreground then
            // all the points that were in that neighbourhood
            // are also foreground since all of those points were taken to calculate that score


            float[][] criteria_xy = new float[image_width][image_height]; // temp (necessary for local neighbourhood)
            for (int i = 0; i < criteria.length; i++) {
                criteria_xy[i%image_width][i/image_width] = criteria[i];
            }



            int circNeighSize = sizeCircularNbhood(radiusCheck);
            float[] circNeighVals = new float[circNeighSize];




            // this is effectively dilatation using real values
            float[] criteria_lmax = new float[criteria.length]; // temp
            for (int i=0; i<criteria.length; i++) {

                int atX = i % image_width;
                int atY = i / image_width;

                extractCircularNbhood(atX, atY, radiusCheck, circNeighVals, criteria_xy);
                criteria_lmax[i] = Toolbox.get_max(circNeighVals);

            }

            for (int i = 0; i < criteria_lmax.length; i++) {
                if (criteria_lmax[i]>globalTh) {
                    mask_xy[i%image_width][i/image_width] = true;
                }
                else {
                    mask_xy[i%image_width][i/image_width] = false;
                }
            }



            criteria_lmax = null;
            for (int i = 0; i < criteria_xy.length; i++) criteria_xy[i] = null;
            criteria_xy=null;

        }

    }

    public static void formRemainingOutputs() // will fill in the remainder of the mapping arrays
    {

        xy2i 	= new int[image_width][image_height];

        int cnt = 0;
        //for (int zz=0; zz<image_length; zz++) {
        for (int xx=0; xx<image_width; xx++) {
            for (int yy=0; yy<image_height; yy++) {

                if (mask_xy[xx][yy]) {
                    xy2i[xx][yy] = cnt;
                    cnt++;
                }
                else {
                    xy2i[xx][yy] = -1;
                }

            }
        }
        //}

//		float perc = (cnt*100f) / (image_width*image_height);
//		IJ.log(String.format("%3.2f %% vol. foreground extracted", perc));

        // component: foreground point list: locations and background estimates
        i2xy = new int[cnt][2];
//		mo.backgroundEst = new byte[cnt];

        cnt =0;
//		for (int zz=0; zz<image_length; zz++) {
        for (int xx=0; xx<image_width; xx++) {
            for (int yy=0; yy<image_height; yy++) {

                if (mask_xy[xx][yy]) {

//						mo.foregroundLocsZXY[cnt][0] = zz;
                    i2xy[cnt][0] = xx;
                    i2xy[cnt][1] = yy;

//						mo.backgroundEst[cnt] = Masker3D.back3[zz][xx][yy];

                    cnt++;

                }

            }
        }
//		}


    }

    public static ImagePlus getCriteria()
    {
        float[][] criteria2d = new float[image_width][image_height];
        for (int i = 0; i < criteria.length; i++) {
            criteria2d[i%image_width][i/image_width] = criteria[i];
        }
        return new ImagePlus("", new FloatProcessor(criteria2d));
    }

}
