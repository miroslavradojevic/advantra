package profile;

import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 7/12/13
 * Time: 11:44 AM
 */
public class Masker extends Thread {

    private int begN, endN; // range of locations to work on

    private static int 		image_width;
    private static int 		image_height;

    public static FloatProcessor    inip;
    public static FloatProcessor    back;    // might be useful when hard deciding whether the point belongs to the background or not
    public static ByteProcessor     maskip;  // output

    private static int              nhood;
    //public static float th;               // maybe use it

    private static int              bgComputationMode;                  // background is lower than mean
    private static boolean          localComputation;
    private static float            currBckg;

	public static float VISIBLE_INTENSITY_DIFF = 3;
	public static int VERY_SMALL_REGION_SIZE = 20;

    public Masker (int n0, int n1) { // complete range would be image_width*image_height
        this.begN = n0;
        this.endN = n1;
    }

    public static void loadTemplate(ImageProcessor inip1, float nhood1, int bgComputationMode1, boolean localComputation1)  // float th1
    {
        /*
        set inip
         */
        inip = new FloatProcessor(inip1.getWidth(), inip1.getHeight());
        for (int i=0; i<inip1.getWidth()*inip1.getHeight(); i++) {
            inip.setf(i, inip1.getf(i));
        }

        /*
        initialize maskip
         */
        maskip = new ByteProcessor(inip1.getWidth(), inip1.getHeight()); // all zeros

        /*
        initialize back
         */
        back = new FloatProcessor(inip1.getWidth(), inip1.getHeight()); // all zeros

        image_height 	= inip.getHeight();
        image_width 	= inip.getWidth();
        nhood = (int) Math.ceil(nhood1);
//        th = th1;
        if (bgComputationMode1==0 || bgComputationMode1==1)
            bgComputationMode = bgComputationMode1;  // 0-mean, 1-median
        else
            bgComputationMode = -1;

        localComputation = localComputation1;

        if (localComputation==false) {
            if (bgComputationMode==0) {
                // mean
                currBckg = 0;
                for (int i=0; i<image_width*image_height; i++) currBckg += inip.getf(i);
                currBckg /= image_height*image_width;
            }
            else if (bgComputationMode==1) {
                // median
                float[] array = new float[image_width*image_height];
                for (int i=0; i<array.length; i++) array[i] = inip.getf(i);
                currBckg = (float) Tools.median_Wirth(array);
            }
            else currBckg = Float.NaN;
        }

    }

    public void 	    run()
    {

        if (localComputation==false) {

            // background is global mean/median
            for (int locIdx=begN; locIdx<endN; locIdx++) {
                int atX = locIdx%image_width;
                int atY = locIdx/image_width;
                back.setf(atX, atY, currBckg);

				if (bgComputationMode==0) {
					// mean
					if (inip.getf(atX, atY) > currBckg + VISIBLE_INTENSITY_DIFF) maskip.set(atX, atY, (byte)255);
				}
				else if (bgComputationMode==1) {
					// median
					if (inip.getf(atX, atY) > currBckg) maskip.set(atX, atY, (byte)255);
				}
				else
					maskip.set(atX, atY, (byte)0);

            }

        }
        else { // localComputation==true

            if (bgComputationMode==0) {

                // local MEAN computation in nhood radius
                for (int locIdx=begN; locIdx<endN; locIdx++) {
                    int atX = locIdx%image_width;
                    int atY = locIdx/image_width;

                    if (atX>nhood && atY>nhood && atX<inip.getWidth()-nhood && atY<inip.getHeight()-nhood) { // is in the image

                        //float[] neigh = new float[(2*nhood+1)*(2*nhood+1)];
                        //int idx = 0;
                        float b = 0;
                        for (int locX = atX-nhood; locX<=atX+nhood; locX++) {
                            for (int locY = atY-nhood; locY<=atY+nhood; locY++) {
                                b += inip.getf(locX, locY);
                                //idx++;
                            }
                        }
                        // take the mean as a bkg estimate
                        b = b / ((2*nhood+1)*(2*nhood+1));//float) Tools.median_Wirth(neigh);
                        back.setf(atX, atY, b);
                        if (inip.getf(atX, atY) > b + VISIBLE_INTENSITY_DIFF) maskip.set(atX, atY, (byte)255);

                    }

                }

            }
            else if (bgComputationMode==1) {

                // local MEDIAN computation in nhood radius
                for (int locIdx=begN; locIdx<endN; locIdx++) {
                    int atX = locIdx%image_width;
                    int atY = locIdx/image_width;

                    if (atX>nhood && atY>nhood && atX<inip.getWidth()-nhood && atY<inip.getHeight()-nhood) { // is in the image

                        float[] neigh = new float[(2*nhood+1)*(2*nhood+1)];
                        int idx = 0;
                        for (int locX = atX-nhood; locX<=atX+nhood; locX++) {
                            for (int locY = atY-nhood; locY<=atY+nhood; locY++) {
                                neigh[idx] = inip.getf(locX, locY);
                                idx++;
                            }
                        }
                        // take the median as a bkg estimate
                        float b = (float) Tools.median_Wirth(neigh);
                        back.setf(atX, atY, b);
                        if (inip.getf(atX, atY)>b) maskip.set(atX, atY, (byte)255);

                    }



                }

            }

        }

    }

}