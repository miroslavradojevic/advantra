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

    public static int 		image_width;
    public static int 		image_height;

    public static FloatProcessor    inip;
    public static FloatProcessor    back;
    public static ByteProcessor     maskip;  // output

    public static int   nhood;
    public static float th;


    public Masker (int n0, int n1) {
        this.begN = n0;
        this.endN = n1;
    }

    public static void loadTemplate(ImageProcessor inip1, int nhood1, float th1)
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
        nhood = nhood1;
        th = th1;

    }

    public void 	    run()
    {
        // loop all locations in defined range
        for (int locIdx=begN; locIdx<endN; locIdx++) {

            int atX = locIdx%image_width;//Profiler.locations[locIdx][0];
            int atY = locIdx/image_width;//Profiler.locations[locIdx][1];

            if (atX>nhood && atY>nhood && atX<inip.getWidth()-nhood && atY<inip.getHeight()-nhood) { // is in the image


                /*
                background using median
                 */

                float[] neigh = new float[(2*nhood+1)*(2*nhood+1)];
                // fill the array in
                int idx = 0;
                for (int locX = atX-nhood; locX<=atX+nhood; locX++) {
                    for (int locY = atY-nhood; locY<=atY+nhood; locY++) {
                        neigh[idx] = inip.getf(locX, locY);
                        idx++;
                    }
                }
                // take the median as a bkg estimate
                float currBkgr = (float) Tools.median_Wirth(neigh);






                back.setf(atX, atY, currBkgr);

                //float currVal  = inip.getf(atX, atY);
                // check the neighbours how they relate to backgraound value

                boolean isForeground =
                        //currVal-currBkgr>th ||

//                        inip.getf(atX-1, atY-1)>currBkgr+th ||
//                        inip.getf(atX-1, atY+0)>currBkgr+th ||
//                        inip.getf(atX-1, atY+1)>currBkgr+th ||
//
//                                inip.getf(atX+0, atY-1)>currBkgr+th ||
                                inip.getf(atX+0, atY+0)>currBkgr+th
//                                ||
//                                inip.getf(atX+0, atY+1)>currBkgr+th ||
//
//                                inip.getf(atX+1, atY-1)>currBkgr+th ||
//                                inip.getf(atX+1, atY+0)>currBkgr+th ||
//                                inip.getf(atX+1, atY+1)>currBkgr+th
                        ;

                if (isForeground) {
                    maskip.set(atX, atY, (byte)255);
                }

            }

        }

    }


}
