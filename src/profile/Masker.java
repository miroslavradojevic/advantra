package profile;

import ij.IJ;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.util.ArrayList;

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

    private static FloatProcessor    inip;

    public static FloatProcessor    backgr;    		// at one point
    public static FloatProcessor    margin;    		// how much above background
	public static ByteProcessor     mask;  		// refers to region around the location

	private static int              	radius;
	private static int					alloc;

//    private static int              bgComputationMode;

	private static   float 			iDiff;
//	public static 	float HYSTERESIS = 1;
//    public static   float DECAY_STD = .25f*iDiff;

    public Masker (int n0, int n1) { // complete range would be image_width*image_height
        this.begN = n0;
        this.endN = n1;
    }

    public static void loadTemplate(ImageProcessor inip1, float radius1, float iDiff1)
    {
        /*
        set inip
         */
        inip = new FloatProcessor(inip1.getWidth(), inip1.getHeight());
        for (int i=0; i<inip1.getWidth()*inip1.getHeight(); i++) {
            inip.setf(i, inip1.getf(i));
        }

        /*
        initialize mask
         */
        mask = new ByteProcessor(inip1.getWidth(), inip1.getHeight());
		//for (int mm=0; mm<inip1.getWidth() * inip1.getHeight(); mm++) mask.setf(mm, -1);

        /*
        initialize local average
         */
		backgr 		= new FloatProcessor(inip1.getWidth(), inip1.getHeight()); // all zeros
		margin 		= new FloatProcessor(inip1.getWidth(), inip1.getHeight()); // all zeros

        image_height 	= inip.getHeight();
        image_width 	= inip.getWidth();
        radius          = (int) Math.ceil(radius1);
		alloc           = circleElements(radius);
		iDiff			= iDiff1;

    }

    public void run()
    {
		for (int locIdx=begN; locIdx<endN; locIdx++) {

			int atX = locIdx%image_width;
            int atY = locIdx/image_width;

            float[] circNeigh = extractCircleVals(atX, atY, radius, alloc);

			float locAvgXY 	= average(circNeigh);

			float locStdXY 	= std(circNeigh, locAvgXY);

			//float locQ3XY 	= quartile3(circNeigh);

            float locMedXY 	= median(circNeigh);   		// background

            float locIdiffXY= locAvgXY + 2 * locStdXY - locMedXY;
            //float locIdiffXY= locQ3XY - locMedXY;

            float locMgnXY 	= (locIdiffXY<=iDiff)? iDiff : 0;// I_DIFF*(float)Math.exp(- (locIdiffXY-I_DIFF)*(locIdiffXY-I_DIFF) / (2*DECAY_STD*DECAY_STD) ) ;

			backgr.setf(atX, atY, locMedXY);
			margin.setf(atX, atY, locMgnXY);

			//create mask (these locations will be considered in detection)
			if (atX>radius && atY>radius && atX<inip.getWidth()-radius && atY<inip.getHeight()-radius) { // is in the image
				float Ixy = (inip.getf(atX, atY) + inip.getf(atX-1, atY) + inip.getf(atX+1, atY) + inip.getf(atX, atY-1) + inip.getf(atX, atY+1)) / 5f;
				if (Ixy > locMedXY + locMgnXY) mask.set(atX, atY, (byte)255);
			}

        }
    }

	public static int getMaskLocationsTotal() // makes sense after run()
	{

		int out = 0;
		for (int a=0; a<mask.getWidth()*mask.getHeight(); a++) {
			if (mask.get(a)==255) {
				out++;
			}
		}
		return out;

	}

	public static float[] extractCircleVals(int x, int y, int r, int alloc)
	{

		float[] out = new float[alloc];

		if (x>r && y>r && x<inip.getWidth()-r && y<inip.getHeight()-r) {
			int idx = 0;
			for (int locX = x-r; locX<=x+r; locX++) {
				for (int locY = y-r; locY<=y+r; locY++) {
					if ((locX-x)*(locX-x)+(locY-y)*(locY-y)<=r*r) {
						out[idx] = inip.getf(locX, locY);
						idx++;
					}
				}
			}
		}

		return out;
	}

	public static int circleElements(int r)
	{
		int cnt = 0;
		for (int a=-r; a<=r; a++) {
			for (int b=-r; b<=r; b++) {
				if (a*a+b*b<=r*r) cnt++;
			}
		}
		return cnt;
	}

	/*
	average
	 */

	public static float[] averageVec(float[] in)
	{
		float[] out = new float[in.length];

		float avg = average(in);

		for (int i=0; i<in.length; i++) {
			out[i] = avg;
		}

		return out;
	}

	public static float average(float[] in)
	{
		float avg = 0;
		for (int i=0; i<in.length; i++) {
			avg += in[i];
		}
		avg /= in.length;
		return avg;
	}

	/*
	median
	 */

	public static float[] medianVec(float[] in)
	{
		float[] out = new float[in.length];

		float avg = median(in);

		for (int i=0; i<in.length; i++) {
			out[i] = avg;
		}

		return out;
	}

	public static float median(float[] a)
	{
		int n = a.length;
		int i, j, l, m, k;
		double x;
		if (n % 2 == 0) k = (n/2)-1;
		else k = (n/2);
		l=0 ; m=n-1 ;
		while (l < m)
		{
			x=a[k] ;
			i = l ;
			j = m ;
			do
			{
				while (a[i] < x) i++ ;
				while (x < a[j]) j-- ;
				if (i <= j) {
					float temp = a[i];
					a[i] = a[j];
					a[j] = temp;
					i++ ; j-- ;
				}
			} while (i <= j) ;
			if (j < k) l = i ;
			if (k < i) m = j ;
		}
		return a[k] ;
	}

	/*
    kth smallest
	 */

	public static float[] quartile3Vec(float[] in)
	{
		float[] out = new float[in.length];

		float q3 = quartile3(in);

		for (int i=0; i<in.length; i++) {
			out[i] = q3;
		}

		return out;

	}

	public static float quartile3(float[] a)
	{
		int n = a.length;
		int i, j, l, m, k;
		double x;

		if (n % 2 == 0) k = ((3*n)/4)-1;
		else k = ((3*n)/4);

		l=0 ; m=n-1 ;
		while (l < m)
		{
			x=a[k] ;
			i = l ;
			j = m ;
			do
			{
				while (a[i] < x) i++ ;
				while (x < a[j]) j-- ;
				if (i <= j) {
					float temp = a[i];
					a[i] = a[j];
					a[j] = temp;
					i++ ; j-- ;
				}
			} while (i <= j) ;
			if (j < k) l = i ;
			if (k < i) m = j ;
		}
		return a[k] ;
	}

	/*
	standard deviation
	 */

	public static float[] stdVec(float[] in, float avg)
	{
		float[] out = new float[in.length];
		float std = std(in, average(in));
		for (int i=0; i<in.length; i++) {
			out[i] = std;
		}
		return out;
	}

	public static float std(float[] in, float avg)
	{
		float std = 0;
		for (int i=0; i<in.length; i++) {
			std += (in[i]-avg)*(in[i]-avg);
		}
		std /= in.length;
		std = (float) Math.sqrt(std);
		return std;
	}

}