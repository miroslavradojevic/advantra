package detection2d;

import detection3d.Stat;
import ij.IJ;
import ij.process.ByteProcessor;

import java.util.Arrays;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 12/16/13
 * Time: 3:24 PM
 * Parallel threaded implementation of masker module - robust foreground point extraction
 */
public class Masker2D extends Thread {

	private int begN, endN;

	public static int 				image_width;
	public static int 				image_height;

	public static float[][]			inimg_xy;

	private static float            radiusCheck;
	private static   float 			iDiff;
	private static int              marginPix;

	/*
	OUTPUT
	 */
	public static byte[][]			back_xy;
	public static boolean[][]		mask_xy;
	public static int[][] 			i2xy;
	public static  int[][]			xy2i;

	public Masker2D (int n0, int n1) {
		this.begN = n0;
		this.endN = n1;
	}

	public static void loadTemplate(float[][] _inimg_xy, int _margin, float _radiusNbhoodCheck, float _iDiff)
	{

		radiusCheck     = _radiusNbhoodCheck;
		iDiff			= _iDiff;

		marginPix       = (int) Math.ceil(radiusCheck);
		marginPix = (_margin>marginPix)? _margin : marginPix ;  // narrow down selection in XY plane

		inimg_xy = _inimg_xy;
		image_height 	= inimg_xy[0].length;
		image_width 	= inimg_xy.length;

		/*
			initialize outputs
		 */
		back_xy = new byte[image_width][image_height];
		mask_xy = new boolean[image_width][image_height];
		i2xy 	= new int[1][1];                       // values known after run()
		xy2i 	= new int[1][1];  // values known after run()

	}

	public void run()
	{

		int circNeighSize = sizeCircularNbhood(radiusCheck);
		float[] circNeigh = new float[circNeighSize];

		for (int locIdx=begN; locIdx<endN; locIdx++) {

			int atX = locIdx%image_width;
			int atY = locIdx/image_width;

			boolean processIt = (atX>=marginPix) && (atY>=marginPix) && //(atZ>=marginLay) &&
								(atX<image_width-marginPix-1) && (atY<image_height-marginPix-1);// && (atZ<image_length-marginLay-1);

			if (processIt) {

				extractCircularNbhood(atX, atY, radiusCheck, circNeigh);
				float m50 	= Stat.median(circNeigh);
				float m95 = Stat.quantile(circNeigh, 15, 20);

				back_xy[atX][atY] = (byte) Math.round(m50);
				if (m95 - m50 > iDiff) {
					mask_xy[atX][atY] = true;
				}

			}

		}
	}

	public static void formRemainingOutputs(){

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

		float perc = (cnt*100f) / (image_width*image_height);

		IJ.log(String.format("%3.2f %% vol. foreground extracted", perc));

		if (perc > 80) {// more than 80 percent is wrong
			System.out.println("warning: a lot of foreground points...");
			//return; // return just the lookup table
		}

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

	public static ByteProcessor getMask(){

		byte[] out = new byte[image_height*image_width];
		for (int i=0; i<out.length; i++) {
			out[i] = mask_xy[i%image_width][i/image_width]? (byte)255 : (byte)0;
		}
		return new ByteProcessor(image_width, image_height, out);

	}

	public static ByteProcessor getBackground(){

		byte[] out = new byte[image_height*image_width];
		for (int i=0; i<out.length; i++) {
			out[i] = back_xy[i%image_width][i/image_width];
		}
		return new ByteProcessor(image_width, image_height, out);

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

	private static void extractCircularNbhood(int atX, int atY, float sphereRadius, float[] values)
	{

		int rPix = Math.round(sphereRadius);
		//int rLay = Math.round(sphereRadius / zDist);

		if (atX-rPix>=0 && atY-rPix>=0 && atX+rPix<image_width && atY+rPix<image_height) {  // && atZ-rLay>=0  // && atZ+rLay < image_length

			int cnt = 0;

			for (int xLoc=atX-rPix; xLoc<=atX+rPix; xLoc++) {
				for (int yLoc=atY-rPix; yLoc<=atY+rPix; yLoc++) {
					//for (int zLoc=atZ-rLay; zLoc<=atZ+rLay; zLoc++) { // loop in layers
					//float c = (zLoc-atZ) * zDist;   // back to pixels
					if ( (xLoc-atX)*(xLoc-atX)+(yLoc-atY)*(yLoc-atY) <= rPix*rPix ) { // +c*c

						values[cnt] = inimg_xy[xLoc][yLoc];// instack3[atZ][xLoc][yLoc];
						cnt++;

					}
					//}
				}
			}
		}
		else {
			Arrays.fill(values, 0f); // set all values to zero
		}

	}


}
