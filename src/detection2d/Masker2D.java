package detection2d;

import aux.Stat;
import ij.IJ;
import ij.ImagePlus;
import ij.io.FileSaver;
import ij.process.ByteProcessor;

import java.io.*;
import java.util.Arrays;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 12/16/13
 * Time: 3:24 PM
 * Parallel threaded implementation of masker module - robust foreground point extraction
 * Separates foreground from the background to reduce the computation by comparing
 * and thresholding the difference between
 * 95th percentile and 50th percentile (median)  of the values within range of radiuses surrounding some location
 * where median is background estimate taken at the smallest radius where the difference above iDiff exists
 */
public class Masker2D extends Thread {

	private int begN, endN;

	public static int 				image_width;
	public static int 				image_height;

	public static float[][]			inimg_xy;

	private static float            radiusCheck;
	private static   float 			iDiff;
	private static int              marginPix;
	private static float[]			rses; // radiuses around the one given as argument
	private static float			alfa = 0.75f;
	private static float 			rmin = 2f;
	private static float			rmax = 20f; // max what you expect the diameter to be


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
		i2xy 	= new int[1][1];  // values known after run()
		xy2i 	= new int[1][1];  // values known after run()

		/*
        multiscale : try different radiuses around the one that was set - to cover more scales
		 */
		rses = new float[3];
		rses[0] = radiusCheck*alfa;
		rses[1] = radiusCheck;
		rses[2] = radiusCheck*(1f/alfa);

		IJ.log("diameters to check:");
		// constrain them
		for (int aa=0; aa<rses.length; aa++) {
			rses[aa] = (rses[aa]<rmin)? rmin : (rses[aa]>rmax)? rmax : rses[aa];
			IJ.log("the neighbourhood diameter is "+rses[aa]);
		}

	}

	public void run()
	{

		//int circNeighSize = sizeCircularNbhood(radiusCheck);
		//float[] circNeigh = new float[circNeighSize];

		for (int locIdx=begN; locIdx<endN; locIdx++) {

			int atX = locIdx%image_width;
			int atY = locIdx/image_width;

			boolean processIt = (atX>=marginPix) && (atY>=marginPix) && //(atZ>=marginLay) &&
								(atX<image_width-marginPix-1) && (atY<image_height-marginPix-1);// && (atZ<image_length-marginLay-1);

			if (processIt) {

				float iDiffMax = Float.NEGATIVE_INFINITY;

				for(int ridx = 0; ridx<rses.length; ridx++) {  // loop rs radiuses

					float[] circNeigh = new float[sizeCircularNbhood(rses[ridx])];
					extractCircularNbhood(atX, atY, rses[ridx], circNeigh);
					float m50 	= Stat.median(circNeigh);
					float m95 = Stat.quantile(circNeigh, 19, 20); // quantile ratio: how much signal there has to be

					if (m95 - m50 > iDiffMax) {
						back_xy[atX][atY] = (byte) Math.round(m50); // background estimate where the difference was the highest
						iDiffMax = m95 - m50;
					}

					if (m95 - m50 > iDiff) {
						mask_xy[atX][atY] = true;
					}

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

    public static void exportI2xyCsv(String file_path) {

        PrintWriter logWriter = null; //initialize writer

        try {
            logWriter = new PrintWriter(file_path);
            logWriter.print("");
            logWriter.close();
        } catch (FileNotFoundException ex) {}   // empty the file before logging...

        try {                                   // initialize detection log file
            logWriter = new PrintWriter(new BufferedWriter(new FileWriter(file_path, true)));
        } catch (IOException e) {}

        for (int ii=0; ii<i2xy.length; ii++) {  //main loop
            for (int jj=0; jj<i2xy[ii].length; jj++) {
                logWriter.print(String.format("%6d", i2xy[ii][jj]));
                if (jj<i2xy[ii].length-1) {
                    logWriter.print(",\t");
                }
                else {
                    logWriter.print("\n");
                }
            }
        }

        logWriter.close(); // close log

    }

	public static void exportXy2iCsv(String file_path) {

		PrintWriter logWriter = null; //initialize writer

		try {
			logWriter = new PrintWriter(file_path);
			logWriter.print("");
			logWriter.close();
		} catch (FileNotFoundException ex) {}   // empty the file before logging...

		try {                                   // initialize detection log file
			logWriter = new PrintWriter(new BufferedWriter(new FileWriter(file_path, true)));
		} catch (IOException e) {}

		for (int ii=0; ii<xy2i.length; ii++) {  //main loop
			for (int jj=0; jj<xy2i[ii].length; jj++) {
				logWriter.print(String.format("%15d", xy2i[ii][jj]));
				if (jj<xy2i[ii].length-1) {
					logWriter.print(",\t");
				}
				else {
					logWriter.print("\n");
				}
			}
		}

		logWriter.close(); // close log

	}

	public static void exportForegroundMask(String file_path) {

		ByteProcessor bp = getMask();
		ImagePlus ip = new ImagePlus("foreground_mask", bp);
		FileSaver fs = new FileSaver(ip);
		fs.saveAsTiff(file_path);

	}

	public static void exportEstBackground(String file_path) {

		ByteProcessor bp = getBackground();
		ImagePlus ip = new ImagePlus("est_background", bp);
		FileSaver fs = new FileSaver(ip);
		fs.saveAsTiff(file_path);

	}

}
