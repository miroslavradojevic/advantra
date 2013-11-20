package detection3d;

import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 11/7/13
 * Time: 2:16 PM
 */
public class Masker3D extends Thread {

	private int begN, endN; 						// range of locations to work on

	public static int 				image_width;
	public static int 				image_height;
	public static int 				image_length;

	public static float[][][]		instack3;
	public static float[][][]		back3;
	public static byte[][][]		mask3;

	private static int              	radius;
	private static float				zDist;
	private static   float 				iDiff;

	public Masker3D (int n0, int n1) {
		// complete range would be image_width*image_height*image_length
		this.begN = n0;
		this.endN = n1;
	}

	public static void loadTemplate(ImageStack inis1, float zDist1, float radius1, float iDiff1)
			// inis1 is a FloatProcessor stack
			// zDist1 is ia z distance between layers in pixels
	{

		image_height 	= inis1.getHeight();
		image_width 	= inis1.getWidth();
		image_length 	= inis1.getSize();
		radius          = (int) Math.ceil(radius1);
		iDiff			= iDiff1;
		zDist			= zDist1;

		/*
			set instack3
		 */
		instack3 = new float[image_length][][];
		for (int ll=0; ll<image_length; ll++) {

			instack3[ll] = new float[image_width][image_height];

			float[] readSlice = (float[]) inis1.getProcessor(ll+1).convertToFloat().getPixels();

			for (int ww=0; ww<image_width; ww++) {
				for (int hh=0; hh<image_height; hh++) {
					instack3[ll][ww][hh] = readSlice[hh*image_width+ww];
				}
			}
		}

		/*
			initialize back3
		 */

		back3 = new float[image_length][image_width][image_height];

		/*
			initialize mask3
		 */
		mask3 = new byte[image_length][image_width][image_height];

	}

	public static void extractCircularNbhood(int atX, int atY, int atZ, float sphereRadius, float[] values)
	{

		int rPix = Math.round(sphereRadius);
		int rLay = Math.round(sphereRadius / zDist);

		if (atX-rPix>=0 && atY-rPix>=0 && atZ-rLay>=0 && atX+rPix<image_width && atY+rPix<image_height && atZ+rLay<image_length) {

			int cnt = 0;

			for (int xLoc=atX-rPix; xLoc<=atX+rPix; xLoc++) {
				for (int yLoc=atY-rPix; yLoc<=atY+rPix; yLoc++) {
					for (int zLoc=atZ-rLay; zLoc<=atZ+rLay; zLoc++) { // loop in layers

						float c = (zLoc-atZ) * zDist;   // back to pixels
						if ( (xLoc-atX)*(xLoc-atX)+(yLoc-atY)*(yLoc-atY)+c*c <= rPix*rPix ) {

							//System.out.println(cnt+" : at (x,y,z): "+xLoc+", "+yLoc+", "+zLoc+" -> " + inis.getProcessor(zLoc+1).getf(xLoc, yLoc));
							values[cnt] = instack3[zLoc][xLoc][yLoc];
							//values[cnt] = inis.getProcessor(zLoc+1).getf(xLoc, yLoc);
							//values[cnt] = instack2[zLoc][xLoc+yLoc*image_width];
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

	public static int sizeCircularNbhood(float sphereRadius)
	{

		// define how many layers will be captured by the radius
		int rLayers 	= Math.round(sphereRadius / zDist);
		int rPix 		= Math.round(sphereRadius);

		int cnt = 0;

		for (int a=-rPix; a<=rPix; a++) {
			for (int b=-rPix; b<=rPix; b++) {
				for (int cLayers=-rLayers; cLayers<=rLayers; cLayers++) { // loop in layers
					float c = cLayers * zDist;   // back to pixels
					if ( a*a+b*b+c*c <= rPix*rPix ) cnt++;
				}
			}
		}

		return cnt;
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

	public static float average(float[] a)
	{
		float meanVal = 0;
		for (int i=0; i<a.length; i++)
			meanVal += a[i];
		return meanVal/a.length;

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

	public void run()
	{
		for (int locIdx=begN; locIdx<endN; locIdx++) {

			int atX = idx2x(locIdx, image_width, image_height);
			int atY = idx2y(locIdx, image_width, image_height);
			int atZ = idx2z(locIdx, image_width, image_height);

			// allocate array where circular neighbourhood values will be stored

			int circNeighSize = sizeCircularNbhood(radius);
			float[] circNeigh = new float[circNeighSize];
			extractCircularNbhood(atX, atY, atZ, radius, circNeigh);

			float locAvgXYZ 	= average(circNeigh);

			float locStdXYZ 	= std(circNeigh, locAvgXYZ);

			float locMedXYZ 	= median(circNeigh);   		// background

			float locIdiffXYZ= locAvgXYZ + 2 * locStdXYZ - locMedXYZ;

			float locMgnXYZ 	= (locIdiffXYZ<=iDiff)? iDiff : 0;

			back3[atZ][atX][atY] = locMedXYZ; //backgr.setf(atX, atY, locMedXY);

			if (atX>0 && atY>0 && atZ>0 && atX<(image_width-radius-1) && atY<(image_height-radius-1) && atZ<(image_length-radius-1)) {

				float Ixy = (instack3[atZ][atX][atY] + instack3[atZ][atX-1][atY] + instack3[atZ][atX+1][atY] + instack3[atZ][atX][atY-1] + instack3[atZ][atX][atY+1]) / 5f;
				if (Ixy > locMedXYZ + locMgnXYZ)
					mask3[atZ][atX][atY] = (byte) 255;

			}

		}
	}

	private static int[] idx2xyz(int idx, int W, int H){
		int[] xyz = new int[3];

		// z = layer
		xyz[2] = idx / (W*H);
		// x ~ width
		xyz[0] = (idx - xyz[2] * H * W) % W;
		// y ~ height
		xyz[1] = (idx - xyz[2] * H * W) / W;

		return xyz;
	}

	private static int idx2x(int idx, int W, int H) {
		return (idx - (idx / (W*H)) * (H*W)) % W;
	}

	private static int idx2y(int idx, int W, int H) {
		return (idx - (idx / (W*H)) * (H*W)) / W;
	}

	private static int idx2z(int idx, int W, int H) {
		return idx / (W*H);
	}

	private static int xyz2ind(int x, int y, int z, int W, int H){

		return H*W*z + y*W + x;

	}

	public static ImageStack getMask()
	{
		ImageStack isOut = new ImageStack(image_width, image_height);
		for (int l=0; l<image_length; l++) {

			byte[] byteArray = new byte[image_width*image_height];

			// align rows
			int cnt = 0;
			for (int y=0; y<image_height; y++) {
				for (int x=0; x<image_width; x++) {
					byteArray[cnt] = mask3[l][x][y];
					cnt++;
				}
			}

			ByteProcessor ipOut = new ByteProcessor(image_width, image_height, byteArray);

			isOut.addSlice(ipOut);


		}
		return isOut;
	}

	public static int getNrLocations() {
		int cnt = 0;
		for (int l=0; l<image_length; l++) {
			for (int x=0; x<image_width; x++) {
				for (int y=0; y<image_height; y++) {
					if (mask3[l][x][y]==(byte)255)
						cnt++;
				}
			}
		}
		return cnt;
	}

	public static ArrayList<int[]> getXyzLocations() {

		// extract list of foreground locations (x,y,z(layer))
		ArrayList<int[]> locs = new ArrayList<int[]>();

		for (int l=0; l<image_length; l++) {
			for (int x=0; x<image_width; x++) {
				for (int y=0; y<image_height; y++) {
					if (mask3[l][x][y]==(byte)255)
						locs.add(new int[]{x, y, l});
				}
			}
		}

		return locs;

	}


}
