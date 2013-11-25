package detection3d;

import ij.ImageStack;
import ij.process.ByteProcessor;

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
	public static int JumpN = 3; 					// calculate every JumpN location (speed up reasons) and vote in between and interpolate background value

	public static int 				image_width;
	public static int 				image_height;
	public static int 				image_length;

	public static float[][][]		instack3;       // input reference

	private static float              	radius;
	private static float				zDist;
	private static   float 				iDiff;

	// outputs
	public static float[][][]		back3;
	public static boolean[][][]		mask3;

	public Masker3D (int n0, int n1) {
		// complete range would be image_width*image_height*image_length
		this.begN = n0;
		this.endN = n1;
	}

	public static void loadTemplate(float[][][] inimg3d_zxy, float zDist1, float radius1, float iDiff1)
	{
		// inis1 is a FloatProcessor stack
		// zDist1 is ia z distance between layers in pixels

		image_height 	= inimg3d_zxy[0][0].length;
		image_width 	= inimg3d_zxy[0].length;
		image_length 	= inimg3d_zxy.length;
		radius          = radius1;//(int) Math.ceil(radius1);
		iDiff			= iDiff1;
		zDist			= zDist1;

		/*
			set instack3
		 */
		instack3 = inimg3d_zxy;

		/*
			initialize back3
		 */

		back3 = new float[image_length][image_width][image_height];

		/*
			initialize mask3
		 */
		mask3 = new boolean[image_length][image_width][image_height];

	}

    public static void setMaskAtPos(int atX, int atY, int atZ, float sphereRadius) {

            ArrayList<int[]> positions = extractCircularNbhoodIndexes(atX, atY, atZ, sphereRadius);
            for (int tt=0; tt<positions.size(); tt++) {
                System.out.println(positions.get(tt)[0]+" , "+positions.get(tt)[1]+" , "+positions.get(tt)[2]);  // zxy
                mask3[positions.get(tt)[0]][positions.get(tt)[1]][positions.get(tt)[2]] = true;
            }


    }

    private static ArrayList<int[]> extractCircularNbhoodIndexes(int atX, int atY, int atZ, float sphereRadius) {

        ArrayList<int[]> out = new ArrayList<int[]>();

        int rPix = Math.round(sphereRadius);
        int rLay = Math.round(sphereRadius / zDist);

        if (atX-rPix>=0 && atY-rPix>=0 && atZ-rLay>=0 && atX+rPix<image_width && atY+rPix<image_height && atZ+rLay<image_length) {

            int cnt = 0;

            for (int xLoc=atX-rPix; xLoc<=atX+rPix; xLoc++) {
                for (int yLoc=atY-rPix; yLoc<=atY+rPix; yLoc++) {
                    for (int zLoc=atZ-rLay; zLoc<=atZ+rLay; zLoc++) { // loop in layers

                        float c = (zLoc-atZ) * zDist;   // back to pixels
                        if ( (xLoc-atX)*(xLoc-atX)+(yLoc-atY)*(yLoc-atY)+c*c <= rPix*rPix ) {

                            out.add(new int[]{zLoc, xLoc, yLoc});

                        }
                    }
                }
            }
        }

        return out;

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

							values[cnt] = instack3[zLoc][xLoc][yLoc];
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

	public static void fill() {

		for (int locIdx=0; locIdx<image_width*image_height*image_length; locIdx++) {

			int atX = idx2x(locIdx, image_width, image_height);
			int atY = idx2y(locIdx, image_width, image_height);
			int atZ = idx2z(locIdx, image_width, image_height);

			if (atX>=JumpN && atY>=JumpN && atX<image_width-JumpN && atY<image_height-JumpN) {

				int modX = atX % JumpN;
				int modY = atY % JumpN;


				if (modX!=0 && modY==0) {

					back3[atZ][atX][atY] = back3[atZ][atX-modX][atY] * ((float)(JumpN-modX)/JumpN) + back3[atZ][atX+(JumpN-modX)][atY] * ((float)(modX)/JumpN);
					mask3[atZ][atX][atY] = mask3[atZ][atX-modX][atY] && mask3[atZ][atX+(JumpN-modX)][atY];

				}
				else if (modX==0 && modY!=0) {

					back3[atZ][atX][atY] = back3[atZ][atX][atY-modY] * ((float)()/) + back3[][][] * ();
					mask3[atZ][atX][atY] = mask3[atZ][atX][atY-modY] && mask3[atZ][atX][atY+(JumpN-modY)];


				}
				else if (modX!=0 && modY!=0) {

					back3[atZ][atX][atY] = ;
					mask3[atZ][atX][atY] = ;

				}



			}






		}
	}

	public void run()
	{

        int circNeighSize = sizeCircularNbhood(radius);
        float[] circNeigh = new float[circNeighSize];                   // allocate array where circular neighbourhood values will be stored

		for (int locIdx=begN; locIdx<endN; locIdx++) {

			int atX = idx2x(locIdx, image_width, image_height);
			int atY = idx2y(locIdx, image_width, image_height);
			int atZ = idx2z(locIdx, image_width, image_height);

			boolean processIt = (atX % JumpN == 0) && (atY % JumpN == 0); // to speed up the calculation

			if (processIt) {

				extractCircularNbhood(atX, atY, atZ, radius, circNeigh); // extract values from circular neighbourhood, will assign zeros to circNeigh if it is out

				float locAvgXYZ 	= average(circNeigh);

				float locStdXYZ 	= std(circNeigh, locAvgXYZ);

				float locMedXYZ 	= median(circNeigh);

				float locMgnXYZ 	= (locAvgXYZ + 2 * locStdXYZ - locMedXYZ > iDiff)? 0 : iDiff;

				back3[atZ][atX][atY] = locMedXYZ;

				float Ixy;
				if (
						(atX-1>=0 && atX+1<image_width) 	&&
						(atY-1>=0 && atY+1<image_height) 	&&
						(atZ-1>=0 && atZ+1<image_length)
				)
				{
					Ixy = (instack3[atZ][atX][atY] + instack3[atZ][atX-1][atY] + instack3[atZ][atX+1][atY] + instack3[atZ][atX][atY-1] + instack3[atZ][atX][atY+1]) / 5f;
				}
				else {
					Ixy = instack3[atZ][atX][atY];
				}

				if (Ixy > locMedXYZ + locMgnXYZ) {
					mask3[atZ][atX][atY] = true;
                }

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
					if (mask3[l][x][y]) {
						byteArray[cnt] = (byte)255;

					}
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
					if (mask3[l][x][y])
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
					if (mask3[l][x][y])
						locs.add(new int[]{x, y, l});
				}
			}
		}

		return locs;

	}

}
