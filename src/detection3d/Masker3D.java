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
	public static int JumpN = 2; 					// calculate every JumpN location (speed up reasons) and vote in between and interpolate background value

	public static int 				image_width;
	public static int 				image_height;
	public static int 				image_length;

	public static float[][][]		instack3;       // input reference

	private static float            radiusCheck;
	private static   float 			iDiff;
    private static int              marginPix;
    private static int              marginLay;

	// outputs
	public static byte[][][]		back3;
	public static boolean[][][]		mask3;

	public Masker3D (int n0, int n1) {
		// complete range would be image_width*image_height*image_length
		this.begN = n0;
		this.endN = n1;
	}

	public static void loadTemplate(float[][][] inimg3d_zxy, int margin1, float radiusNbhoodCheck1, float zDist1, float iDiff1)
	{
		// inis1 is a FloatProcessor stack
		// zDist1 is ia z distance between layers in pixels
        // set radiusCheck at least the outer radius of the sampling points of the sphere, scaled outer radius

		image_height 	= inimg3d_zxy[0][0].length;
		image_width 	= inimg3d_zxy[0].length;
		image_length 	= inimg3d_zxy.length;
		radiusCheck     = radiusNbhoodCheck1;
		iDiff			= iDiff1;

        marginPix       = (int) Math.ceil(radiusCheck);
        marginLay       = (int) Math.ceil(radiusCheck/zDist1);
        // correct margins in case given margin was higher
        marginPix = (margin1>marginPix)? margin1 : marginPix ;  // narrow down selection in XY plane

		/*
			set instack3
		 */
		instack3 = inimg3d_zxy;

		/*
			initialize back3
		 */

		back3 = new byte[image_length][image_width][image_height];

		/*
			initialize mask3
		 */
		mask3 = new boolean[image_length][image_width][image_height];

	}

//    public static void setMaskAtPos(int atX, int atY, int atZ, float sphereRadius) {
//
//        ArrayList<int[]> positions = extractCircularNbhoodIndexes(atX, atY, atZ, sphereRadius);
//            for (int tt=0; tt<positions.size(); tt++) {
//                System.out.println(positions.get(tt)[0]+" , "+positions.get(tt)[1]+" , "+positions.get(tt)[2]);  // zxy
//                mask3[positions.get(tt)[0]][positions.get(tt)[1]][positions.get(tt)[2]] = true;
//            }
//
//    }

//    private static ArrayList<int[]> extractCircularNbhoodIndexes(int atX, int atY, int atZ, float sphereRadius) {
//
//        ArrayList<int[]> out = new ArrayList<int[]>();
//
//        int rPix = Math.round(sphereRadius);
//        int rLay = Math.round(sphereRadius / zDist);
//
//        if (atX-rPix>=0 && atY-rPix>=0 && atZ-rLay>=0 && atX+rPix<image_width && atY+rPix<image_height && atZ+rLay<image_length) {
//
//            int cnt = 0;
//
//            for (int xLoc=atX-rPix; xLoc<=atX+rPix; xLoc++) {
//                for (int yLoc=atY-rPix; yLoc<=atY+rPix; yLoc++) {
//                    for (int zLoc=atZ-rLay; zLoc<=atZ+rLay; zLoc++) { // loop in layers
//
//                        float c = (zLoc-atZ) * zDist;   // back to pixels
//                        if ( (xLoc-atX)*(xLoc-atX)+(yLoc-atY)*(yLoc-atY)+c*c <= rPix*rPix ) {
//
//                            out.add(new int[]{zLoc, xLoc, yLoc});
//
//                        }
//                    }
//                }
//            }
//        }
//
//        return out;
//
//    }

	public static void extractCircularNbhood(int atX, int atY, int atZ, float sphereRadius, float[] values)
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

							values[cnt] = instack3[atZ][xLoc][yLoc];
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

	public static int sizeCircularNbhood(float sphereRadius)
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

    /*
        following statistical measures will be used in ranking the peaks (e.g. median along line)
     */

	public static void fill() {

		for (int locIdx=0; locIdx<image_width*image_height*image_length; locIdx++) {

			int atX = idx2x(locIdx, image_width, image_height);
			int atY = idx2y(locIdx, image_width, image_height);
			int atZ = idx2z(locIdx, image_width, image_height);

			if (atX>=JumpN && atY>=JumpN && atX<image_width-JumpN && atY<image_height-JumpN) {

				int modX = atX % JumpN;
				int modY = atY % JumpN;


				if (modX!=0 && modY==0) {

					back3[atZ][atX][atY] = (byte)Math.round(
							(back3[atZ][atX-modX][atY]      	&0xff) * ((float)(JumpN-modX)/JumpN) +
							(back3[atZ][atX+(JumpN-modX)][atY]	&0xff) * ((float)(modX)      /JumpN)
					);
					mask3[atZ][atX][atY] = mask3[atZ][atX-modX][atY] && mask3[atZ][atX+(JumpN-modX)][atY];

				}
				else if (modX==0 && modY!=0) {

					back3[atZ][atX][atY] = (byte)Math.round(
							(back3[atZ][atX][atY-modY] 		   &0xff) * ((float)(JumpN-modY)/JumpN) +
							(back3[atZ][atX][atY+(JumpN-modY)] &0xff) * ((float)(modY)      /JumpN)
					);
					mask3[atZ][atX][atY] = mask3[atZ][atX][atY-modY] && mask3[atZ][atX][atY+(JumpN-modY)];


				}
				else if (modX!=0 && modY!=0) {

					back3[atZ][atX][atY] = (byte)Math.round(
                            (1f/(JumpN*JumpN)) *
						    (
							(back3[atZ][atX-modX][atY-modY]&0xff) *                    ((float)(JumpN-modX)*(JumpN-modY)) +
							(back3[atZ][atX+(JumpN-modX)][atY-modY]&0xff) *            ((float)(JumpN-modY)*(modX))       +
							(back3[atZ][atX-modX][atY+(JumpN-modY)]&0xff) *            ((float)(JumpN-modX)*(modY))       +
							(back3[atZ][atX+(JumpN-modX)][atY+(JumpN-modY)]&0xff) *    ((float)(modX)*(modY))
                            )
					); // bilinear

					mask3[atZ][atX][atY] =  (mask3[atZ][atX-modX][atY-modY]                 && mask3[atZ][atX+(JumpN-modX)][atY-modY])          ||
                                            (mask3[atZ][atX-modX][atY+(JumpN-modY)]         && mask3[atZ][atX+(JumpN-modX)][atY+(JumpN-modY)])  ||
                                            (mask3[atZ][atX-modX][atY-modY]                 && mask3[atZ][atX-modX][atY+(JumpN-modY)])          ||
                                            (mask3[atZ][atX+(JumpN-modX)][atY-modY]         && mask3[atZ][atX+(JumpN-modX)][atY+(JumpN-modY)]);

				}

			}

		}
	}

	public void run()
	{

        int circNeighSize = sizeCircularNbhood(radiusCheck);
        float[] circNeigh = new float[circNeighSize];                   // allocate array where circular neighbourhood values will be stored

		for (int locIdx=begN; locIdx<endN; locIdx++) {

			int atX = idx2x(locIdx, image_width, image_height);
			int atY = idx2y(locIdx, image_width, image_height);
			int atZ = idx2z(locIdx, image_width, image_height);

			boolean processIt = (atX % JumpN == 0) && (atY % JumpN == 0); // to speed up the calculation
            // exclude those that fit the radius margin
            processIt = processIt && (atX>=marginPix) && (atY>=marginPix) && (atZ>=marginLay) && (atX<image_width-marginPix-1) && (atY<image_height-marginPix-1) && (atZ<image_length-marginLay);


			if (processIt) {

				extractCircularNbhood(atX, atY, atZ, radiusCheck, circNeigh); // extract values from circular neighbourhood, will assign zeros to circNeigh if it is out

				//float locAvgXYZ 	= average(circNeigh);
				//float locStdXYZ 	= std(circNeigh, locAvgXYZ);

				float locMedXYZ 	= Stat.median(circNeigh);

                // calculate 95% median here or take median of smaller circle

				//float locMgnXYZ 	= iDiff; // (locAvgXYZ + 2 * locStdXYZ - locMedXYZ > iDiff)? 0 : iDiff;

				back3[atZ][atX][atY] = (byte)Math.round(locMedXYZ);

				float Ixy = Stat.quantile(circNeigh, 15, 20); // 95%

//				if (
//						(atX-1>=0 && atX+1<image_width) 	&&
//						(atY-1>=0 && atY+1<image_height) 	&&
//						(atZ-1>=0 && atZ+1<image_length)
//				)
//				{
//					Ixy = (instack3[atZ][atX][atY] + instack3[atZ][atX-1][atY] + instack3[atZ][atX+1][atY] + instack3[atZ][atX][atY-1] + instack3[atZ][atX][atY+1]) / 5f;
//				}
//				else {
//					Ixy = instack3[atZ][atX][atY];
//				}

				if (Ixy > locMedXYZ + iDiff) {
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