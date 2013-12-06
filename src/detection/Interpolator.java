package detection;

import ij.IJ;
import ij.ImageStack;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 6/22/13
 * Time: 5:26 PM
 */
public final class Interpolator {

//	public int          H, W, L;
//	public float[][] 	img_array; 			// elements of 2d or 3d image layers x 2d image size

//	public Interpolator(ImageStack img){
//
//		H = img.getHeight();
//		W = img.getWidth();
//		L = img.getSize();
//
//		img_array = new float[L][]; // will store input values
//
//		for (int i = 0; i < L; i++) {
//
//			float[] lay_fl 	= new float[H*W];
//
//			for (int j = 0; j < lay_fl.length; j++) {
//				lay_fl[j] = img.getProcessor(i+1).getf(j);
//			}
//			img_array[i] = lay_fl;
//
//		}
//
//	}

//	public Interpolator(ImageProcessor img){
//
//		H = img.getHeight();
//		W = img.getWidth();
//		L = 1;
//
//		img_array = new float[L][]; // will store input values
//
//		for (int i = 0; i < L; i++) {
//
//			float[] lay_fl 	= new float[H*W];
//
//			for (int j = 0; j < lay_fl.length; j++) {
//				lay_fl[j] = img.getf(j);
//			}
//			img_array[i] = lay_fl;
//
//		}
//
//	}

//	public float 	interpolateAt(float x, float y){
//
//		float value = 0;
//
//		boolean isIn =
//				y>=0 && y<=(H-1) &&
//						x>=0 && x<=(W-1);
//
//		if(isIn){
//
//			int[] p11 = {(int)Math.floor(y),	(int)Math.floor(x)};
//			int[] p12 = {(int)Math.floor(y), 	(int)Math.ceil(x)};
//			int[] p21 = {(int)Math.ceil(y), 	(int)Math.floor(x)};
//			int[] p22 = {(int)Math.ceil(y), 	(int)Math.ceil(x)};
//
//			float I11_1 = img_array[0][p11[0]*W+p11[1]];
//			float I12_1 = img_array[0][p12[0]*W+p12[1]];
//			float I21_1 = img_array[0][p21[0]*W+p21[1]];
//			float I22_1 = img_array[0][p22[0]*W+p22[1]];
//
//			float a = (p12[1]!=p11[1])?(x-p11[1])/(p12[1]-p11[1]) : 0.5f;	//col
//			float b = (p21[0]!=p11[0])?(y-p11[0])/(p21[0]-p11[0]) : 0.5f;	//row
//
//			float I_1 = (1-a)*(1-b)*I11_1 + (1-a)*b*I21_1 + a*(1-b)*I12_1 + a*b*I22_1;
//
//			value = I_1;
//
//		}
//
//		return value;
//
//	}

//	public float 	interpolateAt(float x, float y, float z){
//
//		float value = 0;
//
//		boolean isIn =
//				y>=0 && y<=(H-1) &&
//						x>=0 && x<=(W-1) &&
//						z>=0 && z<=(L-1);
//
//		if(isIn){
//
//			int[] p11 = {(int)Math.floor(y),	(int)Math.floor(x)};
//			int[] p12 = {(int)Math.floor(y), 	(int)Math.ceil(x)};
//			int[] p21 = {(int)Math.ceil(y), 	(int)Math.floor(x)};
//			int[] p22 = {(int)Math.ceil(y), 	(int)Math.ceil(x)};
//
//			int		l1 = (int)Math.floor(	z);
//			int		l2 = (int)Math.ceil(	z);
//
//			float I11_1 = img_array[l1][p11[0]*W+p11[1]];
//			float I12_1 = img_array[l1][p12[0]*W+p12[1]];
//			float I21_1 = img_array[l1][p21[0]*W+p21[1]];
//			float I22_1 = img_array[l1][p22[0]*W+p22[1]];
//
//			float I11_2 = img_array[l2][p11[0]*W+p11[1]];
//			float I12_2 = img_array[l2][p12[0]*W+p12[1]];
//			float I21_2 = img_array[l2][p21[0]*W+p21[1]];
//			float I22_2 = img_array[l2][p22[0]*W+p22[1]];
//
//			float a = (p12[1]!=p11[1])?(x-p11[1])/(p12[1]-p11[1]) : 0.5f;	//col
//			float b = (p21[0]!=p11[0])?(y-p11[0])/(p21[0]-p11[0]) : 0.5f;	//row
//			float c = (l1!=l2)?(z-l1)/(l2-l1) : 0.5f;						//lay
//
//			float I_1 = (1-a)*(1-b)*I11_1 + (1-a)*b*I21_1 + a*(1-b)*I12_1 + a*b*I22_1;
//			float I_2 = (1-a)*(1-b)*I11_2 + (1-a)*b*I21_2 + a*(1-b)*I12_2 + a*b*I22_2;
//
//			value = (1-c)*I_1+c*I_2;
//
//		}
//		return value;
//	}

	/*
		static methods - image given as input argument
	 */
	public static final float	interpolateAt(double x, double y, FloatProcessor inip2d) {

		float value = 0;

		boolean isIn =
				y>=0 && y<=(inip2d.getHeight()-1) &&
						x>=0 && x<=(inip2d.getWidth()-1);

		if(isIn){

			int[] p11 = {(int)Math.floor(y),	(int)Math.floor(x)};
			int[] p12 = {(int)Math.floor(y), 	(int)Math.ceil(x)};
			int[] p21 = {(int)Math.ceil(y), 	(int)Math.floor(x)};
			int[] p22 = {(int)Math.ceil(y), 	(int)Math.ceil(x)};

			float I11_1 = inip2d.getf(p11[0]*inip2d.getWidth()+p11[1]);
			float I12_1 = inip2d.getf(p12[0]*inip2d.getWidth()+p12[1]);
			float I21_1 = inip2d.getf(p21[0]*inip2d.getWidth()+p21[1]);
			float I22_1 = inip2d.getf(p22[0]*inip2d.getWidth()+p22[1]);

			float a = (p12[1]!=p11[1])?((float)x-p11[1])/(p12[1]-p11[1]) : 0.5f;	//col
			float b = (p21[0]!=p11[0])?((float)y-p11[0])/(p21[0]-p11[0]) : 0.5f;	//row

			float I_1 = (1-a)*(1-b)*I11_1 + (1-a)*b*I21_1 + a*(1-b)*I12_1 + a*b*I22_1;

			value = I_1;

		}

		return value;

	}

	public static final float interpolateAt(float atX, float atY, float atZ_layer, float[][][] img3d_zxy) {

//		float value;

        int x1 = (int) atX;//Math.floor(atX);
        int x2 = x1 + 1;//(int) Math.ceil(atX);
        float x_frac = atX - x1;

        int y1 = (int) atY;//Math.floor(atY);
        int y2 = y1 + 1;//(int) Math.ceil(atY);
        float y_frac = atY - y1;

        int z1 = (int) atZ_layer;//Math.floor(atZ_layer);
        int z2 = z1 + 1;//(int) Math.ceil(atZ_layer);
        float z_frac = atZ_layer - z1;


        // technically this check should be there but is skipped to make things faster
//        boolean isIn =
//						y1>=0 && y2<img3d_zxy[0][0].length &&
//						x1>=0 && x1+1<img3d_zxy[0].length &&
//						z1>=0 && z2<img3d_zxy.length;

		//if(isIn){

			// take neighbourhood
			float I11_1 = img3d_zxy[ z1 ][ x1  ][ y1 ];  // upper left
			float I12_1 = img3d_zxy[ z1 ][ x2  ][ y1 ];  // upper right
			float I21_1 = img3d_zxy[ z1 ][ x1  ][ y2  ]; // bottom left
			float I22_1 = img3d_zxy[ z1 ][ x2  ][ y2  ]; // bottom right

			float I11_2 = img3d_zxy[ z2  ][ x1 ][ y1 ]; // upper left
			float I12_2 = img3d_zxy[ z2  ][ x2 ][ y1 ]; // upper right
			float I21_2 = img3d_zxy[ z2  ][ x1 ][ y2 ]; // bottom left
			float I22_2 = img3d_zxy[ z2  ][ x2 ][ y2 ]; // bottom right

//			float a = ( x1 != x2 )?( atX - x1 ) / ( x2 - x1 ) : 0.5f;	// col
//			float b = ( y1 != y2 )?( atY - y1 ) / ( y2 - y1 ) : 0.5f;	//row
//			float c = ( z1 != z2 )?( atZ_layer - z1 )/( z2 - z1 ) : 0.5f; // lay

			float I_1 =

                    (1-z_frac)  * (  (1-y_frac) * ((1-x_frac)*I11_1 + x_frac*I12_1) + (y_frac) * ((1-x_frac)*I21_1 + x_frac*I22_1) )   +
                    z_frac      * (  (1-y_frac) * ((1-x_frac)*I11_2 + x_frac*I12_2) + (y_frac) * ((1-x_frac)*I21_2 + x_frac*I22_2) );


			//a*(1-b)*I12_1 + a*b*I22_1;
//			float I_2 = (1-a)*(1-b)*I11_2 + (1-a)*b*I21_2 + a*(1-b)*I12_2 + a*b*I22_2;

//			value = (1-c)*I_1+c*I_2;

		//}

		return I_1;

	}

}
