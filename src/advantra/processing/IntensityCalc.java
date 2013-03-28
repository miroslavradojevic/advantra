package advantra.processing;

import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.NewImage;
import imagescience.image.Axes;
import imagescience.image.Coordinates;
import imagescience.image.Image;

public class IntensityCalc {
	
	public int H, W, L;
	public float[][] 	img_array; // perhaps make it just 1D instead of 2D for more speed
	double[][] temp;
	
	public IntensityCalc(ImageStack img){
		H = img.getHeight();
		W = img.getWidth();
		L = img.getSize();
		
		img_array = new float[L][]; //H*W
		for (int i = 0; i < L; i++) {
			img_array[i] = (float [])img.getProcessor(i+1).convertToFloat().getPixels(); 
		}
		
		// alternative (just to test)
		temp = new double[H][W];
		ImagePlus  temp_ip = new ImagePlus("ij lib", img);
		temp_ip.show();
		
		Image temp_img =Image.wrap(temp_ip); 
		temp_img.name("eriks lib");
		temp_img.imageplus().show();
		temp_img.axes(Axes.X+Axes.Y);
		temp_img.get(new Coordinates(), temp);
	}
	
	public float 	interpolateAt(float p1, float p2){ // p1~row, p2~col
		
		float value = 0;
		
		boolean isIn = 
				p1>=0 && p1<=(H-1) &&
				p2>=0 && p2<=(W-1);
		
		if(isIn){
			int[] p11 = {(int)Math.floor(p1),	(int)Math.floor(p2)}; 	
			int[] p12 = {(int)Math.floor(p1), 	(int)Math.ceil(p2)}; 		
			int[] p21 = {(int)Math.ceil(p1), 	(int)Math.floor(p2)}; 	
			int[] p22 = {(int)Math.ceil(p1), 	(int)Math.ceil(p2)};
			
			float I11_1 = img_array[0][p11[0]*W+p11[1]]; 
			float I12_1 = img_array[0][p12[0]*W+p12[1]];
			float I21_1 = img_array[0][p21[0]*W+p21[1]]; 
			float I22_1 = img_array[0][p22[0]*W+p22[1]]; 
			
			float a = (p12[1]!=p11[1])?(p2-p11[1])/(p12[1]-p11[1]) : 0.5f;	//col
			float b = (p21[0]!=p11[0])?(p1-p11[0])/(p21[0]-p11[0]) : 0.5f;	//row
			
			float I_1 = (1-a)*(1-b)*I11_1 + (1-a)*b*I21_1 + a*(1-b)*I12_1 + a*b*I22_1;
			
			System.out.format("(%.2f, %.2f), a=%.2f, b=%.2f, (%.2f (%.2f), %.2f (%.2f), %.2f (%.2f), %.2f (%.2f)), I=%.2f\n", 
					p1, p2, a, b, 
					I11_1, temp[p11[0]][p11[1]], 
					I12_1, temp[p12[0]][p12[1]], 
					I21_1, temp[p21[0]][p21[1]],
					I22_1, temp[p22[0]][p22[1]],
					I_1);
			
			value = I_1;
			
			
		}
				
		return value;		
		
	}

	public float 	interpolateAt(float p1, float p2, float p3){// p1, p2, p3 ---> row, col, lay
		
		float value = 0;
		
		boolean isIn = 
				p1>=0 && p1<=(H-1) && 
				p2>=0 && p2<=(W-1) && 
				p3>=0 && p3<=(L-1);
				
		if(isIn){
			int[] p11 = {(int)Math.floor(p1),	(int)Math.floor(p2)}; 	
			int[] p12 = {(int)Math.floor(p1), 	(int)Math.ceil(p2)}; 		
			int[] p21 = {(int)Math.ceil(p1), 	(int)Math.floor(p2)}; 	
			int[] p22 = {(int)Math.ceil(p1), 	(int)Math.ceil(p2)};
			
			int		l1 = (int)Math.floor(	p3);
			int		l2 = (int)Math.ceil(	p3);
			
			float I11_1 = img_array[l1][p11[0]*W+p11[1]]; 
			float I12_1 = img_array[l1][p12[0]*W+p12[1]];
			float I21_1 = img_array[l1][p21[0]*W+p21[1]]; 
			float I22_1 = img_array[l1][p22[0]*W+p22[1]]; 

			float I11_2 = img_array[l2][p11[0]*W+p11[1]]; 
			float I12_2 = img_array[l2][p12[0]*W+p12[1]]; 
			float I21_2 = img_array[l2][p21[0]*W+p21[1]]; 
			float I22_2 = img_array[l2][p22[0]*W+p22[1]];
			
			float a = (p12[1]!=p11[1])?(p2-p11[1])/(p12[1]-p11[1]) : 0.5f;	//col
			float b = (p21[0]!=p11[0])?(p1-p11[0])/(p21[0]-p11[0]) : 0.5f;	//row
			float c = (l1!=l2)?(p3-l1)/(l2-l1) : 0.5f;						//lay
			
			float I_1 = (1-a)*(1-b)*I11_1 + (1-a)*b*I21_1 + a*(1-b)*I12_1 + a*b*I22_1;
			float I_2 = (1-a)*(1-b)*I11_2 + (1-a)*b*I21_2 + a*(1-b)*I12_2 + a*b*I22_2;
			
			value = (1-c)*I_1+c*I_2;
			
		}
		return value;
	}

}