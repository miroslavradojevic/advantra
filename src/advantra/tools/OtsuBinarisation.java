package advantra.tools;

import advantra.general.ImageConversions;
import ij.ImagePlus;
import ij.gui.NewImage;

public class OtsuBinarisation {
	
	/*
	 * based on  N. Otsu
	 * A threshold selection method from gray-level histograms
	 * IEEE Transactions on Systems, Man and Cybernetics, Vol. 9, No. 1. (January 1979), pp. 62-66
	 */

	ImagePlus 	img;
	int 		k; 		// can range from 0 to 255
	
	
	public OtsuBinarisation(ImagePlus img){
		// TODO: check this - maybe conversion on float images can be bad!
		// TODO: make it work with float 
		// make it work with 32bit and 16 bit images!
		
		if(img.getType()!=ImagePlus.GRAY8){
			System.out.println("Will be converted to gray8.. can cause errors if 32-bit... ");
		}
		
		this.img = ImageConversions.ImagePlusToGray8(img);
		this.k = 0;
	}
	
	public ImagePlus run(){ // returns a byte image (that's ok)
		
		int[] 	n 		= new int[256];
		float[]	p		= new float[256];
		float[] omega	= new float[256];
		float[] mi		= new float[256];
		int 	N 	= img.getHeight()*img.getWidth()*img.getStack().getSize();
		
		boolean imgIsStack = img.getStack().getSize()>1;
		
		if(imgIsStack){
			for (int i = 0; i < img.getStack().getSize(); i++) {
				int[] hs = img.getStack().getProcessor(i+1).getHistogram();
				for (int idx = 0; idx < hs.length; idx++) {
					n[idx]+= hs[idx];
				}
			}
		}
		else{
			n = img.getProcessor().getHistogram();
		}
		
		// normalize histogram
		float sum_p = 0;
		float sum_i_p = 0;
		for (int i = 0; i < n.length; i++) {
			p[i] = n[i]/(float)N;
			sum_p = sum_p + p[i];
			sum_i_p = sum_i_p + i*p[i];
			omega[i] = sum_p;
			mi[i] = sum_i_p;
		}
		omega[255] = 1f;
		float mi_T = mi[255];
		
		float sigma_b_2_MAX = Float.MIN_VALUE; // between class variance
		// find optimal k
		for (int i = 0; i < 255; i++) {
			float sigma_b_2 = ((mi_T*omega[i]-mi[i])*(mi_T*omega[i]-mi[i])) / (omega[i]*(1-omega[i]));
			if(sigma_b_2>sigma_b_2_MAX){
				sigma_b_2_MAX = sigma_b_2;
				this.k = i;
			}
		}
		
		// create out binary image
		ImagePlus out = NewImage.createByteImage("after binarisation", img.getWidth(), img.getHeight(), img.getStack().getSize(), NewImage.FILL_BLACK);
		for (int i = 0; i < img.getStack().getSize(); i++) {
			byte[] pix_out 	= (byte[])out.getStack().getProcessor(i+1).getPixels();
			byte[] pix_in 	= (byte[])img.getStack().getProcessor(i+1).getPixels();
			for (int j = 0; j < pix_in.length; j++) {
				boolean C = (int)(pix_in[j]&0xff)<=this.k;
				if(C){
					pix_out[j] = (byte)0;
				}
				else{
					pix_out[j] = (byte)255;
				}
			}
			out.getStack().getProcessor(i+1).setPixels(pix_out);
		}
		
		return out;
		
	}
	
	public int getK() {
		return this.k;
	}
	
}
