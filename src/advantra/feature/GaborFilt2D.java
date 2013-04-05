package advantra.feature;

import java.util.Vector;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.filter.Convolver;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

public class GaborFilt2D {
	
	public static ImagePlus run(ImagePlus input2D, int M, double sigma, double lambda, double bandwidth, double psi, boolean real){
		
		/* This method calculates a bank of Gabor filters over the selected image.
		* 	sigma				-> gaussian envelope, 
		* 	lambda 				-> wavelength 
		* 	(sigma & lambda are related)
		* 	psi 				-> phase offset [in radians] of the cosine factor
		* 	bandwidth			-> spatial frequency bandwidth at half response
		* 	M   				-> number of orientations
		*/
		
		/*
		 *  calculation of the ratio sigma/lambda from bandwidth
		 *  according to Kruizinga and Petkov, 1999 IEEE Trans on Image Processing 8 (10) p.1396
		 *  e.g. BANDWITH = 1 will result in the slratio = 0.56 
		 */
		double slratio = (1/Math.PI) * Math.sqrt( (Math.log(2)/2) ) * ( (Math.pow(2,bandwidth) + 1) / (Math.pow(2,bandwidth) - 1) );
		
		int width 	= input2D.getWidth();
		int height 	= input2D.getHeight();
		
		if(sigma==0){
			sigma = slratio * lambda;
		}
		else if(lambda==0){
			lambda = sigma / slratio;
		}
		
		// Create set of filters
		int filterSizeX = 6 * (int)sigma + 1;
		int filterSizeY = 6 * (int)sigma + 1;
		
		int middleX = (int) Math.round(filterSizeX / 2);
		int middleY = (int) Math.round(filterSizeY / 2);
		
		ImageStack kernels 	= new ImageStack(filterSizeX, filterSizeY);
		 
		double rotationAngle = Math.PI/(double)M;
		
		for (int m=0; m<M; m++){  
				
			double theta = rotationAngle * m;
				
			ImageProcessor filter = new FloatProcessor(filterSizeX, filterSizeY); 
				
			float sumPos = 0;
			float sumNeg = 0;	
			
			// x0 and y0 at patch center (0,0)
		    for (int x=-middleX; x<=middleX; x++){
		        for (int y=-middleY; y<=middleY; y++){           
			        	
		        	double  xr = (double)x * Math.cos(theta) + (double)y * Math.sin(theta);
		        	double  yr = (double)y * Math.cos(theta) - (double)x * Math.sin(theta);
			        	
		        	double env = (1.0/(2.0 * Math.PI * sigma*sigma)) * 
		        			Math.exp(- 0.5 *((xr*xr)/(sigma*sigma) + (yr*yr)/(sigma*sigma)));
			        	
		        	double carr;
		        	if(real){
		            	carr = Math.cos(2 * Math.PI * xr / lambda + psi);
		            }
		            else{
		            	carr = Math.sin(2 * Math.PI * xr / lambda + psi);
		            }
		        	
		        	float coeff = (float)(env*carr);
		        	
		        	filter.setf(x+middleX, y+middleY, coeff);
		        	
		        	if(coeff>=0){
		        		sumPos += coeff;
		        	}
		        	else{
		        		sumNeg += Math.abs(coeff);
		        	}
			        
		        }
		    }
		    
		    // normalize - ensure that the integral of the kernel is 0
		    for (int x=0; x<filterSizeX; x++){
		        for (int y=0; y<filterSizeY; y++){
		        	float val = filter.getPixelValue(x, y);
		        	if(val>=0){
		        		filter.setf(x, y, val/sumPos);
		        	}
		        	else{
		        		filter.setf(x, y, val/sumNeg);
		        	}
		        }
		    }

//		    // to check
//		    float sum = 0;
//		    for (int x=0; x<filterSizeX; x++){
//		        for (int y=0; y<filterSizeY; y++){
//		        	sum += filter.getPixelValue(x, y);
//		        }
//		    }
//		    System.out.format("patch m=%d sigma=%.2f has sum %.2f \n", m, sigma, sum);
		    String slice_name = ((real)?"real":"imag") + "_kernel,theta=" + IJ.d2s(theta,2) + ",sigma=" + IJ.d2s(sigma,2);
		    kernels.addSlice(slice_name, filter);
				
		}
		
		// Show kernels
		if(false){
		ImagePlus ip_kernels = new ImagePlus("kernels", kernels);
		ip_kernels.show();
			for (int i = 0; i < 8; i++) {
				ip_kernels.getCanvas().zoomIn(0, 0);
			}
		}
		
		ImageStack filtered	= new ImageStack(width, height);
		
		// Apply kernels
		for (int m=0; m<M; m++){
		    
				Convolver c = new Convolver();                
				c.setNormalize(false); // important not to normalize (did my own normalization)
		    
				float[] kernel = (float[]) kernels.getProcessor(m+1).getPixels();
				ImageProcessor ip = input2D.getProcessor().convertToFloat().duplicate();
				c.convolveFloat(ip, kernel, filterSizeX, filterSizeY);      
		 
				filtered.addSlice("gabor,m=" + m + ",sigma=" + IJ.d2s(sigma, 2), ip);
				
		}
		
		String result_name = ((real)?"real":"imag")+"_gabor_filt";
		
		ImagePlus result = new ImagePlus (result_name, filtered) ;
		return result;
		
	}

	public static Vector<ImagePlus> runAll(ImagePlus input2D, int M, double[] sigmas, double lambda, double bandwidth, double psi, boolean real){
		
		// do the same just for several sigmas
		Vector<ImagePlus> res = new Vector<ImagePlus>(sigmas.length);
		
		for (int i = 0; i < sigmas.length; i++) {
			ImagePlus img_per_sigma = run(input2D, M, sigmas[i], lambda, bandwidth, psi, real);
			res.add(img_per_sigma);
		}
		
		return res;
		
	}
	
}
//public GaborFilt2D(){}
//public static Vector<Image> calculateArray(Image input2D, double[] sigmas) {
//	Vector<Image> out = new Vector<Image>(sigmas.length);
//	Dimensions dim = input2D.dimensions();
//	double[] aL1 = new double[dim.x];
//	double[] aBness = new double[dim.x];
//	for (int i = 0; i < sigmas.length; i++) {
//		Image Bness = new FloatImage(input2D.dimensions()); Bness.axes(Axes.X);
//		Hessian hs = new Hessian();
//		Vector<Image> L = hs.run(input2D.duplicate(), sigmas[i], false); 
//		L.get(1).axes(Axes.X);
//		Coordinates coords 	= new Coordinates();
//		for (coords.y=0; coords.y<dim.y; ++coords.y) {
//			L.get(1).get(coords,aL1);
//			for (int x=0; x<dim.x; ++x){
//				aBness[x] = Math.abs(aL1[x]);// * aL2[x];
//			}
//			Bness.set(coords, aBness);
//		}
//		out.add(Bness);
//	}
//	return out;
//}