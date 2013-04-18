package advantra.feature;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.filter.Convolver;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

public class GaborFilt2D {
	
	/* calculates bank of Gabor filters over the selected image.
	* 	sigma				-> gaussian envelope, 
	* 	lambda 				-> wavelength 
	* 	(sigma & lambda are related)
	* 	psi 				-> phase offset [in radians] of the cosine factor
	* 	bandwidth			-> spatial frequency bandwidth at half response
	* 	M   				-> number of orientations
	*/
	
	
	/*
	 *  calculation of the ratio sigma/lambda from bandwidth (slratio variable)
	 *  according to Kruizinga and Petkov, 1999 IEEE Trans on Image Processing 8 (10) p.1396
	 *  e.g. BANDWITH = 1 will result in the slratio = 0.56 
	 */
	
	public static int W = 6;
	
	public static ImagePlus run(
			ImagePlus input2D, 
			double theta, 
			double[] sigmas, 
			double[] lambdas, 
			double bandwidth,
			double psi,
			double gamma,
			boolean real){
		
		ImageStack result_stack = new ImageStack(input2D.getWidth(), input2D.getHeight());
		
		for (int i = 0; i < sigmas.length; i++) {
			
			result_stack.addSlice(
					"gabor,theta=" + IJ.d2s(theta, 2) + ",sigma=" + IJ.d2s(sigmas[i], 2),
					run(input2D, theta, sigmas[i], lambdas[i], bandwidth, psi, gamma, real).getProcessor());
			
		}
		
		String result_name = ((real)?"real":"imag")+"_gabor_filt_per_scales";
		
		return new ImagePlus (result_name, result_stack);
		
	}
	
	public static ImagePlus run(
			ImagePlus input2D,
			double[] thetas, 
			double sigma, 
			double lambda, 
			double bandwidth, 
			double psi, 
			double gamma,
			boolean real){

		// if input2D is null just give kernels
		if(input2D==null){
			ImageStack result_kernel = run(null, thetas[0], sigma, lambda, bandwidth, psi, gamma, real).getStack();
			for (int i = 1; i < thetas.length; i++) {
				result_kernel.addSlice(run(null, thetas[i], sigma, lambda, bandwidth, psi, gamma, real).getProcessor());
			}
			return new ImagePlus("gabor_kernels", result_kernel);
		}
		
		ImageStack result_stack = new ImageStack(input2D.getWidth(), input2D.getHeight());
		
		for (int i = 0; i < thetas.length; i++) {
			
			result_stack.addSlice(
					"gabor,theta=" + IJ.d2s(thetas[i], 2) + ",sigma=" + IJ.d2s(sigma, 2),
					run(input2D, thetas[i], sigma, lambda, bandwidth, psi, gamma, real).getProcessor());
		}
		
		String result_name = ((real)?"real":"imag")+"_gabor_filt_per_thetas";
		
		return new ImagePlus (result_name, result_stack);
		
	}
	
	public static ImagePlus run(
			ImagePlus input2D, 
			double theta, 
			double sigma, 
			double lambda, 
			double bandwidth, 
			double psi, 
			double gamma,
			boolean real){

		double slratio = 
				(1/Math.PI) * Math.sqrt((Math.log(2)/2)) * ((Math.pow(2,bandwidth)+1) / (Math.pow(2,bandwidth) - 1) );
		
		if(sigma==0){
			sigma = slratio * lambda;
		}
		else if(lambda==0){
			lambda = sigma / slratio;
		}
		
		double sigma_x = sigma;
		double sigma_y = sigma / gamma;
		double larger_sigma = (sigma_x>sigma_y)? sigma_x : sigma_y ;
		
		// Create set of filters
		int filterSizeX = W * (int)larger_sigma + 1;
		int filterSizeY = W * (int)larger_sigma + 1;
		
		int middleX = (int) Math.round(filterSizeX / 2);
		int middleY = (int) Math.round(filterSizeY / 2);
		
		ImageStack kernels 	= new ImageStack(filterSizeX, filterSizeY);
		ImageProcessor filter = new FloatProcessor(filterSizeX, filterSizeY);  
		
		float sumPos = 0;
		float sumNeg = 0;
		
		// x0 and y0 at patch center (0,0)
	    for (int x=-middleX; x<=middleX; x++){
	        for (int y=-middleY; y<=middleY; y++){           
		        	
	        	double  xr = (double)x * Math.cos(theta) + (double)y * Math.sin(theta);
	        	double  yr = (double)y * Math.cos(theta) - (double)x * Math.sin(theta);
		        	
	        	double env = (1.0/(2.0 * Math.PI * sigma_x*sigma_y)) * 
	        			Math.exp(- 0.5 *((xr*xr)/(sigma_x*sigma_x) + (yr*yr)/(sigma_y*sigma_y)));
		        	
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
		
	    String slice_name = ((real)?"real":"imag") + "_kernel,theta=" + IJ.d2s(theta,2) + ",sigma=" + IJ.d2s(sigma,2);
	    kernels.addSlice(slice_name, filter);
	    
	    if(input2D==null){
	    	return new ImagePlus ("gabor_kernels", kernels);
	    }
	    
	    /*
	     * convolution
	     */
	    
	    Convolver c = new Convolver();                
		c.setNormalize(false); // important not to normalize (did my own normalization)
    
		float[] kernel = (float[]) kernels.getProcessor(1).getPixels();
		ImageProcessor ip = input2D.getProcessor().convertToFloat().duplicate();
		c.convolveFloat(ip, kernel, filterSizeX, filterSizeY);      
 
		int width 	= input2D.getWidth();
		int height 	= input2D.getHeight();
		
		ImageStack filtered	= new ImageStack(width, height);
		filtered.addSlice("gabor,theta=" + IJ.d2s(theta, 2) + ",sigma=" + IJ.d2s(sigma, 2), ip);
		
		String result_name = ((real)?"real":"imag")+"_gabor_filt";
		
		return new ImagePlus (result_name, filtered) ;
		
	}
	
//	public static ImagePlus getKernel(
//			int 		M, 
//			double 		sigma, 
//			double 		lambda, 
//			double 		bandwidth, 
//			double 		psi,
//			double 		gamma,
//			boolean 	real){
//		
//		double slratio = 
//				(1/Math.PI) * Math.sqrt((Math.log(2)/2)) * ((Math.pow(2,bandwidth)+1) / (Math.pow(2,bandwidth) - 1) );
//		
//		if(sigma==0){
//			sigma = slratio * lambda;
//		}
//		else if(lambda==0){
//			lambda = sigma / slratio;
//		}
//		
//		int filterSizeX =  W * (int)sigma + 1;
//		int filterSizeY =  W * (int)sigma + 1;
//		
//		int middleX = (int) Math.round(filterSizeX / 2);
//		int middleY = (int) Math.round(filterSizeY / 2);
//		
//		ImageStack kernels 		= new ImageStack(filterSizeX, filterSizeY);
//		
//		double rotationAngle = Math.PI/(double)M;
//		
//		for (int m=0; m<M; m++){  
//			
//			ImageProcessor filter 	= new FloatProcessor(filterSizeX, filterSizeY);
//			
//			double theta = rotationAngle * m;
//			
//			float sumPos = 0;
//			float sumNeg = 0;
//			
//			// x0 and y0 at patch center (0,0)
//		    for (int x=-middleX; x<=middleX; x++){
//		        for (int y=-middleY; y<=middleY; y++){           
//			        	
//		        	double  xr = (double)x * Math.cos(theta) + (double)y * Math.sin(theta);
//		        	double  yr = (double)y * Math.cos(theta) - (double)x * Math.sin(theta);
//			        	
//		        	double env = (1.0/(2.0 * Math.PI * sigma*sigma)) * 
//		        			Math.exp(- 0.5 *((xr*xr)/(sigma*sigma) + (yr*yr)/(sigma*sigma)));
//			        	
//		        	double carr;
//		        	if(real){
//		            	carr = Math.cos(2 * Math.PI * xr / lambda + psi);
//		            }
//		            else{
//		            	carr = Math.sin(2 * Math.PI * xr / lambda + psi);
//		            }
//		        	
//		        	float coeff = (float)(env*carr);
//		        	
//		        	filter.setf(x+middleX, y+middleY, coeff);
//		        	
//		        	if(coeff>=0){
//		        		sumPos += coeff;
//		        	}
//		        	else{
//		        		sumNeg += Math.abs(coeff);
//		        	}
//			        
//		        }
//		    }
//		    
//		    // normalize
//		    for (int x=0; x<filterSizeX; x++){
//		        for (int y=0; y<filterSizeY; y++){
//		        	float val = filter.getPixelValue(x, y);
//		        	if(val>=0){
//		        		filter.setf(x, y, val/sumPos);
//		        	}
//		        	else{
//		        		filter.setf(x, y, val/sumNeg);
//		        	}
//		        }
//		    }
//		    
//		    String slice_name = 
//		    		((real)?"real":"imag") + "_kernel,theta=" + IJ.d2s(theta,2) + ",sigma=" + IJ.d2s(sigma,2);
//		    kernels.addSlice(slice_name, filter);
//			
//		}
//
//		ImagePlus ip_kernels = new ImagePlus("kernels", kernels);
//		
////		ip_kernels.show();
////		for (int i = 0; i < 8; i++) {
////			ip_kernels.getCanvas().zoomIn(0, 0);
////		}
//		
//		return ip_kernels;
//		
//	}

}