package advantra.feature;

import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.ZProjector;
import ij.plugin.filter.Convolver;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

public class GaborFilt2D {
	
	public static ImagePlus run(ImagePlus input2D, double sigma, double gamma, double psi, double Fx, int nAngles, boolean real){
	
		/*
		* This method calculates a set of Gabor filters over the selected image.
		*
		* Parameters: sigma, gamma, psi, Fx, nAngles
		* 
		* Sigma defining the size of the Gaussian envelope
		* Aspect ratio of the Gaussian curves
		* Phase
		* Frequency of the sinusoidal component
		* Number of diferent orientation angles to use
		*/
		
//		Dimensions dim = input2D.dimensions();
//		Image out = new FloatImage(new Dimensions(dim.x, dim.y, 1));
//		out.axes(Axes.X);
		
		int width 	= input2D.getWidth();//dim.x;
		int height 	= input2D.getHeight();//dim.y;
		
		double sigma_x = sigma;
		double sigma_y = sigma / gamma;
		
		int largerSigma = (sigma_x > sigma_y) ? (int) sigma_x : (int) sigma_y;
		if(largerSigma < 1)
		    largerSigma = 1;
		
		double sigma_x2 = sigma_x * sigma_x;
		double sigma_y2 = sigma_y * sigma_y;
		
		// Create set of filters
		int filterSizeX = 19; //6 * largerSigma + 1;
		int filterSizeY = 19; //6 * largerSigma + 1;
		
		int middleX = (int) Math.round(filterSizeX / 2);
		int middleY = (int) Math.round(filterSizeY / 2);
		
		ImageStack is 		= new ImageStack(width, height);
		ImageStack kernels 	= new ImageStack(filterSizeX, filterSizeY);
		 
		double rotationAngle = Math.PI/(double)nAngles;
		
		for (int i=0; i<nAngles; i++){   
		    double theta = rotationAngle * i;
		    ImageProcessor filter = new FloatProcessor(filterSizeX, filterSizeY);  
		    for (int x=-middleX; x<=middleX; x++){
		        for (int y=-middleY; y<=middleY; y++){           
		        	
		        	double  xPrime = (double)x * Math.cos(theta) + (double)y * Math.sin(theta);
		        	double  yPrime = (double)y * Math.cos(theta) - (double)x * Math.sin(theta);
		                 
		            double a = 1.0 / ( 2.0 * Math.PI * sigma_x * sigma_y ) *
		                            Math.exp(-0.5 * (xPrime*xPrime / sigma_x2 + yPrime*yPrime / sigma_y2) );
		            double c;
		            if(real){
		            	c = Math.cos( 2.0 * Math.PI * (Fx * xPrime) / filterSizeX + psi);
		            }
		            else{
		            	c = Math.sin( 2.0 * Math.PI * (Fx * xPrime) / filterSizeX + psi);
		            }
		             
		            filter.setf(x+middleX, y+middleY, (float)(a*c) );
		        }
		    }
		    kernels.addSlice("kernel angle = " + theta, filter);
		}
		
		// Show kernels
		ImagePlus ip_kernels = new ImagePlus("kernels", kernels);
		ip_kernels.show();
		
		// Apply kernels
		for (int i=0; i<nAngles; i++){
		    
		    Convolver c = new Convolver();                
		    c.setNormalize(false);
		    
		    float[] kernel = (float[]) kernels.getProcessor(i+1).getPixels();
		    
		    ImageProcessor ip = input2D.getProcessor().convertToFloat().duplicate();
		    
		    c.convolveFloat(ip, kernel, filterSizeX, filterSizeY);      
		 
		    is.addSlice("gabor angle = " + i, ip);
		}
		
		ImagePlus projectStack = new ImagePlus("filtered stack",is);
		
		ImageStack resultStack = new ImageStack(width, height);
        
		ZProjector zp = new ZProjector(projectStack);
		zp.setStopSlice(is.getSize());
		for (int i=0;i<=nAngles; i++)
		{
		    zp.setMethod(i);
		    zp.doProjection();
		    resultStack.addSlice("Gabor_" + i 
		            +"_"+sigma+"_" + gamma + "_"+ (int) (psi / (Math.PI/4) ) +"_"+Fx, 
		            zp.getProjection().getChannelProcessor());
		}
		 
		// Display filtered images
		(new ImagePlus("gabor, sigma="+sigma+" gamma="+gamma+ " psi="+psi, is)).show();
		 
		ImagePlus result = new ImagePlus ("Gabor stack projections", resultStack) ;
		
//		double[][] aIn 			= new double[height][width];
//		Coordinates coords 	= new Coordinates();
		
//		for (int i = 0; i < sigmas.length; ++i) {
//			Vector<Image> L = hs.run(input2D, sigmas[i], true);
//			L.get(1).axes(Axes.X); // smaller abs value
//			for (coords.y=0; coords.y<dim.y; ++coords.y) {
//				L.get(1).get(coords,aIn);
//				for (int x = 0; x < dim.x; ++x) {
//					aScale_space[x] = aIn[x];
//				}
//				coords.z = i;
//				out.set(coords, aScale_space);
//				coords.z = 0;
//			}
//		}
		
		return result;
		
	}

}


//public GaborFilt2D(){}

//public static Vector<Image> calculateArray(Image input2D, double[] sigmas) {
//	
//	Vector<Image> out = new Vector<Image>(sigmas.length);
//	
//	Dimensions dim = input2D.dimensions();
//	double[] aL1 = new double[dim.x];
//	double[] aBness = new double[dim.x];
//	
//	for (int i = 0; i < sigmas.length; i++) {
//		
//		Image Bness = new FloatImage(input2D.dimensions()); Bness.axes(Axes.X);
//		
//		Hessian hs = new Hessian();
//		Vector<Image> L = hs.run(input2D.duplicate(), sigmas[i], false); 
//		L.get(1).axes(Axes.X);
//		
//		Coordinates coords 	= new Coordinates();
//		for (coords.y=0; coords.y<dim.y; ++coords.y) {
//			L.get(1).get(coords,aL1);
//			for (int x=0; x<dim.x; ++x){
//				aBness[x] = Math.abs(aL1[x]);// * aL2[x];
//			}
//			Bness.set(coords, aBness);
//		}
//		
//		out.add(Bness);
//		
//	}
//	return out;
//}