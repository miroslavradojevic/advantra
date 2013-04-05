package advantra.feature;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.filter.Convolver;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

public class GaborFilt2D {
	
	public static ImagePlus run(ImagePlus input2D, int N, int M, double sigma, boolean real){
		/*
		 * 
		* This method calculates a bank of Gabor filters over the selected image.
		* Parameters: 
		* sigma				-> gaussian envelope, 
		* psi    			-> phase
		* N		    		-> wavelength lambda = patch_size/N 
		* M   				-> number of orientations
		*/
		
		double psi = 0;
		
		int width 	= input2D.getWidth();
		int height 	= input2D.getHeight();
		
		// correct just in case, lower than 1 doesn't make sense
		sigma = (sigma<1)? 1 : sigma;
		
		// larger sigma_ will define the size of the patch
//		int largerSigma = (int)Math.max(sigma*Math.pow(2, N-1), sigma*Math.pow(2, N-1));
		
//		double a = 1 / sigma_x;
//		double b = 1 / sigma_y;
//		
//		double a2 = a * a;
//		double b2 = b * b;
		
		// Create set of filters
		int filterSizeX = 6 * (int)sigma + 1;
		int filterSizeY = 6 * (int)sigma + 1;
		
		int middleX = (int) Math.round(filterSizeX / 2);
		int middleY = (int) Math.round(filterSizeY / 2);
		
		ImageStack filtered	= new ImageStack(width, height);
		ImageStack kernels 	= new ImageStack(filterSizeX, filterSizeY);
		 
		double rotationAngle = Math.PI/(double)M;
		
		for (int n = 1; n <= N; n++) {
			
			for (int m=0; m<M; m++){  
				
				double theta = rotationAngle * m;
				
				ImageProcessor filter = new FloatProcessor(filterSizeX, filterSizeY); 
				
				// x0 and y0 at patch center (0,0)
			    for (int x=-middleX; x<=middleX; x++){
			        for (int y=-middleY; y<=middleY; y++){           
			        	
			        	double  xr = (double)x * Math.cos(theta) + (double)y * Math.sin(theta);
			        	double  yr = (double)y * Math.cos(theta) - (double)x * Math.sin(theta);
			        	
			        	double env = (1.0/(2.0 * Math.PI * sigma*sigma)) * 
			        			Math.exp(- 0.5 *((xr*xr)/(sigma*sigma) + (yr*yr)/(sigma*sigma)));
			        	
			        	double carr;
			        	if(real){
			            	carr = Math.cos(2 * Math.PI * n * xr / filterSizeX + psi);// lambda is filterSizeX
			            }
			            else{
			            	carr = Math.sin(2 * Math.PI * n * xr / filterSizeX + psi);
			            }
			        	
			        	filter.setf(x+middleX, y+middleY, (float)(env*carr) );
			        	
			        }
			    }
				
			    kernels.addSlice("kernel,theta=" + IJ.d2s(theta,2) + ",lambda=" + IJ.d2s(filterSizeX/(double)n,2), filter);
				
			}
		}
		
		// Show kernels
		ImagePlus ip_kernels = new ImagePlus("kernels", kernels);
		ip_kernels.show();
		
		for (int i = 0; i < 8; i++) {
			ip_kernels.getCanvas().zoomIn(0, 0);
		}
		
//		if(true) return;
		
		
		
		// Apply kernels
		for (int n = 1; n <= N; n++) {
			for (int i=0; i<M; i++){
		    
				Convolver c = new Convolver();                
				c.setNormalize(false); // important not to normalize gabor because they sum to 0 (imagej hack!)
		    
				float[] kernel = (float[]) kernels.getProcessor(i+1).getPixels();
				ImageProcessor ip = input2D.getProcessor().convertToFloat().duplicate();
				c.convolveFloat(ip, kernel, filterSizeX, filterSizeY);      
		 
				filtered.addSlice("gabor,theta=" + i + ",lambda=" + IJ.d2s(filterSizeX/(double)n,2), ip);
				
			}
		}
		
//		ImagePlus projectStack = new ImagePlus("filtered stack",filtered);
//		ImageStack resultStack = new ImageStack(width, height);
		/*
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
		*/
		 
		// Display filtered images
//		(new ImagePlus("gabor, sigma="+sigma+" gamma="+gamma+ " psi="+psi, is)).show();
		 
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
		ImagePlus result = new ImagePlus ("GaborFiltered", filtered) ;
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