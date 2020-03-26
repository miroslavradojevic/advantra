package advantra.feature;

import java.util.Vector;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.ZProjector;
import ij.plugin.filter.Convolver;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

public class GaborFilt2D extends Thread { 
	
	/*
	 * this implementation was chosen since parallelizing per angles was faster and 
	 * there was no loss in information
	 */
	
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
	public static int Ht;
	public static int Wt;
	
	// gabor params
	public static double[]	thetas;
	public static double[]	sigmas;
	public static double[] 	lambdas;
	public static double 	bandwidth;
	public static double 	psi;
	public static double	gamma;
	public static boolean	real;
	
	/*
	 * thetas will be processed in parallel
	 */
	
	public static ImageStack	gabor_directional_responses; // will be stack with lambdas layers		
	public static ImageProcessor input_ip;
	public static Vector<ImageStack> kernels;
	
	private int n0, n1;
	
	public GaborFilt2D(int n0, int n1){
		this.n0 = n0;
		this.n1 = n1;
	}

	public static void load(
			ImagePlus 	input_img99, 
			double[] 	thetas99, 
			double[] 	sigmas99, 
			double[] 	lambdas99, 
			double		bandwidth99,
			double 		psi99,
			double		gamma99,
			boolean		real99){
		
		thetas = new double[thetas99.length];
		for (int i = 0; i < thetas99.length; i++) {
			thetas[i] = thetas99[i];
		}
		
		sigmas = new double[sigmas99.length];
		for (int i = 0; i < sigmas99.length; i++) {
			sigmas[i] = sigmas99[i];
		}
		
		lambdas = new double[lambdas99.length];
		for (int i = 0; i < lambdas99.length; i++) {
			lambdas[i] = lambdas99[i];
		}
		
		bandwidth = bandwidth99;
		psi = psi99;
		gamma = gamma99;
		real = real99;
		
		Ht = input_img99.getHeight();
		Wt = input_img99.getWidth();
		
		input_ip = input_img99.getProcessor().convertToFloat().duplicate();
		
		gabor_directional_responses = new ImageStack(Wt, Ht, thetas.length);
		
		// create filters
		kernels = new Vector<ImageStack>(sigmas.length);

		double slratio = 
				(1/Math.PI) * 
				Math.sqrt((Math.log(2)/2)) * 
				((Math.pow(2,bandwidth)+1)/(Math.pow(2,bandwidth)-1));
		
		for (int scaleIdx = 0; scaleIdx < sigmas.length; scaleIdx++) {
			
			if(sigmas[scaleIdx]==0){
				sigmas[scaleIdx] = slratio * lambdas[scaleIdx];
			}
			else if(lambdas[scaleIdx]==0){
				lambdas[scaleIdx] = sigmas[scaleIdx] / slratio;
			}
			
			double sigma_x = sigmas[scaleIdx];
			double sigma_y = sigmas[scaleIdx] / gamma;
			double larger_sigma = (sigma_x>sigma_y)? sigma_x : sigma_y ;
			
			int filterSizeX = W * (int)larger_sigma + 1;
			int filterSizeY = W * (int)larger_sigma + 1;
			
			int middleX = (int) Math.round(filterSizeX / 2);
			int middleY = (int) Math.round(filterSizeY / 2);
			
			ImageStack kernel_stack = new ImageStack(filterSizeX, filterSizeY);
			
			for (int angIdx = 0; angIdx < thetas.length; angIdx++) {
				
				double curr_angle = thetas[angIdx]; 
				double lambda = lambdas[scaleIdx];
				
				ImageProcessor filter = new FloatProcessor(filterSizeX, filterSizeY);  
				
				float sumPos = 0;
				float sumNeg = 0;
				
				// x0 and y0 at patch center (0,0)
			    for (int x=-middleX; x<=middleX; x++){
			        for (int y=-middleY; y<=middleY; y++){           
				        	
			        	double  xr = (double)x * Math.cos(curr_angle) + (double)y * Math.sin(curr_angle);
			        	double  yr = (double)y * Math.cos(curr_angle) - (double)x * Math.sin(curr_angle);
				        	
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
				
			    kernel_stack.addSlice("", filter);
				
			} // loop thetas
			
			kernels.add(kernel_stack);
			
		}
		
//		for (int i = 0; i < kernels.size(); i++) {
//			new ImagePlus("", kernels.get(i)).show();
//		}
		
	}
	
	public void 	run(){
		
		// run one part of the whole job, one range of patches that make the whole image
		ZProjector zmax1 = new ZProjector();
		for (int i = n0; i < n1; i++) {
			
			ImagePlus g_theta = runTheta(i);

			zmax1.setImage(g_theta);
			zmax1.setStartSlice(1);
			zmax1.setStopSlice(g_theta.getStackSize());
			zmax1.setMethod(ZProjector.MAX_METHOD);
			zmax1.doProjection();
			gabor_directional_responses.setPixels(
					zmax1.getProjection().getChannelProcessor().getPixels(), 
					i+1);

		}
		
	}	
	
	public static ImagePlus runTheta(
			int theta_idx
			){
		
		ImageStack result_stack = new ImageStack(Wt, Ht);
		
		for (int i = 0; i < sigmas.length; i++) {
			
			result_stack.addSlice(
					"gabor,theta=" + IJ.d2s(thetas[theta_idx], 2) + ",sigma=" + IJ.d2s(sigmas[i], 2),
					runOne(theta_idx, i)
					);
			
		}
		
//		String result_name = ((real)?"real":"imag")+"_gabor_filt_per_scales";
		
		return new ImagePlus ("", result_stack);
		
	}
	
	public static ImageProcessor runOne(
			int theta_idx,
			int sigma_idx 
			){

		//convolution
	    Convolver c = new Convolver();                
		c.setNormalize(false); // important not to normalize (did my own normalization)
		
		float[] kernel = (float[]) kernels.get(sigma_idx).getProcessor(theta_idx+1).getPixels();
		
		int filterSizeX = kernels.get(sigma_idx).getWidth();
		int filterSizeY = kernels.get(sigma_idx).getHeight();
		
		ImageProcessor ip = input_ip.duplicate(); // needs to be float stack
		c.convolveFloat(ip, kernel, filterSizeX, filterSizeY);      
		return ip;
		
	}
	
	public static ImagePlus extractDirectionalGabor(
			ImagePlus input, 
			double[] angles_pi, 
			double[] t, 
			double bw, 
			double psi, 
			double gamma,
			boolean isReal){
		
		int w = input.getWidth();
		int h = input.getHeight();
		
		ImageStack 	angul = new ImageStack(w, h); 
		
		ZProjector zmax = new ZProjector();
		
		for (int i = 0; i < angles_pi.length; i++) {
			
			double current_theta = angles_pi[i];
			
			ImagePlus g_theta = runTheta(i);
			
			zmax.setImage(g_theta);
			zmax.setStartSlice(1);
			zmax.setStopSlice(g_theta.getStackSize());
			zmax.setMethod(ZProjector.MAX_METHOD);
			zmax.doProjection();
			angul.addSlice("theta="+IJ.d2s(current_theta, 2), zmax.getProjection().getChannelProcessor());
			
		}
		
		return new ImagePlus("gabor_per_angle_reg", angul);

	}
	
}