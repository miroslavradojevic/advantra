package advantra.feature;

import java.util.Vector;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.ZProjector;
import ij.plugin.filter.Convolver;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import imagescience.image.Axes;
import imagescience.image.Coordinates;
import imagescience.image.Dimensions;
import imagescience.image.FloatImage;
import imagescience.image.Image;

public class GaborFilt2D extends Thread { // 
	
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
	
	// parallel
	public static Vector<ImagePlus> ptches;
	public static int ptch_len;
	public static Vector<ImagePlus> ptches_out;
	public static Vector<Integer> ptch_root_x;
	public static Vector<Integer> ptch_root_y;
	
	private int n0, n1;
	
	public GaborFilt2D(int n0, int n1){
		this.n0 = n0;
		this.n1 = n1;
	}
	
	public void 	run(){
		
		// run one part of the whole job, one range of patches that make the whole image
		
//		System.out.println("calculating range: "+n0+" to "+(n1-1));
		
		for (int i = n0; i < n1; i++) {
			
			ImagePlus to_add = extractDirectionalGabor(ptches.get(i), thetas, sigmas, bandwidth, psi, gamma, real);
			//IJ.log("extracted "+to_add.getHeight()+" x "+to_add.getWidth()+" x "+to_add.getStackSize());
			ptches_out.set(i, to_add);
			//IJ.log("current size: "+ptches_out.size()+" ");
			//IJ.log("added  to "+i);
			
		}
		
	}	
	
	public static void load(
			ImagePlus 	input_img, 
			int 		number_of_patches, 
			double[] 	thetas99, 
			double[] 	sigmas99, 
			double[] 	lambdas99, 
			double		bandwidth99,
			double 		psi99,
			double		gamma99,
			boolean		real99){
		
		
		System.out.println("how much "+thetas99.length);
		
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
		
		Ht = input_img.getHeight();
		Wt = input_img.getWidth();
		
		ptch_len = (int)Math.ceil(Math.sqrt((Ht*Wt)/number_of_patches));
		
		ptches 		= new Vector<ImagePlus>();
		ptches_out 	= new Vector<ImagePlus>();
		ptch_root_x = new Vector<Integer>();
		ptch_root_y = new Vector<Integer>();
		
		// split into patches ptch_len x ptch_len
		// easier wrapped in imagescience
		
		double[][] aPtch = new double[ptch_len][ptch_len];
		
		Image input_image = Image.wrap(input_img);
		input_image.axes(Axes.X+Axes.Y);
		
		Image output_patch = new FloatImage(new Dimensions(ptch_len, ptch_len));
		output_patch.axes(Axes.X+Axes.Y);
		
		Coordinates crd = new Coordinates();
		for (int row = 0; row < Ht; row+=ptch_len) {
			for (int col = 0; col < Wt; col+=ptch_len) {
				
				if(row+ptch_len>Ht){
					aPtch = new double[ptch_len][ptch_len];
				}
				
				if(col+ptch_len>Wt){
					aPtch = new double[ptch_len][ptch_len];
				}
				
				crd.x = col;//Wt-ptch_len+1;
				crd.y = row;//Ht-ptch_len+1;
				input_image.get(crd, aPtch);
				
				crd.x = 0;
				crd.y = 0;
				output_patch.set(crd, aPtch);
				
//				output_patch.imageplus().show();
				
				// add it to Vector
				ptches.add(output_patch.duplicate().imageplus());
				ptches_out.add(new ImagePlus("", new FloatProcessor(ptch_len, ptch_len)));
				ptch_root_x.add(col);
				ptch_root_y.add(row);
				
			}
		}
		
		
		
		// allocate the outputs so that the future 
		//ptches_out.set(0, ptches.get(0));
		
		System.out.println("formed "+ptches.size()+" input patches and allocated "+ptches_out.size()+" output patches");
		
	}

	public static ImagePlus concatenateOutput(){
		
		Image			out_image = new FloatImage(new Dimensions(Wt, Ht, thetas.length));//Image.wrap(out_im);
		out_image.axes(Axes.X+Axes.Y+Axes.Z);
		
		double[][][] aPatch = new double[thetas.length][ptch_len][ptch_len];
		
		Coordinates crd = new Coordinates();

		for (int i = 0; i < ptches_out.size(); i++) {  
			
			Image this_ptch = Image.wrap(ptches_out.get(i));
			this_ptch.axes(Axes.X+Axes.Y+Axes.Z);

			crd.x = 0;
			crd.y = 0;
			
			this_ptch.get(crd, aPatch);
			
			crd.x = ptch_root_x.get(i);
			crd.y = ptch_root_y.get(i);
			
			out_image.set(crd, aPatch);
			
		}
		
		return out_image.imageplus();
		
	}
	
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

		// if input2D is null just give kernels (for visualization)
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
			
			//System.out.println("processing theta = "+i+" / "+ (angles_pi.length-1));
			
			ImagePlus g_theta = run(
					input, current_theta, t, new double[t.length], bw, psi, gamma, isReal);
			
			zmax.setImage(g_theta);
			zmax.setStartSlice(1);
			zmax.setStopSlice(g_theta.getStackSize());
			zmax.setMethod(ZProjector.MAX_METHOD);
			zmax.doProjection();
			angul.addSlice("theta="+IJ.d2s(current_theta, 2), zmax.getProjection().getChannelProcessor());

		}
		
		return new ImagePlus("gabor_per_angle", angul);

	}
	
}