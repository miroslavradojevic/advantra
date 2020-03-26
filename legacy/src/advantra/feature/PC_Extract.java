package advantra.feature;

import java.util.Arrays;

import ij.ImagePlus;
import ij.gui.NewImage;
import ij.io.FileSaver;
import ij.process.ImageProcessor;
import imagescience.fourier.FFT;
import imagescience.image.Axes;
import imagescience.image.FloatImage;
import advantra.filter.FastMedian;
import advantra.filter.LogGabor;
import advantra.filter.Lowpass;
import advantra.filter.Template;
import advantra.filter.Rayleigh;

public class PC_Extract {
	
	//	boolean PC_caculated;
	private static final float SCALE_MEDIAN = (float) Math.sqrt(Math.log(4));
	private static final double PI = Math.PI;
	
	// ALL THE CALCULATION IS WITHIN THE CLASS CONSTRUCTOR
	public static float[][] run(ImageProcessor imp, String im_title, PC_Params params, boolean isCalledFromImageJ, boolean saveIt){
		
		int rows = imp.getHeight();
		int cols = imp.getWidth();
		
		// OUTPUT - EACH ORIENT WILL BE ONE IMAGE
//		FloatImage output 	= new FloatImage(NewImage.createFloatImage("template", 
//				cols, rows, 1, NewImage.FILL_BLACK));
		
		float[][] PC = new float[params.norient][rows*cols];
//		
//		for (int o = 0; o < params.norient; o++) {
//			PC[o] = (float[])output.imageplus().getProcessor().getPixels();
//		}
		
		float[] radius 		= Template.radius(rows, cols);
		float[] theta		= Template.theta(rows, cols);
		
		float epsilon = 1e-4f;

		float[] direction_spread = new float[rows*cols];
		float[] lp = Lowpass.createRow(radius, 0.45f, 15); 
		float[] lg = new float[rows*cols];
		
		//// FOURIER TF ////
		FloatImage real  		= new FloatImage(new ImagePlus("REAL FFT", imp));
		FloatImage imag 		= new FloatImage(NewImage.createFloatImage("IMAG FFT", cols, rows, 1, NewImage.FILL_BLACK));
		
		// IMPORTANT! 
		float[] real_array 	= (float[])real.imageplus().getProcessor().getPixels();
		float[] imag_array 	= (float[])imag.imageplus().getProcessor().getPixels();
		
		FFT fft = new FFT();
		
		fft.forward(real, imag, new Axes(true, true)); // real & imag will contain FFT coefficients
		// store FFT in float[] format
		float[] real_fft 	= new float[rows*cols];
		float[] imag_fft 	= new float[rows*cols];
		
		for (int i = 0; i < real_fft.length; i++) {
			real_fft[i] = real_array[i];
			imag_fft[i] = imag_array[i];
		}
		
		// allocate space to store the results
		float[][][] EO_real = new float[params.nscale][params.norient][rows*cols];
		float[][][] EO_imag = new float[params.nscale][params.norient][rows*cols];	
		
		// form weighting that penalizes frequency distributions that are particularly narrow
		float[] width 	= new float[rows*cols];
		float[] weight 	= new float[rows*cols];
		
		
		// for this orientation
		float[] sum_E 	= new float[rows*cols];
		float[] sum_O 	= new float[rows*cols];
		float[] sum_An 	= new float[rows*cols]; 
		float[] An 		= new float[rows*cols];
		float[] max_An	= new float[rows*cols];
		float[] Energy	= new float[rows*cols];
		
		// weighted mean filter response
		// weighted mean phase angle
		float[] Mean_E 	= new float[rows*cols];
		float[] Mean_O 	= new float[rows*cols];

		
		
		// MAIN LOOP
		for (int o = 0; o < params.norient; o++) {
			
			//long start_orientation = System.currentTimeMillis();
			
			double orient = o * Math.PI / params.norient;
			
			//// DIRECTION SPREAD ////
			// angular distance from the specified orientation
			for (int i = 0; i < direction_spread.length; i++) {
				// set the values for the direction
				double ds		= Math.sin(theta[i]) * Math.cos(orient) - Math.cos(theta[i]) * Math.sin(orient); 
				double dc 		= Math.cos(theta[i]) * Math.cos(orient) + Math.sin(theta[i]) * Math.sin(orient); 
				double dtheta 	= Math.abs(Math.atan2(ds, dc));
				dtheta 	= ( dtheta > PI/2   )? PI-dtheta  : dtheta ;// keep the mirror side of the filter
				dtheta 	=   dtheta * (double)params.norient / 2; 	// scale it 
				dtheta 	= ( dtheta > PI     )? PI         : dtheta ;// clamp to PI
				direction_spread[i] = (float) ((Math.cos(dtheta) + 1) / 2);
			}

			float tau = 0;
			
			Arrays.fill(sum_E, 0f);
			Arrays.fill(sum_O, 0f);
			Arrays.fill(sum_An, 0f);
			Arrays.fill(An, 0f);
			Arrays.fill(max_An, 0f);
			Arrays.fill(Energy, 0f);
			
			for (int s = 0; s < params.nscale; s++) { 
				
				// for each scale convolve the image with oriented filter returning the result in EO
				
				//// LOWPASS FILTER & LOG GABOR ////
				LogGabor.createRow(radius, params.minWvl, params.mult, params.sigmaOnf, s, lg);
				
				// FILTERING
				for (int i = 0; i < rows*cols; i++) {
					real_array[i] =  real_fft[i] * lp[i] * lg[i] * direction_spread[i];
					imag_array[i] =  imag_fft[i] * lp[i] * lg[i] * direction_spread[i];
				}
				
				//// INVERSE FOURIER TF ////
				fft.inverse(real, imag, new Axes(true, true)); // real_array and imag_array will still hold the values
					
				for (int i = 0; i < rows*cols; i++) {
					
					// sum them for each scale
					An[i]		= (float) Math.sqrt(real_array[i]*real_array[i] + imag_array[i]*imag_array[i]);
					sum_An[i] 	+= An[i];
					sum_E[i]	+= real_array[i];
					sum_O[i]	+= imag_array[i];
					
					if(An[i]>max_An[i]){
						
						max_An[i] = An[i];
					
					}
				}
				
				// store real_array and imag_array (IFFT result) into EO_real and EO_imag
						for (int i = 0; i < rows*cols; i++) {
							
							EO_real[s][o][i] = real_array[i];
							EO_imag[s][o][i] = imag_array[i];
						
						}
						
				//if the scale is smallest (higher frequencies) - estimate noise characteristics
				if(s==0){
					if(params.noiseMethod==-1){
						tau 		= (float) (FastMedian.median_Wirth(sum_An) / SCALE_MEDIAN);
					}
					else if (params.noiseMethod==-2){
						
						double[] sum_An_converted = new double[sum_An.length];
						
						for (int i = 0; i < sum_An_converted.length; i++) {
							
							sum_An_converted[i] = (double)sum_An[i];
							
						}
						
						tau = Rayleigh.mode(sum_An_converted, 50);
					
					}
				}
				
			} // 1 : nscales
			
			// System.out.println("elapsed time for orientation loop : "+(System.currentTimeMillis()-start_orientation));
			
			// long start_the_rest = System.currentTimeMillis(); 
			
			Arrays.fill(Mean_E, 0f);
			Arrays.fill(Mean_O, 0f);

			for (int i = 0; i < rows*cols; i++) {
				Mean_E[i] = (float) (sum_E[i] / ( Math.sqrt(sum_E[i]*sum_E[i]+sum_O[i]*sum_O[i]) + epsilon ));
				Mean_O[i] = (float) (sum_O[i] / ( Math.sqrt(sum_E[i]*sum_E[i]+sum_O[i]*sum_O[i]) + epsilon ));
			}
			
			
			// calculate energy
			for (int s = 0; s < params.nscale; s++) {
				
				for (int i = 0; i < rows*cols; i++) {
					
					Energy[i] += EO_real[s][o][i] * Mean_E[i] + EO_imag[s][o][i] * Mean_O[i] - Math.abs(EO_real[s][o][i]*Mean_O[i] - EO_imag[s][o][i]*Mean_E[i]);
				
				}
				
			}
			
			// calculate noise threshold
			float T;
			if(params.noiseMethod>=0){
				
				T = (float) params.noiseMethod;
				
			}else {
				
				T = (float) (tau * ((1-Math.pow(1/params.mult, params.nscale)) / (1-1/params.mult)) * (Math.sqrt(Math.PI / 2) + params.k * Math.sqrt((4-Math.PI)/2)));  
			
			}
			
			// apply noise threshold
			for (int i = 0; i < rows*cols; i++) {
				
				Energy[i] = (Energy[i]-T < 0)? 0 : (Energy[i]-T) ;
			
			}
			
			for (int i = 0; i < rows*cols; i++) {
				
				width[i] 	= (sum_An[i]/(max_An[i]+epsilon)-1)/(params.nscale-1);
				// calculate weight for this orientation
				weight[i] 	= (float) (1.0 / ( 1.0 + Math.exp( (params.cutOff - width[i])*params.g ) ));
				
			}
			
			//weight the energy and calculate the phase congruency
			for (int i = 0; i < rows*cols; i++) {
				PC[o][i] = weight[i]*Energy[i]/sum_An[i];
			}
			
			// PC[orientation] calculated !
			
		} // 1 : norients
		
		if(isCalledFromImageJ || saveIt){
			
			// ImagePlus objects are needed, as many as there are orientations
			ImagePlus[] images = new ImagePlus[params.norient];
			
			for (int i = 0; i < params.norient; i++) {
				
				String image_name = im_title+"_pc_"+Double.toString(i*180/params.norient);
				
				images[i] = NewImage.createFloatImage (
						image_name, cols, rows, 1, NewImage.FILL_BLACK);
				
				images[i].getProcessor().setPixels(PC[i]);
			
			}
			
			if(isCalledFromImageJ){
				
				// show them!
				for (int i = 0; i < params.norient; i++) {
					
					images[i].show();
					
				}
				
				
			}
			
			if(saveIt){
				
				// save it to HD... 
					
				boolean savedTiff; 
				//FileSaver[] file_saver 	= new FileSaver[params.norient];//(output.imageplus());
				
				for (int i = 0; i < params.norient; i++) {
					
					String output_title		= im_title+"_pc_"+Integer.toString(i)+".png";
					savedTiff = (new FileSaver(images[i])).saveAsPng(output_title); 
					
					if(!savedTiff){
						
			  			System.out.println("Couldn't save the png.");
			  		
					}
			  		else{
			  			
			  			System.out.println(output_title+" exported.");
			  		
			  		}
				}
			}
			
		}

		return PC;
		
	}//PC_Extract constructor

}