package advantra.filter;

import java.util.ArrayList;
import java.util.Vector;

import advantra.general.ImageConversions;

import ij.ImagePlus;
import ij.gui.NewImage;
import imagescience.feature.Differentiator;
import imagescience.feature.Hessian;
import imagescience.image.FloatImage;
import imagescience.image.Image;

public class Frangi {
	
	/* 
	 * calculates frangi's features (vesselness, eig.vals) 
	 * for an image stack (3d)
	 */
	
	private ImagePlus 	img;
	private double[]	sigma;
	
	private double 		alpha;
	private double		beta;				 
	private double		C;
	
	// 			[nr_scales][rows][cols][size]
	private float[][][][]	vess;		// vesselness per scale
	private float[][][][]	lambda_1;	// |lambda_1|<=|lambda_2|<=|lambda_3|
	private float[][][][]	lambda_2;
	private float[][][][]	lambda_3;
	
	public Frangi(
			ImagePlus 	img, 
			double 		start_s, 
			double 		step_s, 
			double 		end_s
			){
		
		if(img.getType()!=ImagePlus.GRAY8 || img.getStack().getSize()==1){
			return;
		}
		
		this.img 		= img;
		
		int nr_scales = 0;
		for (double sigma_i = start_s; sigma_i <= end_s; sigma_i+=step_s) {
			nr_scales++;
		}
		this.sigma      = new double[nr_scales];
		int idx = 0;
		for (double sigma_i = start_s; sigma_i <= end_s; sigma_i+=step_s) {
			this.sigma[idx] = sigma_i;
			idx++;
		}
		
		// default values for parameters
		this.alpha 	= 0.5;
		this.beta 	= 0.5;
		this.C 		= 2.0;
		
		this.vess 		= null;
		this.lambda_1	= null;
		this.lambda_2	= null;
		this.lambda_3	= null;
		
	}
	
	public void run(){
		
		vess   		= new float[sigma.length][img.getHeight()][img.getWidth()][img.getStack().getSize()];
		lambda_1	= new float[sigma.length][img.getHeight()][img.getWidth()][img.getStack().getSize()];
		lambda_2	= new float[sigma.length][img.getHeight()][img.getWidth()][img.getStack().getSize()];
		lambda_3	= new float[sigma.length][img.getHeight()][img.getWidth()][img.getStack().getSize()];
		
		// go through all the scales & extract the values
		Hessian 		hess 			= new Hessian();
		Vector<Image>	hessian_imgs 	= new Vector<Image>();
		
		for (int i = 0; i < sigma.length; i++) {
			System.out.println("\n---\nhessian at scale "+sigma[i]+"  :");
			long t1 = System.currentTimeMillis();
			hessian_imgs = hess.run(Image.wrap(ImageConversions.ImagePlusToFloat(img)), sigma[i], true);
			System.out.println("elapsed "+(System.currentTimeMillis()-t1)/1000f+" s, "+hessian_imgs.size()+" images");
			
			// save them as images in form: hess_eig_X_sigma_X_.tif
			
		}
		
	}
	
}
