package advantra.plugins;

import java.io.File;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import imagescience.feature.Differentiator;
import imagescience.image.Axes;
import imagescience.image.Coordinates;
import imagescience.image.Dimensions;
import imagescience.image.FloatImage;
import imagescience.image.Image;

public class ExtractBifurcationFeatures implements PlugInFilter {

	ImagePlus img;
	
	double 	s_1 	= 1.0;
	double 	s_2 	= 1.0;
	int		nr 		= 1;
	
	Image Lx, 	Ly, 	Lxx, 	Lyy, 	Lxy; 
	Image Lxyy, Lxxy, 	Lxxx, 	Lyyy;
	
	/*
	 * ImageJ methods
	 * 
	 */
	
	public void run(ImageProcessor arg0) {
		
		// reset calibration before going further
		Calibration cal = new Calibration(img);
		cal.pixelWidth = cal.pixelHeight = cal.pixelDepth = 1.0;
		cal.setUnit("pixel");
		img.setCalibration(cal);
		
		// constrain on 2d images
		if(img.getStackSize()>1){
			IJ.error("cannot process image stacks");
			return;
		}
		
		// dialog to enter input values
		GenericDialog gd = new GenericDialog("Differential Structure", IJ.getInstance());

		gd.addMessage("Choose scales:");
		
		gd.addNumericField( "s start:", s_1, 0, 5, "pix" );
		gd.addNumericField( "s end  :", s_2, 0, 5, "pix" );	
		gd.addNumericField( "number of s scales : ", nr,  0, 5, "");
		
		gd.addMessage("Show descriptors:");
		
		gd.addCheckbox("gradient", 				false);//true); // 1
		gd.addCheckbox("laplacian", 			false);//true);
		gd.addCheckbox("ridge detection", 		false);//true);
		gd.addCheckbox("isophote curvature",	false);//true);
		gd.addCheckbox("flowline curvature",	false);//true);
		gd.addCheckbox("isophote density", 		false);//true);
		gd.addCheckbox("corner detector", 		false);//true);
		gd.addCheckbox("shape index", 			false);//true);
		gd.addCheckbox("curvedness", 			false);//true);
		gd.addCheckbox("hessian determinant",	false);//true); 
		gd.addCheckbox("mean curvature",		false);//true);
		gd.addCheckbox("gaussian extremality",	false);//true);
		gd.addCheckbox("t junction likeliness",	false);//true);
		
		gd.addMessage("* * *");
		
		gd.showDialog();
		if (gd.wasCanceled()) return;
		
		s_1 			= (double)gd.getNextNumber();
		s_2 			= (double)gd.getNextNumber();
		nr 				= (int)gd.getNextNumber();

			
	}

	public int setup(String arg0, ImagePlus im) {
		this.img = im;
		return DOES_8G+NO_CHANGES;
	}
	
	/*
	 * main() for terminal use
	 * java -cp JAR_LIBRARIES:advantra_.jar advantra.plugins.ExtractDifferentialStructure ARG1 ARG2 ARG3 ARG4
	 */
	
	public static void main(String[] args){
		
		double 	sigma_1 = 1;
		double 	sigma_2 = 1;
		int 	nr 		= 1;
		String image_path = "";
		
		if(args.length!=4){
			System.out.println(
					"/-----------------------------------------------------------------\n" +
					"Plugin extracts differential structure images of the input image  \n" +
					"/-----------------------------------------------------------------\n" +
					"usage: SIGMA1 SIGMA2 NR_SCALES IMAGE_PATH");
			return;
		}
		else{
			sigma_1 = Double.parseDouble(args[0]);
			sigma_2 = Double.parseDouble(args[1]);
			nr 		= Integer.parseInt(args[2]);
			image_path = args[3];
		}
		
		image_path = (new File(image_path)).getAbsolutePath();
		
		ImagePlus img; 
		if((new File(image_path)).exists()){
			img = new ImagePlus(image_path);
		}
		else{
			System.out.println("File does not exist!");
			return;
		}
		
		DifferentialStructure d_struct = new DifferentialStructure(img, sigma_1, sigma_2, nr);
		
		// saving outputs to disk...
		String export_path = "";
		System.out.println("exporting results... ");
		
		export_path = img.getTitle()+"_gradient.tif";
		System.out.println(((new File(export_path)).getAbsolutePath()));
		if(nr>1) 	(new FileSaver(d_struct.getGradientMagnitude())).saveAsTiffStack(export_path);
		else 		(new FileSaver(d_struct.getGradientMagnitude())).saveAsTiff(export_path);
		
		
		export_path = img.getTitle()+"_laplacian.tif";
		System.out.println(((new File(export_path)).getAbsolutePath()));
		if(nr>1) 	(new FileSaver(d_struct.getLaplacian())).saveAsTiffStack(export_path);
		else 		(new FileSaver(d_struct.getLaplacian())).saveAsTiff(export_path);
	
		export_path = img.getTitle()+"_ridge_detector.tif";
		System.out.println(((new File(export_path)).getAbsolutePath()));
		if(nr>1) 	(new FileSaver(d_struct.getRidgeDet())).saveAsTiffStack(export_path);
		else 		(new FileSaver(d_struct.getRidgeDet())).saveAsTiff(export_path);
		
		export_path = img.getTitle()+"_isophote_curv.tif";
		System.out.println(((new File(export_path)).getAbsolutePath()));
		if(nr>1) 	(new FileSaver(d_struct.getIsophoteCurvature())).saveAsTiffStack(export_path);
		else 		(new FileSaver(d_struct.getIsophoteCurvature())).saveAsTiff(export_path);
		
		export_path = img.getTitle()+"_flowline_curv.tif";
		System.out.println(((new File(export_path)).getAbsolutePath()));
		if(nr>1) 	(new FileSaver(d_struct.getFlowlineCurv())).saveAsTiffStack(export_path);
		else 		(new FileSaver(d_struct.getFlowlineCurv())).saveAsTiff(export_path);
		
		export_path = img.getTitle()+"_isophote_den.tif";
		System.out.println(((new File(export_path)).getAbsolutePath()));
		if(nr>1) 	(new FileSaver(d_struct.getIsophoteDensity())).saveAsTiffStack(export_path);
		else 		(new FileSaver(d_struct.getIsophoteDensity())).saveAsTiff(export_path);
		
		export_path = img.getTitle()+"_corner_det.tif";
		System.out.println(((new File(export_path)).getAbsolutePath()));
		if(nr>1) 	(new FileSaver(d_struct.getCornerDetector())).saveAsTiffStack(export_path);
		else 		(new FileSaver(d_struct.getCornerDetector())).saveAsTiff(export_path);	
		
		export_path = img.getTitle()+"_shape_index.tif";
		System.out.println(((new File(export_path)).getAbsolutePath()));
		if(nr>1) 	(new FileSaver(d_struct.getShapeIndex())).saveAsTiffStack(export_path);
		else 		(new FileSaver(d_struct.getShapeIndex())).saveAsTiff(export_path);	
		
		export_path = img.getTitle()+"_curvedness.tif";
		System.out.println(((new File(export_path)).getAbsolutePath()));
		if(nr>1) 	(new FileSaver(d_struct.getCurvedness())).saveAsTiffStack(export_path);
		else 		(new FileSaver(d_struct.getCurvedness())).saveAsTiff(export_path);
		
		export_path = img.getTitle()+"_hessian_determinant.tif";
		System.out.println(((new File(export_path)).getAbsolutePath()));
		if(nr>1) 	(new FileSaver(d_struct.getHessianDeterminant())).saveAsTiffStack(export_path);
		else 		(new FileSaver(d_struct.getHessianDeterminant())).saveAsTiff(export_path);
		
		export_path = img.getTitle()+"_mean_curvature.tif";
		System.out.println(((new File(export_path)).getAbsolutePath()));
		if(nr>1) 	(new FileSaver(d_struct.getMeanCurvature())).saveAsTiffStack(export_path);
		else 		(new FileSaver(d_struct.getMeanCurvature())).saveAsTiff(export_path);
		
		export_path = img.getTitle()+"_gaussian_extremality.tif";
		System.out.println(((new File(export_path)).getAbsolutePath()));
		if(nr>1) 	(new FileSaver(d_struct.getGaussianExtremality())).saveAsTiffStack(export_path);
		else 		(new FileSaver(d_struct.getGaussianExtremality())).saveAsTiff(export_path);
		
		export_path = img.getTitle()+"_t_junction_likeliness.tif";
		System.out.println(((new File(export_path)).getAbsolutePath()));
		if(nr>1) 	(new FileSaver(d_struct.getTJunctionLikeliness())).saveAsTiffStack(export_path);
		else 		(new FileSaver(d_struct.getTJunctionLikeliness())).saveAsTiff(export_path);
		
	} // main()
	
}

/*
 * classes used by both terminal call and ImageJ
 */

class ScaleSpace{
	
	/*
	 * from Msc thesis "Vessels & Bifurcations in Scale Space"
	 * and 				"Vascular Bifurcation Detection in Scale-Space"
	 * by Daniel-Marian Baboiu and Ghassan Hamarneh
	 */
	ScaleSpace(ImagePlus img, double s1, double s2, double nr){
		
	}
	
}