package advantra.plugins;

import java.io.File;

import advantra.feature.DifferentialStructure;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

import imagescience.image.Image;

public class ExtractDifferentialStructure implements PlugInFilter {

	ImagePlus img;
	
	double 	sigma_1 = 1.0;
	double 	sigma_2 = 1.0;
	int		nr 		= 1;
	
	boolean		gradient_show, 
				laplacian_show, 
				ridge_det_show, 
				isophote_curv_show, 
				flowline_curv_show, 
				isophote_density_show,
				corner_det_show,
				shape_index_show,
				curvedness_show,
				hessian_det_show,
				mean_curvature_show,
				gaussian_extremality_show,
				t_junction_likeliness_show;
	
	Image Lx, 	Ly, 	Lxx, 	Lyy, 	Lxy; 
	Image Lxyy, Lxxy, 	Lxxx, 	Lyyy;
	
	/*
	 * ImageJ methods
	 * 
	 */
	
	public int setup(String arg0, ImagePlus im) {
		this.img = im;
		return DOES_8G+NO_CHANGES;
	}
	
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
		
		gd.addNumericField( "sigma start:", sigma_1, 0, 5, "pix" );
		gd.addNumericField( "sigma end  :", sigma_2, 0, 5, "pix" );	
		gd.addNumericField( "number of scales : ", nr,  0, 5, "");
		
		gd.addMessage("Show descriptors:");
		
		gd.addCheckbox("gradient", 				false);//true); // 1
		gd.addCheckbox("laplacian", 			false);//true);
		gd.addCheckbox("ridge detection", 		false);//true);
		gd.addCheckbox("isophote curvature",	false);//true);
		gd.addCheckbox("flowline curvature",	false);//true); // 5
		gd.addCheckbox("isophote density", 		false);//true);
		gd.addCheckbox("corner detector", 		false);//true);
		gd.addCheckbox("shape index", 			false);//true);
		gd.addCheckbox("curvedness", 			false);//true);
		gd.addCheckbox("hessian determinant",	false);//true); // 10
		gd.addCheckbox("mean curvature",		false);//true);
		gd.addCheckbox("gaussian extremality",	false);//true);
		gd.addCheckbox("t junction likeliness",	false);//true);
		
		gd.addMessage("* * *");
		
		gd.showDialog();
		if (gd.wasCanceled()) return;
		
		sigma_1 		= (double)gd.getNextNumber();
		sigma_2 		= (double)gd.getNextNumber();
		nr 				= (int)gd.getNextNumber();
		
		gradient_show 			= (boolean)gd.getNextBoolean();
		laplacian_show 			= (boolean)gd.getNextBoolean();
		ridge_det_show 			= (boolean)gd.getNextBoolean();
		isophote_curv_show 		= (boolean)gd.getNextBoolean();
		flowline_curv_show 		= (boolean)gd.getNextBoolean();
		isophote_density_show 	= (boolean)gd.getNextBoolean();
		corner_det_show			= (boolean)gd.getNextBoolean();
		shape_index_show		= (boolean)gd.getNextBoolean();
		curvedness_show			= (boolean)gd.getNextBoolean();
		hessian_det_show		= (boolean)gd.getNextBoolean();
		mean_curvature_show		= (boolean)gd.getNextBoolean();
		gaussian_extremality_show 	= (boolean)gd.getNextBoolean();
		t_junction_likeliness_show 	= (boolean)gd.getNextBoolean();
		
		DifferentialStructure d_struct = new DifferentialStructure(img, sigma_1, sigma_2, nr);
		
		if (gradient_show)		d_struct.getGradientMagnitude().show();
		if (laplacian_show) 	d_struct.getLaplacian().show();
		if (ridge_det_show) 	d_struct.getRidgeDet().show();
		if(isophote_curv_show) 	d_struct.getIsophoteCurvature().show();
		if(flowline_curv_show)	d_struct.getFlowlineCurv().show();
		if(isophote_density_show)d_struct.getIsophoteDensity().show();
		if(corner_det_show)		d_struct.getCornerDetector().show();
		if(shape_index_show)	d_struct.getShapeIndex().show();
		if(curvedness_show)		d_struct.getCurvedness().show();
		if(hessian_det_show)	d_struct.getHessianDeterminant().show();
		if(mean_curvature_show) d_struct.getMeanCurvature().show();
		if(gaussian_extremality_show)d_struct.getGaussianExtremality().show();
		if(t_junction_likeliness_show)d_struct.getTJunctionLikeliness().show();
		
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