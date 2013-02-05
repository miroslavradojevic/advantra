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
	
	public void run(ImageProcessor arg0) {
		
		// reset calibration before going further
		Calibration cal = new Calibration(img);
		cal.pixelWidth = 1.0;
		cal.pixelHeight = 1.0;
		cal.pixelDepth = 1.0;
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

//class used by both terminal call and ImageJ
class DifferentialStructure {
	
	Image image; // input
	
	Image Lx, 	Ly;
	Image Lxx, 	Lyy, 	Lxy; 
	Image Lxyy, Lxxy, 	Lxxx, 	Lyyy;
	
	Image gradient_image, 			laplacian_image, 			ridge_det_image, 			isophote_curv_image, 	flowline_curv_image;
	Image isophote_density_image, 	corner_det_image, 			shape_index_image, 			curvedness_image, 		hessian_det_image;
	Image mean_curvature_image, 	gaussian_extremality_image, t_junction_likeliness_image;
	
	double[] sc;
	
	DifferentialStructure(ImagePlus img, double sigma_1, double sigma_2, int nr){
		
		sc = new double[nr];
		for (int i = 0; i < nr; i++) {
			sc[i] = (i==0)? sigma_1 : sigma_1+i*((sigma_2-sigma_1)/(nr-1));
		}
		
		image = new FloatImage(Image.wrap(img));
		final Dimensions dims 			= image.dimensions();
		final Dimensions filter_dims	= new Dimensions(dims.x, dims.y, nr); // as many layers as there are scales
		final Differentiator df 		= new Differentiator();
		
		// calculate derivatives
		Lx 		= new FloatImage(dims); 
		Ly 		= new FloatImage(dims); 
		Lxx 	= new FloatImage(dims); 
		Lyy		= new FloatImage(dims); 
		Lxy 	= new FloatImage(dims);
		Lxyy	= new FloatImage(dims);
		Lxxy	= new FloatImage(dims);
		Lxxx	= new FloatImage(dims);
		Lyyy	= new FloatImage(dims);
		
		gradient_image				= new FloatImage(filter_dims); gradient_image.axes(Axes.X);
		laplacian_image				= new FloatImage(filter_dims); laplacian_image.axes(Axes.X);
		ridge_det_image				= new FloatImage(filter_dims); ridge_det_image.axes(Axes.X);
		isophote_curv_image			= new FloatImage(filter_dims); isophote_curv_image.axes(Axes.X);
		flowline_curv_image			= new FloatImage(filter_dims); flowline_curv_image.axes(Axes.X);
		isophote_density_image		= new FloatImage(filter_dims); isophote_density_image.axes(Axes.X);
		corner_det_image			= new FloatImage(filter_dims); corner_det_image.axes(Axes.X);
		shape_index_image			= new FloatImage(filter_dims); shape_index_image.axes(Axes.X);
		curvedness_image			= new FloatImage(filter_dims); curvedness_image.axes(Axes.X);
		hessian_det_image			= new FloatImage(filter_dims); hessian_det_image.axes(Axes.X);
		mean_curvature_image		= new FloatImage(filter_dims); mean_curvature_image.axes(Axes.X);
		gaussian_extremality_image	= new FloatImage(filter_dims); gaussian_extremality_image.axes(Axes.X);
		t_junction_likeliness_image	= new FloatImage(filter_dims); t_junction_likeliness_image.axes(Axes.X);
		
		// arrays 
		double[] aLx 	= new double[dims.x];
		double[] aLy 	= new double[dims.x];
		double[] aLxx 	= new double[dims.x];
		double[] aLyy 	= new double[dims.x];
		double[] aLxy 	= new double[dims.x];
		double[] aLxyy 	= new double[dims.x];
		double[] aLxxy 	= new double[dims.x];
		double[] aLxxx 	= new double[dims.x];
		double[] aLyyy 	= new double[dims.x];
					
		double[] aGradient 	= new double[dims.x];
		double[] aLaplacian = new double[dims.x];
		double[] aRidgeDet 	= new double[dims.x];
		double[] aIsophoteCurv = new double[dims.x];
		double[] aFlowlineCurv = new double[dims.x];
		double[] aIsophoteDensity = new double[dims.x];
		double[] aCornerDet = new double[dims.x];
		double[] aShapeIndex = new double[dims.x];
		double[] aCurvedness = new double[dims.x];
		double[] aHessianDet = new double[dims.x];
		double[] aMeanCurvature = new double[dims.x];
		double[] aGaussianExtremality = new double[dims.x];
		double[] aTJunctionLikeliness = new double[dims.x];
		
		for (int i = 0; i < nr; i++) { // run through scales
			
			System.out.println("processing scale "+sc[i]+" ...");
			
			Lx 		= df.run(image.duplicate(), sc[i], 1, 0, 0); Lx.axes(Axes.X);
			Ly	 	= df.run(image.duplicate(), sc[i], 0, 1, 0); Ly.axes(Axes.X);
			Lxx	 	= df.run(image.duplicate(), sc[i], 2, 0, 0); Lxx.axes(Axes.X);
			Lyy 	= df.run(image.duplicate(), sc[i], 0, 2, 0); Lyy.axes(Axes.X);
			Lxy 	= df.run(image.duplicate(), sc[i], 1, 1, 0); Lxy.axes(Axes.X);
			Lxyy 	= df.run(image.duplicate(), sc[i], 1, 2, 0); Lxyy.axes(Axes.X);
			Lxxy 	= df.run(image.duplicate(), sc[i], 2, 1, 0); Lxxy.axes(Axes.X);
			Lxxx 	= df.run(image.duplicate(), sc[i], 3, 0, 0); Lxxx.axes(Axes.X);
			Lyyy 	= df.run(image.duplicate(), sc[i], 0, 3, 0); Lyyy.axes(Axes.X);
			
			Coordinates coords 	= new Coordinates();
			for (coords.y=0; coords.y<dims.y; ++coords.y) {
				
				Lx.get(coords,aLx); 
				Ly.get(coords,aLy);
				Lxx.get(coords,aLxx);
				Lyy.get(coords,aLyy);
				Lxy.get(coords,aLxy);
				Lxyy.get(coords,aLxyy);
				Lxxy.get(coords,aLxxy);
				Lxxx.get(coords,aLxxx);
				Lyyy.get(coords,aLyyy);
				
				// operations with arrays
				for (int x=0; x<dims.x; ++x){
					aGradient[x] 	= Math.sqrt(aLx[x]*aLx[x]+aLy[x]*aLy[x]);
					
					aLaplacian[x] 	= aLxx[x]+aLyy[x]; 
					
					aRidgeDet[x]	= (aGradient[x]>0)? ((-2*aLx[x]*aLxy[x]*aLy[x]+aLxx[x]*aLy[x]*aLy[x]+aLx[x]*aLx[x]*aLyy[x])/(aGradient[x]*aGradient[x])):0;
					
					aIsophoteCurv[x] = (aGradient[x]>0)?(-aRidgeDet[x]/aGradient[x]):0;
					
					aFlowlineCurv[x] = (aGradient[x]>0)?(-(-aLx[x]*aLx[x]*aLxy[x]+aLxy[x]*aLy[x]*aLy[x]+aLx[x]*aLy[x]*(aLxx[x]-aLyy[x]))/(aGradient[x]*aGradient[x])):0;		
					
					aIsophoteDensity[x] = (aGradient[x]>0)?(-(aLx[x]*aLx[x]*aLxx[x]+2*aLx[x]*aLxy[x]*aLy[x]+aLy[x]*aLy[x]*aLyy[x])/(aGradient[x]*aGradient[x]*aGradient[x])):0;
					
					aCornerDet[x] = 2*aLx[x]*aLxy[x]*aLy[x]-aLxx[x]*aLy[x]*aLy[x]-aLx[x]*aLx[x]*aLyy[x];
					
					aShapeIndex[x] = (2/Math.PI)*Math.atan((-aLxx[x]-aLyy[x])/Math.sqrt(4*aLxy[x]*aLxy[x]+(aLxx[x]-aLyy[x])*(aLxx[x]-aLyy[x])));
					
					aCurvedness[x] = 0.5*Math.sqrt(aLxx[x]*aLxx[x]+2*aLxy[x]*aLxy[x]+aLyy[x]*aLyy[x]);
					
					aHessianDet[x] = -aLxy[x]*aLxy[x]+aLxx[x]*aLyy[x];
					
					aMeanCurvature[x] = 0.5*(aLxx[x]+aLyy[x]);
					
					aGaussianExtremality[x] = 
									-(aLxy[x]*aLxy[x]*(aLxxx[x]*aLxxx[x]+2*aLxxx[x]*aLxyy[x]-3*(aLxxy[x]*aLxxy[x]+aLxyy[x]*aLxyy[x]))+
									aLxxx[x]*aLxyy[x]*(aLxx[x]-aLyy[x])*(aLxx[x]-aLyy[x])+2*aLxxx[x]*aLxxy[x]*aLxy[x]*(-aLxx[x]+aLyy[x])+
									(aLxx[x]*aLxx[x]*aLxxy[x]-2*aLxy[x]*aLxyy[x]*aLyy[x]+2*aLxx[x]*(aLxy[x]*aLxyy[x]-aLxxy[x]*aLyy[x])+aLxxy[x]*(2*aLxy[x]*aLxy[x]+aLyy[x]*aLyy[x]))*aLyyy[x]+
									aLxy[x]*aLxy[x]*aLyyy[x]*aLyyy[x])/(aLxx[x]*aLxx[x]+4*aLxy[x]*aLxy[x]-2*aLxx[x]*aLyy[x]+aLyy[x]*aLyy[x]);
					
					aTJunctionLikeliness[x] = 
							Math.pow(aLx[x],5)*aLxyy[x]+Math.pow(aLy[x], 4)*(-2*aLxy[x]*aLxy[x]+aLxxy[x]*aLy[x]-aLxx[x]*aLyy[x])+
							Math.pow(aLx[x], 3)*aLy[x]*(6*aLxx[x]*aLxy[x]+aLxxx[x]*aLy[x]-aLxyy[x]*aLy[x]-6*aLxy[x]*aLyy[x])+
							aLx[x]*Math.pow(aLy[x], 3)*(-6*aLxx[x]*aLxy[x]+aLxxx[x]*aLy[x]-2*aLxyy[x]*aLy[x]+6*aLxy[x]*aLyy[x])-
							Math.pow(aLx[x], 4)*(2*aLxy[x]*aLxy[x]+2*aLxxy[x]*aLy[x]+aLxx[x]*aLyy[x]-aLy[x]*aLyyy[x])+
							aLx[x]*aLx[x]*aLy[x]*aLy[x]*(-3*aLxx[x]*aLxx[x]+8*aLxy[x]*aLxy[x]-aLxxy[x]*aLy[x]+4*aLxx[x]*aLyy[x]-3*aLyy[x]*aLyy[x]+aLy[x]*aLyyy[x]);
				}
				
				coords.z = i; // scale will index the layer
				
				gradient_image.set(coords,aGradient);
				laplacian_image.set(coords,aLaplacian);
				ridge_det_image.set(coords, aRidgeDet);
				isophote_curv_image.set(coords, aIsophoteCurv);
				flowline_curv_image.set(coords, aFlowlineCurv);
				isophote_density_image.set(coords, aIsophoteDensity);
				corner_det_image.set(coords, aCornerDet);
				shape_index_image.set(coords, aShapeIndex);
				curvedness_image.set(coords, aCurvedness);
				hessian_det_image.set(coords, aHessianDet);
				mean_curvature_image.set(coords, aMeanCurvature);
				gaussian_extremality_image.set(coords, aGaussianExtremality);
				t_junction_likeliness_image.set(coords, aTJunctionLikeliness);
				
				coords.z = 0;
				
			}
					
		}//scales
		
		System.out.println("done.");
		
	}
	
	public ImagePlus getGradientMagnitude(){
		ImagePlus gradient_image_ip = gradient_image.imageplus();
		gradient_image_ip.setTitle("Gradient magnitude");
		return gradient_image_ip;
	}
	
	public ImagePlus getLaplacian(){
		ImagePlus laplacian_image_ip = laplacian_image.imageplus();
		laplacian_image_ip.setTitle("Laplacian");
		return laplacian_image_ip;
	}
	
	public ImagePlus getRidgeDet(){
		ImagePlus ridge_det_image_ip = ridge_det_image.imageplus();
		ridge_det_image_ip.setTitle("Ridge detection");
		return ridge_det_image_ip;
	}
	
	public ImagePlus getIsophoteCurvature(){
		ImagePlus isophote_curv_image_ip = isophote_curv_image.imageplus();
		isophote_curv_image_ip.setTitle("Isophote curvature");
		return isophote_curv_image_ip;
	}
	
	public ImagePlus getFlowlineCurv(){
		ImagePlus flowline_curv_image_ip = flowline_curv_image.imageplus();
		flowline_curv_image_ip.setTitle("Flowline curvature");
		return flowline_curv_image_ip;
	}
	
	public ImagePlus getIsophoteDensity(){
		ImagePlus isophote_density_image_ip = isophote_density_image.imageplus();
		isophote_density_image_ip.setTitle("Isophote density");
		return isophote_density_image_ip;
	}
	
	public ImagePlus getCornerDetector(){
		ImagePlus corner_det_image_ip = corner_det_image.imageplus();
		corner_det_image_ip.setTitle("Affine invariant corner detector");
		return corner_det_image_ip;
	}
	
	public ImagePlus getShapeIndex(){
		ImagePlus shape_index_image_ip = shape_index_image.imageplus();
		shape_index_image_ip.setTitle("Shape index");
		return shape_index_image_ip;
	}
	
	public ImagePlus getCurvedness(){
		ImagePlus curvedness_image_ip = curvedness_image.imageplus();
		curvedness_image_ip.setTitle("Hessian determinant (gaussian curvature)");
		return curvedness_image_ip; 
	}
	
	public ImagePlus getHessianDeterminant(){
		ImagePlus hessian_det_image_ip = hessian_det_image.imageplus();
		hessian_det_image_ip.setTitle("Hessian determinant (gaussian curvature)");
		return hessian_det_image_ip; 
	}
	
	public ImagePlus getMeanCurvature(){
		ImagePlus mean_curvature_image_ip = mean_curvature_image.imageplus();
		mean_curvature_image_ip.setTitle("Mean curvature");
		return mean_curvature_image_ip;
	}
	
	public ImagePlus getGaussianExtremality(){
		ImagePlus gaussian_extremality_image_ip = gaussian_extremality_image.imageplus();
		gaussian_extremality_image_ip.setTitle("Extremality");
		return gaussian_extremality_image_ip;
	}
	
	public ImagePlus getTJunctionLikeliness(){
		ImagePlus t_junction_likeliness_image_ip = t_junction_likeliness_image.imageplus();
		t_junction_likeliness_image_ip.setTitle("T-junction likeliness");
		return t_junction_likeliness_image_ip;
	}
}