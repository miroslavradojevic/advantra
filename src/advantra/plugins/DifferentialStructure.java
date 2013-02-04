package advantra.plugins;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import imagescience.feature.Differentiator;
import imagescience.image.Axes;
import imagescience.image.Coordinates;
import imagescience.image.Dimensions;
import imagescience.image.FloatImage;
import imagescience.image.Image;

public class DifferentialStructure implements PlugInFilter {

	/*
	 * terminal use
	 */
	
	public static void main(String[] args){
		
		double 	sigma_1 = 1;
		double 	sigma_2 = 1;
		int 	nr 		= 1;
		String image_path = "";
		
		if(args.length!=4){
			System.out.println("usage: arg0=SIGMA1 arg1=SIGMA2 arg2=NR_SCALES arg3=IMAGE_PATH");
		}
		else{
			sigma_1 = Double.parseDouble(args[0]);
			sigma_2 = Double.parseDouble(args[1]);
			nr 		= Integer.parseInt(args[2]);
			image_path = args[3];
		}
		
		Image Lx, 	Ly, 	Lxx, 	Lyy, 	Lxy; 
		Image Lxyy, Lxxy, 	Lxxx, 	Lyyy;
		
		double[] sc = new double[nr];
		for (int i = 0; i < nr; i++) {
			sc[i] = (i==0)? sigma_1 : sigma_1+i*((sigma_2-sigma_1)/(nr-1));
		}
		
		// TODO: finish it, pretty much copy what's in imagej part
		
	}
	
	
	/*
	 * imagej part
	 * 
	 */
	
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
	// 3D image
	Image Lz, 	Lzz, 	Lxz, 	Lyz; 
	
	public void run(ImageProcessor arg0) {
		
		// dialog to enter input values
		GenericDialog gd = new GenericDialog("DifferentialStructure", IJ.getInstance());

		gd.addMessage("Choose scales:");
		
		gd.addNumericField( "sigma start:", sigma_1, 0, 5, "pix" );
		gd.addNumericField( "sigma end  :", sigma_2, 0, 5, "pix" );	
		gd.addNumericField( "number of scales : ", nr,  0, 5, "");
		
		gd.addMessage("Show descriptors:");
		
		gd.addCheckbox("gradient", 				true); // 1
		gd.addCheckbox("laplacian", 			true);
		gd.addCheckbox("ridge detection", 		true);
		gd.addCheckbox("isophote curvature",	true);
		gd.addCheckbox("flowline curvature",	true); // 5
		gd.addCheckbox("isophote density", 		true);
		gd.addCheckbox("corner detector", 		true);
		gd.addCheckbox("shape index", 			true);
		gd.addCheckbox("curvedness", 			true);
		gd.addCheckbox("hessian determinant",	true); // 10
		gd.addCheckbox("mean curvature",		true);
		gd.addCheckbox("gaussian extremality",	true);
		gd.addCheckbox("t junction likeliness",	true);
		
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
		
		double[] sc = new double[nr];
		for (int i = 0; i < nr; i++) {
			sc[i] = (i==0)? sigma_1 : sigma_1+i*((sigma_2-sigma_1)/(nr-1));
		}
		
		Image image = new FloatImage(Image.wrap(img));
		final Dimensions dims 		= image.dimensions();
		final Dimensions filter_dims= new Dimensions(dims.x, dims.y, nr); // as many layers as there are scales
		final Differentiator df 	= new Differentiator();
		
		if (dims.z == 1) { // 2D case
			
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
			
			Image gradient_image	= new FloatImage(filter_dims); gradient_image.axes(Axes.X);
			Image laplacian_image	= new FloatImage(filter_dims); laplacian_image.axes(Axes.X);
			Image ridge_det_image	= new FloatImage(filter_dims); ridge_det_image.axes(Axes.X);
			Image isophote_curv_image	= new FloatImage(filter_dims); isophote_curv_image.axes(Axes.X);
			Image flowline_curv_image	= new FloatImage(filter_dims); flowline_curv_image.axes(Axes.X);
			Image isophote_density_image	= new FloatImage(filter_dims); isophote_density_image.axes(Axes.X);
			Image corner_det_image	= new FloatImage(filter_dims); corner_det_image.axes(Axes.X);
			Image shape_index_image	= new FloatImage(filter_dims); shape_index_image.axes(Axes.X);
			Image curvedness_image	= new FloatImage(filter_dims); curvedness_image.axes(Axes.X);
			Image hessian_det_image	= new FloatImage(filter_dims); hessian_det_image.axes(Axes.X);
			Image mean_curvature_image	= new FloatImage(filter_dims); mean_curvature_image.axes(Axes.X);
			Image gaussian_extremality_image	= new FloatImage(filter_dims); gaussian_extremality_image.axes(Axes.X);
			Image t_junction_likeliness_image	= new FloatImage(filter_dims); t_junction_likeliness_image.axes(Axes.X);
			
			// gradients
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
					
//					if(i!=0){
//						Lw_image.get(coords,aLw);
//					}
//					for (int x=0; x<dims.x; ++x){
//						double grad = Math.sqrt(aLx[x]*aLx[x] + aLy[x]*aLy[x]);
//						if(i==0){
//							aLw[x] = grad;
//						}
//						else{
//							if(grad>aLw[x]) aLw[x] = grad;
//						}
//					}
//					Lw_image.set(coords,aLw);
					
				}
						
			}//scales
			
			if (gradient_show){
				ImagePlus gradient_image_ip = gradient_image.imageplus();
				gradient_image_ip.setTitle("Gradient sqrt(Lx^2+Ly^2)");
				gradient_image_ip.show();
			}
			if (laplacian_show){
				ImagePlus laplacian_image_ip = laplacian_image.imageplus();
				laplacian_image_ip.setTitle("Laplacian (Lxx+Lyy)");
				laplacian_image_ip.show();
			}
			if (ridge_det_show){
				ImagePlus ridge_det_image_ip = ridge_det_image.imageplus();
				ridge_det_image_ip.setTitle("Ridge detection (Lvv)");
				ridge_det_image_ip.show();
			}
			if(isophote_curv_show){
				ImagePlus isophote_curv_image_ip = isophote_curv_image.imageplus();
				isophote_curv_image_ip.setTitle("Isophote curvature (-Lvv/Lw)");
				isophote_curv_image_ip.show();
			}
			if(flowline_curv_show){
				ImagePlus flowline_curv_image_ip = flowline_curv_image.imageplus();
				flowline_curv_image_ip.setTitle("Flowline curvature (-Lvw/Lw)");
				flowline_curv_image_ip.show();
			}
			if(isophote_density_show){
				ImagePlus isophote_density_image_ip = isophote_curv_image.imageplus();
				isophote_density_image_ip.setTitle("Isophote density (-Lww/Lw)");
				isophote_density_image_ip.show();
			}
			if(corner_det_show){
				ImagePlus corner_det_image_ip = corner_det_image.imageplus();
				corner_det_image_ip.setTitle("Affine invariant corner detector (LvvLw^2)");
				corner_det_image_ip.show();
			}
			if(shape_index_show){
				ImagePlus shape_index_image_ip = shape_index_image.imageplus();
				shape_index_image_ip.setTitle("Shape index");
				shape_index_image_ip.show();
			}
			if(curvedness_show){
				ImagePlus curvedness_image_ip = curvedness_image.imageplus();
				curvedness_image_ip.setTitle("Curvedness");
				curvedness_image_ip.show();
			}
			if(hessian_det_show){
				ImagePlus hessian_det_image_ip = hessian_det_image.imageplus();
				hessian_det_image_ip.setTitle("Hessian determinant (gaussian curvature)");
				hessian_det_image_ip.show();
			}
			if(mean_curvature_show){
				ImagePlus mean_curvature_image_ip = mean_curvature_image.imageplus();
				mean_curvature_image_ip.setTitle("Mean curvature (k1+k2)/2");
				mean_curvature_image_ip.show();
			}
			if(gaussian_extremality_show){
				ImagePlus gaussian_extremality_image_ip = gaussian_extremality_image.imageplus();
				gaussian_extremality_image_ip.setTitle("Extremality");
				gaussian_extremality_image_ip.show();
			}
			if(t_junction_likeliness_show){
				ImagePlus t_junction_likeliness_image_ip = t_junction_likeliness_image.imageplus();
				t_junction_likeliness_image_ip.setTitle("T-junction likeliness");
				t_junction_likeliness_image_ip.show();
			}
		}
		else{ //3D case
			// TODO
		}
		
		System.out.println("done calculating filters.");
		
	}

	public int setup(String arg0, ImagePlus im) {
		this.img = im;
		return DOES_8G+NO_CHANGES;
	}
	
}