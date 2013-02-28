package advantra.feature;

import java.util.Vector;

import ij.ImagePlus;
import imagescience.feature.Differentiator;
import imagescience.feature.Hessian;
import imagescience.image.Axes;
import imagescience.image.Coordinates;
import imagescience.image.Dimensions;
import imagescience.image.FloatImage;
import imagescience.image.Image;

public class DifferentialFeatures {

	/*
	 * “Front-End Vision & Multi-Scale Image Analysis”, chapter 6: Differential structure of images
	 * and
	 * "Vascular Bifurcation Detection in Scale-Space", 
	 */

	/*
	 * class members
	 */
	
	// Input
	Image image; 
	// Auxiliary
	Image Lx, 	Ly;
	Image Lxx, 	Lyy, 	Lxy; 
	Image Lxyy, Lxxy, 	Lxxx, 	Lyyy;
	Image Lambda1, 		Lambda2;
	// Outputs (14 currently)
	public static int FEATS_NR = 14;
	Image gradient_image, 			laplacian_image, 			ridge_det_image, 			isophote_curv_image, 	flowline_curv_image;
	Image isophote_density_image, 	corner_det_image, 			shape_index_image, 			curvedness_image;
	Image mean_curvature_image, 	gaussian_extremality_image, t_junction_likeliness_image;
	Image ballness_filter,			DoH_filter; // |Lambda1|
	// Scales
	double[] sc;

	/*
	 * constructor
	 */
	
	public DifferentialFeatures(ImagePlus img, double sigma_1, double sigma_2, int nr){
		
		// set image
		if(img!=null){
			image = new FloatImage(Image.wrap(img));
			if(img.getStackSize()>1){
				System.out.println("Cannot work with 3d!");
				return;
			}
		}
			
		else
			image = null;
		
		// scales
		sc = new double[nr];
		for (int i = 0; i < nr; i++) {
			sc[i] = (i==0)? sigma_1 : sigma_1+i*((sigma_2-sigma_1)/(nr-1));
		}
		
		// the rest is set to null at the initialization
		Lx = Ly = null;
		Lxx = Lyy = Lxy = null; 
		Lxyy = Lxxy = Lxxx = Lyyy = null;
		Lambda1 = Lambda2 = null;
		// Outputs (14 currently)
		gradient_image = laplacian_image = ridge_det_image = isophote_curv_image = flowline_curv_image = null;
		isophote_density_image = corner_det_image = shape_index_image = curvedness_image = null;
		mean_curvature_image = gaussian_extremality_image = t_junction_likeliness_image = null;
		ballness_filter = DoH_filter = null; 
		
	}
	
	public void calculateFeatures(){
		
		System.out.println("calculating differential features...");
		
		final Dimensions dims 			= image.dimensions();
		final Dimensions filter_dims	= new Dimensions(dims.x, dims.y, sc.length); 
		// as many layers as there are scales
		final Differentiator 	df 		= new Differentiator();
		final Hessian			hs		= new Hessian();
		
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
		Lambda1	= new FloatImage(dims);
		Lambda2 = new FloatImage(dims);
		
		gradient_image				= new FloatImage(filter_dims); gradient_image.axes(Axes.X);
		laplacian_image				= new FloatImage(filter_dims); laplacian_image.axes(Axes.X);
		ridge_det_image				= new FloatImage(filter_dims); ridge_det_image.axes(Axes.X);
		isophote_curv_image			= new FloatImage(filter_dims); isophote_curv_image.axes(Axes.X);
		flowline_curv_image			= new FloatImage(filter_dims); flowline_curv_image.axes(Axes.X);
		isophote_density_image		= new FloatImage(filter_dims); isophote_density_image.axes(Axes.X);
		corner_det_image			= new FloatImage(filter_dims); corner_det_image.axes(Axes.X);
		shape_index_image			= new FloatImage(filter_dims); shape_index_image.axes(Axes.X);
		curvedness_image			= new FloatImage(filter_dims); curvedness_image.axes(Axes.X);
		//hessian_det_image			= new FloatImage(filter_dims); hessian_det_image.axes(Axes.X);
		mean_curvature_image		= new FloatImage(filter_dims); mean_curvature_image.axes(Axes.X);
		gaussian_extremality_image	= new FloatImage(filter_dims); gaussian_extremality_image.axes(Axes.X);
		t_junction_likeliness_image	= new FloatImage(filter_dims); t_junction_likeliness_image.axes(Axes.X);
		ballness_filter				= new FloatImage(filter_dims); ballness_filter.axes(Axes.X);
		DoH_filter					= new FloatImage(filter_dims); DoH_filter.axes(Axes.X);
		
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
		double[] aLambda1 = new double[dims.x];
		double[] aLambda2 = new double[dims.x];
					
		double[] aGradient 	= new double[dims.x];
		double[] aLaplacian = new double[dims.x];
		double[] aRidgeDet 	= new double[dims.x];
		double[] aIsophoteCurv = new double[dims.x];
		double[] aFlowlineCurv = new double[dims.x];
		double[] aIsophoteDensity = new double[dims.x];
		double[] aCornerDet = new double[dims.x];
		double[] aShapeIndex = new double[dims.x];
		double[] aCurvedness = new double[dims.x];
		//double[] aHessianDet = new double[dims.x];
		double[] aMeanCurvature = new double[dims.x];
		double[] aGaussianExtremality = new double[dims.x];
		double[] aTJunctionLikeliness = new double[dims.x];
		double[] aBallness = new double[dims.x];
		double[] aDoH = new double[dims.x];
		
		for (int i = 0; i < sc.length; i++) { // run through scales
			
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
			
			Vector<Image> l = hs.run(image.duplicate(), sc[i], true);
			Lambda2 = l.elementAt(0);
			Lambda1 = l.elementAt(1);
			
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
				
				Lambda1.get(coords, aLambda1);
				Lambda2.get(coords, aLambda2);
				
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
					
					//aHessianDet[x] = -aLxy[x]*aLxy[x]+aLxx[x]*aLyy[x];
					
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
					
					aBallness[x] 	= Math.abs(aLambda1[x]);
					
					aDoH[x]			= aLambda1[x]*aLambda2[x];	
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
				//hessian_det_image.set(coords, aHessianDet);
				mean_curvature_image.set(coords, aMeanCurvature);
				gaussian_extremality_image.set(coords, aGaussianExtremality);
				t_junction_likeliness_image.set(coords, aTJunctionLikeliness);
				
				ballness_filter.set(coords, aBallness);
				DoH_filter.set(coords, aDoH);
				
				coords.z = 0;
				
			}
					
		}//scales
		
		System.out.println("done.");
		
	}

	/*
	 * methods
	 */
	// gradient mag.
	public ImagePlus 	getGradientMagnitude(){
		ImagePlus gradient_image_ip;
		if(gradient_image!=null) 
			gradient_image_ip = gradient_image.imageplus();
		else {
			System.out.println("Differential features are not calculated.");
			return null;
		}
		gradient_image_ip.setTitle("Gradient magnitude");
		return gradient_image_ip;
	}
	
	public double 		getGradientMagnitude(int[] at_pos){  // at_pos contains (row,col)
		Coordinates coords = new Coordinates(at_pos[0], at_pos[1], at_pos[2]); // coords contains (col,row,lay)
		
		if(gradient_image!=null)
			return gradient_image.get(coords);
		else {
			System.out.println("Differential features are not calculated.");
 			return Double.NaN;
		}
			
	}
	
	// laplacian
 	public ImagePlus 	getLaplacian(){
		ImagePlus laplacian_image_ip;
		if(laplacian_image!=null)
			laplacian_image_ip = laplacian_image.imageplus();
		else{
			System.out.println("Differential features are not calculated.");
			return null;
		}
		laplacian_image_ip.setTitle("Laplacian");
		return laplacian_image_ip;
	}
 	
 	public double 		getLaplacian(int[] at_pos){
 		Coordinates coords = new Coordinates(at_pos[0], at_pos[1], at_pos[2]);
 		if(laplacian_image!=null)
 			return laplacian_image.get(coords);
 		else{
 			System.out.println("Differential features are not calculated.");
 			return Double.NaN;
 		}
 	}
	
 	// ridge detection
	public ImagePlus 	getRidgeDet(){
		ImagePlus ridge_det_image_ip;
		if(ridge_det_image!=null) 
			ridge_det_image_ip = ridge_det_image.imageplus();
		else {
			System.out.println("Differential features are not calculated.");
			return null;
		}
		ridge_det_image_ip.setTitle("Ridge detection");
		return ridge_det_image_ip;
	}
	
	public double  		getRidgeDet(int[] at_pos){ 
		Coordinates coords = new Coordinates(at_pos[0], at_pos[1], at_pos[2]);
 		if(ridge_det_image!=null)
 			return ridge_det_image.get(coords);
 		else{
 			System.out.println("Differential features are not calculated.");
 			return Double.NaN;
 		}
	}
	
	// isophote curvature
	public ImagePlus getIsophoteCurvature(){
		ImagePlus isophote_curv_image_ip;
		if(isophote_curv_image!=null) 
			isophote_curv_image_ip = isophote_curv_image.imageplus();
		else {
			System.out.println("Differential features are not calculated.");
			return null;
		}
		isophote_curv_image_ip.setTitle("Isophote curvature");
		return isophote_curv_image_ip;
	}
	
	public double  	getIsophoteCurvature(int[] at_pos){
		Coordinates coords = new Coordinates(at_pos[0], at_pos[1], at_pos[2]);
 		if(isophote_curv_image!=null)
 			return isophote_curv_image.get(coords);
 		else{
 			System.out.println("Differential features are not calculated.");
 			return Double.NaN;
 		}
	}
	
	// flowline curvature
	public ImagePlus 	getFlowlineCurv(){
		ImagePlus flowline_curv_image_ip;
		if(flowline_curv_image!=null) 
			flowline_curv_image_ip = flowline_curv_image.imageplus();
		else {
			System.out.println("Differential features are not calculated.");
			return null;
		}
		flowline_curv_image_ip.setTitle("Flowline curvature");
		return flowline_curv_image_ip;
	}
	
	public double 		getFlowlineCurv(int[] at_pos){
		Coordinates coords = new Coordinates(at_pos[0], at_pos[1], at_pos[2]);
 		if(flowline_curv_image!=null) 
 			return flowline_curv_image.get(coords);
 		else{
 			System.out.println("Differential features are not calculated.");
 			return Double.NaN;
 		}
	}
	
	// isophote density
	public ImagePlus 	getIsophoteDensity(){
		ImagePlus isophote_density_image_ip;
		if(isophote_density_image!=null) 
			isophote_density_image_ip = isophote_density_image.imageplus();
		else {
			System.out.println("Differential features are not calculated.");
			return null;
		}
		isophote_density_image_ip.setTitle("Isophote density");
		return isophote_density_image_ip;
	}
	
	public 	double		getIsophoteDensity(int[] at_pos){
		Coordinates coords = new Coordinates(at_pos[0], at_pos[1], at_pos[2]);
 		if(isophote_density_image!=null)
 			return isophote_density_image.get(coords);
 		else{
 			System.out.println("Differential features are not calculated.");
 			return Double.NaN;
 		}
	}
	
	// corner detector
	public ImagePlus getCornerDetector(){
		ImagePlus corner_det_image_ip;
		if(corner_det_image!=null) 
			corner_det_image_ip = corner_det_image.imageplus();
		else {
			System.out.println("Differential features are not calculated.");
			return null;
		}
		corner_det_image_ip.setTitle("Affine invariant corner detector");
		return corner_det_image_ip;
	}
	
	public double	getCornerDetector(int[] at_pos){
		Coordinates coords = new Coordinates(at_pos[0], at_pos[1], at_pos[2]);
 		if(corner_det_image!=null) 
 			return corner_det_image.get(coords);
 		else{
 			System.out.println("Differential features are not calculated.");
 			return Double.NaN;
 		}
	}
	
	// shape index
	public ImagePlus getShapeIndex(){
		ImagePlus shape_index_image_ip;
		if(shape_index_image!=null) 
			shape_index_image_ip = shape_index_image.imageplus();
		else {
			System.out.println("Differential features are not calculated.");
			return null;
		}
		shape_index_image_ip.setTitle("Shape index");
		return shape_index_image_ip;
	}
	
	public double 	getShapeIndex(int[] at_pos){
		Coordinates coords = new Coordinates(at_pos[0], at_pos[1], at_pos[2]);
 		if(shape_index_image!=null) 
 			return shape_index_image.get(coords);
 		else{
 			System.out.println("Differential features are not calculated.");
 			return Double.NaN;
 		}
	}
	
	// curvedness
	public ImagePlus getCurvedness(){
		ImagePlus curvedness_image_ip;
		if(curvedness_image!=null) 
			curvedness_image_ip = curvedness_image.imageplus();
		else {
			System.out.println("Differential features are not calculated.");
			return null;
		}
		curvedness_image_ip.setTitle("Hessian determinant (gaussian curvature)");
		return curvedness_image_ip; 
	}
	
	public double	getCurvedness(int[] at_pos){
		Coordinates coords = new Coordinates(at_pos[0], at_pos[1], at_pos[2]);
 		if(curvedness_image!=null) 
 			return curvedness_image.get(coords);
 		else{
 			System.out.println("Differential features are not calculated.");
 			return Double.NaN;
 		}
	}
	
	// DoH
//	public ImagePlus getHessianDeterminant(){
//		ImagePlus hessian_det_image_ip = hessian_det_image.imageplus();
//		hessian_det_image_ip.setTitle("Hessian determinant (gaussian curvature)");
//		return hessian_det_image_ip; 
//	}
//	
//	public double 	getHessianDeterminant(int[] at_pos){
//		Coordinates coords = new Coordinates(at_pos[0], at_pos[1], at_pos[2]);
// 		return hessian_det_image.get(coords);
//	}
	
	// mean curv.
	public ImagePlus getMeanCurvature(){
		ImagePlus mean_curvature_image_ip;
		if(mean_curvature_image!=null) 
			mean_curvature_image_ip = mean_curvature_image.imageplus();
		else {
			System.out.println("Differential features are not calculated.");
			return null;
		}
		mean_curvature_image_ip.setTitle("Mean curvature");
		return mean_curvature_image_ip;
	}
	
	public double getMeanCurvature(int[] at_pos){
		
		Coordinates coords = new Coordinates(at_pos[0], at_pos[1], at_pos[2]);
 		if(mean_curvature_image!=null)
 			return mean_curvature_image.get(coords);
 		else{
 			System.out.println("Differential features are not calculated.");
 			return Double.NaN;
 		}
	}
	
	// gaussian extremality
	public ImagePlus getGaussianExtremality(){
		ImagePlus gaussian_extremality_image_ip;
		if(gaussian_extremality_image!=null) 
			gaussian_extremality_image_ip = gaussian_extremality_image.imageplus();
		else {
			System.out.println("Differential features are not calculated.");
			return null;
		}
		gaussian_extremality_image_ip.setTitle("Extremality");
		return gaussian_extremality_image_ip;
	}
	
	public double getGaussianExtremality(int[] at_pos){
		Coordinates coords = new Coordinates(at_pos[0], at_pos[1], at_pos[2]);
 		if(gaussian_extremality_image!=null) 
 			return gaussian_extremality_image.get(coords);
 		else{
 			System.out.println("Differential features are not calculated.");
 			return Double.NaN;
 		}
	}
	
	// t-junct.
	public ImagePlus getTJunctionLikeliness(){
		ImagePlus t_junction_likeliness_image_ip;
		if(t_junction_likeliness_image!=null) 
			t_junction_likeliness_image_ip = t_junction_likeliness_image.imageplus();
		else {
			System.out.println("Differential features are not calculated.");
			return null;
		}
		t_junction_likeliness_image_ip.setTitle("T-junction likeliness");
		return t_junction_likeliness_image_ip;
	}
	
	public double 	getTJunctionLikeliness(int[] at_pos){
		Coordinates coords = new Coordinates(at_pos[0], at_pos[1], at_pos[2]);
 		if(t_junction_likeliness_image!=null)
 			return t_junction_likeliness_image.get(coords);
 		else{
 			System.out.println("Differential features are not calculated.");
 			return Double.NaN;
 		}
	}
	
	// ballness
	public ImagePlus getBallness(){
		ImagePlus ballness_filter_ip;
		if(ballness_filter!=null) 
			ballness_filter_ip = ballness_filter.imageplus();
		else {
			System.out.println("Differential features are not calculated.");
			return null;
		}
		ballness_filter_ip.setTitle("Ballness");
		return ballness_filter_ip;
	}
	
	public double 	getBallness(int[] at_pos){
		Coordinates coords = new Coordinates(at_pos[0], at_pos[1], at_pos[2]);
 		if(ballness_filter!=null) 
 			return ballness_filter.get(coords);
 		else{
 			System.out.println("Differential features are not calculated.");
 			return Double.NaN;
 		}
	}
	
	// DoH (just to compare)
	public ImagePlus getDoH(){
		ImagePlus DoH_filter_ip;
		if(DoH_filter!=null) 
			DoH_filter_ip = DoH_filter.imageplus();
		else {
			System.out.println("Differential features are not calculated.");
			return null;
		}
		DoH_filter_ip.setTitle("DoH");
		return DoH_filter_ip;
	}
	
	public double 	getDoH(int[] at_pos){
		Coordinates coords = new Coordinates(at_pos[0], at_pos[1], at_pos[2]);
 		if(DoH_filter!=null) 
 			return DoH_filter.get(coords);
 		else {
 			System.out.println("Differential features are not calculated.");
 			return Double.NaN;
 		}
	}	
	
	/*
	 * export features for different scales:
	 */
	
	public double[][] 		exportFeatures(double[][] locations, boolean[] which_features){// locations contains row, col obtained from each line
		
		if(which_features.length!=FEATS_NR){
			System.out.println("DifferentialFeatures:exportFeatures(): Boolean array describing the number of features has to have length "+FEATS_NR);
			return null;
		}
		
		int nr_features = 0;
		for (int i = 0; i < which_features.length; i++) {
			if(which_features[i]) nr_features++;
		}
		
		int nr_scales 	= sc.length;
		int[] pos 		= new int[3];
		
		double[][] feat = new double[locations.length][nr_scales*nr_features];
		
		for (int feat_row = 0; feat_row < locations.length; feat_row++) {
			
			int feat_col = 0;
			
			/*
			 *  switching cols & rows here below at pos[0] and pos[1] because imagescience 
			 *  considers Coords.x as column position and Coords.y as row position
			 */
			
			pos[0] = (int)Math.round(locations[feat_row][0]); // col (location[0]) will be assigned to Coords.x
			pos[1] = (int)Math.round(locations[feat_row][1]); // row (location[1]) will be assigned to Coords.y
			//System.out.println("pos(0,1,2): "+pos[0]+", "+pos[1]+", "+pos[2]);
			//System.out.println("here:");
			//for (int i = 0; i < which_features.length; i++) {
			//	if(which_features[i]) System.out.print(i+", ");
			//}
			//System.out.println();
			
			for (int i = 0; i < nr_scales; i++) { 
				
				pos[2] = i; // scale will define the layer in the images
				
				if(which_features[0]) feat[feat_row][feat_col++] = getGradientMagnitude(pos);
				if(which_features[1]) feat[feat_row][feat_col++] = getLaplacian(pos); 
				if(which_features[2]) feat[feat_row][feat_col++] = getRidgeDet(pos); 
				if(which_features[3]) feat[feat_row][feat_col++] = getIsophoteCurvature(pos); 
				if(which_features[4]) feat[feat_row][feat_col++] = getFlowlineCurv(pos); 
				if(which_features[5]) feat[feat_row][feat_col++] = getIsophoteDensity(pos);
				if(which_features[6]) feat[feat_row][feat_col++] = getCornerDetector(pos); 
				if(which_features[7]) feat[feat_row][feat_col++] = getShapeIndex(pos); 
				if(which_features[8]) feat[feat_row][feat_col++] = getCurvedness(pos); 
				if(which_features[9]) feat[feat_row][feat_col++] = getDoH(pos);
				if(which_features[10]) feat[feat_row][feat_col++] = getMeanCurvature(pos); 
				if(which_features[11]) feat[feat_row][feat_col++] = getGaussianExtremality(pos); 
				if(which_features[12]) feat[feat_row][feat_col++] = getTJunctionLikeliness(pos); 
				if(which_features[13]) feat[feat_row][feat_col++] = getBallness(pos); 
				
			}
			
		}
		return feat;
	}
	
	public static String[]	exportAllLabels(){
		
		String[] all_feat_names = new String[]{
				"GradientMag",
				"Laplacian",
				"RidgeDetector",
				"IsophoteCurvature",
				"FlowlineCurvature",
				"IsophoteDensity",
				"CornerDetector",
				"ShapeIndex",
				"Curvedness",
				"DoH",
				"MeanCurvature",
				"GaussianExtremality",
				"T-JunctionLikeliness",
				"Ballness"
		};
		return all_feat_names;
		
	}
	
	public String[] 		exportFeatureLabels(boolean with_scales, boolean[] which_features){// locations contains row, col obtained from each line
		
		if(which_features.length!=FEATS_NR){
			System.out.println("DifferentialFeatures:exportFeatures(): Boolean array describing the number of features has to have length "+FEATS_NR);
			return null;
		}
		
		int nr_features = 0;
		for (int i = 0; i < FEATS_NR; i++) {
			if(which_features[i]) nr_features++;
		}
		
		String[] feat_names;
		int feat_col = 0;
		if(with_scales){
			int nr_scales = sc.length;
			feat_names = new String[nr_scales*nr_features];
			for (int i = 0; i < nr_scales; i++) { 
				if(which_features[0]) feat_names[feat_col++] = String.format("GradientMag_s%.2f", 		sc[i]);
				if(which_features[1]) feat_names[feat_col++] = String.format("Laplacian_s%.2f", 		sc[i]);
				if(which_features[2]) feat_names[feat_col++] = String.format("RidgeDetector_s%.2f", 	sc[i]);
				if(which_features[3]) feat_names[feat_col++] = String.format("IsophoteCurvature_s%.2f", sc[i]);
				if(which_features[4]) feat_names[feat_col++] = String.format("FlowlineCurvature_s%.2f", sc[i]);
				if(which_features[5]) feat_names[feat_col++] = String.format("IsophoteDensity_s%.2f", 	sc[i]);
				if(which_features[6]) feat_names[feat_col++] = String.format("CornerDetector_s%.2f", 	sc[i]);
				if(which_features[7]) feat_names[feat_col++] = String.format("ShapeIndex_s%.2f", 		sc[i]);
				if(which_features[8]) feat_names[feat_col++] = String.format("Curvedness_s%.2f", 		sc[i]);
				if(which_features[9]) feat_names[feat_col++] = String.format("DoH_s%.2f", 				sc[i]);
				if(which_features[10]) feat_names[feat_col++] = String.format("MeanCurvature_s%.2f",	sc[i]);
				if(which_features[11]) feat_names[feat_col++] = String.format("GaussianExtremality_s%.2f",sc[i]);
				if(which_features[12]) feat_names[feat_col++] = String.format("T-JunctionLikeliness_s%.2f",sc[i]);
				if(which_features[13]) feat_names[feat_col++] = String.format("Ballness_s%.2f",			sc[i]);
			}
		}
		else{
			feat_names = new String[nr_features];
			if(which_features[0]) feat_names[feat_col++] 	= "GradientMag";
			if(which_features[1]) feat_names[feat_col++] 	= "Laplacian";
			if(which_features[2]) feat_names[feat_col++] 	= "RidgeDetector";
			if(which_features[3]) feat_names[feat_col++] 	= "IsophoteCurvature";
			if(which_features[4]) feat_names[feat_col++] 	= "FlowlineCurvature";
			if(which_features[5]) feat_names[feat_col++] 	= "IsophoteDensity";
			if(which_features[6]) feat_names[feat_col++] 	= "CornerDetector";
			if(which_features[7]) feat_names[feat_col++] 	= "ShapeIndex";
			if(which_features[8]) feat_names[feat_col++] 	= "Curvedness";
			if(which_features[9]) feat_names[feat_col++] 	= "DoH";
			if(which_features[10]) feat_names[feat_col++] 	= "MeanCurvature";
			if(which_features[11]) feat_names[feat_col++] 	= "GaussianExtremality";
			if(which_features[12]) feat_names[feat_col++] 	= "T-JunctionLikeliness";
			if(which_features[13]) feat_names[feat_col++] 	= "Ballness";
		}
		return feat_names;
	}

}