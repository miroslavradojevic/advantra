package advantra.feature;

import java.util.Vector;

import advantra.general.ArrayHandling;

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
		
		if(img.getStackSize()>1){
			System.out.println("DifferentialFeatures:DifferentialFeatures() works with 2D images!");
			return;
		}
		
		sc = new double[nr];
		for (int i = 0; i < nr; i++) {
			sc[i] = (i==0)? sigma_1 : sigma_1+i*((sigma_2-sigma_1)/(nr-1));
		}
		
		image = new FloatImage(Image.wrap(img));
		final Dimensions dims 			= image.dimensions();
		final Dimensions filter_dims	= new Dimensions(dims.x, dims.y, nr); 
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
		System.out.println("size dims(x,y,z): "+dims.x+" , "+dims.y+" , "+dims.z);
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
		ImagePlus gradient_image_ip = gradient_image.imageplus();
		gradient_image_ip.setTitle("Gradient magnitude");
		return gradient_image_ip;
	}
	
	public double 		getGradientMagnitude(int[] at_pos){  // at_pos contains (row,col)
		Coordinates coords = new Coordinates(at_pos[0], at_pos[1], at_pos[2]); // coords contains (col,row,lay)
		return gradient_image.get(coords);
	}
	
	// laplacian
 	public ImagePlus 	getLaplacian(){
		ImagePlus laplacian_image_ip = laplacian_image.imageplus();
		laplacian_image_ip.setTitle("Laplacian");
		return laplacian_image_ip;
	}
 	
 	public double 		getLaplacian(int[] at_pos){
 		Coordinates coords = new Coordinates(at_pos[0], at_pos[1], at_pos[2]);
 		return laplacian_image.get(coords);
 	}
	
 	// ridge detection
	public ImagePlus 	getRidgeDet(){
		ImagePlus ridge_det_image_ip = ridge_det_image.imageplus();
		ridge_det_image_ip.setTitle("Ridge detection");
		return ridge_det_image_ip;
	}
	
	public double  		getRidgeDet(int[] at_pos){ 
		Coordinates coords = new Coordinates(at_pos[0], at_pos[1], at_pos[2]);
 		return ridge_det_image.get(coords);
	}
	
	// isophote curvature
	public ImagePlus getIsophoteCurvature(){
		ImagePlus isophote_curv_image_ip = isophote_curv_image.imageplus();
		isophote_curv_image_ip.setTitle("Isophote curvature");
		return isophote_curv_image_ip;
	}
	
	public double  	getIsophoteCurvature(int[] at_pos){
		Coordinates coords = new Coordinates(at_pos[0], at_pos[1], at_pos[2]);
 		return isophote_curv_image.get(coords);
	}
	
	// flowline curvature
	public ImagePlus 	getFlowlineCurv(){
		ImagePlus flowline_curv_image_ip = flowline_curv_image.imageplus();
		flowline_curv_image_ip.setTitle("Flowline curvature");
		return flowline_curv_image_ip;
	}
	
	public double 		getFlowlineCurv(int[] at_pos){
		Coordinates coords = new Coordinates(at_pos[0], at_pos[1], at_pos[2]);
 		return flowline_curv_image.get(coords);
	}
	
	// isophote density
	public ImagePlus 	getIsophoteDensity(){
		ImagePlus isophote_density_image_ip = isophote_density_image.imageplus();
		isophote_density_image_ip.setTitle("Isophote density");
		return isophote_density_image_ip;
	}
	
	public 	double		getIsophoteDensity(int[] at_pos){
		Coordinates coords = new Coordinates(at_pos[0], at_pos[1], at_pos[2]);
 		return isophote_density_image.get(coords);
	}
	
	// corner detector
	public ImagePlus getCornerDetector(){
		ImagePlus corner_det_image_ip = corner_det_image.imageplus();
		corner_det_image_ip.setTitle("Affine invariant corner detector");
		return corner_det_image_ip;
	}
	
	public double	getCornerDetector(int[] at_pos){
		Coordinates coords = new Coordinates(at_pos[0], at_pos[1], at_pos[2]);
 		return corner_det_image.get(coords);
	}
	
	// shape index
	public ImagePlus getShapeIndex(){
		ImagePlus shape_index_image_ip = shape_index_image.imageplus();
		shape_index_image_ip.setTitle("Shape index");
		return shape_index_image_ip;
	}
	
	public double 	getShapeIndex(int[] at_pos){
		Coordinates coords = new Coordinates(at_pos[0], at_pos[1], at_pos[2]);
 		return shape_index_image.get(coords);
	}
	
	// curvedness
	public ImagePlus getCurvedness(){
		ImagePlus curvedness_image_ip = curvedness_image.imageplus();
		curvedness_image_ip.setTitle("Hessian determinant (gaussian curvature)");
		return curvedness_image_ip; 
	}
	
	public double	getCurvedness(int[] at_pos){
		Coordinates coords = new Coordinates(at_pos[0], at_pos[1], at_pos[2]);
 		return curvedness_image.get(coords);
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
		ImagePlus mean_curvature_image_ip = mean_curvature_image.imageplus();
		mean_curvature_image_ip.setTitle("Mean curvature");
		return mean_curvature_image_ip;
	}
	
	public double getMeanCurvature(int[] at_pos){
		Coordinates coords = new Coordinates(at_pos[0], at_pos[1], at_pos[2]);
 		return mean_curvature_image.get(coords);
	}
	
	// gaussian extremality
	public ImagePlus getGaussianExtremality(){
		ImagePlus gaussian_extremality_image_ip = gaussian_extremality_image.imageplus();
		gaussian_extremality_image_ip.setTitle("Extremality");
		return gaussian_extremality_image_ip;
	}
	
	public double getGaussianExtremality(int[] at_pos){
		Coordinates coords = new Coordinates(at_pos[0], at_pos[1], at_pos[2]);
 		return gaussian_extremality_image.get(coords);
	}
	
	// t-junct.
	public ImagePlus getTJunctionLikeliness(){
		ImagePlus t_junction_likeliness_image_ip = t_junction_likeliness_image.imageplus();
		t_junction_likeliness_image_ip.setTitle("T-junction likeliness");
		return t_junction_likeliness_image_ip;
	}
	
	public double 	getTJunctionLikeliness(int[] at_pos){
		Coordinates coords = new Coordinates(at_pos[0], at_pos[1], at_pos[2]);
 		return t_junction_likeliness_image.get(coords);
	}
	
	// ballness
	public ImagePlus getBallness(){
		ImagePlus ballness_filter_ip = ballness_filter.imageplus();
		ballness_filter_ip.setTitle("Ballness");
		return ballness_filter_ip;
	}
	
	public double 	getBallness(int[] at_pos){
		Coordinates coords = new Coordinates(at_pos[0], at_pos[1], at_pos[2]);
 		return ballness_filter.get(coords);
	}
	
	// DoH (just to compare)
	public ImagePlus getDoH(){
		ImagePlus DoH_filter_ip = DoH_filter.imageplus();
		DoH_filter_ip.setTitle("DoH");
		return DoH_filter_ip;
	}
	
	public double 	getDoH(int[] at_pos){
		Coordinates coords = new Coordinates(at_pos[0], at_pos[1], at_pos[2]);
 		return DoH_filter.get(coords);
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
		
		//ArrayHandling.print2DArray(locations);
		
		for (int feat_row = 0; feat_row < locations.length; feat_row++) {
			
			//System.out.println("loc: "+locations[feat_row][0]+" , "+locations[feat_row][1]+"  ");
			
			int feat_col = 0;
			
			/*
			 *  switching cols & rows here below at pos[0] and pos[1] because imagescience 
			 *  considers Coords.x as column position and Coords.y as row position
			 */
			
			pos[1] = (int)Math.round(locations[feat_row][0]); // row location will be assigned to Coords.y
			pos[0] = (int)Math.round(locations[feat_row][1]); // col location will be assigned to Coords.x
			
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
	
	public static String[] 	exportFeatureLabels(int nr_scales, boolean[] which_features){// locations contains row, col obtained from each line
		
		if(which_features.length!=FEATS_NR){
			System.out.println("DifferentialFeatures:exportFeatures(): Boolean array describing the number of features has to have length "+FEATS_NR);
			return null;
		}
		
		int nr_features = 0;
		for (int i = 0; i < FEATS_NR; i++) {
			if(which_features[i]) nr_features++;
		}
		
		String[] feat_names = new String[nr_scales*nr_features];
		int feat_col = 0;
		
		for (int i = 0; i < nr_scales; i++) { // scale will define the layer
		
			if(which_features[0]) feat_names[feat_col++] = String.format("df_ft_01s%d", i);
			if(which_features[1]) feat_names[feat_col++] = String.format("df_ft_02s%d", i);
			if(which_features[2]) feat_names[feat_col++] = String.format("df_ft_03s%d", i);
			if(which_features[3]) feat_names[feat_col++] = String.format("df_ft_04s%d", i);
			if(which_features[4]) feat_names[feat_col++] = String.format("df_ft_05s%d", i);
			if(which_features[5]) feat_names[feat_col++] = String.format("df_ft_06s%d", i);
			if(which_features[6]) feat_names[feat_col++] = String.format("df_ft_07s%d", i);
			if(which_features[7]) feat_names[feat_col++] = String.format("df_ft_08s%d", i);
			if(which_features[8]) feat_names[feat_col++] = String.format("df_ft_09s%d", i);
			if(which_features[9]) feat_names[feat_col++] = String.format("df_ft_10s%d", i);
			if(which_features[10]) feat_names[feat_col++] = String.format("df_ft_11s%d", i);
			if(which_features[11]) feat_names[feat_col++] = String.format("df_ft_12s%d", i);
			if(which_features[12]) feat_names[feat_col++] = String.format("df_ft_13s%d", i);
			if(which_features[13]) feat_names[feat_col++] = String.format("df_ft_14s%d", i);
			
		}
		return feat_names;
	}

}