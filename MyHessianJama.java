package advantra.feature;

import java.util.Vector;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;


import imagescience.feature.Hessian;
import imagescience.image.Aspects;
import imagescience.image.Axes;
import imagescience.image.Coordinates;
import imagescience.image.Dimensions;
import imagescience.image.FloatImage;
import imagescience.image.Image;

public class MyHessianJama extends Hessian {
	
	public MyHessianJama(){
		super();
	}
	
	public double[] eigL(final double[][] mat){
		
		Matrix aa = new Matrix(mat);
		EigenvalueDecomposition ed = new EigenvalueDecomposition(aa);
		Matrix readD = ed.getD();
		double l1 = readD.get(0, 0);
		double l2 = readD.get(1, 1);
		double l3 = readD.get(2, 2);
		return new double[]{l1, l2, l3};
		
	}
	
	public double[][] eigV(final double[][] mat){
		
		Matrix aa = new Matrix(mat);
		EigenvalueDecomposition ed = new EigenvalueDecomposition(aa);
		Matrix bb = ed.getV();
		return bb.getArrayCopy();
		
	}
	
	public Vector<Image> eigs(final Image image, final double scale, final boolean absolute) {
		
		if (scale <= 0) throw new IllegalArgumentException("Smoothing scale less than or equal to 0");
		
		final Dimensions dims = image.dimensions();
		
		final Aspects asps = image.aspects();
		if (asps.x <= 0) throw new IllegalStateException("Aspect-ratio value in x-dimension less than or equal to 0");
		if (asps.y <= 0) throw new IllegalStateException("Aspect-ratio value in y-dimension less than or equal to 0");
		if (asps.z <= 0) throw new IllegalStateException("Aspect-ratio value in z-dimension less than or equal to 0");
		
		final Image smoothImage = (image instanceof FloatImage) ? image : new FloatImage(image);
		Vector<Image> eigenimages = null;
		
		if (dims.z == 1) { // 2D case
			
			final Image Hxx = differentiator.run(smoothImage.duplicate(),scale,2,0,0); 	// lambda2 will be here (the higher one)
			final Image Hxy = differentiator.run(smoothImage.duplicate(),scale,1,1,0); 	// lambda1 will be here (smaller one)
			final Image Hyy = differentiator.run(smoothImage,scale,0,2,0);
			final Image V1  = new FloatImage(smoothImage.dimensions());					// v1(1) will be here (first coord)
			final Image V2  = new FloatImage(smoothImage.dimensions());					// v1(2) will be here (second coord)
			
			Hxx.axes(Axes.X); Hxy.axes(Axes.X); Hyy.axes(Axes.X);
			V1.axes(Axes.X);  V2.axes(Axes.X);	
			
			// arrays to facilitate the value filling-up of the Image instances
			final double[] ahxx = new double[dims.x];
			final double[] ahxy = new double[dims.x];
			final double[] ahyy = new double[dims.x];
			
			final double[] av1  = new double[dims.x];
			final double[] av2  = new double[dims.x];
			
			
			double[][] a_hess_mat = new double[2][2];
			
			final Coordinates coords = new Coordinates();
			
			double[][] 	eig_vecs;
			
			if (absolute) {
				for (coords.c=0; coords.c<dims.c; ++coords.c)
					for (coords.t=0; coords.t<dims.t; ++coords.t)
						for (coords.y=0; coords.y<dims.y; ++coords.y) {
							
							Hxx.get(coords,ahxx);
							Hxy.get(coords,ahxy);
							Hyy.get(coords,ahyy);
							
							for (int x=0; x<dims.x; ++x) {
								
								a_hess_mat[0][0] = ahxx[x];
								a_hess_mat[0][1] = ahxy[x];
								a_hess_mat[1][0] = ahxy[x];
								a_hess_mat[1][1] = ahyy[x];
								
								EigenvalueDecomposition ed = new EigenvalueDecomposition(Matrix.constructWithCopy(a_hess_mat));
								
								eig_vecs = ed.getV().getArrayCopy();
								
								double absh1, absh2;
								Matrix readD = ed.getD();
								absh1 = Math.abs(readD.get(0, 0));
								absh2 = Math.abs(readD.get(1, 1));
								
								
								if (absh1 > absh2) {
									ahxx[x] = absh1;
									ahyy[x] = absh2;
									av1[x]  = eig_vecs[0][1];
									av2[x]  = eig_vecs[1][1];
								} else {
									ahxx[x] = absh2;
									ahyy[x] = absh1;
									av1[x]  = eig_vecs[0][0];
									av2[x]  = eig_vecs[1][0];
								}
							}
							
							Hxx.set(coords,ahxx);
							Hyy.set(coords,ahyy);
							V1.set(coords,av1);
							V2.set(coords,av2);
						}
			} else {
				for (coords.c=0; coords.c<dims.c; ++coords.c)
					for (coords.t=0; coords.t<dims.t; ++coords.t)
						for (coords.y=0; coords.y<dims.y; ++coords.y) {
							
							Hxx.get(coords,ahxx);
							Hxy.get(coords,ahxy);
							Hyy.get(coords,ahyy);
							
							for (int x=0; x<dims.x; ++x) {
								
								a_hess_mat[0][0] = ahxx[x];
								a_hess_mat[0][1] = ahxy[x];
								a_hess_mat[1][0] = ahxy[x];
								a_hess_mat[1][1] = ahyy[x];
								
								EigenvalueDecomposition ed = new EigenvalueDecomposition(Matrix.constructWithCopy(a_hess_mat));
								
								eig_vecs = ed.getV().getArrayCopy();
								
								double h1, h2;
								Matrix readD = ed.getD();
								h1 = readD.get(0, 0);
								h2 = readD.get(1, 1);
								
								if (h1 > h2) {
									ahxx[x] = h1;
									ahyy[x] = h2;
									av1[x]  = eig_vecs[0][1];
									av2[x]  = eig_vecs[1][1];
								} else {
									ahxx[x] = h2;
									ahyy[x] = h1;
									av1[x]  = eig_vecs[0][0];
									av2[x]  = eig_vecs[1][0];
								}
							}
							Hxx.set(coords,ahxx);
							Hyy.set(coords,ahyy);
							V1.set(coords,av1);
							V2.set(coords,av2);
						}
			}
			
			Hxx.name("L2(JAMA)");
			Hyy.name("L1(JAMA)");
			V1.name("V11(JAMA)");
			V2.name("V12(JAMA)");
			
			Hxx.aspects(asps.duplicate());
			Hyy.aspects(asps.duplicate());
			V1.aspects(asps.duplicate());
			V2.aspects(asps.duplicate());
			
			eigenimages = new Vector<Image>(4); // lambda1,2 and v1(1,2)
			eigenimages.add(Hxx); 
			eigenimages.add(Hyy);
			eigenimages.add(V1);
			eigenimages.add(V2);
			
		} else { // 3D case
			
			// Compute Hessian components:
			final Image Hxx = differentiator.run(smoothImage.duplicate(),scale,2,0,0);
			final Image Hxy = differentiator.run(smoothImage.duplicate(),scale,1,1,0);
			final Image Hxz = differentiator.run(smoothImage.duplicate(),scale,1,0,1);
			final Image Hyy = differentiator.run(smoothImage.duplicate(),scale,0,2,0);
			final Image Hyz = differentiator.run(smoothImage.duplicate(),scale,0,1,1);
			final Image Hzz = differentiator.run(smoothImage,scale,0,0,2);
			final Image V1  = new FloatImage(smoothImage.dimensions());					// v1(1) will be here (first coord)
			final Image V2  = new FloatImage(smoothImage.dimensions());					// v1(2) will be here (second coord)
			final Image V3  = new FloatImage(smoothImage.dimensions());					// v1(2) will be here (third coord)
			
			// Compute eigenimages (Hxx, Hyy, Hzz are reused to save memory):
			Hxx.axes(Axes.X); Hxy.axes(Axes.X); Hxz.axes(Axes.X);
			Hyy.axes(Axes.X); Hyz.axes(Axes.X); Hzz.axes(Axes.X);
			V1.axes(Axes.X);  V2.axes(Axes.X);  V3.axes(Axes.X); 
			
			final double[] ahxx = new double[dims.x];
			final double[] ahxy = new double[dims.x];
			final double[] ahxz = new double[dims.x];
			final double[] ahyy = new double[dims.x];
			final double[] ahyz = new double[dims.x];
			final double[] ahzz = new double[dims.x];
			final double[] av1  = new double[dims.x];
			final double[] av2  = new double[dims.x];
			final double[] av3  = new double[dims.x];
			
			
			double[][] a_hess_mat = new double[3][3];
			
			final Coordinates coords = new Coordinates();
			
			double[][] 	eig_vecs;
			
			if (absolute) {
				for (coords.c=0; coords.c<dims.c; ++coords.c)
					for (coords.t=0; coords.t<dims.t; ++coords.t)
						for (coords.z=0; coords.z<dims.z; ++coords.z)
							for (coords.y=0; coords.y<dims.y; ++coords.y) {
								
								Hxx.get(coords,ahxx);
								Hxy.get(coords,ahxy);
								Hxz.get(coords,ahxz);
								Hyy.get(coords,ahyy);
								Hyz.get(coords,ahyz);
								Hzz.get(coords,ahzz);
								
								for (int x=0; x<dims.x; ++x) {
									
									a_hess_mat[0][0] = ahxx[x];
									a_hess_mat[0][1] = ahxy[x];
									a_hess_mat[0][2] = ahxz[x];
									a_hess_mat[1][0] = ahxy[x];
									a_hess_mat[1][1] = ahyy[x];
									a_hess_mat[1][2] = ahyz[x];
									a_hess_mat[1][0] = ahxz[x];
									a_hess_mat[1][1] = ahyz[x];
									a_hess_mat[1][2] = ahzz[x];
									
									EigenvalueDecomposition ed = new EigenvalueDecomposition(Matrix.constructWithCopy(a_hess_mat));
									
									eig_vecs = ed.getV().getArrayCopy();
									
									double absh1, absh2, absh3;
									Matrix readD = ed.getD();
									absh1 = Math.abs(readD.get(0, 0));
									absh2 = Math.abs(readD.get(1, 1));
									absh3 = Math.abs(readD.get(2, 2));
									
									if (absh2 < absh3) {  
										double tmp = absh2; absh2 = absh3; absh3 = tmp; // swap 2,3
										tmp = eig_vecs[0][1]; eig_vecs[0][1] = eig_vecs[0][2]; eig_vecs[0][2] = tmp; 		
										tmp = eig_vecs[1][1]; eig_vecs[1][1] = eig_vecs[1][2]; eig_vecs[1][2] = tmp;
										tmp = eig_vecs[2][1]; eig_vecs[2][1] = eig_vecs[2][2]; eig_vecs[2][2] = tmp;
									}
									if (absh1 < absh2) { 
										double tmp = absh1; absh1 = absh2; absh2 = tmp; // swap 1,2
										tmp = eig_vecs[0][0]; eig_vecs[0][0] = eig_vecs[0][1]; eig_vecs[0][1] = tmp; 
										tmp = eig_vecs[1][0]; eig_vecs[1][0] = eig_vecs[1][1]; eig_vecs[1][1] = tmp;
										tmp = eig_vecs[2][0]; eig_vecs[2][0] = eig_vecs[2][1]; eig_vecs[2][1] = tmp;
										
										if (absh2 < absh3) { 
											double tmp1 = absh2; absh2 = absh3; absh3 = tmp1; // swap 2,3
											tmp1 = eig_vecs[0][1]; eig_vecs[0][1] = eig_vecs[0][2]; eig_vecs[0][2] = tmp1;
											tmp1 = eig_vecs[1][1]; eig_vecs[1][1] = eig_vecs[1][2]; eig_vecs[1][2] = tmp1;
											tmp1 = eig_vecs[2][1]; eig_vecs[2][1] = eig_vecs[2][2]; eig_vecs[2][2] = tmp1;
										}
									}
									
									ahxx[x] = absh1;
									ahyy[x] = absh2;
									ahzz[x] = absh3;
									av1[x] = eig_vecs[0][2];
									av2[x] = eig_vecs[1][2];
									av3[x] = eig_vecs[2][2];
								}
								
								Hxx.set(coords,ahxx);
								Hyy.set(coords,ahyy);
								Hzz.set(coords,ahzz);
								V1.set(coords, av1);
								V2.set(coords, av2);
								V3.set(coords, av3);
							}
			} else {
				for (coords.c=0; coords.c<dims.c; ++coords.c)
					for (coords.t=0; coords.t<dims.t; ++coords.t)
						for (coords.z=0; coords.z<dims.z; ++coords.z)
							for (coords.y=0; coords.y<dims.y; ++coords.y) {
								
								Hxx.get(coords,ahxx);
								Hxy.get(coords,ahxy);
								Hxz.get(coords,ahxz);
								Hyy.get(coords,ahyy);
								Hyz.get(coords,ahyz);
								Hzz.get(coords,ahzz);
								
								for (int x=0; x<dims.x; ++x) {
									
									a_hess_mat[0][0] = ahxx[x];
									a_hess_mat[0][1] = ahxy[x];
									a_hess_mat[0][2] = ahxz[x];
									a_hess_mat[1][0] = ahxy[x];
									a_hess_mat[1][1] = ahyy[x];
									a_hess_mat[1][2] = ahyz[x];
									a_hess_mat[1][0] = ahxz[x];
									a_hess_mat[1][1] = ahyz[x];
									a_hess_mat[1][2] = ahzz[x];
									
									EigenvalueDecomposition ed = new EigenvalueDecomposition(Matrix.constructWithCopy(a_hess_mat));
									
									eig_vecs = ed.getV().getArrayCopy();
									
									double h1, h2, h3;
									Matrix readD = ed.getD();
									h1 = readD.get(0, 0);
									h2 = readD.get(1, 1);
									h3 = readD.get(2, 2);
									
									//sort them
									if (h2 < h3) { 
										double tmp = h2; h2 = h3; h3 = tmp; 
										tmp = eig_vecs[0][1]; eig_vecs[0][1] = eig_vecs[0][2]; eig_vecs[0][2] = tmp; 		
										tmp = eig_vecs[1][1]; eig_vecs[1][1] = eig_vecs[1][2]; eig_vecs[1][2] = tmp;
										tmp = eig_vecs[2][1]; eig_vecs[2][1] = eig_vecs[2][2]; eig_vecs[2][2] = tmp;
									}
									if (h1 < h2) { 
										double tmp1 = h1; h1 = h2; h2 = tmp1;
										tmp1 = eig_vecs[0][0]; eig_vecs[0][0] = eig_vecs[0][1]; eig_vecs[0][1] = tmp1; 		
										tmp1 = eig_vecs[1][0]; eig_vecs[1][0] = eig_vecs[1][1]; eig_vecs[1][1] = tmp1;
										tmp1 = eig_vecs[2][0]; eig_vecs[2][0] = eig_vecs[2][1]; eig_vecs[2][1] = tmp1;
										if (h2 < h3) { 
											double tmp2 = h2; h2 = h3; h3 = tmp2; 
											tmp2 = eig_vecs[0][1]; eig_vecs[0][1] = eig_vecs[0][2]; eig_vecs[0][2] = tmp2; 
											tmp2 = eig_vecs[1][1]; eig_vecs[1][1] = eig_vecs[1][2]; eig_vecs[1][2] = tmp2;
											tmp2 = eig_vecs[2][1]; eig_vecs[2][1] = eig_vecs[2][2]; eig_vecs[2][2] = tmp2;
										}
									}
									
									ahxx[x] = h1;
									ahyy[x] = h2;
									ahzz[x] = h3;
									av1[x] = eig_vecs[0][2];
									av2[x] = eig_vecs[1][2];
									av3[x] = eig_vecs[2][2];
								}
								Hxx.set(coords,ahxx);
								Hyy.set(coords,ahyy);
								Hzz.set(coords,ahzz);
								V1.set(coords, av1);
								V2.set(coords, av2);
								V3.set(coords, av3);
							}
			}
			
			Hxx.name("L3(JAMA)");
			Hyy.name("L2(JAMA)");
			Hzz.name("L1(JAMA)");
			V1.name("V11(JAMA)");
			V2.name("V12(JAMA)");
			V3.name("V13(JAMA)");
			
			Hxx.aspects(asps.duplicate());
			Hyy.aspects(asps.duplicate());
			Hzz.aspects(asps.duplicate());
			V1.aspects(asps.duplicate());
			V2.aspects(asps.duplicate());
			V3.aspects(asps.duplicate());
			
			eigenimages = new Vector<Image>(6);
			eigenimages.add(Hxx);
			eigenimages.add(Hyy);
			eigenimages.add(Hzz);
			eigenimages.add(V1);
			eigenimages.add(V2);
			eigenimages.add(V3);
		}

		return eigenimages;
	}
	
}
/*
private double[] refineStartPoint(double[] point3d, int range){
	
	// expand Sphere around the seed
	Sphere startSphere	= new Sphere(point3d[0], point3d[1], point3d[2], range);
	
	// take voxel values & locations
	int[][] roi_coord = new int[3][	Sphere.numberOfVoxInSphere(range)];	
	int[]   roi_vals  = new int[	Sphere.numberOfVoxInSphere(range)];
	
	int count_sphere_voxels = startSphere.extractVox(traced_img, roi_coord, roi_vals);
	
	double[] momts = new double[9]; // allocate space to store moments from extracted roi		
	// extract moments: momts[0], momts[1], momts[2] define CENTROID
	double sum_of_intensities = Moments.extract_moments_3D(roi_coord, roi_vals, count_sphere_voxels, momts);		
	
	double[] refined_point3d = new double[3];
	if(sum_of_intensities>0){
		refined_point3d[0] = momts[0];
		refined_point3d[1] = momts[1];
		refined_point3d[2] = momts[2];
	}
	else{
		System.out.println("Point was not refined... sum of the intensities was "+sum_of_intensities+" ...");
		refined_point3d[0] = point3d[0];
		refined_point3d[1] = point3d[1];
		refined_point3d[2] = point3d[2];
	}

	
	return refined_point3d;

}
*/
