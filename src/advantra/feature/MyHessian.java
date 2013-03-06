package advantra.feature;

import java.util.Vector;

import flanagan.math.Matrix;

import imagescience.feature.Hessian;
import imagescience.image.Aspects;
import imagescience.image.Axes;
import imagescience.image.Coordinates;
import imagescience.image.Dimensions;
import imagescience.image.FloatImage;
import imagescience.image.Image;

public class MyHessian extends Hessian {
	
	public MyHessian(){
		super();
	}
	
	public double[] eigL(final double[][] mat){
		
		Matrix aa = new Matrix(mat);
		return aa.getEigenValues();
		
	}
	
	public double[][] eigV(final double[][] mat){
		
		Matrix aa = new Matrix(mat);
		return aa.getEigenVectorsAsRows();
		
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
			
			Matrix hess_mat = new Matrix(new double[2][2]);
			
			final Coordinates coords = new Coordinates();
			
			double[] 	eig_vals;	//1x2
			double[][] 	eig_vecs; 	//2x2
			
			if (absolute) {
				for (coords.c=0; coords.c<dims.c; ++coords.c)
					for (coords.t=0; coords.t<dims.t; ++coords.t)
						for (coords.y=0; coords.y<dims.y; ++coords.y) {
							
							Hxx.get(coords,ahxx);
							Hxy.get(coords,ahxy);
							Hyy.get(coords,ahyy);
							
							for (int x=0; x<dims.x; ++x) {
								
								hess_mat.setElement(0, 0, ahxx[x]);
								hess_mat.setElement(0, 1, ahxy[x]);
								hess_mat.setElement(1, 0, ahxy[x]);
								hess_mat.setElement(1, 1, ahyy[x]);
								
								eig_vals = hess_mat.getEigenValues();
								eig_vecs = hess_mat.getEigenVectorsAsRows();
								
								double absh1, absh2;
								absh1 = Math.abs(eig_vals[0]);
								absh2 = Math.abs(eig_vals[1]);
								
								if (absh1 > absh2) {
									ahxx[x] = absh1;
									ahyy[x] = absh2;
									av1[x]  = eig_vecs[1][0];
									av2[x]  = eig_vecs[1][1];
								} else {
									ahxx[x] = absh2;
									ahyy[x] = absh1;
									av1[x]  = eig_vecs[0][0];
									av2[x]  = eig_vecs[0][1];
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
								
								hess_mat.setElement(0, 0, ahxx[x]);
								hess_mat.setElement(0, 1, ahxy[x]);
								hess_mat.setElement(1, 0, ahxy[x]);
								hess_mat.setElement(1, 1, ahyy[x]);
								
								eig_vals = hess_mat.getEigenValues();
								eig_vecs = hess_mat.getEigenVectorsAsRows();
								
								double h1, h2;
								
								h1 = eig_vals[0];
								h2 = eig_vals[1];
								
								if (h1 > h2) {
									ahxx[x] = h1;
									ahyy[x] = h2;
									av1[x]  = eig_vecs[1][0];
									av2[x]  = eig_vecs[1][1];
								} else {
									ahxx[x] = h2;
									ahyy[x] = h1;
									av1[x]  = eig_vecs[0][0];
									av2[x]  = eig_vecs[0][1];
								}
							}
							Hxx.set(coords,ahxx);
							Hyy.set(coords,ahyy);
							V1.set(coords,av1);
							V2.set(coords,av2);
						}
			}
			
			Hxx.name("L2(FLAN.)");
			Hyy.name("L1(FLAN.)");
			V1.name("V11(FLAN.)");
			V2.name("V12(FLAN.)");
			
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
			
			Matrix hess_mat = new Matrix(new double[3][3]);
			
			final Coordinates coords = new Coordinates();
			
			double[] 	eig_vals;	
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
									
									hess_mat.setElement(0, 0, ahxx[x]);
									hess_mat.setElement(0, 1, ahxy[x]);
									hess_mat.setElement(0, 2, ahxz[x]);
									hess_mat.setElement(1, 0, ahxy[x]);
									hess_mat.setElement(1, 1, ahyy[x]);
									hess_mat.setElement(1, 2, ahyz[x]);
									hess_mat.setElement(2, 0, ahxz[x]);
									hess_mat.setElement(2, 1, ahyz[x]);
									hess_mat.setElement(2, 2, ahzz[x]);
									
									eig_vals = hess_mat.getEigenValues();
									eig_vecs = hess_mat.getEigenVectorsAsRows();
									
									double absh1, absh2, absh3;
									absh1 = Math.abs(eig_vals[0]);
									absh2 = Math.abs(eig_vals[1]);
									absh3 = Math.abs(eig_vals[2]);
									double[] v1 = eig_vecs[0];
									double[] v2 = eig_vecs[1];
									double[] v3 = eig_vecs[2];
									
									if (absh2 < absh3) {  
										double tmp = absh2; absh2 = absh3; absh3 = tmp; // swap 2,3
										tmp = v2[0]; v2[0] = v3[0]; v3[0] = tmp; 		// eigenvectors
										tmp = v2[1]; v2[1] = v3[1]; v3[1] = tmp;
										tmp = v2[2]; v2[2] = v3[2]; v3[2] = tmp;
									}
									if (absh1 < absh2) { 
										double tmp = absh1; absh1 = absh2; absh2 = tmp; // swap 1,2
										tmp = v1[0]; v1[0] = v2[0]; v2[0] = tmp; 		// eigenvectors
										tmp = v1[1]; v1[1] = v2[1]; v2[1] = tmp;
										tmp = v1[2]; v1[2] = v2[2]; v2[2] = tmp;
										
										if (absh2 < absh3) { 
											double tmp1 = absh2; absh2 = absh3; absh3 = tmp1; // swap 2,3
											tmp1 = v2[0]; v2[0] = v3[0]; v3[0] = tmp1; 		  // eigenvectors
											tmp1 = v2[1]; v2[1] = v3[1]; v3[1] = tmp1;
											tmp1 = v2[2]; v2[2] = v3[2]; v3[2] = tmp1;
										}
									}
									
									ahxx[x] = absh1;
									ahyy[x] = absh2;
									ahzz[x] = absh3;
									av1[x] = v3[0];
									av2[x] = v3[1];
									av3[x] = v3[2];
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
									
									hess_mat.setElement(0, 0, ahxx[x]);
									hess_mat.setElement(0, 1, ahxy[x]);
									hess_mat.setElement(0, 2, ahxz[x]);
									hess_mat.setElement(1, 0, ahxy[x]);
									hess_mat.setElement(1, 1, ahyy[x]);
									hess_mat.setElement(1, 2, ahyz[x]);
									hess_mat.setElement(2, 0, ahxz[x]);
									hess_mat.setElement(2, 1, ahyz[x]);
									hess_mat.setElement(2, 2, ahzz[x]);
									
									eig_vals = hess_mat.getEigenValues();
									eig_vecs = hess_mat.getEigenVectorsAsRows();
									
									double h1, h2, h3;
									h1 = eig_vals[0];
									h2 = eig_vals[1];
									h3 = eig_vals[2];
									double[] v1 = eig_vecs[0];
									double[] v2 = eig_vecs[1];
									double[] v3 = eig_vecs[2];
									
									//sort them
									if (h2 < h3) { 
										double tmp = h2; h2 = h3; h3 = tmp; 
										tmp = v2[0]; v2[0] = v3[0]; v3[0] = tmp; 		
										tmp = v2[1]; v2[1] = v3[1]; v3[1] = tmp;
										tmp = v2[2]; v2[2] = v3[2]; v3[2] = tmp;
									}
									if (h1 < h2) { 
										double tmp1 = h1; h1 = h2; h2 = tmp1;
										tmp1 = v1[0]; v1[0] = v2[0]; v2[0] = tmp1; 		
										tmp1 = v1[1]; v1[1] = v2[1]; v2[1] = tmp1;
										tmp1 = v1[2]; v1[2] = v2[2]; v2[2] = tmp1;
										if (h2 < h3) { 
											double tmp2 = h2; h2 = h3; h3 = tmp2; 
											tmp2 = v2[0]; v2[0] = v3[0]; v3[0] = tmp2; 
											tmp2 = v2[1]; v2[1] = v3[1]; v3[1] = tmp2;
											tmp2 = v2[2]; v2[2] = v3[2]; v3[2] = tmp2;
										}
									}
									
									ahxx[x] = h1;
									ahyy[x] = h2;
									ahzz[x] = h3;
									av1[x] = v3[0];
									av2[x] = v3[1];
									av3[x] = v3[2];
								}
								Hxx.set(coords,ahxx);
								Hyy.set(coords,ahyy);
								Hzz.set(coords,ahzz);
								V1.set(coords, av1);
								V2.set(coords, av2);
								V3.set(coords, av3);
							}
			}
			
			Hxx.name("L3(FLAN.)");
			Hyy.name("L2(FLAN.)");
			Hzz.name("L1(FLAN.)");
			V1.name("V11(FLAN.)");
			V2.name("V12(FLAN.)");
			V3.name("V13(FLAN.)");
			
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
