package advantra.feature;

import java.util.Vector;

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
	
	// duplicate the code of the private members of the super-class
	// "A subclass does not inherit the private members of its parent class"
	private static final double TWOPI = 2*Math.PI;
	
//	private void logstatus(final String s) {
//		
//		messenger.log(s);
//		messenger.status(s+"...");
//	}
	
	// add modified run() method as a new class method
	public Vector<Image> eig(final Image image, final double scale, final boolean absolute) {
		
		//messenger.log(ImageScience.prelude()+"Hessian");
		
		//final Timer timer = new Timer();
		//timer.messenger.log(messenger.log());
		//timer.start();
		
		// Initialize:
		//messenger.log("Checking arguments");
		if (scale <= 0) throw new IllegalArgumentException("Smoothing scale less than or equal to 0");
		
		final Dimensions dims = image.dimensions();
		//messenger.log("Input image dimensions: (x,y,z,t,c) = ("+dims.x+","+dims.y+","+dims.z+","+dims.t+","+dims.c+")");
		
		final Aspects asps = image.aspects();
		//messenger.log("Element aspect-ratios: ("+asps.x+","+asps.y+","+asps.z+","+asps.t+","+asps.c+")");
		if (asps.x <= 0) throw new IllegalStateException("Aspect-ratio value in x-dimension less than or equal to 0");
		if (asps.y <= 0) throw new IllegalStateException("Aspect-ratio value in y-dimension less than or equal to 0");
		if (asps.z <= 0) throw new IllegalStateException("Aspect-ratio value in z-dimension less than or equal to 0");
		
		final Image smoothImage = (image instanceof FloatImage) ? image : new FloatImage(image);
		Vector<Image> eigenimages = null;
		final String name = image.name();
		
		//differentiator.messenger.log(messenger.log());
		//differentiator.progressor.parent(progressor);
		
		// Compute Hessian matrix and eigenimages:
		if (dims.z == 1) { // 2D case
			
			// final double[] pls = {0, 0.32, 0.64, 0.96, 1}; int pl = 0;
			
			// Compute Hessian components:
			// logstatus("Computing Hxx"); progressor.range(pls[pl],pls[++pl]);
			final Image Hxx = differentiator.run(smoothImage.duplicate(),scale,2,0,0);
			// logstatus("Computing Hxy"); progressor.range(pls[pl],pls[++pl]);
			final Image Hxy = differentiator.run(smoothImage.duplicate(),scale,1,1,0);
			// logstatus("Computing Hyy"); progressor.range(pls[pl],pls[++pl]);
			final Image Hyy = differentiator.run(smoothImage,scale,0,2,0);
			
			// Compute eigenimages (Hxx and Hyy are reused to save memory):
			// logstatus("Computing eigenimages");
			//progressor.steps(dims.c*dims.t*dims.y);
			//progressor.range(pls[pl],pls[++pl]);
			Hxx.axes(Axes.X); Hxy.axes(Axes.X); Hyy.axes(Axes.X);
			final double[] ahxx = new double[dims.x];
			final double[] ahxy = new double[dims.x];
			final double[] ahyy = new double[dims.x];
			final Coordinates coords = new Coordinates();
			
			//progressor.start();
			if (absolute) {
				//messenger.log("Comparing and storing absolute eigenvalues");
				for (coords.c=0; coords.c<dims.c; ++coords.c)
					for (coords.t=0; coords.t<dims.t; ++coords.t)
						for (coords.y=0; coords.y<dims.y; ++coords.y) {
							Hxx.get(coords,ahxx);
							Hxy.get(coords,ahxy);
							Hyy.get(coords,ahyy);
							for (int x=0; x<dims.x; ++x) {
								final double b = -(ahxx[x] + ahyy[x]);
								final double c = ahxx[x]*ahyy[x] - ahxy[x]*ahxy[x];
								final double q = -0.5*(b + (b < 0 ? -1 : 1)*Math.sqrt(b*b - 4*c));
								double absh1, absh2;
								if (q == 0) {
									absh1 = 0;
									absh2 = 0;
								} else {
									absh1 = Math.abs(q);
									absh2 = Math.abs(c/q);
								}
								if (absh1 > absh2) {
									ahxx[x] = absh1;
									ahyy[x] = absh2;
								} else {
									ahxx[x] = absh2;
									ahyy[x] = absh1;
								}
							}
							Hxx.set(coords,ahxx);
							Hyy.set(coords,ahyy);
							//progressor.step();
						}
			} else {
				//messenger.log("Comparing and storing actual eigenvalues");
				for (coords.c=0; coords.c<dims.c; ++coords.c)
					for (coords.t=0; coords.t<dims.t; ++coords.t)
						for (coords.y=0; coords.y<dims.y; ++coords.y) {
							Hxx.get(coords,ahxx);
							Hxy.get(coords,ahxy);
							Hyy.get(coords,ahyy);
							for (int x=0; x<dims.x; ++x) {
								final double b = -(ahxx[x] + ahyy[x]);
								final double c = ahxx[x]*ahyy[x] - ahxy[x]*ahxy[x];
								final double q = -0.5*(b + (b < 0 ? -1 : 1)*Math.sqrt(b*b - 4*c));
								double h1, h2;
								if (q == 0) {
									h1 = 0;
									h2 = 0;
								} else {
									h1 = q;
									h2 = c/q;
								}
								if (h1 > h2) {
									ahxx[x] = h1;
									ahyy[x] = h2;
								} else {
									ahxx[x] = h2;
									ahyy[x] = h1;
								}
							}
							Hxx.set(coords,ahxx);
							Hyy.set(coords,ahyy);
							//progressor.step();
						}
			}
			//progressor.stop();
			
			Hxx.name(name+" largest Hessian eigenvalues");
			Hyy.name(name+" smallest Hessian eigenvalues");
			
			Hxx.aspects(asps.duplicate());
			Hyy.aspects(asps.duplicate());
			
			eigenimages = new Vector<Image>(6); // lambda1,2 and v1,2
			eigenimages.add(Hxx); 
			eigenimages.add(Hyy);
			// add vector v1 : 2 layers
			
			// add vector v2 : 2 layers
			
			
		} else { // 3D case
			
			//final double[] pls = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1}; int pl = 0;
			
			// Compute Hessian components:
			//logstatus("Computing Hxx"); progressor.range(pls[pl],pls[++pl]);
			final Image Hxx = differentiator.run(smoothImage.duplicate(),scale,2,0,0);
			//logstatus("Computing Hxy"); progressor.range(pls[pl],pls[++pl]);
			final Image Hxy = differentiator.run(smoothImage.duplicate(),scale,1,1,0);
			//logstatus("Computing Hxz"); progressor.range(pls[pl],pls[++pl]);
			final Image Hxz = differentiator.run(smoothImage.duplicate(),scale,1,0,1);
			//logstatus("Computing Hyy"); progressor.range(pls[pl],pls[++pl]);
			final Image Hyy = differentiator.run(smoothImage.duplicate(),scale,0,2,0);
			//logstatus("Computing Hyz"); progressor.range(pls[pl],pls[++pl]);
			final Image Hyz = differentiator.run(smoothImage.duplicate(),scale,0,1,1);
			//logstatus("Computing Hzz"); progressor.range(pls[pl],pls[++pl]);
			final Image Hzz = differentiator.run(smoothImage,scale,0,0,2);
			
			// Compute eigenimages (Hxx, Hyy, Hzz are reused to save memory):
			//logstatus("Computing eigenimages");
			//progressor.steps(dims.c*dims.t*dims.z*dims.y);
			//progressor.range(pls[pl],pls[++pl]);
			Hxx.axes(Axes.X); Hxy.axes(Axes.X); Hxz.axes(Axes.X);
			Hyy.axes(Axes.X); Hyz.axes(Axes.X); Hzz.axes(Axes.X);
			final double[] ahxx = new double[dims.x];
			final double[] ahxy = new double[dims.x];
			final double[] ahxz = new double[dims.x];
			final double[] ahyy = new double[dims.x];
			final double[] ahyz = new double[dims.x];
			final double[] ahzz = new double[dims.x];
			final Coordinates coords = new Coordinates();
			
			//progressor.start();
			if (absolute) {
				messenger.log("Comparing and storing absolute eigenvalues");
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
									final double fhxx = ahxx[x];
									final double fhxy = ahxy[x];
									final double fhxz = ahxz[x];
									final double fhyy = ahyy[x];
									final double fhyz = ahyz[x];
									final double fhzz = ahzz[x];
									final double a = -(fhxx + fhyy + fhzz);
									final double b = fhxx*fhyy + fhxx*fhzz + fhyy*fhzz - fhxy*fhxy - fhxz*fhxz - fhyz*fhyz;
									final double c = fhxx*(fhyz*fhyz - fhyy*fhzz) + fhyy*fhxz*fhxz + fhzz*fhxy*fhxy - 2*fhxy*fhxz*fhyz;
									final double q = (a*a - 3*b)/9;
									final double r = (a*a*a - 4.5*a*b + 13.5*c)/27;
									final double sqrtq = (q > 0) ? Math.sqrt(q) : 0;
									final double sqrtq3 = sqrtq*sqrtq*sqrtq;
									double absh1, absh2, absh3;
									if (sqrtq3 == 0) {
										absh1 = 0;
										absh2 = 0;
										absh3 = 0;
									} else {
										final double rsqq3 = r/sqrtq3;
										final double angle = (rsqq3*rsqq3 <= 1) ? Math.acos(rsqq3) : Math.acos(rsqq3 < 0 ? -1 : 1);
										absh1 = Math.abs(-2*sqrtq*Math.cos(angle/3) - a/3);
										absh2 = Math.abs(-2*sqrtq*Math.cos((angle + TWOPI)/3) - a/3);
										absh3 = Math.abs(-2*sqrtq*Math.cos((angle - TWOPI)/3) - a/3);
									}
									if (absh2 < absh3) { final double tmp = absh2; absh2 = absh3; absh3 = tmp; }
									if (absh1 < absh2) { final double tmp1 = absh1; absh1 = absh2; absh2 = tmp1;
									if (absh2 < absh3) { final double tmp2 = absh2; absh2 = absh3; absh3 = tmp2; }}
									ahxx[x] = absh1;
									ahyy[x] = absh2;
									ahzz[x] = absh3;
								}
								Hxx.set(coords,ahxx);
								Hyy.set(coords,ahyy);
								Hzz.set(coords,ahzz);
								//progressor.step();
							}
			} else {
				messenger.log("Comparing and storing actual eigenvalues");
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
									final double fhxx = ahxx[x];
									final double fhxy = ahxy[x];
									final double fhxz = ahxz[x];
									final double fhyy = ahyy[x];
									final double fhyz = ahyz[x];
									final double fhzz = ahzz[x];
									final double a = -(fhxx + fhyy + fhzz);
									final double b = fhxx*fhyy + fhxx*fhzz + fhyy*fhzz - fhxy*fhxy - fhxz*fhxz - fhyz*fhyz;
									final double c = fhxx*(fhyz*fhyz - fhyy*fhzz) + fhyy*fhxz*fhxz + fhzz*fhxy*fhxy - 2*fhxy*fhxz*fhyz;
									final double q = (a*a - 3*b)/9;
									final double r = (a*a*a - 4.5*a*b + 13.5*c)/27;
									final double sqrtq = (q > 0) ? Math.sqrt(q) : 0;
									final double sqrtq3 = sqrtq*sqrtq*sqrtq;
									double h1, h2, h3;
									if (sqrtq3 == 0) {
										h1 = 0;
										h2 = 0;
										h3 = 0;
									} else {
										final double rsqq3 = r/sqrtq3;
										final double angle = (rsqq3*rsqq3 <= 1) ? Math.acos(rsqq3) : Math.acos(rsqq3 < 0 ? -1 : 1);
										h1 = -2*sqrtq*Math.cos(angle/3) - a/3;
										h2 = -2*sqrtq*Math.cos((angle + TWOPI)/3) - a/3;
										h3 = -2*sqrtq*Math.cos((angle - TWOPI)/3) - a/3;
									}
									if (h2 < h3) { final double tmp = h2; h2 = h3; h3 = tmp; }
									if (h1 < h2) { final double tmp1 = h1; h1 = h2; h2 = tmp1;
									if (h2 < h3) { final double tmp2 = h2; h2 = h3; h3 = tmp2; }}
									ahxx[x] = h1;
									ahyy[x] = h2;
									ahzz[x] = h3;
								}
								Hxx.set(coords,ahxx);
								Hyy.set(coords,ahyy);
								Hzz.set(coords,ahzz);
								//progressor.step();
							}
			}
			//progressor.stop();
			
			Hxx.name(name+" largest Hessian eigenvalues");
			Hyy.name(name+" middle Hessian eigenvalues");
			Hzz.name(name+" smallest Hessian eigenvalues");
			
			Hxx.aspects(asps.duplicate());
			Hyy.aspects(asps.duplicate());
			Hzz.aspects(asps.duplicate());
			
			eigenimages = new Vector<Image>(3);
			eigenimages.add(Hxx);
			eigenimages.add(Hyy);
			eigenimages.add(Hzz);
		}
		
		messenger.status("");
		
		//timer.stop();
		
		return eigenimages;
	}
	
	
}
