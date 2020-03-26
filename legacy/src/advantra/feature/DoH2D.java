package advantra.feature;

import imagescience.feature.Hessian;
import imagescience.image.Axes;
import imagescience.image.Coordinates;
import imagescience.image.Dimensions;
import imagescience.image.FloatImage;
import imagescience.image.Image;

import java.util.Vector;

public class DoH2D {
	
	public static Vector<Image> calculateArray(Image input2D, double[] sigmas) {
		
		Vector<Image> out = new Vector<Image>(sigmas.length);
		
		Dimensions dim = input2D.dimensions();
		double[] aL1 = new double[dim.x];
		double[] aL2 = new double[dim.x];
		double[] aDoh = new double[dim.x];
		
		for (int i = 0; i < sigmas.length; i++) {
			
			Image Doh = new FloatImage(input2D.dimensions()); Doh.axes(Axes.X);
			
			Hessian hs = new Hessian();
			Vector<Image> L = hs.run(input2D.duplicate(), sigmas[i], false); 
			L.get(0).axes(Axes.X);
			L.get(1).axes(Axes.X);
			
			Coordinates coords 	= new Coordinates();
			for (coords.y=0; coords.y<dim.y; ++coords.y) {
				L.get(0).get(coords,aL1); // largest
				L.get(1).get(coords,aL2);
				for (int x=0; x<dim.x; ++x){
					aDoh[x] = aL1[x] * aL2[x];
				}
				Doh.set(coords, aDoh);
			}
			
			out.add(Doh);
			
		}
		return out;
	}
	
	public static Image calculateImg(Image input2D, double[] sigmas) {
		
		Dimensions dim = input2D.dimensions();
		
		Image out = new FloatImage(new Dimensions(dim.x, dim.y, sigmas.length));
		out.axes(Axes.X);
		
		double[] aL1 = new double[dim.x];
		double[] aL2 = new double[dim.x];
		double[] aDoh = new double[dim.x];
		
		Hessian hs = new Hessian();
		
		Coordinates coords 	= new Coordinates();
		
		for (int i = 0; i < sigmas.length; i++) {
			
			Image Doh = new FloatImage(input2D.dimensions()); Doh.axes(Axes.X);
			
			Vector<Image> L = hs.run(input2D.duplicate(), sigmas[i], false); 
			L.get(0).axes(Axes.X); // largest
			L.get(1).axes(Axes.X);
			
			for (coords.y=0; coords.y<dim.y; ++coords.y) {
				L.get(0).get(coords,aL1); // largest
				L.get(1).get(coords,aL2);
				for (int x=0; x<dim.x; ++x){
					aDoh[x] = aL1[x] * aL2[x];
				}
				
				coords.z = i;
				out.set(coords, aDoh);
				coords.z = 0;
				
			}
			
		}
		
		out.name("DoH_s_");
		
		return out;
		
	}

}
