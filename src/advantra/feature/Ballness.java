package advantra.feature;

import imagescience.feature.Hessian;
import imagescience.image.Axes;
import imagescience.image.Coordinates;
import imagescience.image.Dimensions;
import imagescience.image.FloatImage;
import imagescience.image.Image;

import java.util.Vector;

public class Ballness {
	
	public Ballness(){}
	
	public Vector<Image> run(Image input2D, double[] sigmas) {
		
		Vector<Image> out = new Vector<Image>(sigmas.length);
		
		Dimensions dim = input2D.dimensions();
		double[] aL1 = new double[dim.x];
		//double[] aL2 = new double[dim.x];
		double[] aBness = new double[dim.x];
		
		for (int i = 0; i < sigmas.length; i++) {
			
			Image Bness = new FloatImage(input2D.dimensions()); Bness.axes(Axes.X);
			
			Hessian hs = new Hessian();
			Vector<Image> L = hs.run(input2D.duplicate(), sigmas[i], false); 
			//L.get(0).axes(Axes.X);
			L.get(1).axes(Axes.X);
			
			Coordinates coords 	= new Coordinates();
			for (coords.y=0; coords.y<dim.y; ++coords.y) {
				//L.get(0).get(coords,aL1);
				L.get(1).get(coords,aL1);
				for (int x=0; x<dim.x; ++x){
					aBness[x] = Math.abs(aL1[x]);// * aL2[x];
				}
				Bness.set(coords, aBness);
			}
			
			out.add(Bness);
			
		}
		return out;
	}

}
