package advantra.feature;

import java.util.Vector;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.PlugIn;
import imagescience.image.FloatImage;
import imagescience.image.Image;
import weka.core.Utils;

public class Test implements PlugIn {

	public static void main(String[] args){
		
		
		
		
	}

	public void run(String arg0) {
		
		double[][] mat = new double[][]
				{{ 5, 	-4, 	 2},
				 {-4, 	 5, 	 2},
				 { 2,  	 2, 	-1}
				};
			
			MyHessianFlanagan 	myhf = new MyHessianFlanagan();
			MyHessianJama		myhj = new MyHessianJama();
			
			System.out.println("FLANAGAN SAYS:");
			
			double[] a = myhf.eigL(mat);
			System.out.println(Utils.arrayToString(a));
			
			double[][] aa = myhf.eigV(mat);
			System.out.println(Utils.arrayToString(aa));
			
			System.out.println("JAMA SAYS:");
			
			a = myhf.eigL(mat);
			System.out.println(Utils.arrayToString(a));
			
			aa = myhf.eigV(mat);
			System.out.println(Utils.arrayToString(aa));
			
			ImagePlus ip = IJ.openImage();
			ImageStack is = new ImageStack(ip.getWidth(), ip.getHeight());
			for (int i = 0; i < ip.getStackSize(); i++) {
				is.addSlice(ip.getStack().getProcessor(i+1).convertToFloat());
			}
			
			ImagePlus ip_float = new ImagePlus("converted", is);
			ip_float.show();
			Image img = new FloatImage(ip_float);
			
			long t1, t2;
			
			int dimensionality = (img.dimensions().z==1)? 1 : 3 ; 
			
			t1 = System.currentTimeMillis();
			Vector<Image> outputsFL = myhf.eigs(img.duplicate(), 2.0, true);
			t2 = System.currentTimeMillis();
			System.out.println("FLANAGAN! \ncalculated in: "+((t2-t1)/1000f)+" sec.");
			
			for (int i = 0; i < 2*dimensionality; i++) {
				outputsFL.get(i).imageplus().show();
			}
			
			t1 = System.currentTimeMillis();
			Vector<Image> outputsJM = myhj.eigs(img, 2.0, true);
			t2 = System.currentTimeMillis();
			System.out.println("JAMA! \ncalculated in: "+((t2-t1)/1000f)+" sec.");
			
			for (int i = 0; i < 2*dimensionality; i++) {
				outputsJM.get(i).imageplus().show();
			}
		
			// save them as csv files if 2d
			
			if(dimensionality==1){
				
				for (int i = 0; i < outputsFL.size(); i++) {
					System.out.println("save as csv..." );
					//outputsFL.get(i).imageplus().show();
				}
				
			}
			
			
	}
	
}
