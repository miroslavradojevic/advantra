package advantra.critpoint;

import java.awt.Color;

import advantra.feature.CircHaarFeat;
import advantra.feature.CircularFilterSet;
import advantra.general.Sort;
import ij.IJ;
import ij.gui.Plot;
import ij.plugin.PlugIn;

public class DemoFeatures implements PlugIn {

	public void run(String arg0) {
	
		boolean exec = false;
		
//	if(exec){
//
//		CircHaarFeat circhaar = new CircHaarFeat(16);
//		circhaar.createFeatures();
//		circhaar.showFeatures();
//		System.out.println("nr fts. " + circhaar.nrFeats());
//
//		double[] profile = new double[16];
//		profile[0] = 10;
//		for (int i = 1; i < profile.length; i++) {
//			profile[i] = 1;
//		}
//
//		long t11 = System.currentTimeMillis();
//		float[] res = circhaar.allFeatScore(profile, 1);
//		System.out.println(""+IJ.d2s((System.currentTimeMillis()-t11)/1000f, 2)+" s.");
//
//
//		Plot p = new Plot("profile response", "feat #", "intens.");
//		p.setLimits(0, circhaar.nrFeats(), Sort.findMin(res), Sort.findMax(res));
//		p.setColor(Color.RED);
//		float[] indx = new float[res.length];
//		for (int i = 0; i < res.length; i++) {
//			indx[i] = i;
//		}
//		p.addPoints(indx, res, Plot.LINE);
//		p.draw();
//		p.show();
//	}

		/*
		 * test profile features
		 */

		int[] angleRes = new int[]{40, 60};
		CircularFilterSet cft = new CircularFilterSet(angleRes);
		cft.showConfigs();
		IJ.log("nr. configurations: "+ cft.filts.size());
		
	}


}
