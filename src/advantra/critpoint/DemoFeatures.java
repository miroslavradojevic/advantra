package advantra.critpoint;

import java.awt.Color;

import advantra.feature.CircHaarFeat;
import advantra.feature.ProfileFeatures;
import advantra.general.Sort;
import ij.IJ;
import ij.gui.Plot;
import ij.plugin.PlugIn;

public class DemoFeatures implements PlugIn {

	public void run(String arg0) {
		
if(false){		
		
		CircHaarFeat circhaar = new CircHaarFeat(16);
		circhaar.createFeatures();
		circhaar.showFeatures();
		System.out.println("nr fts. " + circhaar.nrFeats());
		
		double[] profile = new double[16];
		profile[0] = 10;
		for (int i = 1; i < profile.length; i++) {
			profile[i] = 1;
		}
		
		long t11 = System.currentTimeMillis();
		float[] res = circhaar.allFeatScore(profile, 1);
		System.out.println(""+IJ.d2s((System.currentTimeMillis()-t11)/1000f, 2)+" s.");

		
		Plot p = new Plot("profile response", "feat #", "intens.");
		p.setLimits(0, circhaar.nrFeats(), Sort.findMin(res), Sort.findMax(res));
		p.setColor(Color.RED);
		float[] indx = new float[res.length];
		for (int i = 0; i < res.length; i++) {
			indx[i] = i;
		}
		p.addPoints(indx, res, Plot.LINE);
		p.draw();
		p.show();
}

		/*
		 * test profile features
		 */
		
		// create angular profile
		double radius		= 10;
		double dr 			= 1;
		double darc 		= 1;
		double rratio 		= 0.2;
		
		// count the number of points
		int cnt = 0;
		for (double r = radius; r >= radius*rratio; r-=dr) {
			for (double arc = 0; arc < 2*r*Math.PI; arc+=darc) {
				cnt++;
			}
		}
		
		double[] angProfile 	= new double[cnt];
		
		cnt = 0;
		for (double r = radius; r >= radius*rratio; r-=dr) {
			for (double arc = 0; arc < 2*r*Math.PI; arc+=darc) {
				angProfile[cnt++] = arc/r;
				
			}
		}
		
		System.out.println("created " + cnt + " angles");
		int nrFeats = 128;
		double angStep = Math.PI/8;
		ProfileFeatures pft = new ProfileFeatures(nrFeats, angStep);
		pft.createFeatures();
		pft.showFeatures();
	    
	}


}
