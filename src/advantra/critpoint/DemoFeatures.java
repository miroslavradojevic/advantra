package advantra.critpoint;

import java.awt.Color;

import advantra.feature.CircHaarFeat;
import advantra.feature.FilterSet;
import advantra.general.Sort;
import ij.IJ;
import ij.gui.GenericDialog;
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

        int         a1=40, a2=40, as=20;
        double      r1=0.2, r2=0.4;
        int         rn=2;

        GenericDialog gd = new GenericDialog("Visualize Features");

        gd.addMessage("angular scale:");
        gd.addNumericField("min scale", a1, 0, 5, "deg");
        gd.addNumericField("max scale", a2, 0, 5, "deg");
        gd.addNumericField("ang. step", as, 0, 5, "deg");

        gd.addMessage("radial scale:");
        gd.addNumericField("min scale", r1, 1);
        gd.addNumericField("max scale", r2, 1);
        gd.addNumericField("ang. step", rn, 0, 5, "#");

        gd.showDialog();
        if (gd.wasCanceled()) return;

        a1 	= 		    (int)gd.getNextNumber();
        a2	= 		    (int)gd.getNextNumber();
        as	= 		    (int)gd.getNextNumber();

        r1 	= 		    gd.getNextNumber();
        r2	= 		    gd.getNextNumber();
        rn	= 		    (int)gd.getNextNumber();

        int aSclNr = 0;
        for (int i = a1; i <= a2; i+=as) aSclNr++;

		int[] aScl = new int[aSclNr];
        for (int i = 0; i < aSclNr; i++) aScl[i] = a1+i*as;

        double[] rScl = new double[rn];
        for (int i = 0; i < rn; i ++){
            rScl[i] = (i==0)? r1 : r1+i*((r2-r1)/(rn-1));
        }

        for (int i = 0; i < aSclNr; i++) System.out.println(i + " : " + aScl[i]);
        for (int i = 0; i < rn; i++) System.out.println(i + " : " + rScl[i]);



		FilterSet fs = new FilterSet(aScl, rScl);
        fs.print();
		fs.showConfigs();
		System.out.println("total nr. configurations: "+ (fs.circConfs.size()+fs.radlConfs.size()));
		
	}


}
