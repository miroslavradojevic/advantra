package profile;

import ij.IJ;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.plugin.PlugIn;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 9/1/13
 * Time: 4:44 PM
 */
public class FuzzyDemo implements PlugIn {


	public void run(String s) {

		float iDiff, ratio_iDiff;    // iBackgr
		float Q_veryhigh_A, Q_veryhigh_B;
		float Q_high_A, Q_high_B;
		float Q_moderate_A, Q_moderate_B;

		//iBackgr       		= (float)	Prefs.get("advantra.critpoint.fuzzy.iBackgr", 10);
		iDiff				= (float) 	Prefs.get("advantra.critpoint.fuzzy.iDiff", 10);
		ratio_iDiff         = (float) 	Prefs.get("advantra.critpoint.fuzzy.ratio_iDiff", 0.3);

//		Q_veryhigh_A       	= (float)	Prefs.get("advantra.critpoint.fuzzy.Q_veryhigh_A", 0.5);
//		Q_veryhigh_B       	= (float)	Prefs.get("advantra.critpoint.fuzzy.Q_veryhigh_B", 0.8);
//
//		Q_high_A       		= (float)	Prefs.get("advantra.critpoint.fuzzy.Q_high_A", 0.6);
//		Q_high_B       		= (float)	Prefs.get("advantra.critpoint.fuzzy.Q_high_B", 0.9);
//
//		Q_moderate_A       	= (float)	Prefs.get("advantra.critpoint.fuzzy.Q_moderate_A", 0.6);
//		Q_moderate_B       	= (float)	Prefs.get("advantra.critpoint.fuzzy.Q_moderate_B", 0.9);

		GenericDialog gd = new GenericDialog("FUZZY");

		//gd.addNumericField("iBackground ",  iBackgr, 	0, 5, "[0,255]");
		gd.addNumericField("iDiff ",  		iDiff, 		0, 5, "[0,255]");
		gd.addNumericField("ratio iDiff ", 	ratio_iDiff,1, 5, "");
//		gd.addMessage("VERY HIGH");
//		gd.addNumericField("A",  Q_veryhigh_A, 	1, 5, "");
//		gd.addNumericField("B",  Q_veryhigh_B, 	1, 5, "");
//		gd.addMessage("HIGH");
//		gd.addNumericField("A",  Q_high_A, 	1, 5, "");
//		gd.addNumericField("B",  Q_high_B, 	1, 5, "");
//		gd.addMessage("MODERATE");
//		gd.addNumericField("A",  Q_moderate_A, 	1, 5, "");
//		gd.addNumericField("B",  Q_moderate_B, 	1, 5, "");

		gd.showDialog();
		if (gd.wasCanceled()) return;

		//iBackgr       = (float) gd.getNextNumber();
		//Prefs.set("advantra.critpoint.fuzzy.iBackgr", iBackgr);

		iDiff         = (float) gd.getNextNumber();
		Prefs.set("advantra.critpoint.fuzzy.iDiff",   iDiff);

		ratio_iDiff   = (float) gd.getNextNumber();
		Prefs.set("advantra.critpoint.fuzzy.ratio_iDiff",   ratio_iDiff);

//		Q_veryhigh_A       = (float) gd.getNextNumber();
//		Prefs.set("advantra.critpoint.fuzzy.Q_veryhigh_A", Q_veryhigh_A);
//		Q_veryhigh_B       = (float) gd.getNextNumber();
//		Prefs.set("advantra.critpoint.fuzzy.Q_veryhigh_B", Q_veryhigh_B);
//
//		Q_high_A       = (float) gd.getNextNumber();
//		Prefs.set("advantra.critpoint.fuzzy.Q_high_A", Q_high_A);
//		Q_high_B       = (float) gd.getNextNumber();
//		Prefs.set("advantra.critpoint.fuzzy.Q_high_B", Q_high_B);
//
//		Q_moderate_A       = (float) gd.getNextNumber();
//		Prefs.set("advantra.critpoint.fuzzy.Q_moderate_A", Q_moderate_A);
//		Q_moderate_B       = (float) gd.getNextNumber();
//		Prefs.set("advantra.critpoint.fuzzy.Q_moderate_B", Q_moderate_B);

		/*
		*****************************************
		 */

		Fuzzy f = new Fuzzy(iDiff, ratio_iDiff);
		f.demo();


	}
}
