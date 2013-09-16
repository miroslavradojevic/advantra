package detection;

import ij.Prefs;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 9/1/13
 * Time: 4:44 PM
 */
public class FuzzyDemo implements PlugIn {


	public void run(String s) {

		float iDiff;

		iDiff				= (float) 	Prefs.get("advantra.critpoint.fuzzy.iDiff", 10);
		GenericDialog gd = new GenericDialog("FUZZY");
		gd.addNumericField("iDiff ",  		iDiff, 		0, 5, "[0,255]");

		gd.showDialog();
		if (gd.wasCanceled()) return;

		iDiff         = (float) gd.getNextNumber();
		Prefs.set("advantra.critpoint.fuzzy.iDiff",   iDiff);

		/*
		*****************************************
		 */

		Fuzzy f = new Fuzzy(iDiff);
		f.demo();

	}
}
