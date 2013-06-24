import ij.IJ;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;

import java.io.File;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 6/23/13
 * Time: 5:38 PM
 */
public class BifDetectML implements PlugIn {

	int neuronDiamMin, neuronDiamMax;
	double scaleMin, scaleMax;
	String train_folder, test_folder;

	public void run(String s) {

		neuronDiamMin     	= (int) Prefs.get("advantra.critpoint.neuron_diam_min", 2);
		neuronDiamMax     	= (int) Prefs.get("advantra.critpoint.neuron_diam_max", 4);

		scaleMin     	= Prefs.get("advantra.critpoint.neuron_diam_min", 2);
		scaleMax     	= Prefs.get("advantra.critpoint.neuron_diam_max", 4);

		train_folder    = (String)Prefs.get("advantra.critpoint.train_folder", (System.getProperty("user.home")+ File.separator));
		test_folder     = (String)Prefs.get("advantra.critpoint.test_folder", (System.getProperty("user.home")+File.separator));

		GenericDialog gd = new GenericDialog("Bif. Detection");
		gd.addNumericField("n. diam. (min)", neuronDiamMin, 0, 5, "");
		gd.addNumericField("n. diam. (min)", neuronDiamMax, 0, 5, "");

		gd.addStringField("train folder : ", train_folder, 	40);
		gd.addStringField("test  folder : ", test_folder, 	40);
		gd.addCheckbox("images got same dimensions", false);
		gd.addCheckbox("equal # (+) and (-)", false);

		gd.showDialog();
		if (gd.wasCanceled()) return;

		neuronDiamMin =  (int)gd.getNextNumber();
		Prefs.set("advantra.critpoint.neuron_diam_min", 	neuronDiamMin);

		train_folder = 	gd.getNextString();
		test_folder = 	gd.getNextString();
		if(! new File(train_folder).exists())	{IJ.showMessage("train folder does not exist!"); 	return;}
		if(! new File(test_folder).exists()) 	{IJ.showMessage("test folder does not exist!"); 	return;}
		Prefs.set("advantra.critpoint.train_folder",    train_folder);
		Prefs.set("advantra.critpoint.test_folder",     test_folder);

		boolean sameSize = gd.getNextBoolean();
		boolean equal = gd.getNextBoolean();



	}

}


