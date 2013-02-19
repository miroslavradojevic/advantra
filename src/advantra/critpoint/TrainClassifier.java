package advantra.critpoint;

import java.awt.Checkbox;
import java.io.File;

import javax.swing.JFileChooser;

import ij.IJ;
import ij.gui.GenericDialog;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;

public class TrainClassifier  implements PlugIn {

	/*
	 * choose classifier and the train dataset, train it and save!
	 */
	
	String 	train_dataset_path    = "";
	File 	train_dataset_file;

	String[] cls_labels = new String[]
			{"Tree", "NB"};
	boolean[] cls_enable = new boolean[]
			{true, true};
	
	String 		out_dir    		= System.getProperty("user.home");
	
	public void run(String arg0) {
		
		IJ.showMessage("HINT: Browse TRAIN_FEATS_*****.ARFF file with datasets for training!");
		
		/*
		 * based on File_Opener.java ij source code 
		 * this sequence enables to browse any file 
		 * and get the path to it
		 */
		File dir=null;
		JFileChooser fc = null;
		try {fc = new JFileChooser();}
		catch (Throwable e) {IJ.error("This plugin requires Java 2 or Swing."); return;}
		//fc.setMultiSelectionEnabled(true);
		if (dir==null) {
			String sdir = OpenDialog.getDefaultDirectory();
			if (sdir!=null)
				dir = new File(sdir);
		}
		if (dir!=null)
			fc.setCurrentDirectory(dir);
		int returnVal = fc.showOpenDialog(IJ.getInstance());
		if (returnVal!=JFileChooser.APPROVE_OPTION)
			return;
		train_dataset_file = fc.getSelectedFile();
		train_dataset_path = fc.getSelectedFile().getAbsolutePath();
		if(train_dataset_file.exists())IJ.showMessage("opened: "+train_dataset_path);
		else {
			IJ.showMessage("file "+train_dataset_path+" does not exist");
			return;
		}
		/*
		 * 
		 */
		
		GenericDialog gd = new GenericDialog("Critical Point features");
		
		final int menu_width = 50;
		
		gd.addMessage("Choose classifier(s)...");
		//gd.addStringField("folder with annotations", folder_name, menu_width);
		gd.addCheckboxGroup(3, 5, cls_labels, cls_enable);
		
		gd.addMessage("Classification params...");
		
		gd.addMessage("Save dataset...");
		gd.addStringField("destination directory", out_dir, menu_width);
		
		gd.showDialog();
		if (gd.wasCanceled()) return;
		
		for (int i = 0; i < cls_enable.length; i++) {
			cls_enable[i] = ((Checkbox)gd.getCheckboxes().get(i)).getState();
		}
		
		out_dir			= gd.getNextString();
		
	}

}
