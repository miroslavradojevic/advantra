package advantra.critpoint;

import java.awt.Checkbox;
import java.io.File;
import java.io.FileOutputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStream;

import javax.swing.JFileChooser;

import weka.classifiers.Classifier;
import weka.classifiers.trees.J48;
import weka.core.Instances;
import weka.core.converters.ArffLoader;

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
			{"J48", "NB"};
	boolean[] cls_enable = new boolean[]
			{false, false};
	
//	String 		out_dir    		= System.getProperty("user.home");
	
	public void run(String arg0) {
		
		//IJ.showMessage("HINT: Browse FEATS_*****.ARFF file with datasets for training!");
		
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
		if(!train_dataset_file.exists()){
			IJ.showMessage("file "+train_dataset_path+" does not exist");
			return;
		}
		
		/*
		 * generic dialog
		 */
		//final int menu_width = 50;
		GenericDialog gd = new GenericDialog("Train Classifier");
		gd.addMessage("Choose classifier(s)...");
		gd.addCheckboxGroup(3, 5, cls_labels, cls_enable);
//		gd.addMessage("Classification params...");
//		gd.addMessage("Save dataset...");
//		gd.addStringField("destination directory", out_dir, menu_width);
		gd.showDialog();
		if (gd.wasCanceled()) return;
		for (int i = 0; i < cls_enable.length; i++) {
			cls_enable[i] = ((Checkbox)gd.getCheckboxes().get(i)).getState();
		}
//		out_dir			= gd.getNextString();
		/*
		 * generic dialog
		 */
		
		try {
			
			/*
			 * main 
			 */
			
			Instances train_dataset = loadTrainDataset(train_dataset_file);
			train_dataset.setClassIndex(0);
			
			/*
			 * train & save classifier
			 */
			
			if(cls_enable[0]){
				IJ.log("training J48 classifier...");
				// J48 train & save
				Classifier classifier = new J48();
				classifier.buildClassifier(train_dataset);
				// save classifier
				String name = train_dataset_file.getName();//.substring(0, endIndex);
				name = name.substring(0, name.length()-5);
				String file = train_dataset_file.getParent()+File.separator+name+"_J48.arff";
				OutputStream os = new FileOutputStream(file);
				ObjectOutputStream objectOutputStream = new ObjectOutputStream(os);
				objectOutputStream.writeObject(classifier);
				objectOutputStream.close();
				IJ.log("done, exported "+file);
			}
			if(cls_enable[1]){
				IJ.log("training NB classifier");
				IJ.log("done, exported "+"nothing");
			}
			
		} 
		catch (Exception e) {
			e.printStackTrace();
		}
		
	}
	
	private Instances loadTrainDataset(File train_dataset_file) throws Exception {
		// read dataset from file
		ArffLoader loader = new ArffLoader();
		loader.setFile(train_dataset_file);
		Instances dataset = loader.getDataSet();
		return dataset;		
	}

}
