package advantra.critpoint;

import java.awt.Checkbox;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.FilenameFilter;

import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.converters.ArffSaver;

import advantra.feature.DifferentialFeatures;
import advantra.file.AnalyzeCSV;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;

public class ExtractFeatures implements PlugIn {

	String 		train_dir 	= System.getProperty("user.dir")+File.separator+"train";
	String 		test_dir 	= System.getProperty("user.dir")+File.separator+"test";
	
	File[] csv_train 		= null; // train
	File[] tif_train 		= null;
	File[] tif_test 		= null;
	
	
	
	String[] 	diff_feature_labels = new String[DifferentialFeatures.FEATS_NR];
	boolean[] 	diff_feature_enable = new boolean[DifferentialFeatures.FEATS_NR];
	
	/*
	 * some other features labels, enable
	 * ????_feature_labels
	 * ????_feature_enable
	 */
	
	ImagePlus 	train_img;
	
	double[][] 	all_feats 		= null;
	int[] 		all_cls 		= null;
	String[] 	all_labels 		= null;
	
	String 		out_dir    		= System.getProperty("user.home");
	String 		out_file   		= "";
	
	
	
	/*
	 *  extraction parameters (standard deviations of the Gaussian used for smoothing/scaling & number of scales)
	 */

	double 		sigma_1 		= 2.0;
	double 		sigma_2 		= 3.0;
	int			nr 				= 2;
	
	public void run(String arg0) {
		
		IJ.log("Legend: \n" +
				"f01 -> gradient magnitude \n" +
				"f02 -> laplacian \n" +
				"f03 -> ridge detection \n" +
				"f04 -> isophote curvature \n" +
				"f05 -> flowline curvature \n" +
				"f06 -> isophote density \n" +
				"f07 -> corner detector \n" +
				"f08 -> shape index \n" +
				"f09 -> curvedness \n" +
				"f10 -> DoH \n" +
				"f11 -> Mean Curvature \n" +
				"f12 -> Gaussian Extremality \n" + 
				"f13 -> T-junction likeliness \n" +
				"f14 -> Ballness \n");
		
		for (int i = 0; i < diff_feature_labels.length; i++) {
			diff_feature_labels[i] = String.format("diff_feature_%02d", (i+1));
		}
		
		GenericDialog gd = new GenericDialog("Critical Point features");
		
		final int menu_width = 50;
		
		gd.addMessage("Locate folder...");
		gd.addStringField("folder with annotations", folder_name, menu_width);
		
		gd.addMessage("Choose diff. features...");
		gd.addCheckboxGroup(3, 5, diff_feature_labels, diff_feature_enable);
		
		gd.addMessage("Choose some other features...");
		
		gd.addMessage(		"Extraction params...");
		gd.addNumericField( "sigma start:", sigma_1, 	0, 5, "" );
		gd.addNumericField( "sigma end  :", sigma_2, 	0, 5, "" );	
		gd.addNumericField( "number of scales : ", nr,  0, 5, "");
		
		gd.addMessage("Save dataset...");
		gd.addStringField("destination directory", out_dir, menu_width);
		
		gd.showDialog();
		if (gd.wasCanceled()) return;
		
		folder_name 	= gd.getNextString();
		
		for (int i = 0; i < diff_feature_enable.length; i++) {
			diff_feature_enable[i] = ((Checkbox)gd.getCheckboxes().get(i)).getState();
		}
		
		/*
		 * read other features here...
		 */
		
		sigma_1 		= (double)gd.getNextNumber();
		sigma_2 		= (double)gd.getNextNumber();
		nr 				= (int)gd.getNextNumber();
		
		out_dir			= gd.getNextString();
		
		out_file		= "f_"; // prefix
		
		for (int i = 1; i <= diff_feature_enable.length; i++) {
			if(diff_feature_enable[i-1]) out_file += Integer.toString(i)+",";
		}
		
		out_file += String.format("s_%.2f,%.2f,%d", sigma_1, sigma_2, nr);
		out_dir = out_dir+File.separator+"training"+File.separator;
		//CreateDirectory.createOneDir(out_dir);
		
		File dir = new File(folder_name);
		folder_name = dir.getAbsolutePath();
		
		if(!dir.isDirectory()){
			IJ.error(folder_name+" is not a directory!");
			return;
		}
		
		csv_files = dir.listFiles(new FilenameFilter() {
		    	public boolean accept(File dir, String name) {
		    		return name.toLowerCase().endsWith(".csv");
		    	}
			}
		);
		
		File[] tif_files = dir.listFiles(
				new FilenameFilter() {
					public boolean accept(File dir, String name) {
						return name.toLowerCase().endsWith(".tif");
				}
			}
		);
		
		if(csv_files.length<=0){
			IJ.log("\nError: There was no csv annotation files.\n");
			return;
		}
		
		for (int i = 0; i < csv_files.length; i++) {
			String 	current_csv_name 	= csv_files[i].getName();
			System.out.print("extracting "+current_csv_name+" ... ");
			boolean found = false;
			int idx_found = 0;
			// is there a .tif pair, remember it's index
			for (int j = 0; j < tif_files.length; j++) {
				String current_tif_name = tif_files[j].getName();
				if(removeExt(current_csv_name).equals(removeExt(current_tif_name))){
					found = true;
					idx_found = j;
					break;
				}
			}
			
			if(!found){
				System.out.println("FAILED");
				continue; // try other .csv-s
			}
			
			System.out.println();
			System.out.println("processing "+tif_files[idx_found].getName()+"...");
			
			
			// match was found, extract image, locations and class value
			train_img = new ImagePlus(tif_files[idx_found].getAbsolutePath());
			
			AnalyzeCSV reader_csv = new AnalyzeCSV(csv_files[i].getAbsolutePath());
			double[][] train_loc	= reader_csv.readLn(2);
			int[]      train_cls 	= reader_csv.readLastCol();
			
			all_cls 	= concatenate(all_cls, train_cls);
			
			
			/*
			 * actual feature extraction 
			 */
			
			DifferentialFeatures df = new DifferentialFeatures(train_img, sigma_1, sigma_2, nr);
			double[][] diff_feats  = df.exportFeatures(train_loc, diff_feature_enable);
			
			/*
			 * add here if there's more features, careful when concatenating (below)
			 */
			
			all_labels = df.exportFeatureLabels(diff_feature_enable);// diff_feature_labels; // concatenate feature_labels if new feature types are added
			all_feats = concatenateRows(all_feats, diff_feats); // will be necessary to concatenate on diff_feats for new feature types
			
			
			
		}
		
		try {
			saveDataset(all_feats, all_labels, all_cls);
			IJ.log("saved trainset!");
		} 
		catch (Exception e) {
			e.printStackTrace();
		}
		
	}
	
	private void saveDataset(double[][] feats, String[] feats_labels, int[] cls) throws Exception {
		
		FastVector attrs = new FastVector();
		
		FastVector class_vals = new FastVector();
		class_vals.addElement("yes");
		class_vals.addElement("no");
		Attribute class_attr = new Attribute("critpoint", class_vals); // nominal
		attrs.addElement(class_attr); 
		for (int i = 0; i < feats_labels.length; i++) {
			Attribute attr = new Attribute(feats_labels[i]); //numeric
			attrs.addElement(attr);
		}
		
		Instances dataset = new Instances("my_dataset", attrs, 0); // capacity is zero
		
		for (int i = 0; i < feats.length; i++) {
			
			// add instances
			double[] attValues = new double[dataset.numAttributes()];
			attValues[0] = (cls[i]==1)?class_attr.indexOfValue("yes"):class_attr.indexOfValue("no");
			for (int j = 1; j <= feats[0].length; j++) {
				attValues[j] = feats[i][j-1];
			}
			dataset.add(new Instance(1.0, attValues));
		}
		dataset.setClassIndex(0); // first column is class
		
		// save dataset
		ArffSaver saver = new ArffSaver();
		File dataset_file 	= new File(out_dir+File.separator+out_file+"_trainset.arff");
		saver.setInstances(dataset);
		saver.setFile(dataset_file);
		saver.writeBatch();
		IJ.log("train dataset saved to: "+dataset_file.getAbsolutePath());
		
		// save parameters to a file in the same folder
		String params_path 	= out_dir+File.separator+out_file+"_params.txt";
		BufferedWriter bw 	= new BufferedWriter(new FileWriter(params_path), 32768);
		bw.write("sigma_1 = "+Double.toString(sigma_1)+"\n"+"sigma_2 = "+Double.toString(sigma_2)+"\n"+"nr = "+Integer.toString(nr)+"\n");
		bw.write("images used for training: \n");
		for (int i = 0; i < csv_files.length; i++) {
			bw.write(csv_files[i].getName()+"\n");
		}
		bw.close();
		IJ.log("feature extraction params saved to: "+dataset_file.getAbsolutePath());
		
	}
	
 	private String removeExt(String in_name){
		int name_len = in_name.length();
		if(name_len>4){
			name_len -= 4;
			return in_name.substring(0, name_len);
		}
		return "";
	}
	
	private int[] concatenate(int[] in1, int[] in2){
		
		if(in1==null){
			int[] out = new int[in2.length];
			for (int i = 0; i < in2.length; i++) {
				 out[i] = in2[i];
			 }
			return out;
		}
		
		int[] out = new int[in1.length+in2.length];
		 
		for (int i = 0; i < in1.length; i++) {
			out[i] = in1[i];
		}
		 
		for (int i = 0; i < in2.length; i++) {
			out[i+in1.length] = in2[i];
		}
		 
		return out;
	}
	
	private double[][] concatenateRows(double[][] in11, double[][] in21){
		
		if(in11==null){
			double[][] out = new double[in21.length][in21[0].length];
			for (int i = 0; i < in21.length; i++) {
				 for (int j = 0; j < in21[0].length; j++) {
					 out[i][j] = in21[i][j];
				 }
			}
			return out;
		}
		
		double[][] out = new double[in11.length+in21.length][in11[0].length];
		
		 for (int i = 0; i < in11.length; i++) {
			for (int j = 0; j < in11[0].length; j++) {
				out[i][j] = in11[i][j];
			}
		 }
		 
		 for (int i = 0; i < in21.length; i++) {
			 for (int j = 0; j < in21[0].length; j++) {
				 out[i+in11.length][j] = in21[i][j];
			 }
		 }
		
		return out;
		
	}
	
}