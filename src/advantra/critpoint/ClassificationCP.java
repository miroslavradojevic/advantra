package advantra.critpoint;

import java.awt.Checkbox;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;

import weka.classifiers.Classifier;
import weka.classifiers.trees.J48;
import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.converters.CSVLoader;
import weka.filters.Filter;
import weka.filters.unsupervised.attribute.Remove;
import advantra.feature.DifferentialFeatures;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;

public class ClassificationCP implements PlugIn {

	String 		train_dir 	= System.getProperty("user.home")+File.separator+"train";
	String 		test_dir 	= System.getProperty("user.home")+File.separator+"test";
	String 		out_dir    	= System.getProperty("user.home")+File.separator+"res";
	
	File[] csv_train 		= null;
	File[] tif_train 		= null;
	File[] tif_test 		= null;
	File[] csv_test			= null;
	
	int 		total_feat_nr 	= DifferentialFeatures.FEATS_NR;//+some other features
	String[] 	feature_labels 	= new String[total_feat_nr]; // read from the generic dialog
	boolean[] 	feature_enable 	= new boolean[total_feat_nr];
	
	// extraction parameters (standard deviations of the Gaussian used for smoothing/scaling & number of scales)
	FEparam params;
	
	public void run(String arg0) {

		GenericDialog gd = new GenericDialog("Critical Point Classification");
		
		final int menu_width = 50;
		
		gd.addStringField("folder with annotations", train_dir, menu_width);
		gd.addMessage("Choose features");
		gd.addMessage("f(orig.)");
		gd.addCheckboxGroup(3, 5, FE.exportAllLabels(), feature_enable);
		gd.addMessage("f(DoH.)");
		gd.addMessage("f(|lambda1|.)");
		
		gd.addNumericField( "sigma start:", 			2, 	0, 5, "" );
		gd.addNumericField( "sigma end  :", 			3, 	0, 5, "" );	
		gd.addNumericField( "number of scales : ", 		2,  0, 5, "");
		
		gd.addStringField("test data directory", test_dir, menu_width);
		//gd.addCheckboxGroup(1, 3, cls_labels, cls_enable);
		
		//gd.addMessage("***OUTPUT***");
		//gd.addStringField("destination directory", out_dir, menu_width);
		
		gd.showDialog();
		if (gd.wasCanceled()) return;
		
		train_dir 	= gd.getNextString();
		
		for (int i = 0; i < feature_enable.length; i++) {
			feature_enable[i] = ((Checkbox)gd.getCheckboxes().get(i)).getState();
		}
		
		double 		sigma_1 		=  (double)gd.getNextNumber();
		double 		sigma_2 		=  (double)gd.getNextNumber();
		int 		nr				=  (int)gd.getNextNumber();
		
		test_dir		= gd.getNextString();
		
		/*
		 * dialog read done
		 */
		
		params = new FEparam(sigma_1, sigma_2, nr, feature_enable);
		if(params.getNrFeatures()<=0) {
			IJ.showMessage("no features were selected...");
			return;
		}
		
		File dir1 = new File(train_dir);
		train_dir = dir1.getAbsolutePath();
		File dir2 = new File(test_dir);
		if(!dir1.isDirectory() || !dir2.isDirectory()){
			IJ.error("Wrong train/test directory!");
			return;
		}
		
		// train dir csv files
		csv_train = dir1.listFiles(
				new FilenameFilter() {
					public boolean accept(File dir, String name) {
						return name.toLowerCase().endsWith(".csv");
		    	}
			}
		);
		// train dir tif files
		tif_train = dir1.listFiles(
				new FilenameFilter() {
					public boolean accept(File dir, String name) {
						return name.toLowerCase().endsWith(".tif");
				}
			}
		);
		// test dir tif files
		tif_test = dir2.listFiles(
				new FilenameFilter() {
					public boolean accept(File dir, String name) {
						return name.toLowerCase().endsWith(".tif");
				}
			}
		);
		// test dir csv files
		csv_test = dir2.listFiles(
				new FilenameFilter() {
					public boolean accept(File dir, String name) {
						return name.toLowerCase().endsWith(".csv");
				}
			}
		);
		
		if(csv_train.length<=0){
			IJ.showMessage("\nThere was no csv files in train dir.\n");
			return;
		}
		
		//loop through train set
		ArrayList<Instances> 		clss 		= new ArrayList<Instances>();
		ArrayList<Instances> 		data12 		= new ArrayList<Instances>();
		ArrayList<ImagePlus> 		train_imgs 	= new ArrayList<ImagePlus>();
		
		for (int i = 0; i < csv_train.length; i++) { // for each .csv file
			
			String 	current_csv_name 	= csv_train[i].getName();
			System.out.print(current_csv_name+" ... ");
			boolean found = false;
			int idx_found = 0; // is there a .tif pair, remember it's index
			
			for (int j = 0; j < tif_train.length; j++) {
				String current_tif_name = tif_train[j].getName();
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
			
			// match was found, add to the extraction list
			System.out.println("OK");
			
			Instances data = loadFromCsvFile(csv_train[i]);
			// data12
			Remove keep12 = new Remove();
			keep12.setInvertSelection(true); 		// attribs are kept
			int[] first_two = new int[]{0, 1};
			keep12.setAttributeIndicesArray(first_two);
			try {
				keep12.setInputFormat(data);
				data12.add(Filter.useFilter(data, keep12));
			} catch (Exception e) {
				System.out.println("There was a problem adding locations data.");
				e.printStackTrace();
			}
			// clss
			Remove keepLast = new Remove();
			keepLast.setInvertSelection(true);
			int[] last_idx = new int[]{(data.numAttributes()-1)};
			keepLast.setAttributeIndicesArray(last_idx);
			try {
				keepLast.setInputFormat(data);
				clss.add(Filter.useFilter(data, keepLast));
			}
			catch (Exception e) {
				System.out.println("There was a problem adding class data.");
				e.printStackTrace();
			}
			
			train_imgs.add(new ImagePlus(tif_train[idx_found].getAbsolutePath()));
			
		}
		
		if(train_imgs.size()<=0 || clss.size()<=0 || data12.size()<=0){
			System.out.println("No images/locations to extract features from.");
			return;
		}
		
		// concatenate the clss 
		Instances classes = new Instances(clss.get(0));
		
		for (int k = 1; k < clss.size(); k++) {
			for (int l = 0; l < clss.get(k).numInstances(); l++) {
				classes.add(clss.get(k).instance(l));
			}
		}
		
		FE fe 				= new FE(params);
		double[][] feats 	= fe.extractFeats(train_imgs, data12); // each instance corr. to one location
		String[] labels 	= fe.extractFeatLabels(true); // with scales
		
		FastVector fv = new FastVector(labels.length);
		for (int i = 0; i < labels.length; i++) {
			Attribute att = new Attribute(labels[i]); // numeric att.
			fv.addElement(att);
		}
		
		Instances feat_set = new Instances("features for training", fv, feats.length);
		
		for (int i = 0; i < feats.length; i++) {
			feat_set.add(new Instance(1.0, feats[i]));
		}
		
		Instances train_set = Instances.mergeInstances(feat_set, classes);
		train_set.setClassIndex((train_set.numAttributes()-1));
		
		IJ.showMessage("Train set formed!\n" +
				""+train_set.numInstances()+" instances, "+train_set.numAttributes()+" attributes.\n" +
						"saved to:\n" +
						"");
		
		System.out.println(train_set);
		
		// define classifiers
		//ArrayList<Classifier> classifiers = new ArrayList<Classifier>();
		//ArrayList<String> classifiers_id = new ArrayList<String>();
		
		Classifier tree = new J48();
		
		//classifiers.add(new J48());
		//classifiers_id.add("J48tree");
		
		// build classifiers
		try {
			
			tree.setOptions(new String[]{"-U"});
			tree.buildClassifier(train_set);
			
		} 
		catch (Exception e) {
			e.printStackTrace();
		}

		IJ.showMessage("Classifiers built!");
		System.out.println(tree);
		// cross-validation
		
		
		if(tif_test.length>0){
		 	
		for (int i = 0; i < tif_test.length; i++) {
			
			String 	current_tif_name 	= tif_test[i].getName();
			
			System.out.println("loading "+current_tif_name+" ... ");
			
			ImagePlus test_img = new ImagePlus(tif_test[i].getAbsolutePath());
			test_img.show();
			// will use the same parameters
			
			feats 	= fe.extractFeats(test_img); // each instance one location
			labels 	= fe.extractFeatLabels(true);
			
			classes = new Instances(clss.get(0), feats.length); // empty dataset with classes, reuse
			
			fv = new FastVector(labels.length);
			for (int k = 0; k < labels.length; k++) {
				Attribute att = new Attribute(labels[k]); // numeric att.
				fv.addElement(att);
			}
			
			feat_set = new Instances("features for classification", fv, feats.length); //reuse
			
			for (int k = 0; k < feats.length; k++) {
				classes.add(new Instance(1)); // add instance with 1 attrib, undefined value
				feat_set.add(new Instance(1.0, feats[k]));
			}
			
			Instances test_set = Instances.mergeInstances(feat_set, classes); // each instance one location
			test_set.setClassIndex((test_set.numAttributes()-1));
			
			int sz = test_img.getHeight()*test_img.getWidth();
			try{
				for (int el = 0; el < sz; el++) {
					// classify and set the value
					double cls_label = tree.classifyInstance(test_set.instance(el));
					test_set.instance(el).setClassValue(cls_label);
				}
			}
			catch(Exception e) {
				e.printStackTrace();
			}
			
			// classify the test dataset, basically classifying all the pixels of an image
			double[] classification_result = new double[sz];
			// take the dataset where each instance is one location 	
			for (int j = 0; j < test_set.numInstances(); j++) {
				// value predicted by classifier	
				classification_result[j] = test_set.instance(j).classValue();
			}
			
			int wd = test_img.getWidth();
			int ht = test_img.getHeight();	
				
			ImageStack cls_stk = new ImageStack(wd, ht);
			cls_stk.addSlice(new FloatProcessor(wd, ht, classification_result));
			String output_name = "classified_"+tif_test[i].getName();
			(new ImagePlus(output_name, cls_stk)).show();
			
		}
		
		}
		
		IJ.showMessage("finished");
		
	}
	
 	private String removeExt(String in_name){
		int name_len = in_name.length();
		if(name_len>4){
			name_len -= 4;
			return in_name.substring(0, name_len);
		}
		return "";
	}
	
 	private Instances loadFromCsvFile(File f)  {
 		CSVLoader loader 	= new CSVLoader();
 		Instances data;
 		
 		try {
 			loader.setSource(f);
			data = loader.getDataSet();
			return data;
		} 
 		catch (IOException e) {
			System.out.println("There was a problem loading "+f.getName()+" file. Returning null.");
			e.printStackTrace();
			return null;
		}

 	}
	

}


/*
private void saveDataset(double[][] feats, String[] feats_labels, int[] cls) {
	
	try {
	
	FastVector attrs 		= new FastVector();
	FastVector class_vals 	= new FastVector();
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
	File dataset_file 	= new File(train_dataset_file);//(out_dir+File.separator+out_file+"_trainset.arff");
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
	catch (Exception e) {
		e.printStackTrace();
	}
	
}
*/

//private int[] concatenate(int[] in1, int[] in2){
//	
//	if(in1==null){
//		int[] out = new int[in2.length];
//		for (int i = 0; i < in2.length; i++) {
//			 out[i] = in2[i];
//		 }
//		return out;
//	}
//	
//	int[] out = new int[in1.length+in2.length];
//	 
//	for (int i = 0; i < in1.length; i++) {
//		out[i] = in1[i];
//	}
//	 
//	for (int i = 0; i < in2.length; i++) {
//		out[i+in1.length] = in2[i];
//	}
//	 
//	return out;
//}
//private static double[][] extractCols(Instances data, int[] attribs_to_keep) throws Exception {
//
//Remove rm = new Remove();
//rm.setInvertSelection(true); // columns are kept
//rm.setAttributeIndicesArray(attribs_to_keep);//""+data.numAttributes()+"");
//rm.setInputFormat(data);
//Instances data_no_last_att = Filter.useFilter(data, rm);
//
//// convert to double[][]
//double[][] out = new double[data_no_last_att.numInstances()][data_no_last_att.numAttributes()];
//
//for (int i = 0; i < data_no_last_att.numInstances(); i++) {
//	out[i] = data_no_last_att.instance(i).toDoubleArray();
//}
//
//return out;
//
//}
//// save dataset snippet
//ArffSaver saver = new ArffSaver();
//File dataset_file 	= new File(train_dataset_file);//(out_dir+File.separator+out_file+"_trainset.arff");
//saver.setInstances(dataset);
//saver.setFile(dataset_file);
//saver.writeBatch();
//IJ.log("train dataset saved to: "+dataset_file.getAbsolutePath());
//// save classifier snippet
//String name = train_dataset_file.getName();
//name = name.substring(0, name.length()-5);
//String file = train_dataset_file.getParent()+File.separator+name+"_J48.arff";
//OutputStream os = new FileOutputStream(file);
//ObjectOutputStream objectOutputStream = new ObjectOutputStream(os);
//objectOutputStream.writeObject(classifier);
//objectOutputStream.close();