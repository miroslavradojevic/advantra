package advantra.critpoint;

import java.awt.Checkbox;
import java.io.File;
import java.io.FilenameFilter;
import java.util.ArrayList;

import weka.classifiers.Classifier;
import weka.classifiers.trees.J48;
import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
import advantra.feature.DifferentialFeatures;
import advantra.file.AnalyzeCSV;
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
	String train_dataset_file = "";
	
	File[] csv_train 		= null; // train
	File[] tif_train 		= null;
	File[] tif_test 		= null;
	
	int 		total_feat_nr 	= DifferentialFeatures.FEATS_NR;//+some other features
	String[] 	feature_labels 	= new String[total_feat_nr]; // read from the generic dialog
	boolean[] 	feature_enable 	= new boolean[total_feat_nr];
	
	/*
	 *  extraction parameters (standard deviations of the Gaussian used for smoothing/scaling & number of scales)
	 */
	FEparam params;// = new FEparam();

	String[] cls_labels = new String[]
			{"J48", "NB", "NN"};
	boolean[] cls_enable = new boolean[]
			{true, false, false};

	
	public void run(String arg0) {
		
		IJ.log("Legend: \n" +
				"f01 -> gradient magnitude " +
				"f02 -> laplacian \n" +
				"f03 -> ridge detection " +
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
		
		for (int i = 0; i < DifferentialFeatures.FEATS_NR; i++) {
			feature_labels[i] = String.format("df.ft %02d", (i+1));
		}
		
		GenericDialog gd = new GenericDialog("Critical Point Classification");
		
		final int menu_width = 50;
		
		gd.addMessage("***TRAIN***");
		gd.addStringField("folder with annotations", train_dir, menu_width);
		gd.addMessage("Choose features");
		gd.addCheckboxGroup(1, 14, feature_labels, feature_enable);
		
		gd.addMessage(		"Extraction params...");
		gd.addNumericField( "sigma start:", 			2, 	0, 5, "" );
		gd.addNumericField( "sigma end  :", 			3, 	0, 5, "" );	
		gd.addNumericField( "number of scales : ", 		2,  0, 5, "");
		
		gd.addMessage("***TEST***");
		gd.addStringField("test data directory", test_dir, menu_width);
		gd.addCheckboxGroup(1, 3, cls_labels, cls_enable);
		
		gd.addMessage("***OUTPUT***");
		gd.addStringField("destination directory", out_dir, menu_width);
		
		gd.showDialog();
		if (gd.wasCanceled()) return;
		
		train_dir 	= gd.getNextString();
		
		for (int i = 0; i < feature_enable.length; i++) {
			feature_enable[i] = ((Checkbox)gd.getCheckboxes().get(i)).getState();
			System.out.println("feature "+i+" : "+feature_enable[i]);
		}
		
		double 		sigma_1 		=  (double)gd.getNextNumber();
		double 		sigma_2 		=  (double)gd.getNextNumber();
		int 		nr				=  (int)gd.getNextNumber();
		
		test_dir		= gd.getNextString();
		// continue reading checkboxes
		for (int i = feature_enable.length; i < feature_enable.length+cls_enable.length; i++) {
			cls_enable[i-feature_enable.length] = ((Checkbox)gd.getCheckboxes().get(i)).getState();
		}
		
		out_dir			= gd.getNextString();
		
		params = new FEparam(sigma_1, sigma_2, nr, feature_enable);
		if(params.getNrFeatures()<=0) {
			IJ.log("no features were selected... closing");
			return;
		}
		
		train_dataset_file		= "feats_"; // prefix
		for (int i = 1; i <= feature_enable.length; i++) {
			if(feature_enable[i-1]) train_dataset_file += Integer.toString(i)+",";
		}
		//train_dataset_file += String.format("s_%.2f,%.2f,%d", sigma_1, sigma_2, nr);
		train_dataset_file = out_dir+File.separator+train_dataset_file+"_trainset.arff";
		
		File dir = new File(train_dir);
		train_dir = dir.getAbsolutePath();
		
		if(!dir.isDirectory()){
			IJ.error(train_dir+" is not a directory!");
			return;
		}
		
		csv_train = dir.listFiles(
				new FilenameFilter() {
					public boolean accept(File dir, String name) {
						return name.toLowerCase().endsWith(".csv");
		    	}
			}
		);
		
		tif_train = dir.listFiles(
				new FilenameFilter() {
					public boolean accept(File dir, String name) {
						return name.toLowerCase().endsWith(".tif");
				}
			}
		);
		
		if(csv_train.length<=0){
			IJ.log("\nError: There was no csv annotation files.\n");
			return;
		}
		
		int[] 		all_cls 		= null;
		
		//loop through train set
		ArrayList<ImagePlus> train_imgs 	= new ArrayList<ImagePlus>();
		ArrayList<double[][]> locs2D 		= new ArrayList<double[][]>();
		
		for (int i = 0; i < csv_train.length; i++) {
			
			String 	current_csv_name 	= csv_train[i].getName();
			System.out.print(current_csv_name+" ... ");
			boolean found = false;
			int idx_found = 0;
			// is there a .tif pair, remember it's index
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
			
			System.out.println("found match! "+tif_train[idx_found].getName()+"...");
			
			// match was found, add to the extraction list
			
			train_imgs.add(new ImagePlus(tif_train[idx_found].getAbsolutePath()));
			
			AnalyzeCSV reader_csv = new AnalyzeCSV(csv_train[i].getAbsolutePath());
			
			locs2D.add(reader_csv.readLn(2));//double[][] train_loc	= ;
			all_cls 	= concatenate(all_cls, reader_csv.readLastCol());
			
		}
		
		System.out.println("done: "+train_imgs.size()+" train images, "+params.getNrFeatures()+" features (each with "+params.nr_sigmas+" scales) used for training!");
			
		/*
		 *  feature extraction (fe considers all available features)
		 */
			
		FE fe 				= new FE(params);
		double[][] feats 	= fe.extractFeats(train_imgs, locs2D);
		String[] labels 	= fe.extractFeatLabels();
		
		System.out.println("features extracted!");

		/*
		 * make weka
		 */
		FastVector 	class_vals 	= new FastVector();
		class_vals.addElement("yes");
		class_vals.addElement("no");
		Attribute 	class_attr = new Attribute("critpoint", class_vals); // nominal
		
		FastVector 	attrs 		= new FastVector();
		attrs.addElement(class_attr); 
		for (int i = 0; i < labels.length; i++) {
			Attribute attr = new Attribute(labels[i]); //numeric
			attrs.addElement(attr);
		}
		
		Instances train_dataset = new Instances("train_dataset", attrs, 0); // capacity is zero
		
		for (int i = 0; i < feats.length; i++) {
			double[] attValues = new double[train_dataset.numAttributes()];
			attValues[0] = (all_cls[i]==1)?class_attr.indexOfValue("yes"):class_attr.indexOfValue("no");
			for (int j = 1; j <= feats[0].length; j++) {
				attValues[j] = feats[i][j-1];
			}
			train_dataset.add(new Instance(1.0, attValues));
		}
		train_dataset.setClassIndex(0); // first column is class
		
		/*
		 * train classifier(s)
		 */
		Classifier classifier;
		//if(cls_enable[0]){
			IJ.log("training J48 classifier...");
			// J48 train
			classifier = new J48();
			try {
				classifier.buildClassifier(train_dataset);
			} catch (Exception e) {
				e.printStackTrace();
			}
			IJ.log("done.");
//		}
//		if(cls_enable[1]){
//			IJ.log("training NB classifier");
//			IJ.log("empty now...");
//		}
//		if(cls_enable[2]){
//			IJ.log("training NN classifier");
//		}
		
		/*
		 * classification (fe+classify)
		 */
		dir = new File(test_dir);
		test_dir = dir.getAbsolutePath();
		
		if(!dir.isDirectory()){
			IJ.error(test_dir+" is not a directory!");
			return;
		}
		
		tif_test = dir.listFiles(
				new FilenameFilter() {
					public boolean accept(File dir, String name) {
						return name.toLowerCase().endsWith(".tif");
		    	}
			}
		);
		
		if(tif_test.length<=0){
			IJ.log("\nError: There was no tiff files in test set.\n");
			return;
		}
		
		//loop through test set
		ImagePlus test_img;
		ImageStack cls_stk;
		for (int i = 0; i < tif_test.length; i++) {
			
			String 	current_tif_name 	= tif_test[i].getName();
			System.out.print("loading "+current_tif_name+" ... ");
			
			test_img = new ImagePlus(tif_test[i].getAbsolutePath());
			
			feats 	= fe.extractFeats(test_img);
			labels 	= fe.extractFeatLabels();

			/*
			 * make weka
			 */
			class_vals 	= new FastVector();
			class_vals.addElement("yes");
			class_vals.addElement("no");
			class_attr = new Attribute("critpoint", class_vals); // nominal
			
			attrs 		= new FastVector();
			attrs.addElement(class_attr);
			for (int i1 = 0; i1 < labels.length; i1++) {
				Attribute attr = new Attribute(labels[i1]); //numeric
				attrs.addElement(attr);
			}
			
			Instances test_dataset = new Instances("test_dataset", attrs, 0); // capacity is zero
			
			for (int i1 = 0; i1 < feats.length; i1++) {
				double[] attValues = new double[test_dataset.numAttributes()]; // add instances to the dataset
				attValues[0] = class_attr.indexOfValue("no");
				for (int j = 1; j <= feats[0].length; j++) {
					attValues[j] = feats[i1][j-1];
				}
				test_dataset.add(new Instance(1.0, attValues));
			}
			test_dataset.setClassIndex(0); // first column is class
			
			// classify the test dataset, basically classifying all the pixels of an image
			int sz = test_img.getHeight()*test_img.getWidth();
			int wd = test_img.getWidth();
			int ht = test_img.getHeight();
			
			try{
				for (int el = 0; el < sz; el++) {
					classifier.classifyInstance(test_dataset.instance(el));
				}
			}
			catch(Exception e) {
				e.printStackTrace();
			}
			
			double[] classification_result = new double[sz];
			
			for (int j = 0; j < test_dataset.numInstances(); j++) {
				
				classification_result[j] = test_dataset.instance(j).classValue();
				
				// show what you got
				//int[] row_col = FE.index2location(j, wd);
				//System.out.format("instance %d (%d,%d) -> %s \n", j, row_col[0], row_col[1], test_dataset.attribute(test_dataset.classIndex()).value((int)test_dataset.instance(j).classValue()));
			}
			
			cls_stk = new ImageStack(wd, ht);
			cls_stk.addSlice(new FloatProcessor(wd, ht, classification_result));
			String output_name = "classified_"+tif_test[i].getName();
			(new ImagePlus(output_name, cls_stk)).show();
			
		}
		System.out.println("done.");
		
		/*
		 * classify every pixel of the image
		 */

		
		
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
	
//	private double[][] concatenateRows(double[][] in11, double[][] in21){
//		
//		if(in11==null){
//			double[][] out = new double[in21.length][in21[0].length];
//			for (int i = 0; i < in21.length; i++) {
//				 for (int j = 0; j < in21[0].length; j++) {
//					 out[i][j] = in21[i][j];
//				 }
//			}
//			return out;
//		}
//		
//		double[][] out = new double[in11.length+in21.length][in11[0].length];
//		
//		 for (int i = 0; i < in11.length; i++) {
//			for (int j = 0; j < in11[0].length; j++) {
//				out[i][j] = in11[i][j];
//			}
//		 }
//		 
//		 for (int i = 0; i < in21.length; i++) {
//			 for (int j = 0; j < in21[0].length; j++) {
//				 out[i+in11.length][j] = in21[i][j];
//			 }
//		 }
//		
//		return out;
//		
//	}
	
}

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