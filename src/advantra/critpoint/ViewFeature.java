package advantra.critpoint;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.io.OpenDialog;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import imagescience.image.Coordinates;
import imagescience.image.FloatImage;
import imagescience.image.Image;

import java.io.File;
import java.io.IOException;

import javax.swing.JFileChooser;

import advantra.feature.Ballness2D;
import advantra.feature.DoH2D;
import advantra.feature.Laplacian2D;
import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.converters.CSVLoader;
import weka.core.converters.ConverterUtils.DataSink;
import weka.filters.Filter;
import weka.filters.unsupervised.attribute.Remove;

public class ViewFeature implements PlugIn {

//	String 	dir_path 		= System.getProperty("user.home")+File.separator+"train";
//	File[] csv_files 		= null;
//	File[] tif_files 		= null;
	
	public void run(String arg0) {
		
		/*
		 * open image
		 */
		
		OpenDialog open_image 	= new OpenDialog("Select image for CP feature extraction", null);
		String image_dir = open_image.getDirectory();
		if(image_dir!=null){
			if (!image_dir.endsWith("/")) image_dir += "/";
		}
		else{
			System.out.println("image was null!");
			return;
		}
		
		ImagePlus imp = new ImagePlus(image_dir+open_image.getFileName());
		
		// reset calibration before going further
		Calibration cal = new Calibration(imp);
		cal.pixelWidth = cal.pixelHeight = cal.pixelDepth = 1.0;
		cal.setUnit("pixel");
		imp.setCalibration(cal);
		
		Image im = new FloatImage(Image.wrap(imp));

		imp.show();
		
		/*
		 * select filter and scales
		 */
		
		String[] choices = new String[3];
		choices[0] = "laplacian";
		choices[1] = "DoH";
		choices[2] = "abs(L1)";
		
		GenericDialog gd = new GenericDialog("Critical Point feature analysis");
		gd.addMessage("Choose features");
		gd.addChoice("proba", choices, choices[0]);
		
		gd.addNumericField( "sigma start:", 		2, 	 0, 5, "");
		gd.addNumericField( "sigma end  :", 		6, 	 0, 5, "");	
		gd.addNumericField( "# scales   :", 		10,  0, 5, "");
		
		gd.showDialog();
		if (gd.wasCanceled()) return;
		
		int idx =  gd.getNextChoiceIndex();
		double 		sigma_1 		=  (double)gd.getNextNumber();
		double 		sigma_2 		=  (double)gd.getNextNumber();
		int 		nr				=  (int)gd.getNextNumber();

		double[] sigma = new double[nr];
		for (int i = 0; i < nr; i++) {
			sigma[i] = (i==0)? sigma_1 : sigma_1+i*((sigma_2-sigma_1)/(nr-1));
		}
		
		/*
		 * extract filter at scales
		 */
		
		Image feats = new FloatImage(im.dimensions());// later on it is allocated with full range of scales, layers are added
		
		switch (idx) {
		case 0:	// laplacian
			feats = Laplacian2D.calculateImg(im, sigma);
			feats.imageplus().show();
			break;
		case 1: // DoH
			feats = DoH2D.calculateImg(im, sigma);
			feats.imageplus().show();
			break;
		case 2: // abs(L1)
			feats = Ballness2D.calculateImg(im, sigma);
			feats.imageplus().show();
			break;
		default:
			
			break;
		}
		
		/*
		 * 
		 */
		
		OpenDialog open_csv = new OpenDialog("Select csv file with CP annotations (optional)", null);
		String csv_dir = open_csv.getDirectory();
		if(csv_dir!=null){
			if (!csv_dir.endsWith("/")) csv_dir += "/";
			System.out.println("csv was loaded!");
		}
		else{
			System.out.println("csv was null! won't be extracting the features");
			return;
		}
		
		/*
		 * extract features
		 */
		
		File csv_file = new File(csv_dir+open_csv.getFileName()); 
		Instances data = loadFromCsvFile(csv_file);
		
		Instances 		locs;
		Instances		clss;
		Instances 		train;
		
		Remove keep12 = new Remove();
		keep12.setInvertSelection(true); 	
		int[] first_two = new int[]{0, 1};
		keep12.setAttributeIndicesArray(first_two);
		
		Remove keepLast = new Remove();
		keepLast.setInvertSelection(true);
		int[] last_idx = new int[]{(data.numAttributes()-1)};
		keepLast.setAttributeIndicesArray(last_idx);
		
		FastVector attributes = new FastVector();
		for (int i = 0; i < nr; i++) {
			String label = String.format(choices[idx]+"s_%.2f", sigma[i]);
			attributes.addElement(new Attribute(label));
		}
		
		try {
			
			keepLast.setInputFormat(data);
			clss = Filter.useFilter(data, keepLast);
			attributes.addElement(clss.attribute(0));
			
			keep12.setInputFormat(data);
			locs = Filter.useFilter(data, keep12);
			
			System.out.println("extracted dataset with locations & assigned class for each.");
			
			train = new Instances("train", attributes, locs.numInstances());
			for (int i = 0; i < locs.numInstances(); i++) {
				train.add(new Instance(nr+1)); // fill them with missing values
			}
			train.setClassIndex(train.numAttributes()-1);
			
			int cnt = 0;
			for (int loc_idx = 0; loc_idx < locs.numInstances(); loc_idx++) {
			
				for (int scale_idx = 0; scale_idx < nr; scale_idx++) {
				
					int col = (int)Math.round( locs.instance(loc_idx).value(0) );
					int row = (int)Math.round( locs.instance(loc_idx).value(1) );
					Coordinates at_pos = new Coordinates(col, row, scale_idx);
					double value = feats.get(at_pos);
					train.instance(cnt).setValue(scale_idx, value);
					
				}
				
				train.instance(cnt).setValue(nr, clss.instance(loc_idx).value(0));
				cnt++;
				
			}
			
			String train_path = System.getProperty("user.home");
			if (!train_path.endsWith("/")) train_path += "/";
			train_path += "train.csv";
			
			System.out.println("wants to save to: "+train_path);
			
			DataSink.write(train_path, train);
			
			System.out.println("done, saved to "+train_path);
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}
	
 	private String 		removeExt(String in_name){
		int name_len = in_name.length();
		if(name_len>4){
			name_len -= 4;
			return in_name.substring(0, name_len);
		}
		return "";
	}
 	
 	private Instances 	loadFromCsvFile(File f)  {
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

	public String openLocation(){
		File dir=null;
		JFileChooser fc = null;
		try {fc = new JFileChooser();}
		catch (Throwable e) {return null;}
		if (dir==null) {
			String sdir = OpenDialog.getDefaultDirectory();
			if (sdir!=null)
				dir = new File(sdir);
		}
		if (dir!=null)
			fc.setCurrentDirectory(dir);
		int returnVal = fc.showOpenDialog(IJ.getInstance());
		if (returnVal!=JFileChooser.APPROVE_OPTION)
			return null;
		
		System.out.println("The");
		
		File 	train_dataset_file = fc.getSelectedFile();
		String 	train_dataset_path = fc.getSelectedFile().getAbsolutePath();
		if(!train_dataset_file.exists()){
			return train_dataset_path;
		}
		else{
			return null;
		}
	}
 	
}

//dir_path = (new File(dir_path)).getAbsolutePath();
//File dir = new File(dir_path);
//if(!dir.isDirectory()){
//	IJ.error("Wrong source directory!");
//	return;
//}
//// csv files
//csv_files = dir.listFiles(
//		new FilenameFilter() {
//			public boolean accept(File dir, String name) {
//				return name.toLowerCase().endsWith(".csv");
//    	}
//	}
//);
//// tif files
//tif_files = dir.listFiles(
//		new FilenameFilter() {
//			public boolean accept(File dir, String name) {
//				return name.toLowerCase().endsWith(".tif");
//		}
//	}
//);
//
//if(csv_files.length<=0){
//	IJ.showMessage("\nThere was no csv files in source.\n");
//	return;
//}
////loop through train set
//ArrayList<Instances> 		clss 		= new ArrayList<Instances>();
//ArrayList<Instances> 		locs 		= new ArrayList<Instances>();
//ArrayList<Image> 			imgs 		= new ArrayList<Image>();
//
//for (int i = 0; i < csv_files.length; i++) { // for each .csv file
//	
//	String 	current_csv_name 	= csv_files[i].getName();
//	System.out.print(current_csv_name+" ... ");
//	boolean found = false;
//	int idx_found = 0; // is there a .tif pair, remember it's index
//	
//	for (int j = 0; j < tif_files.length; j++) {
//		String current_tif_name = tif_files[j].getName();
//		if(removeExt(current_csv_name).equals(removeExt(current_tif_name))){
//			found = true;
//			idx_found = j;
//			break;
//		}
//	}
//	if(!found){
//		System.out.println("FAILED");
//		continue; // try other .csv-s
//	}
//	// match was found, add to the extraction list
//	ImagePlus ip 	= new ImagePlus(tif_files[idx_found].getAbsolutePath());

//	Image im = new FloatImage(Image.wrap(ip));
//	imgs.add(im);
//	System.out.println("stored locations, classes, and image"); // imgs, clss, locs
//}
//if(imgs.size()<=0 || clss.size()<=0 || locs.size()<=0){
//	System.out.println("No images/locations to extract features from.");
//	return;
//}
//
//int nr_locs = 0;
//for (int j = 0; j < locs.size(); j++) {
//	nr_locs += locs.get(j).numInstances();
//}
//// fill it with values 
//int fill_up_idx = 0;
//
//switch (idx) {
//
//case 0:
//	/* 
//	 * laplacian
//	 */
//	Laplacian lp = new Laplacian();
//	
//	for (int image_idx = 0; image_idx < imgs.size(); image_idx++) {
//		
//		System.out.println("extracting laplacians for image "+image_idx+"/"+(imgs.size()-1));
//		
//		Vector<Image> laplacians = lp.run(imgs.get(image_idx), sigma);
//		
//		for (int loc_idx = 0; loc_idx < locs.get(image_idx).numInstances(); loc_idx++) {
//		
//			for (int scale_idx = 0; scale_idx < laplacians.size(); scale_idx++) {
//			
//				int col = (int)Math.round( locs.get(image_idx).instance(loc_idx).value(0) );
//				int row = (int)Math.round( locs.get(image_idx).instance(loc_idx).value(1) );
//				Coordinates at_pos = new Coordinates(col, row);
//				double value = laplacians.get(scale_idx).get(at_pos);
//				train.instance(fill_up_idx).setValue(scale_idx, value);
//				
//			}
//			
//			train.instance(fill_up_idx).setValue(laplacians.size(), clss.get(image_idx).instance(loc_idx).value(0));
//			fill_up_idx++;
//			
//		}
//		
//	}
//	
//	try {
//		DataSink.write(train_path, train);
//		System.out.println("done, saved to "+train_path);
//	} catch (Exception e) {
//		e.printStackTrace();
//	}
//	
//
//	break;
//case 1:
//	/*
//	 *  DoH = l1*l2
//	 */
//	DoH doh = new DoH();
//	
//	for (int image_idx = 0; image_idx < imgs.size(); image_idx++) {
//		
//		System.out.println("extracting DoH for image "+image_idx+"/"+(imgs.size()-1));
//		
//		Vector<Image> dohs = doh.run(imgs.get(image_idx), sigma);
//		
//		for (int loc_idx = 0; loc_idx < locs.get(image_idx).numInstances(); loc_idx++) {
//			
//			for (int scale_idx = 0; scale_idx < dohs.size(); scale_idx++) {
//				
//				int col = (int)Math.round( locs.get(image_idx).instance(loc_idx).value(0) );
//				int row = (int)Math.round( locs.get(image_idx).instance(loc_idx).value(1) );
//				Coordinates at_pos = new Coordinates(col, row);
//				double value = dohs.get(scale_idx).get(at_pos);
//				train.instance(fill_up_idx).setValue(scale_idx, value);
//			
//			}
//			
//			train.instance(fill_up_idx).setValue(dohs.size(), clss.get(image_idx).instance(loc_idx).value(0));
//			fill_up_idx++;
//		
//		}
//		
//	}
//	
//	try {
//		DataSink.write(train_path, train);
//		System.out.println("done, saved to "+train_path);
//	} catch (Exception e) {
//		e.printStackTrace();
//	}
//	
//	break;
//case 2:
//	/*
//	 *  |l1|, Ballness
//	 */
//	Ballness bness = new Ballness();
//	
//	for (int image_idx = 0; image_idx < imgs.size(); image_idx++) {
//		
//		System.out.println("extracting |l1| for image "+image_idx+"/"+(imgs.size()-1));
//		
//		Vector<Image> bnesses = bness.run(imgs.get(image_idx), sigma);
//		
//		for (int loc_idx = 0; loc_idx < locs.get(image_idx).numInstances(); loc_idx++) {
//			
//			for (int scale_idx = 0; scale_idx < bnesses.size(); scale_idx++) {
//				
//				int col = (int)Math.round( locs.get(image_idx).instance(loc_idx).value(0) );
//				int row = (int)Math.round( locs.get(image_idx).instance(loc_idx).value(1) );
//				Coordinates at_pos = new Coordinates(col, row);
//				double value = bnesses.get(scale_idx).get(at_pos);
//				train.instance(fill_up_idx).setValue(scale_idx, value);
//			
//			}
//			
//			train.instance(fill_up_idx).setValue(bnesses.size(), clss.get(image_idx).instance(loc_idx).value(0));
//			fill_up_idx++;
//		
//		}
//		
//	}


//public static void main(String[] args) throws Exception {
//
//File f 				= new File(args[0]);
//double[][] f_data 	= extractCols(f, new int[]{0, 1});
//System.out.println(Utils.arrayToString(f_data));
//}	

//CSVLoader loader = new CSVLoader();
//loader.setSource(new File(args[0]));
//Instances gnd_tth = loader.getDataSet();
//
//// remove last attribute 
//Remove rm = new Remove();
//rm.setAttributeIndices(""+gnd_tth.numAttributes());
//rm.setInputFormat(gnd_tth);
//Instances gnd_tth_locs = Filter.useFilter(gnd_tth, rm);
//double[][] a = new double[gnd_tth_locs.numInstances()][gnd_tth_locs.numAttributes()];
//for (int i = 0; i < gnd_tth_locs.numInstances(); i++) {
//	a[i] = gnd_tth_locs.instance(i).toDoubleArray();
//}
//Instances data0 = DataSource.read( (new File(args[0])).getAbsolutePath() );
//Instances data1 = DataSource.read( (new File(args[1])).getAbsolutePath() );
//System.out.println("data0:\n"+data0.numInstances()+" x "+data0.numAttributes());
//System.out.println("data1:\n"+data1.numInstances()+" x "+data1.numAttributes());
//String[] options = Utils.splitOptions("help");//new String[]{"append"};
//Instances.main(options);
//Instances data  = Instances.mergeInstances(data0, data1);
//System.out.println(" *** AFTER MERGING: *** \n\n"+data);
//private static double[][] extractCols(File csv_file, int[] attribs_to_keep) throws Exception {
//
//CSVLoader loader = new CSVLoader();
//loader.setSource(csv_file);
//Instances data = loader.getDataSet();
//
//// remove last attribute 
//Remove rm = new Remove();
//rm.setInvertSelection(true); // columns are kept
//rm.setAttributeIndicesArray(attribs_to_keep);
//rm.setInputFormat(data);
//Instances data_no_last_att = Filter.useFilter(data, rm);
//
//// convert to double[][]
//
//double[][] out = new double[data_no_last_att.numInstances()][data_no_last_att.numAttributes()];
//
//for (int i = 0; i < data_no_last_att.numInstances(); i++) {
//	out[i] = data_no_last_att.instance(i).toDoubleArray();
//}
//
//return out;
//
//}
