package advantra.critpoint;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import imagescience.image.Coordinates;
import imagescience.image.FloatImage;
import imagescience.image.Image;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Vector;

import advantra.feature.Ballness;
import advantra.feature.DoH;
import advantra.feature.Laplacian;

import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.converters.CSVLoader;
import weka.core.converters.ConverterUtils.DataSink;
import weka.filters.Filter;
import weka.filters.unsupervised.attribute.Remove;

public class ExtractFeatures implements PlugIn {

	String 	dir_path 		= System.getProperty("user.home")+File.separator+"train";
	File[] csv_files 		= null;
	File[] tif_files 		= null;
	
	public void run(String arg0) {
		
		GenericDialog gd = new GenericDialog("Critical Point Feats.");
		gd.addStringField("folder with annotations", dir_path, 50);
		gd.addMessage("Choose features");
		String[] choices = new String[3];
		choices[0] = "laplacian";
		choices[1] = "DoH";
		choices[2] = "|lambda1|";
		
		gd.addChoice("Feat.", choices, choices[0]);
		gd.addNumericField( "sigma start:", 			2, 	 0, 5, "" );
		gd.addNumericField( "sigma end  :", 			6, 	 0, 5, "" );	
		gd.addNumericField( "number of scales : ", 		10,  0, 5, "");
		
		gd.showDialog();
		if (gd.wasCanceled()) return;
		
		dir_path 					=  gd.getNextString();
		int idx 					=  gd.getNextChoiceIndex();
		double 		sigma_1 		=  (double)gd.getNextNumber();
		double 		sigma_2 		=  (double)gd.getNextNumber();
		int 		nr				=  (int)gd.getNextNumber();
		
		// scales
		double[] sigma = new double[nr];
		for (int i = 0; i < nr; i++) {
			sigma[i] = (i==0)? sigma_1 : sigma_1+i*((sigma_2-sigma_1)/(nr-1));
		}
		
		dir_path = (new File(dir_path)).getAbsolutePath();
		File dir = new File(dir_path);
		if(!dir.isDirectory()){
			IJ.error("Wrong source directory!");
			return;
		}
		
		// csv files
		csv_files = dir.listFiles(
				new FilenameFilter() {
					public boolean accept(File dir, String name) {
						return name.toLowerCase().endsWith(".csv");
		    	}
			}
		);
		// tif files
		tif_files = dir.listFiles(
				new FilenameFilter() {
					public boolean accept(File dir, String name) {
						return name.toLowerCase().endsWith(".tif");
				}
			}
		);
		
		if(csv_files.length<=0){
			IJ.showMessage("\nThere was no csv files in source.\n");
			return;
		}

		//loop through train set
		ArrayList<Instances> 		clss 		= new ArrayList<Instances>();
		ArrayList<Instances> 		locs 		= new ArrayList<Instances>();
		ArrayList<Image> 			imgs 		= new ArrayList<Image>();
		
		for (int i = 0; i < csv_files.length; i++) { // for each .csv file
			
			String 	current_csv_name 	= csv_files[i].getName();
			System.out.print(current_csv_name+" ... ");
			boolean found = false;
			int idx_found = 0; // is there a .tif pair, remember it's index
			
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
			
			// match was found, add to the extraction list
			
			Instances data = loadFromCsvFile(csv_files[i]);
			// data12
			Remove keep12 = new Remove();
			keep12.setInvertSelection(true); 		// attribs are kept
			int[] first_two = new int[]{0, 1};
			keep12.setAttributeIndicesArray(first_two);
			try {
				keep12.setInputFormat(data);
				locs.add(Filter.useFilter(data, keep12));
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
			
			ImagePlus ip 	= new ImagePlus(tif_files[idx_found].getAbsolutePath());
			// reset calibration before going further
			Calibration cal = new Calibration(ip);
			cal.pixelWidth = cal.pixelHeight = cal.pixelDepth = 1.0;
			cal.setUnit("pixel");
			ip.setCalibration(cal);
			Image im = new FloatImage(Image.wrap(ip));
			
			imgs.add(im);
			
			System.out.println("stored locations, classes, and image"); // imgs, clss, locs
			
		}
		
		
		
		if(imgs.size()<=0 || clss.size()<=0 || locs.size()<=0){
			System.out.println("No images/locations to extract features from.");
			return;
		}

		int nr_locs = 0;
		for (int j = 0; j < locs.size(); j++) {
			nr_locs += locs.get(j).numInstances();
		}
		
		/*
		 * prepare the dataset
		 */
		FastVector attributes = new FastVector();
		
		for (int i = 0; i < nr; i++) {
			String label = String.format("scale_%.2f", sigma[i]);
			attributes.addElement(new Attribute(label));
		}
		
		// add the class attribute
		attributes.addElement(clss.get(0).attribute(0));
		
		Instances train = new Instances("train", attributes, nr_locs);
		for (int i = 0; i < nr_locs; i++) {
			train.add(new Instance(nr+1)); // fill them with missing values
		}
		train.setClassIndex(train.numAttributes()-1);

		String train_path = System.getProperty("user.home")+File.separator+"train.csv";
		
		// fill it with values 
		int fill_up_idx = 0;
		
		switch (idx) {
		
		case 0:
			/* 
			 * laplacian
			 */
			Laplacian lp = new Laplacian();
			
			for (int image_idx = 0; image_idx < imgs.size(); image_idx++) {
				
				System.out.println("extracting laplacians for image "+image_idx+"/"+(imgs.size()-1));
				
				Vector<Image> laplacians = lp.run(imgs.get(image_idx), sigma);
				
				for (int loc_idx = 0; loc_idx < locs.get(image_idx).numInstances(); loc_idx++) {
				
					for (int scale_idx = 0; scale_idx < laplacians.size(); scale_idx++) {
					
						int col = (int)Math.round( locs.get(image_idx).instance(loc_idx).value(0) );
						int row = (int)Math.round( locs.get(image_idx).instance(loc_idx).value(1) );
						Coordinates at_pos = new Coordinates(col, row);
						double value = laplacians.get(scale_idx).get(at_pos);
						train.instance(fill_up_idx).setValue(scale_idx, value);
						
					}
					
					train.instance(fill_up_idx).setValue(laplacians.size(), clss.get(image_idx).instance(loc_idx).value(0));
					fill_up_idx++;
					
				}
				
			}
			
			try {
				DataSink.write(train_path, train);
				System.out.println("done, saved to "+train_path);
			} catch (Exception e) {
				e.printStackTrace();
			}
			

			break;
		case 1:
			/*
			 *  DoH = l1*l2
			 */
			DoH doh = new DoH();
			
			for (int image_idx = 0; image_idx < imgs.size(); image_idx++) {
				
				System.out.println("extracting DoH for image "+image_idx+"/"+(imgs.size()-1));
				
				Vector<Image> dohs = doh.run(imgs.get(image_idx), sigma);
				
				for (int loc_idx = 0; loc_idx < locs.get(image_idx).numInstances(); loc_idx++) {
					
					for (int scale_idx = 0; scale_idx < dohs.size(); scale_idx++) {
						
						int col = (int)Math.round( locs.get(image_idx).instance(loc_idx).value(0) );
						int row = (int)Math.round( locs.get(image_idx).instance(loc_idx).value(1) );
						Coordinates at_pos = new Coordinates(col, row);
						double value = dohs.get(scale_idx).get(at_pos);
						train.instance(fill_up_idx).setValue(scale_idx, value);
					
					}
					
					train.instance(fill_up_idx).setValue(dohs.size(), clss.get(image_idx).instance(loc_idx).value(0));
					fill_up_idx++;
				
				}
				
			}
			
			try {
				DataSink.write(train_path, train);
				System.out.println("done, saved to "+train_path);
			} catch (Exception e) {
				e.printStackTrace();
			}
			
			break;
		case 2:
			/*
			 *  |l1|, Ballness
			 */
			Ballness bness = new Ballness();
			
			for (int image_idx = 0; image_idx < imgs.size(); image_idx++) {
				
				System.out.println("extracting |l1| for image "+image_idx+"/"+(imgs.size()-1));
				
				Vector<Image> bnesses = bness.run(imgs.get(image_idx), sigma);
				
				for (int loc_idx = 0; loc_idx < locs.get(image_idx).numInstances(); loc_idx++) {
					
					for (int scale_idx = 0; scale_idx < bnesses.size(); scale_idx++) {
						
						int col = (int)Math.round( locs.get(image_idx).instance(loc_idx).value(0) );
						int row = (int)Math.round( locs.get(image_idx).instance(loc_idx).value(1) );
						Coordinates at_pos = new Coordinates(col, row);
						double value = bnesses.get(scale_idx).get(at_pos);
						train.instance(fill_up_idx).setValue(scale_idx, value);
					
					}
					
					train.instance(fill_up_idx).setValue(bnesses.size(), clss.get(image_idx).instance(loc_idx).value(0));
					fill_up_idx++;
				
				}
				
			}
			
			try {
				DataSink.write(train_path, train);
				System.out.println("done, saved to "+train_path);
			} catch (Exception e) {
				e.printStackTrace();
			}
			
			break;
		default:
			System.out.println("nothing...");
			break;
		}
		
		/// just debug
		//System.out.println(locs.get(0));
		
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
