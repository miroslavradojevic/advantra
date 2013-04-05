package advantra.critpoint;

import java.util.Vector;

import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.measure.Calibration;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import advantra.feature.GaborFilt2D;
import advantra.general.Transf;

public class TestGabor implements PlugInFilter {

	ImagePlus img;
	
	public void run(ImageProcessor arg0) {
		
		double 	s1, s2;
		int 	sn;
		double 	lambda 		= 0; // will be determined by sigma
		double 	bandwidth 	= 1; // to correlate lambda & sigma
		int 	M;
		double 	psi 		= 0;
		
		s1 		= Prefs.get("advantra.critpoint.start_scale", 	2.0);
		s2    	= Prefs.get("advantra.critpoint.end_scale", 	5.0);
		sn   	= (int)Prefs.get("advantra.critpoint.nr_scales", 4);
		M 		= (int)Prefs.get("advantra.critpoint.nr_angles", 8);
		
		GenericDialog gd = new GenericDialog("Gabor demo");
		gd.addNumericField( "sigma start:", s1, 	 	2, 5, "");
		gd.addNumericField( "sigma end  :", s2, 	 	2, 5, "");
		gd.addNumericField( "sigma nr   :", sn, 	 	0, 5, "");
		gd.addNumericField( "M:    ", 		M, 	 		0, 5, "");
		gd.showDialog();
		if (gd.wasCanceled()) return;
		s1 		= gd.getNextNumber();
		s2 		= gd.getNextNumber();
		sn 		= (int)gd.getNextNumber();
		M 		= (int)gd.getNextNumber();
		
		Prefs.set("advantra.critpoint.start_scale", s1);
		Prefs.set("advantra.critpoint.end_scale", s2);
		Prefs.set("advantra.critpoint.nr_scales", sn);
		Prefs.set("advantra.critpoint.nr_angles", M);
		
		// reset calibration before going further
		Calibration cal = img.getCalibration();
		cal.pixelWidth = cal.pixelHeight = cal.pixelDepth = 1.0;
		cal.setUnit("pixel");
		img.setCalibration(cal);
		
		// form scales vector
		double[] s = new double[sn];
		for (int i = 0; i < sn; i++) {
			s[i] = (i==0)? s1 : s1+i*((s2-s1)/(sn-1));
		}
		
	    if(!img.getProcessor().isDefaultLut()) return;
		img.setProcessor("testImage", img.getProcessor().convertToFloat().duplicate());
	    
		Vector<ImagePlus> g = GaborFilt2D.runAll(img, M, s, lambda, bandwidth, psi, true);
		
		for (int i = 0; i < g.size(); i++) {
			g.get(i).show();
		}
		
		// now design filter 
		// filter params
		double 		angStep = (10/180f)*Math.PI;
		int 		Nang 	= 36;
		double[] 	ro 		= new double[]{0, 3, 6};
		
		double[][] rofi = new double[(ro.length-1)*Nang+1][];
		rofi[0] = new double[]{ro[0], 0};
		
		int idx = 1;
		
		for (int i = 1; i < ro.length; i++) {
			
			for (int j = 0; j < Nang; j++) {
				
				double ang = j*(2*Math.PI/Nang);
				rofi[idx] = new double[]{ro[i], ang};
				System.out.println("rofi["+idx+"]: "+ro[i]+" , "+ang);
				idx++;
				
			}
			
		}
		
		
		
		// create ro/fi pairs
//		
		
//		for (int i = 0; i < fi.length; i++) {
//			fi[i] = i * (2*Math.PI) / M1;
//		}
		// configuration (extracts tuples, considers surrounding locations)
		
		
//		ImagePlus outIm = GaborFilt2D.run(img, M, sigma, lambda, bandwidth, psi, false);
//		outRe.show();
//		outIm.show();

		
		
//		ImagePlus projectStack = new ImagePlus("filtered stack",filtered);
//		ImageStack resultStack = new ImageStack(width, height);
		/*
		ZProjector zp = new ZProjector(projectStack);
		zp.setStopSlice(is.getSize());
		for (int i=0;i<=nAngles; i++)
		{
		    zp.setMethod(i);
		    zp.doProjection();
		    resultStack.addSlice("Gabor_" + i 
		            +"_"+sigma+"_" + gamma + "_"+ (int) (psi / (Math.PI/4) ) +"_"+Fx, 
		            zp.getProjection().getChannelProcessor());
		}
		*/
		
	}

	public int setup(String arg0, ImagePlus arg1) {
		
		img = arg1;
		
		return DOES_8G+NO_CHANGES;
	}
	
}
//	ImageProcessor ip = img.getProcessor();
//	float[] pix = (float[])ip_float.getPixels();
//	float[] axis_x = new float[pix.length];
//	for (int i = 0; i < axis_x.length; i++) {
//		axis_x[i] = i;
//	}
//	
//	IJ.log("float "+ip_float.getMax()+"     ...     "+ip_float.maxValue());
//	for (int i = 0; i < 10; i++) {
//	IJ.log(i+ ": " +pix[i]+" , "+ip.getPixelValue(i, 0));
//	}
//	Plot p = new Plot("see it", "sample", "value", new double[0], new double[0]);
//	p.setLimits(0, pix.length, 0, 255);
//	
//	p.addPoints(axis_x, pix, Plot.LINE);
//	p.show();
	


//OpenDialog open_image 	= new OpenDialog("Select image for CP feature extraction", null);
//String image_dir = open_image.getDirectory();
//if(image_dir!=null){
//	if (!image_dir.endsWith("/")) image_dir += "/";
//}
//else{
//	System.out.println("image was null!");
//	return;
//}
//ImagePlus imp = new ImagePlus(image_dir+open_image.getFileName());
//Image im = new FloatImage(Image.wrap(imp));
//imp.show();
///*
// * select filter and scales
// */
//String[] choices = new String[3];
//choices[0] = "laplacian";
//choices[1] = "DoH";
//choices[2] = "abs(L1)";
//
//double[] sigma = new double[nr];
//for (int i = 0; i < nr; i++) {
//	sigma[i] = (i==0)? sigma_1 : sigma_1+i*((sigma_2-sigma_1)/(nr-1));
//}
//
///*
// * extract filter at scales
// */
//
//Image feats = new FloatImage(im.dimensions());// later on it is allocated with full range of scales, layers are added
//
//switch (idx) {
//case 0:	// laplacian
//	feats = Laplacian2D.calculateImg(im, sigma);
//	feats.imageplus().show();
//	break;
//case 1: // DoH
//	feats = DoH2D.calculateImg(im, sigma);
//	feats.imageplus().show();
//	break;
//default:
//	
//	break;
//}
//
///*
// * 
// */
//
//OpenDialog open_csv = new OpenDialog("Select csv file with CP annotations (optional)", null);
//String csv_dir = open_csv.getDirectory();
//if(csv_dir!=null){
//	if (!csv_dir.endsWith("/")) csv_dir += "/";
//	System.out.println("csv was loaded!");
//}
//else{
//	System.out.println("csv was null! won't be extracting the features");
//	return;
//}
//
///*
// * extract features
// */
//
//File csv_file = new File(csv_dir+open_csv.getFileName()); 
//Instances data = loadFromCsvFile(csv_file);
//System.out.println("data class formed!");
//
//Instances 	x;
//Instances	y;
//
//// filters
//Remove keep12 = new Remove();
//Remove keepLast = new Remove();
//
//try {
//	keep12.setInvertSelection(true); 	
//	int[] first_two = new int[]{0, 1};
//	keep12.setAttributeIndicesArray(first_two);
//	keep12.setInputFormat(data);
//
//	keepLast.setInvertSelection(true);
//	int[] last_idx = new int[]{(data.numAttributes()-1)};
//	keepLast.setAttributeIndicesArray(last_idx);
//	keepLast.setInputFormat(data);
//
//	x = Filter.useFilter(data, keep12);
//	y = Filter.useFilter(data, keepLast);
//	
//	ArrayList<Instances> 		fx = new ArrayList<Instances>(nr);
//
//	for (int i = 0; i < nr; i++) {
//	
//		// fx(i) form
//		FastVector attributes = new FastVector();
//		String label = String.format(choices[idx]+"_s_%.2f", sigma[i]);
//		attributes.addElement(new Attribute(label));
//		attributes.addElement(y.attribute(0));
//		Instances train = new Instances(label, attributes, x.numInstances());
//
//		int cnt = 0;
//		for (int loc_idx = 0; loc_idx < x.numInstances(); loc_idx++) {
//		
//			int col = (int)Math.round( x.instance(loc_idx).value(0) );
//			int row = (int)Math.round( x.instance(loc_idx).value(1) );
//			Coordinates at_pos = new Coordinates(col, row, i);
//			double value = feats.get(at_pos);
//			train.instance(cnt).setValue(0, value);
//			train.instance(cnt).setValue(nr, y.instance(loc_idx).value(0));
//			cnt++;
//		
//		}
//	
//	fx.set(i, train);
//	System.out.println("fx("+i+") formed!");
//	
//	}
//	
//} catch (Exception e) {
//	e.printStackTrace();
//}
//FastVector attributes = new FastVector();
//for (int i = 0; i < nr; i++) {
//	String label = String.format(choices[idx]+"_s_%.2f", sigma[i]);
//	attributes.addElement(new Attribute(label));
//	attributes.addElement(Filter.useFilter(data, keepLast).attribute(0));
//	fx.set(i, element);
//}
//try {
//	for (int i = 0; i < locs.numInstances(); i++) {
//		train.add(new Instance(nr+1)); // fill them with missing values
//	}
//	train.setClassIndex(train.numAttributes()-1);
//	
//	int cnt = 0;
//	for (int loc_idx = 0; loc_idx < locs.numInstances(); loc_idx++) {
//	
//		for (int scale_idx = 0; scale_idx < nr; scale_idx++) {
//		
//			int col = (int)Math.round( locs.instance(loc_idx).value(0) );
//			int row = (int)Math.round( locs.instance(loc_idx).value(1) );
//			Coordinates at_pos = new Coordinates(col, row, scale_idx);
//			double value = feats.get(at_pos);
//			train.instance(cnt).setValue(scale_idx, value);
//			
//		}
//		
//		train.instance(cnt).setValue(nr, clss.instance(loc_idx).value(0));
//		cnt++;
//		
//	}
//	
//	String train_path = System.getProperty("user.home");
//	if (!train_path.endsWith("/")) train_path += "/";
//	train_path += "train.csv";
//	
//	System.out.println("wants to save to: "+train_path);
//	
//	DataSink.write(train_path, train);
//	
//	System.out.println("done, saved to "+train_path);
//	
//} catch (Exception e) {
//	e.printStackTrace();
//}
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
