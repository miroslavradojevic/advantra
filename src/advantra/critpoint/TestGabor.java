package advantra.critpoint;

import java.awt.Color;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

import flanagan.analysis.Stat;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.measure.Calibration;
import ij.plugin.ZProjector;
import ij.plugin.filter.PlugInFilter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import imagescience.image.Axes;
import imagescience.image.Coordinates;
import imagescience.image.Dimensions;
import imagescience.image.FloatImage;
import imagescience.image.Image;
import advantra.feature.GaborFilt2D;
import advantra.general.Sort;

public class TestGabor implements PlugInFilter, MouseListener {

	ImagePlus 	img;
	ImagePlus 	circ_img;
	int 		M;
	float 		circ_img_max, circ_img_min;
	
	public void run(ImageProcessor arg0) {
		
		double 	s1, s2;
		int 	sn;
//		double 	lambda 		= 0; // will be determined by sigma
		double 	bandwidth 	= 1; // to correlate lambda & sigma
		double 	psi 		= 0;
		
		s1 		= Prefs.get("advantra.critpoint.start_scale", 	2.0);
		s2    	= Prefs.get("advantra.critpoint.end_scale", 	5.0);
		sn   	= (int)Prefs.get("advantra.critpoint.nr_scales", 4);
		M 		= (int)Prefs.get("advantra.critpoint.nr_angles", 8);
		
		GenericDialog gd = new GenericDialog("Gabor demo");
		gd.addNumericField( "sigma start	:", s1, 	 	2, 5, "");
		gd.addNumericField( "sigma end  	:", s2, 	 	2, 5, "");
		gd.addNumericField( "sigma nr   	:", sn, 	 	0, 5, "");
		gd.addNumericField( "angles per pi	:", M, 	 		0, 5, "");
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
		
		// sigmas
		double[] s = new double[sn];
		for (int i = 0; i < sn; i++) {
			s[i] = (i==0)? s1 : s1+i*((s2-s1)/(sn-1));
		}
		
		double[] thetas_pi = new double[M];
		double[] thetas_2pi = new double[2*M];
		
		for (int i = 0; i < 2*M; i++) {
			if(i<M){
				thetas_2pi[i] = thetas_pi[i] = i*Math.PI/M;
			}
			else{
				thetas_2pi[i] = thetas_2pi[i-M]+Math.PI;
			}
		}
		
	    if(!img.getProcessor().isDefaultLut()) return;
		img.setProcessor("testImage", img.getProcessor().convertToFloat().duplicate());
	    
		// calculate g_theta_lambda for every theta 0-pi
		// take g_ value by maximizing over scales
		
		ImageStack gst = new ImageStack(img.getWidth(), img.getHeight()); 
		
		ZProjector zp = new ZProjector();
		for (int i = 0; i < M; i++) {
			double current_theta = thetas_pi[i];
			ImagePlus g_theta = GaborFilt2D.run(img, current_theta, s, new double[s.length], bandwidth, psi, true);
			zp.setImage(g_theta);
			zp.setStartSlice(1);
			zp.setStopSlice(g_theta.getStackSize());
			zp.setMethod(ZProjector.MAX_METHOD);
			zp.doProjection();
			gst.addSlice("theta="+IJ.d2s(current_theta, 2), zp.getProjection().getChannelProcessor());
		}
		
		ImagePlus gim = new ImagePlus("gabor_responses_per_scale", gst);
		gim.show();
		
		// now design entropy filter 
		
		double		ro 		= 1.5*s2;
		int 		margin 	= (int)Math.ceil(ro);

		ImageStack	circ_stack 	= new ImageStack(img.getWidth(), img.getHeight());
		
		// define a shift for every position on the circle
		for (int i = 0; i < 2*M; i++) {
			
			double dx = ro*Math.sin(thetas_2pi[i]);
			double dy = -ro*Math.cos(thetas_2pi[i]);
			
//			System.out.print("\n"+IJ.d2s(dx)+", "+IJ.d2s(dy));
			
			ImageProcessor 	circ_vals 	= new FloatProcessor(img.getWidth(), img.getHeight());
			
			// set values of circ_vals
			for (int x = margin; x < img.getWidth()-margin; x++) {
				for (int y = margin; y < img.getHeight()-margin; y++) {
					int src_x = (int)Math.round(x+dx);
					int src_y = (int)Math.round(y+dy);
					
					if(i<M){
						circ_vals.setf(x, y, gim.getStack().getProcessor(i+1).getPixelValue(src_x, src_y)); // * gim.getStack().getProcessor(i+1).getPixelValue(x, y)
					}
					else{
						circ_vals.setf(x, y, gim.getStack().getProcessor(i-M+1).getPixelValue(src_x, src_y)); // * gim.getStack().getProcessor(i-M+1).getPixelValue(x, y)
					}
					
				}
			}
			
			// add it to the output stack
			circ_stack.addSlice(circ_vals);
//			System.out.print("done.");
			
		}
		
		circ_img = new ImagePlus("circular", circ_stack);
		circ_img.show();
		circ_img.getCanvas().addMouseListener(this);
		img.getCanvas().addMouseListener(this);
		
		// find min/max of the total circular response
		circ_img_min = Float.MAX_VALUE;
		circ_img_max = Float.MIN_VALUE;
		
		for (int z = 0; z < circ_img.getStackSize(); z++) {
			for (int x = 0; x < circ_img.getWidth(); x++) {
				for (int y = 0; y < circ_img.getHeight(); y++) {
					
					float take_val = circ_img.getStack().getProcessor(z+1).getPixelValue(x, y);
					if(take_val<circ_img_min){
						circ_img_min = take_val;
					}
					if(take_val>circ_img_max){
						circ_img_max = take_val;
					}
					
				}
			}
		}
		
//		// Suppress those lower than some ratio of the max
//		boolean suppress = false;
//		if(suppress){
//			for (int z = 0; z < circ_img.getStackSize(); z++) {
//				for (int x = 0; x < circ_img.getWidth(); x++) {
//					for (int y = 0; y < circ_img.getHeight(); y++) {
//						
//						float take_val = circ_img.getStack().getProcessor(z+1).getPixelValue(x, y);
//						if(take_val<0.1*circ_img_max){
//							circ_img.getStack().getProcessor(z+1).setf(x, y, 0);
//						}
//						
//					}
//				}
//			}
//		}
		
		// extract the circular features out the circ_img with angular responses
		Image inimg = Image.wrap(circ_img);
		inimg.axes(Axes.Z);
		Dimensions outd = new Dimensions(inimg.dimensions().x, inimg.dimensions().y);
		Image outimg_entropy = new FloatImage(outd);
		outimg_entropy.axes(Axes.X+Axes.Y);
		
		double[] circ_vals = new double[inimg.dimensions().z];
		Coordinates coord = new Coordinates();
		
		for (coord.x = 0; coord.x < inimg.dimensions().x; coord.x++) {
			for (coord.y = 0; coord.y < inimg.dimensions().y; coord.y++) {
				inimg.get(coord, circ_vals);
				
				double circ_vals_min = Sort.findMin(circ_vals);
				double circ_vals_max = Sort.findMax(circ_vals);
				double entropy = 0;
				
				if(Math.abs(circ_vals_max-circ_vals_min)>0.0001){
				
				double[][] distr = Stat.histogramBins(
			    		circ_vals, 
			    		(circ_vals_max-circ_vals_min)/10, 
			    		circ_vals_min, 
			    		circ_vals_max);
				
			    double hist_norm = Sort.sum(distr[1]);
			    for (int i = 0; i < distr[1].length; i++) {
					double p = distr[1][i]/hist_norm;
			    	entropy += p * Math.log(p);
				}
			    
				}
				
			    outimg_entropy.set(coord, -entropy);
			    
			}
		}
		outimg_entropy.imageplus().show();
		
	}

	public int setup(String arg0, ImagePlus arg1) {
		
		img = arg1;
		return DOES_8G+NO_CHANGES;
		
	}

	public void mouseClicked(MouseEvent e) {

		// location
		int mouseX = 	img.getWindow().getCanvas().offScreenX(e.getX());
		int mouseY = 	img.getWindow().getCanvas().offScreenY(e.getY());
		
		// angles
		double[] phis = new double[2*M];
		
		for (int i = 0; i < 2*M; i++) {
			phis[i] = i*(2*Math.PI)/(2*M);
		}
		
		Image read_img = Image.wrap(circ_img);
		read_img.axes(Axes.Z);
		double[] circ_vals = new double[2*M];
		Coordinates coord = new Coordinates(mouseX, mouseY);
		read_img.get(coord, circ_vals);
		
		double circ_vals_min = Sort.findMin(circ_vals);
		double circ_vals_max = Sort.findMax(circ_vals);
		double[] circ_vals_sum_pos_neg = Sort.sum_pos_neg(circ_vals);
		
		// min-max normalize
		double[] mn_mx = new double[circ_vals.length];
		for (int i = 0; i < circ_vals.length; i++) {
			mn_mx[i] = (circ_vals[i]-circ_vals_min)/circ_vals_max;
		}
		// normalize
		double[] norm = new double[circ_vals.length];
		for (int i = 0; i < circ_vals.length; i++) {
			if(circ_vals[i]>=0){
				norm[i] = circ_vals[i]/circ_vals_sum_pos_neg[0];
			}
			else{
				norm[i] = circ_vals[i]/circ_vals_sum_pos_neg[1];
			}
		}
		
		// plot
	    ImagePlus 	show_gabor_img = new ImagePlus();
	    ImageStack 	show_gabor_stk = new ImageStack(528, 255);
	    
	    Plot plot = new Plot("At,X="+mouseX+",Y="+mouseY, "angle [rad]", "gabor filter output");
	    plot.setLimits(0, 2*Math.PI, circ_img_min, circ_img_max); 
	    plot.setColor(Color.RED);
	    plot.addPoints(phis, circ_vals, Plot.LINE);
	    show_gabor_stk.addSlice("orig", plot.getProcessor());
	    
	    plot = new Plot("At,X="+mouseX+",Y="+mouseY, "angle [rad]", "gabor filter output");
	    plot.setLimits(0, 2*Math.PI, 0, 1); 
	    plot.setColor(Color.BLUE);
	    plot.addPoints(phis, mn_mx, Plot.LINE);
	    show_gabor_stk.addSlice("min max normalized", plot.getProcessor());
	    
	    plot = new Plot("At,X="+mouseX+",Y="+mouseY, "angle [rad]", "gabor filter output");
	    plot.setLimits(0, 2*Math.PI, -1, 1); 
	    plot.setColor(Color.GREEN);
	    plot.addPoints(phis, norm, Plot.LINE);
	    show_gabor_stk.addSlice("normalized", plot.getProcessor());
	    
	    double[][] distr = Stat.histogramBins(
	    		circ_vals, 
	    		(circ_vals_max-circ_vals_min)/10, 
	    		circ_vals_min, 
	    		circ_vals_max);
	    
	    plot = new Plot("At,X="+mouseX+",Y="+mouseY, "angle [rad]", "gabor filter output");
	    plot.setLimits(Sort.findMin(distr[0]), Sort.findMax(distr[0]), Sort.findMin(distr[1]), Sort.findMax(distr[1])); 
	    plot.setColor(Color.YELLOW);
	    plot.addPoints(distr[0], distr[1], Plot.LINE);
	    show_gabor_stk.addSlice("histogram", plot.getProcessor());

	    show_gabor_img.setStack("gabor data", show_gabor_stk);
	    show_gabor_img.show();
	    
	}

	public void mousePressed(MouseEvent e) {}

	public void mouseReleased(MouseEvent e) {}

	public void mouseEntered(MouseEvent e) {}

	public void mouseExited(MouseEvent e) {}
	
}

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
