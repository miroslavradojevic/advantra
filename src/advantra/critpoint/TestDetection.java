package advantra.critpoint;

import java.awt.Color;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.File;
import java.util.ArrayList;

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
import ij.process.ImageProcessor;
import imagescience.feature.Differentiator;
import imagescience.image.Axes;
import imagescience.image.ByteImage;
import imagescience.image.Coordinates;
import imagescience.image.Dimensions;
import imagescience.image.FloatImage;
import imagescience.image.Image;
import imagescience.transform.Translate;
import advantra.feature.CircHaarFeat;
import advantra.feature.GaborFilt2D;
import advantra.file.AnalyzeCSV;
import advantra.general.Sort;
import advantra.plugins.MyOpener;
import advantra.tools.AdaBoost;

public class TestDetection implements PlugInFilter, MouseListener {

	ImagePlus 	img;
	ImagePlus	scales;
	ImagePlus 	profile;
	
	double[] 	t; // scales (gauss var.)
	double[]	s; // gauss std.
	double		r;
	int 		M;
	
	double[] thetas_pi;
	double[] thetas_2pi;
	
	ImagePlus 	gabor_circ;
	ImagePlus 	gabor_angle;
	ImagePlus 	orig_circ;
	
	float 		circ_img_max, circ_img_min;
	
	float 		gabor_angle_min, gabor_angle_max;
	
	ArrayList<float[]> 	pos_fts;
	ArrayList<float[]>	neg_fts;
	
	public void run(ImageProcessor arg0) {
		
		/*
		 *  reset calibration
		 */
		Calibration cal = img.getCalibration();
		cal.pixelWidth = cal.pixelHeight = cal.pixelDepth = 1.0;
		cal.setUnit("pixel");
		img.setCalibration(cal);
		/*
		 * reset calibration
		 */
		
		double 	t1, t2;
		int 	tn;
		double 	bandwidth 	= 1; // to correlate lambda & sigma
		double 	psi 		= 0;
		
		t1 		= Prefs.get("advantra.critpoint.start_scale", 	2.0);
		t2    	= Prefs.get("advantra.critpoint.end_scale", 	5.0);
		tn   	= (int)Prefs.get("advantra.critpoint.nr_scales", 4);
		M 		= (int)Prefs.get("advantra.critpoint.nr_angles", 8);
		
		GenericDialog gd = new GenericDialog("Gabor demo");
		gd.addNumericField( "sigma start(t1):", t1, 	 	2, 5, "");
		gd.addNumericField( "sigma end  (t2):", t2, 	 	2, 5, "");
		gd.addNumericField( "sigma nr   	:", tn, 	 	0, 5, "");
		gd.addNumericField( "angles(per pi)	:", M, 	 		0, 5, "");
		gd.showDialog();
		if (gd.wasCanceled()) return;
		t1 		= gd.getNextNumber();
		t2 		= gd.getNextNumber();
		tn 		= (int)gd.getNextNumber();
		M 		= (int)gd.getNextNumber();
		
		Prefs.set("advantra.critpoint.start_scale", t1);
		Prefs.set("advantra.critpoint.end_scale", 	t2);
		Prefs.set("advantra.critpoint.nr_scales", 	tn);
		Prefs.set("advantra.critpoint.nr_angles", 	M);
		
		String pos_path = MyOpener.open("Open Positives");
		pos_path = new File(pos_path).getAbsolutePath();
		String neg_path = MyOpener.open("Open Negatives");
		neg_path = new File(neg_path).getAbsolutePath();
		
		AnalyzeCSV loader;
		double[][] pos_locs;
		double[][] neg_locs;
		
		loader = new AnalyzeCSV(pos_path);
		pos_locs = loader.readLn(2);
		
		loader = new AnalyzeCSV(neg_path);
		neg_locs = loader.readLn(2);
		
		// scales of the scale-space
		t = new double[tn];
		for (int i = 0; i < tn; i++) {
			t[i] = (i==0)? t1 : t1+i*((t2-t1)/(tn-1));
		}
		s = new double[tn];
		for (int i = 0; i < t.length; i++) {
			s[i] = Math.sqrt(t[i]);
		}
		
		thetas_pi = new double[M];
		thetas_2pi = new double[2*M];
		
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
	    
		img.getCanvas().addMouseListener(this);
		
		int w = img.getWidth();
		int h = img.getHeight();

		/*
		 * main
		 */
		// scale-space analysis
		Image inimg = Image.wrap(img);
		Differentiator 	df 		= new Differentiator();
		Coordinates 	coord 	= new Coordinates();
		Dimensions 		dim 	= inimg.dimensions();
		Dimensions		outdim	= new Dimensions(dim.x, dim.y, tn);
		Image 			Lnorm 	= new FloatImage(outdim);
		
		for (int i = 0; i < tn; i++) {
			
			Image Lxx = df.run(inimg.duplicate(), s[i], 2, 0, 0);
			Image Lyy = df.run(inimg.duplicate(), s[i], 0, 2, 0);
			Image Lxy = df.run(inimg.duplicate(), s[i], 1, 1, 0);
			
			for (coord.x = 0; coord.x < outdim.x; coord.x++) {
				for (coord.y = 0; coord.y < outdim.y; coord.y++) {
					
					coord.z = 0;
					double val = 
							Math.pow(t[i], 3) * 
							Math.pow((Lyy.get(coord) + Lxx.get(coord)),2) * 
							(Math.pow((Lxx.get(coord) - Lyy.get(coord)),2) + 4 * Math.pow(Lxy.get(coord), 2));
					coord.z = i;
					Lnorm.set(coord, val);
					
				}
			}
			
		}
		
		Lnorm.name("Normalized_Ridgeness");
		Lnorm.imageplus().show();
		
		/*
		 * show gabor kernels
		 */
		
		for (int i = 0; i < t.length; i++) {
			ImagePlus krn = GaborFilt2D.getKernel(M, t[i], 0, bandwidth, psi, true);
			krn.setTitle("gabor_kernel_scale"+IJ.d2s(t[i],2));
			krn.show();
			for (int j = 0; j < 5; j++) {
				krn.getCanvas().zoomIn(0, 0);
			}
		}
		
		ImageStack 	angul = new ImageStack(w, h); 
		
		ZProjector zmax = new ZProjector();
		
		for (int i = 0; i < M; i++) {
			
			double current_theta = thetas_pi[i];
			
			System.out.println("processing theta = "+i+" / "+ (M-1));
			
			ImagePlus g_theta = GaborFilt2D.run(
					img, current_theta, t, new double[t.length], bandwidth, psi, true);
			
			zmax.setImage(g_theta);
			zmax.setStartSlice(1);
			zmax.setStopSlice(g_theta.getStackSize());
			zmax.setMethod(ZProjector.MAX_METHOD);
			zmax.doProjection();
			angul.addSlice("gabor,theta="+IJ.d2s(current_theta, 2), zmax.getProjection().getChannelProcessor());

		}
		
		gabor_angle = new ImagePlus("gabor_responses_per_angle", angul);
		gabor_angle.show();
		
		// find max of the angular responses
		ZProjector zp = new ZProjector(gabor_angle);
		zp.setStartSlice(1);
		zp.setStopSlice(gabor_angle.getStackSize());
		zp.setMethod(ZProjector.MAX_METHOD);
		zp.doProjection();
		ImagePlus gabor_angle_proj = new ImagePlus("max_gab_ang_sc", zp.getProjection().getChannelProcessor()); 
		gabor_angle_proj.show();
		
		// extract profiles
		double surr = t[t.length-1];
		
		Translate tr = new Translate();
		
		ImageStack profile_template = new ImageStack(w, h);
		
		for (int j = 0; j < 2*M; j++) {
			
			double dx = surr * Math.sin(thetas_2pi[j]);
			double dy = surr * Math.cos(thetas_2pi[j]);
			
			int stack_idx = (j>=M)?(j-M):j;
			ImagePlus img_to_shift = new ImagePlus("", gabor_angle.getStack().getProcessor(stack_idx+1));
			Image to_shift = Image.wrap(img_to_shift); 
			
			profile_template.addSlice(tr.run(to_shift, dx, dy, 0, Translate.CUBIC).imageplus().getProcessor());
		}
		
		profile = new ImagePlus("profile", profile_template);
		profile.show();
		
		float[] min_max = calculateStackMinMax(profile.getStack());
		circ_img_min = min_max[0];
		circ_img_max = min_max[1];
		
		CircHaarFeat circft = new CircHaarFeat(2*M);
		circft.createFeatures();
		circft.showFeatures();
		
		/*
		 * extract features for pos. locations
		 */
		
		System.out.print("extract features for positives ");
		
		Image profile_img = Image.wrap(profile);
		profile_img.axes(Axes.Z);
		
		int nrPos = pos_locs.length;
		double[][] profilePos 	= new double[nrPos][];
		float[][] featsPos 		= new float[nrPos][];
		double[] vals = new double[profile.getStackSize()];
		coord = new Coordinates();
		
		for (int i = 0; i < nrPos; i++) {
			
			coord.x = (int)Math.round(pos_locs[i][0]);
			coord.y = (int)Math.round(pos_locs[i][1]);
			profile_img.get(coord, vals);
			
			profilePos[i] = vals;
			
			// extract features for the profile
			featsPos[i] = circft.allFeatScore(profilePos[i], 1);
			
		}
		System.out.println("done.");
		
		/* 
		 * extract features for neg. locations
		 */
		
		System.out.print("extract features for negatives ");
		
		int nrNeg = neg_locs.length;
		double[][] profileNeg 	= new double[nrNeg][];
		float[][] featsNeg 		= new float[nrNeg][];
		
		for (int i = 0; i < nrNeg; i++) {
			
			coord.x = (int)Math.round(neg_locs[i][0]);
			coord.y = (int)Math.round(neg_locs[i][1]);
			profile_img.get(coord, vals);
			
			profileNeg[i] = vals;
			
			// extract features for the profile
			featsNeg[i] = circft.allFeatScore(profileNeg[i], 1);
			
		}
		System.out.println("done.");
		
		/*
		 * train
		 */
		
		AdaBoost ab = new AdaBoost(featsPos, featsNeg, 5);
		double[][] adaboost = ab.run();
		// show selected features
		int[] selected_f = new int[adaboost.length];
		for (int i = 0; i < adaboost.length; i++) {
			selected_f[i] = (int)adaboost[i][0];
		}
		circft.showFeatures(selected_f).show();
		
//		if(true) return;
		
		String test_path = MyOpener.open("Open Test Image");
		test_path = new File(test_path).getAbsolutePath();
		
		/*
		 * test
		 */
		ImagePlus test_img = new ImagePlus(test_path);
		if(!test_img.getProcessor().isDefaultLut()) return;
		test_img.setProcessor("testImage", test_img.getProcessor().convertToFloat().duplicate());
		test_img.show();
		
		// calculate g_theta_lambda for every theta 0-pi
		// take g_ value by maximizing response over scales at different angles
		
		ImageStack test_stk = new ImageStack(test_img.getWidth(), test_img.getHeight()); 
		zp  = new ZProjector();
		for (int i = 0; i < M; i++) {
			double current_theta = thetas_pi[i];
			System.out.println("gabor for theta = "+current_theta+" / PI");
			ImagePlus g_theta = GaborFilt2D.run(test_img, current_theta, s, new double[s.length], bandwidth, psi, true);
			zp.setImage(g_theta);
			zp.setStartSlice(1);
			zp.setStopSlice(g_theta.getStackSize());
			zp.setMethod(ZProjector.MAX_METHOD);
			zp.doProjection();
			test_stk.addSlice("theta="+IJ.d2s(current_theta, 2), zp.getProjection().getChannelProcessor());
		}
		
		Image test_profiles = Image.wrap(new ImagePlus("", test_stk));
		test_profiles.axes(Axes.Z);
		
		test_profiles.imageplus().show();
		
		if(true)return;
		
		Dimensions  din = test_profiles.dimensions();
		Dimensions	dout = new Dimensions(din.x, din.y);
		Image 		outimg = new ByteImage(dout);
		double[] 	profileTst 		= new double[din.z];
		float[] 	featsTst;// 		= new float[din.z];
		
		Coordinates cin = new Coordinates();
		for (cin.x = 0; cin.x < din.x; cin.x+=2) { 		 
			for (cin.y = 0; cin.y < din.y; cin.y+=2) { 
				
				test_profiles.get(cin, profileTst);
				
				featsTst = circft.allFeatScore(profileTst, 1);
				
				int outcome = ab.apply(featsTst);
				
                if (outcome == 1) {
                	outimg.set(cin, 255);
                } 
                
//                System.out.print(".");
                
			}
		}
		
		outimg.imageplus().show();
		
		System.out.println("DONE");
		
	}

	public int setup(String arg0, ImagePlus arg1) {
		
		img = arg1;
		return DOES_8G+NO_CHANGES;
		
	}

	private float[] calculateStackMinMax(ImageStack instack){
		
		// find min/max of the total image stack
		float[] mnmx = new float[2];
		mnmx[0] = Float.MAX_VALUE;
		mnmx[1] = Float.MIN_VALUE;
				
		for (int z = 0; z < instack.getSize(); z++) {
			for (int x = 0; x < instack.getWidth(); x++) {
				for (int y = 0; y < instack.getHeight(); y++) {

					float take_val = instack.getProcessor(z+1).getPixelValue(x, y);
						
					if(take_val<mnmx[0]){
						mnmx[0] = take_val;
					}
					
					if(take_val>mnmx[1]){
						mnmx[1] = take_val;
					}
							
				}
			}
		}
		
		return mnmx;
		
	}
	
	public void mouseClicked(MouseEvent e) {

		// location
		int mouseX = 	img.getWindow().getCanvas().offScreenX(e.getX());
		int mouseY = 	img.getWindow().getCanvas().offScreenY(e.getY());
		
		Image read_img = Image.wrap(profile);
		read_img.axes(Axes.Z);
		Dimensions dims = read_img.dimensions();
		
		Coordinates coord = new Coordinates();
		coord.x = mouseX;
		coord.y = mouseY;
		
		double[] xax = new double[dims.z];
		for (int i = 0; i < xax.length; i++) {
			xax[i] = i;
		}
		
		double[] read = new double[dims.z];
		read_img.get(coord, read);
		
		Plot p = new Plot("", "", "");
		p.setLimits(0, xax.length, circ_img_min, circ_img_max);
		p.addPoints(xax, read, Plot.LINE);
		p.draw();
		p.show();
		
	}

	public void mousePressed(MouseEvent e) {}

	public void mouseReleased(MouseEvent e) {}

	public void mouseEntered(MouseEvent e) {}

	public void mouseExited(MouseEvent e) {}
	
}
