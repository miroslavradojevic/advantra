package advantra.critpoint;

import java.awt.Color;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.util.Vector;

import advantra.feature.GaborFilt2D;
import advantra.feature.MyHessian;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.Line;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.measure.Calibration;
import ij.plugin.ZProjector;
import ij.plugin.filter.PlugInFilter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import imagescience.image.Axes;
import imagescience.image.Image;

public class VizFeatures implements PlugInFilter, MouseListener {

	ImagePlus 	img, gab, gabAll;
	ImagePlus	weighted;
	ImagePlus	vector_field;
	ImagePlus	neuriteness;
	ImagePlus 	Vx, Vy;
	
	// multi-scale
	double 		t1, t2; 
	int 		tn; 
	double[] 	s, t;
	
	// gabor filter 
	int			M; // number of angles per pi
	double[] 	theta_2pi;
	double[]  	theta_pi;
	double 		bandwidth = 1;
	double 		psi = 0;
	double 		gamma = 1.0;
	boolean 	isReal = true;
	
	// features
	double 		radius; // calculated wrt. hughest scale but can be fixed
	double 		dr, darc, rratio;
	int			surr=3;
	
	int 		W, H;
	
	int			nr_proc;
	
	public void run(ImageProcessor arg0) {
		
		t1  	= Prefs.get("advantra.critpoint.start_scale", 		3.0);
		t2    	= Prefs.get("advantra.critpoint.end_scale", 		7.0);
		tn 		= (int)Prefs.get("advantra.critpoint.nr_scales", 	3);
		M		= (int)Prefs.get("advantra.critpoint.nr_angles", 	8);
		dr    	= Prefs.get("advantra.critpoint.dr", 				1.0);
		darc    = Prefs.get("advantra.critpoint.darc", 				1.0);
		rratio	= Prefs.get("advantra.critpoint.rratio", 			0.2);
		nr_proc = (int)Prefs.get("advantra.critpoint.nr_proc", 		4);
		
		GenericDialog gd = new GenericDialog("VizFeatures");
		
		gd.addNumericField("start scale", 		t1, 1);
		gd.addNumericField("end   scale", 		t2, 1);
		gd.addNumericField("nr   scales", 		tn, 0, 5, "");
		gd.addNumericField("angles(per 180 deg)",M,	0, 5, "");
		
		gd.addMessage("circular extraction parameters");
		gd.addNumericField("x(lagest scale std)",	surr, 	0);
		gd.addNumericField("radius step", 			dr, 	1);
		gd.addNumericField("arc    step", 			darc, 	1);
		gd.addNumericField("rratio step", 			rratio, 1);
		
		gd.addMessage("parallelization");
		gd.addNumericField("CPU #",					nr_proc, 0);
		
		gd.showDialog();
		if (gd.wasCanceled()) return;
		
		t1 	= 		gd.getNextNumber();
		t2	= 		gd.getNextNumber();
		tn	= (int)	gd.getNextNumber();
		M 	= (int)	gd.getNextNumber();
		
		surr = (int)gd.getNextNumber();
		dr  = 		gd.getNextNumber();
		darc  = 		gd.getNextNumber();
		rratio  = 		gd.getNextNumber();
		
		nr_proc = (int)gd.getNextNumber();
		
		Prefs.set("advantra.critpoint.start_scale", t1);
		Prefs.set("advantra.critpoint.end_scale", 	t2);
		Prefs.set("advantra.critpoint.nr_scales", 	tn);
		Prefs.set("advantra.critpoint.nr_angles", 	M);
		Prefs.set("advantra.critpoint.dr", 			dr);
		Prefs.set("advantra.critpoint.darc", 		darc);
		Prefs.set("advantra.critpoint.rratio", 		rratio);
		Prefs.set("advantra.critpoint.nr_proc", 	nr_proc);
		
		// scales
		t = new double[tn];
		s = new double[tn];
		for (int i = 0; i < tn; i++) {
			t[i] = (i==0)? t1 : t1+i*((t2-t1)/(tn-1));
			s[i] = Math.sqrt(t[i]);
		}
		
		radius = surr*Math.sqrt(t[t.length-1]); // x GaussianStd.
		
		// angles theta angle it makes with x axis (angle it makes with the first row)
		theta_pi 	= new double[M];
		theta_2pi	= new double[2*M];
		for (int i = 0; i < 2*M; i++) {
			theta_2pi[i] = i * (Math.PI/(double)(M));
		}
		for (int i = 0; i < M; i++) {
			theta_pi[i] = i * (Math.PI/(double)M);
		}
		
		/*
		 *  convert to float , prepare image
		 */
		// reset the calibration
		Calibration c = img.getCalibration();
		c.pixelWidth 	= c.pixelHeight = c.pixelDepth = 1;
		c.setUnit("pixel");
		img.setCalibration(c);
		
		if(!img.getProcessor().isDefaultLut()) return;
		img.setProcessor("testImage", img.getProcessor().convertToFloat().duplicate());
		img.getWindow().getCanvas().addMouseListener(this);
		
		H = img.getHeight();
		W = img.getWidth();
		
		/*
		 * extract directional gabor responses
		 */
		
		System.out.print("extracting gabor, "+t.length+" scales, "+theta_pi.length+" angles... ");
		
		int N 				= theta_pi.length; // parallelize
		
		GaborFilt2D.load(
				img, 
				theta_pi,
				t,
				new double[t.length],
				bandwidth,
				psi,
				gamma,
				isReal);

		long t0, t1;
		t0 = System.currentTimeMillis();
		GaborFilt2D gab_jobs[] = new GaborFilt2D[nr_proc];
		for (int i = 0; i < gab_jobs.length; i++) {
			
			int start_interval 	= i*N/nr_proc;
			int end_interval	= (i+1)*N/nr_proc;
			
			gab_jobs[i] = new GaborFilt2D(start_interval,  end_interval);
			gab_jobs[i].start();
		}
		for (int i = 0; i < gab_jobs.length; i++) {
			try {
				gab_jobs[i].join();
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		gab = new ImagePlus("", GaborFilt2D.gabor_directional_responses);
		gab.setTitle("gabor filter, directional responses");
		t1 = System.currentTimeMillis();
		System.out.println("done, "+((t1-t0)/1000f)+" seconds");
		gab.show();
		
		/*
		 * extract maximal directional gabor response in every point
		 */

		ZProjector zmax = new ZProjector();
		zmax.setImage(gab);
		zmax.setStartSlice(1);	
		zmax.setStopSlice(gab.getStackSize());
		zmax.setMethod(ZProjector.MAX_METHOD);
		zmax.doProjection();
		gabAll = new ImagePlus("gabor filter, directional responses, max", zmax.getProjection().getChannelProcessor());
		gabAll.show();
		
		boolean normalize = false;
		if(normalize){
			weighted = new ImagePlus("gabor filter, directional responses, max, min-max normalized", 
					normalizeStackMinMax(gabAll.getStack()));
		}
		else{
			weighted = new ImagePlus("gabor filter, directional responses, max, min-max normalized", 
					gabAll.getStack());
		}
		
		weighted.show();
		
		/*
		 *  extract neuriteness & eigen vecs
		 */

		long t11 = System.currentTimeMillis();
		Vector<ImagePlus> nness = extractNeuritenessAndEigenVec(img, s);
		System.out.println("to get nness: "+((System.currentTimeMillis()-t11)/1000f)+" seconds");
		
		neuriteness = nness.get(0);
		neuriteness.setTitle("neuriteness");
		neuriteness.show();
		
		Vx = nness.get(1);
		Vx.setTitle("Vx");
		//Vx.show();
		
		Vy = nness.get(2);
		Vy.setTitle("Vy");
		//Vy.show();

		/*
		 * eigenvec on original
		 */
		vector_field = new ImagePlus("vector_field", img.getProcessor());
		Overlay eigenVecs = new Overlay();
		for (int row = 0; row < H; row++) {
			for (int col = 0; col < W; col++) {
				double vx 	= Vx.getProcessor().getf(col, row);
				double vy 	= Vy.getProcessor().getf(col, row);
				if(vx*vx+vy*vy>0)
					eigenVecs.add(new Line(col+0.5, row+0.5, col+0.5+vx, row+0.5+vy));
			}
		}
		vector_field.setOverlay(eigenVecs);
		vector_field.show();

		System.out.println("click now to see features at the  point!");
		
//		for (int i = 0; i < a.size(); i++) {
//			new ImagePlus("source"+i, plotValues(x_val, a.get(i))).show();
//		}
		
		System.out.println("done");
		
	}

	public int setup(String arg0, ImagePlus arg1) {
		img = arg1;
		return DOES_8G+NO_CHANGES;
	}

	public ImagePlus imageIntensityMetric(ImagePlus in, double in_max, double a, double lmbda){
		
		int w = in.getWidth();
		int h = in.getHeight();
		
		ImageStack is = new ImageStack(w, h);
		
		for (int l = 0; l < in.getStackSize(); l++) {
			
			ImageProcessor ip = new FloatProcessor(w, h);
			
			float[] temp = (float[])in.getStack().getProcessor(l+1).getPixels(); // temp is just the pointer
			
			for (int i = 0; i < temp.length; i++) {
				
				// exponential - stretching
				//ip.setf(i, 2 -(float)Math.exp(lmbda*(1-temp[i]/in_max)*(1-temp[i]/in_max)));
				
				// exponential - stretching
				//ip.setf(i, (float)Math.exp(lmbda*((temp[i]/in_max)-1)));
				
				// sigmoid
//				double a = 0.5f;
				double ImRatio = a*in_max;
				double t = ((temp[i]-ImRatio)/ImRatio);
				float val = (float)( 1 / (1+Math.exp(-lmbda*t)) );
				ip.setf(i, val);
			}
			
			is.addSlice(ip);
			
		}
		ImagePlus out = new ImagePlus("intensityMetric", is);
		
		return out;
		
	}
	
	public static float[] calculateStackMinMax(ImageStack instack){
		
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
	
	public static ImageStack normalizeStackMinMax(ImageStack instack){
		
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
		
		// create normalized stack
		ImageStack outstack = new ImageStack(instack.getWidth(), instack.getHeight());
		for (int i = 0; i < instack.getSize(); i++) {
			outstack.addSlice(new FloatProcessor(instack.getWidth(), instack.getHeight()));
		}
		
		if(Math.abs(mnmx[1]-mnmx[0])>0.0001){
			for (int z = 0; z < instack.getSize(); z++) {
						
				float[] take_val = (float[])instack.getProcessor(z+1).getPixels(); 
				for (int i = 0; i < take_val.length; i++) {
					take_val[i] = (take_val[i]-mnmx[0])/(mnmx[1]-mnmx[0]);
				}
				outstack.setPixels(take_val, z+1);
								
			}
		
		}	
			
		return outstack;
		
	}
	
	public static Vector<ImagePlus> extractNeuritenessAndEigenVec(ImagePlus in, double[] s){
		
		int H = in.getHeight();
		int W = in.getWidth();
		
		Image inimg = Image.wrap(in); inimg.axes(Axes.X+Axes.Y);
		
		Vector<ImagePlus> out = new Vector<ImagePlus>(3);
		
		ImagePlus neuriteness 	= new ImagePlus("neuriteness", new FloatProcessor(W, H));
		ImagePlus Vx 			= new ImagePlus("neuriteness", new FloatProcessor(W, H));
		ImagePlus Vy 			= new ImagePlus("neuriteness", new FloatProcessor(W, H));
		
		MyHessian myhess = new MyHessian();
		
		for (int i = 0; i < s.length; i++) {
			
			System.out.println("eigen analysis at scale = " + IJ.d2s(s[i]*s[i], 2) + ", with Gaussian sigma = " + IJ.d2s(s[i], 2));
			
			Vector<Image> h1 = myhess.eigs(inimg.duplicate(), s[i], false);
			Vector<Image> h2 = myhess.eigs(inimg.duplicate(), s[i], true);
			
			ImagePlus L2 	= h1.get(0).imageplus();
			ImagePlus L1 	= h1.get(1).imageplus();
			
			ImagePlus V1 	= h2.get(4).imageplus();
			ImagePlus V2 	= h2.get(5).imageplus();
			
			// smallest lambda over pixels
			float[] mn_mx_L = calculateStackMinMax(L1.getStack());
			float minL1 = mn_mx_L[0];
			
			for (int j = 0; j < W*H; j++) {
				
				float ll1 = L1.getProcessor().getf(j);
				float ll2 = L2.getProcessor().getf(j);
				float ll = (Math.abs(ll1)>Math.abs(ll2))?ll1:ll2; // larger in magnitude over 2 values
				if(ll<0){
					
					float ro = ll/minL1;
					if(ro>neuriteness.getProcessor().getf(j)){
						
						neuriteness.getProcessor().setf(j, ro);
						Vx.getProcessor().setf(j, V1.getProcessor().getf(j));
						Vy.getProcessor().setf(j, V2.getProcessor().getf(j));
						
					}
						
				}
				else{
					neuriteness.getProcessor().setf(j, 0);
				}
				
			}
		
		}
		
		out.add(neuriteness);
		out.add(Vx);
		out.add(Vy);
		
		return out;
		
	}
	
	public void mouseClicked(MouseEvent e) {
		
		int mouseX = 	img.getWindow().getCanvas().offScreenX(e.getX());
		int mouseY = 	img.getWindow().getCanvas().offScreenY(e.getY());
		
		Overlay ov1 = new Overlay();
		Overlay ov2 = new Overlay();
		Overlay ov3 = new Overlay();
		
		/*
		int startX 	= mouseX-15;
		startX = startX>=0? startX : 0;
		int stopX	= mouseX+15;
		stopX = stopX>=img.getWidth()? (img.getWidth()-1) : stopX;
		int startY 	= mouseY-15;
		startY = startY>=0? startY : 0;
		int stopY	= mouseY+15;
		stopY = stopY>=img.getHeight()? (img.getHeight()-1) : stopY;
		
		for (int i = startX; i <= stopX; i++) {
			for (int j = startY; j <= stopY; j++) {
				
				float mult1 = weighted.getProcessor().getPixelValue(i, j);
				float mult2 = neuriteness.getProcessor().getPixelValue(i, j);
				float mult3 = img.getProcessor().getPixelValue(i, j); 
				
				float v1 = Vx.getProcessor().getPixelValue(i, j);
				float v2 = Vy.getProcessor().getPixelValue(i, j);
				
				ov1.add(new Line(i+0.5, j+0.5, i+0.5+mult1*v1, j+0.5+mult1*v2));
				ov2.add(new Line(i+0.5, j+0.5, i+0.5+mult2*v1, j+0.5+mult2*v2));
				ov3.add(new Line(i+0.5, j+0.5, i+0.5+mult3*v1, j+0.5+mult3*v2));
				
			}
		}
		*/
		
		int nrang 		= (int)Math.ceil((Math.PI*2)/(darc/radius));
		
		System.out.print("extracting " + nrang + " directional responses... ");
		
		double[] 	orts			= new double[nrang];
		
		double[] 	profile1;//  		= new double[nrang];
		//int[]		profile1_cnt	= new int[nrang];
		double 		profile1_max 	= Double.MIN_VALUE;
		double 		profile1_min 	= Double.MAX_VALUE;
		double[] 	coeffs1;//  		= new double[nrang];
		
		double[] 	profile2;//  		= new double[nrang];
		//int[]		profile2_cnt	= new int[nrang];
		double 		profile2_max 	= Double.MIN_VALUE;
		double 		profile2_min 	= Double.MAX_VALUE;
		double[] 	coeffs2;//  		= new double[nrang];
		
		double[] 	profile3;//  		= new double[nrang];
		//int[]		profile3_cnt	= new int[nrang];
		double 		profile3_max 	= Double.MIN_VALUE;
		double 		profile3_min 	= Double.MAX_VALUE;
		double[] 	coeffs3;//  		= new double[nrang];
		
		profile1 = extractProfile(weighted, Vx, Vy, mouseX, mouseY, radius, dr, darc, rratio, nrang);
		coeffs1 = extractProfile(weighted, mouseX, mouseY, radius, dr, darc, rratio, nrang);
		
		profile2 = extractProfile(neuriteness, Vx, Vy, mouseX, mouseY, radius, dr, darc, rratio, nrang);
		coeffs2 = extractProfile(neuriteness, mouseX, mouseY, radius, dr, darc, rratio, nrang);
		
		profile3 = extractProfile(img, Vx, Vy, mouseX, mouseY, radius, dr, darc, rratio, nrang);
		coeffs3 = extractProfile(img, mouseX, mouseY, radius, dr, darc, rratio, nrang);
		
		for (int i = 0; i < nrang; i++) {
			orts[i] = i*(2*Math.PI/nrang);
		}
		
		ImageStack show_profiles1 = new ImageStack(800, 400);
		show_profiles1.addSlice(plotValues(orts, profile1, Plot.LINE));
		show_profiles1.addSlice(plotValues(orts, coeffs1, Plot.CROSS));
		show_profiles1.addSlice(plotValues(orts, profile2, Plot.LINE));
		show_profiles1.addSlice(plotValues(orts, coeffs2, Plot.CROSS));
		show_profiles1.addSlice(plotValues(orts, profile3, Plot.LINE));
		show_profiles1.addSlice(plotValues(orts, coeffs3, Plot.CROSS));
		new ImagePlus("profiles", show_profiles1).show();
		
		for (double r = radius; r >= radius*rratio; r-=dr) {
			for (double arc = 0; arc < 2*r*Math.PI; arc+=darc) {
				
				double ang = arc/r;
				
				double d1 = Math.sin(ang);
				double d2 = -Math.cos(ang);
				
				double x2 = mouseX + r * d1;
				double y2 = mouseY + r * d2;
				
				int x_loc = (int)Math.round(x2);
				int y_loc = (int)Math.round(y2);
				
				float mult1 = weighted.getProcessor().getPixelValue(x_loc, y_loc);
				float mult2 = neuriteness.getProcessor().getPixelValue(x_loc, y_loc);
				float mult3 = img.getProcessor().getPixelValue(x_loc, y_loc); 
				
				float v1 = Vx.getProcessor().getPixelValue(x_loc, y_loc);
				float v2 = Vy.getProcessor().getPixelValue(x_loc, y_loc);
				
				ov1.add(new Line(x2+0.5, y2+0.5, x2+0.5+mult1*v1, y2+0.5+mult1*v2));
				ov2.add(new Line(x2+0.5, y2+0.5, x2+0.5+mult2*v1, y2+0.5+mult2*v2));
				ov3.add(new Line(x2+0.5, y2+0.5, x2+0.5+mult3*v1, y2+0.5+mult3*v2));
				
			}
			
		}
		
		ImageStack show_profiles = new ImageStack(800, 400);  

		for (int k = 0; k < nrang; k++) {
			
			if(coeffs1[k]>profile1_max) profile1_max = coeffs1[k];
			if(coeffs1[k]<profile1_min) profile1_min = coeffs1[k];
			
			if(coeffs2[k]>profile2_max) profile2_max = coeffs2[k];
			if(coeffs2[k]<profile2_min) profile2_min = coeffs2[k];
			
			if(coeffs3[k]>profile3_max) profile3_max = coeffs3[k];
			if(coeffs3[k]<profile3_min) profile3_min = coeffs3[k];
			
		}
		
		Plot p1 = new Plot("", "angle", "gaborAng");
		p1.setLimits(0, 2*Math.PI, profile1_min-0.1*Math.abs(profile1_min), profile1_max+0.1*Math.abs(profile1_max));
		p1.setSize(800, 400);
		p1.setColor(Color.RED);
		p1.addPoints(orts, profile1, Plot.LINE);
		p1.setColor(Color.BLACK);
		p1.addPoints(orts, coeffs1, Plot.X);
		show_profiles.addSlice(p1.getProcessor());
		
		Plot p2 = new Plot("", "angle", "neurness");
		p2.setLimits(0, 2*Math.PI, profile2_min-0.1*Math.abs(profile2_min), profile2_max+0.1*Math.abs(profile2_max));
		p2.setSize(800, 400);
		p2.setColor(Color.RED);
		p2.addPoints(orts, profile2, Plot.LINE);
		p2.setColor(Color.BLACK);
		p2.addPoints(orts, coeffs2, Plot.X);
		show_profiles.addSlice(p2.getProcessor());
		
		Plot p3 = new Plot("", "angle", "pix.val");
		p3.setLimits(0, 2*Math.PI, 0-0.1*Math.abs(profile3_min), profile3_max+0.1*Math.abs(profile3_max));
		p3.setSize(800, 400);
		p3.setColor(Color.RED);
		p3.addPoints(orts, profile3, Plot.LINE);
		p3.setColor(Color.BLACK);
		p3.addPoints(orts, coeffs3, Plot.X);
		show_profiles.addSlice(p3.getProcessor());
		
		new ImagePlus("resp", show_profiles).show();
		
		weighted.setOverlay(ov1);
		neuriteness.setOverlay(ov2);
		img.setOverlay(ov3);
		
		Vector<ImagePlus> imgs = new Vector<ImagePlus>();
		imgs.add(weighted);
		imgs.add(neuriteness);
		imgs.add(img);
		
		// use the function to extract the profiles (just for the check here - should be used with automatic extraction)
		Vector<double[]> a = extractProfiles(imgs, Vx, Vy, mouseX, mouseY, radius, dr, darc, rratio);
		double[] x_val = new double[a.get(0).length];
		for (int i = 0; i < x_val.length; i++) x_val[i] = i;
		new ImagePlus("", plotValues(x_val, a)).show();
		
	}
	
	public ImageProcessor plotValues(double[] x_val, double[] y_val, int shape){
		
		// range
		double x_val_min = Double.MAX_VALUE;
		double x_val_max = Double.MIN_VALUE;
		double y_val_min = Double.MAX_VALUE;
		double y_val_max = Double.MIN_VALUE;
		
		for (int i = 0; i < y_val.length; i++) {
			if(x_val[i]<x_val_min) x_val_min = x_val[i];
			if(x_val[i]>x_val_max) x_val_max = x_val[i];
			if(y_val[i]<y_val_min) y_val_min = y_val[i];
			if(y_val[i]>y_val_max) y_val_max = y_val[i];
		}
		
		Plot p = new Plot("", "", "");
		p.setLimits(x_val_min, x_val_max, y_val_min-0.1*Math.abs(y_val_min), y_val_max+0.1*Math.abs(y_val_max)); // dangerous if limits are same
		p.setSize(800, 400);
		p.addPoints(x_val, y_val, shape);
		
		return p.getProcessor();
		
	}
	
	public ImageStack plotValues(double[] x_val, Vector<double[]> y_val){
		
		ImageStack show_profiles = new ImageStack(800, 400);
		
		for (int i = 0; i < y_val.size(); i++) {
			
			// range
			double x_val_min = Double.MAX_VALUE;
			double x_val_max = Double.MIN_VALUE;
			double y_val_min = Double.MAX_VALUE;
			double y_val_max = Double.MIN_VALUE;
			
			for (int j = 0; j < y_val.get(i).length; j++) {
				if(x_val[j]<x_val_min) x_val_min = x_val[j];
				if(x_val[j]>x_val_max) x_val_max = x_val[j];
				if(y_val.get(i)[j]<y_val_min) y_val_min = y_val.get(i)[j];
				if(y_val.get(i)[j]>y_val_max) y_val_max = y_val.get(i)[j];
			}
			
			Plot p = new Plot("", "", "");
			p.setLimits(x_val_min, x_val_max, y_val_min-0.1*Math.abs(y_val_min), y_val_max+0.1*Math.abs(y_val_max)); // dangerous if limits are same
			p.setSize(800, 400);
			p.addPoints(x_val, y_val.get(i), Plot.LINE);
			show_profiles.addSlice(p.getProcessor());
			
		}
		
		return show_profiles;
		
	}
	
	public static Vector<double[]> extractProfiles(Vector<ImagePlus> imgs, ImagePlus Vx, ImagePlus Vy, int atX, int atY, double radius, double dr, double darc, double rratio){
		
		int nrang 		= (int)Math.ceil((Math.PI*2)/(darc/radius));
		
		Vector<double[]> 	profile 		= new Vector<double[]>(imgs.size());
		Vector<int[]> 		profile_cnt 	= new Vector<int[]>(imgs.size());
		//Vector<Double> 		profile_min 	= new Vector<Double>(imgs.size());
		//Vector<Double> 		profile_max 	= new Vector<Double>(imgs.size());
		
		for (int i = 0; i < imgs.size(); i++) {
			
			profile.add(new double[nrang]);
			profile_cnt.add(new int[nrang]);
			//profile_min.add(Double.MAX_VALUE);
			//profile_max.add(Double.MIN_VALUE);
			
		}
		
		for (double r = radius; r >= radius*rratio; r-=dr) {
			for (double arc = 0; arc < 2*r*Math.PI; arc+=darc) {
				
				double ang = arc/r;
				
				double d1 = Math.sin(ang);
				double d2 = -Math.cos(ang);
				
				double x2 = atX + r * d1;
				double y2 = atY + r * d2;
				
				int x_loc = (int)Math.round(x2);
				int y_loc = (int)Math.round(y2);
				
				float v1 = Vx.getProcessor().getPixelValue(x_loc, y_loc);
				float v2 = Vy.getProcessor().getPixelValue(x_loc, y_loc);
				
				int idxAng = (int)Math.floor(ang/(2*Math.PI/nrang));
				
				for (int i = 0; i < imgs.size(); i++) {
					
					float mult = imgs.get(i).getProcessor().getPixelValue(x_loc, y_loc);
					profile.get(i)[idxAng] 		+= Math.abs(mult*v1*d1+mult*v2*d2);
					profile_cnt.get(i)[idxAng]	+= 1;
					
				}
				
			}
		}
		
		for (int i = 0; i < profile.size(); i++) {
			for (int k = 0; k < nrang; k++) {
				if(profile_cnt.get(i)[k]>1){
					profile.get(i)[k] = profile.get(i)[k] / profile_cnt.get(i)[k];
				}
			}
		}
		
		return profile;
		
	}
	
	public static double[] extractProfile(ImagePlus img, ImagePlus Vx, ImagePlus Vy, int atX, int atY, double radius, double dr, double darc, double rratio, int nrang){
		
		int angular_resolution = nrang;
		
		double[] out_profile = new double[angular_resolution];
		int[] out_profile_cnt = new int[angular_resolution];
		
		for (double r = radius; r >= radius*rratio; r-=dr) {
			for (double arc = 0; arc < 2*r*Math.PI; arc+=darc) {
				
				double ang = arc/r;
				
				double d1 = Math.sin(ang);
				double d2 = -Math.cos(ang);
				
				double x2 = atX + r * d1;
				double y2 = atY + r * d2;
				
				int x_loc = (int)Math.round(x2);
				int y_loc = (int)Math.round(y2);
				
				float v1 = Vx.getProcessor().getPixelValue(x_loc, y_loc);
				float v2 = Vy.getProcessor().getPixelValue(x_loc, y_loc);
				
				int idxAng = (int)Math.floor(ang/(2*Math.PI/angular_resolution));
				
				float mult = img.getProcessor().getPixelValue(x_loc, y_loc);
				out_profile[idxAng] 		+= Math.abs(mult*v1*d1+mult*v2*d2);
				out_profile_cnt[idxAng]		+= 1;
				
			}
		}
		
		for (int k = 0; k < angular_resolution; k++) {
			if(out_profile_cnt[k]>1){
				out_profile[k] = out_profile[k] / out_profile_cnt[k];
			}
		}
		
		return out_profile;
	}
	
	public static double[] extractProfile(ImagePlus img, int atX, int atY, double radius, double dr, double darc, double rratio, int nrang){
		
		int angular_resolution = nrang;
		
		double[] out_profile = new double[angular_resolution];
		int[] out_profile_cnt = new int[angular_resolution];
		
		for (double r = radius; r >= radius*rratio; r-=dr) {
			for (double arc = 0; arc < 2*r*Math.PI; arc+=darc) {
				
				double ang = arc/r;
				
				double d1 = Math.sin(ang);
				double d2 = -Math.cos(ang);
				
				double x2 = atX + r * d1;
				double y2 = atY + r * d2;
				
				int x_loc = (int)Math.round(x2);
				int y_loc = (int)Math.round(y2);
				
//				float v1 = Vx.getProcessor().getPixelValue(x_loc, y_loc);
//				float v2 = Vy.getProcessor().getPixelValue(x_loc, y_loc);
				
				int idxAng = (int)Math.floor(ang/(2*Math.PI/angular_resolution));
				
				float mult = img.getProcessor().getPixelValue(x_loc, y_loc);
				out_profile[idxAng] 		+= mult;//Math.abs(mult*v1*d1+mult*v2*d2);
				out_profile_cnt[idxAng]		+= 1;
				
			}
		}
		
		for (int k = 0; k < angular_resolution; k++) {
			if(out_profile_cnt[k]>1){
				out_profile[k] = out_profile[k] / out_profile_cnt[k];
			}
		}
		
		return out_profile;
	}
	
	public void mousePressed(MouseEvent e) {}

	public void mouseReleased(MouseEvent e) {}

	public void mouseEntered(MouseEvent e) {}

	public void mouseExited(MouseEvent e) {}
	
}

//double[] profile1 = new double[theta_2pi.length];
//double[] profile2 = new double[theta_2pi.length];
//double[] profile3 = new double[theta_2pi.length];
//// use circular_profile to obtain line profile4
//double[] profile4 = new double[circular_profile3.getWidth()];
//double[] profile5 = new double[circular_profile3.getWidth()];
//int capture = (2*M)/4;
//for (int i = 0; i < circular_profile3.getWidth(); i++) {
////	IJ.log("--->"+i+" : (capture "+capture);
//	double 	sum = 0;
//	int 	cnt = 0;
//	double 	prod = 1;
//	for (int j = 0; j < circular_profile3.getHeight(); j++) {
//		sum 	+= circular_profile3.getPixelValue(i, j);
//		prod 	*= circular_profile3.getPixelValue(i, j);
//		cnt++;
//		
//		/*
//		int range = (int)Math.floor(-capture*  (j/(double)circular_profile3.getHeight())  +  capture);
////		IJ.log("radial distance "+j+" : +/-  "+range);
//		for (int k = i-range; k <= i+range; k++) {
//			int s = (k<0)?(k+(2*M)):(k>=(2*M))?(k-(2*M)):k;
////			IJ.log(s+",");
//			sum += circular_profile3.getPixelValue(s, j);
//			cnt++;
//		}
//		*/
//	}
//	IJ.log(i+" :  sum "+sum+", prod = "+prod+" (cnt="+cnt+")");
//	profile4[i] = sum;
//	profile5[i] = prod;
//}
//PointRoi center = ;
//center.setColor(Color.GREEN);
//ov.add(new PointRoi(mouseX+0.5, mouseY+0.5));
//	int k1 = (k>=theta_pi.length)? k-theta_pi.length : k ;
//	profile1[k] = gab.getStack().getProcessor(k1+1).getInterpolatedValue(x2, y2);
//	profile2[k] = gabAll.getProcessor().getInterpolatedValue(x2, y2);
//	profile3[k] = img.getProcessor().getInterpolatedValue(x2, y2);