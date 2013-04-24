package advantra.critpoint;

import java.awt.Color;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.util.Vector;

import advantra.feature.GaborFilt2D;
import advantra.feature.MyHessian;
import advantra.tools.MeanShift3DSphere;
import advantra.trace.Tracing;

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
	Image		inimg; // imagescience version
	
	// gabor filter extractions params
	double 		t1, t2; 
	int 		tn; 
	double[] 	s;
	double[] 	t;
	int			M; // number of angles per pi
	double[] 	theta_2pi;
	double[]  	theta_pi;
	double 		bandwidth = 1;
	double 		psi = 0;
	double 		gamma = 1.0;
	boolean 	isReal = true;
	
	// circ features
	double 		radius; // calculated wrt. hughest scale but can be fixed
	double 		dr, darc, dratio;
	
	int 		W, H;
	
	public void run(ImageProcessor arg0) {
		
		t1  	= Prefs.get("advantra.critpoint.start_scale", 		3.0);
		t2    	= Prefs.get("advantra.critpoint.end_scale", 		11.0);
		tn 		= (int)Prefs.get("advantra.critpoint.nr_scales", 	5);
		M		= (int)Prefs.getInt("advantra.critpoint.nr_angles", 8);
		dr    	= Prefs.get("advantra.critpoint.dr", 				1.0);
		darc    = Prefs.get("advantra.critpoint.darc", 				1.0);
		dratio	= Prefs.get("advantra.critpoint.dratio", 			0.2);
		
		GenericDialog gd = new GenericDialog("VizFeatures");
		
		gd.addNumericField("start scale", 		t1, 2);
		gd.addNumericField("end   scale", 		t2, 4);
		gd.addNumericField("nr   scales", 		tn, 0, 5, "");
		gd.addNumericField("angles(per 180 deg)",M,	0, 5, "");
		
		gd.addMessage("circular extraction parameters");// NumericField("", gamma, 2);
		
		gd.addNumericField("radius step", 		dr, 1);
		gd.addNumericField("arc    step", 		darc, 1);
		gd.addNumericField("rratio step", 		dratio, 1);
		
		gd.showDialog();
		if (gd.wasCanceled()) return;
		
		t1 	= 		gd.getNextNumber();
		t2	= 		gd.getNextNumber();
		tn	= (int)	gd.getNextNumber();
		M 	= (int)	gd.getNextNumber();
		dr  = 		gd.getNextNumber();
		darc  = 		gd.getNextNumber();
		dratio  = 		gd.getNextNumber();
		
		Prefs.set("advantra.critpoint.start_scale", t1);
		Prefs.set("advantra.critpoint.end_scale", 	t2);
		Prefs.set("advantra.critpoint.nr_scales", 	tn);
		Prefs.set("advantra.critpoint.nr_angles", 	M);
		Prefs.set("advantra.critpoint.dr", 			dr);
		Prefs.set("advantra.critpoint.darc", 		darc);
		Prefs.set("advantra.critpoint.dratio", 		dratio);
		
		// scales
		t = new double[tn];
		s = new double[tn];
		for (int i = 0; i < tn; i++) {
			t[i] = (i==0)? t1 : t1+i*((t2-t1)/(tn-1));
			s[i] = Math.sqrt(t[i]);
		}
		
		radius = 3*Math.sqrt(t[t.length-1]); // xGaussianStd.
		System.out.println("surrounding radius is " + radius);
		
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
		
		float[] mn_mx = calculateStackMinMax(img.getStack());
		
		/*
		 * show gabor kernels for selected scales
		 */
		if(false){
		for (int i = 0; i < t.length; i++) {
			
			ImagePlus krn = GaborFilt2D.run(null, theta_pi, t[i], 0, bandwidth, psi, gamma, isReal);
			krn.setTitle("gabor_kernels_scale"+IJ.d2s(t[i],2));
			krn.show();
			for (int j = 0; j < 2; j++) {
				krn.getCanvas().zoomIn(0, 0);
			}
			
		}
		}
		
		
		/*
		 * extract directional gabor responses
		 */
		

		
		int nr_proc = 4;
		int nr_splits = 4;
		
		GaborFilt2D.load(
				img, 
				nr_splits, 
				theta_pi,
				t,
				new double[t.length],
				bandwidth,
				psi,
				gamma,
				isReal);
		
		long t11 = System.currentTimeMillis();

		int how_many_patches = GaborFilt2D.ptches.size();
		
		GaborFilt2D gab_jobs[] = new GaborFilt2D[nr_proc];
		for (int i = 0; i < gab_jobs.length; i++) {
			
			int start_interval 	= i*how_many_patches/nr_proc;
			int end_interval	= (i+1)*how_many_patches/nr_proc;
			
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
		
//		for (int i = 0; i < GaborFilt2D.ptches_out.size(); i++) {
//			System.out.println(
//					i+" "+
//					GaborFilt2D.ptches_out.get(i).getHeight()+" x "+
//					GaborFilt2D.ptches_out.get(i).getWidth()+" x "+
//					GaborFilt2D.ptches_out.get(i).getStackSize()+" at "+
//					GaborFilt2D.ptch_root_x.get(i)+" , "+
//					GaborFilt2D.ptch_root_y.get(i)
//					);
//			//GaborFilt2D.ptches_out.get(i).show();
//		}
		
		
		ImagePlus atlast = GaborFilt2D.concatenateOutput();
		
		long t12 = System.currentTimeMillis();
		
		System.out.println("parallel threading took: "+((t12-t11)/1000f)+" s.");
		
		atlast.show();
		
		long t21 = System.currentTimeMillis();
		gab = GaborFilt2D.extractDirectionalGabor(img, theta_pi, t, bandwidth, psi, gamma, isReal); 
		long t22 = System.currentTimeMillis();
		
		System.out.println("without 	threading took: "+((t22-t21)/1000f)+" s.");
		
		gab.show();
		
		if(true) return;
		
		mn_mx = calculateStackMinMax(gab.getStack());
		System.out.println("gabor outputs: "+mn_mx[0]+" / "+mn_mx[1]);
		
		/*
		 * extract maximal directional gabor response in every point
		 */
if(true){
		IJ.log("maximal directional response in every point...");
		ZProjector zmax = new ZProjector();
		zmax.setImage(gab);
		zmax.setStartSlice(1);
		zmax.setStopSlice(gab.getStackSize());
		zmax.setMethod(ZProjector.MAX_METHOD);
		zmax.doProjection();
		gabAll = new ImagePlus("allThetas", zmax.getProjection().getChannelProcessor());
		gabAll.show();
		
		mn_mx = calculateStackMinMax(gabAll.getStack());
		IJ.log("redefining gabor extrema : "+mn_mx[0]+" / "+mn_mx[1]);
		
		weighted = imageIntensityMetric(gabAll, mn_mx[1], 0.4, 2);
		weighted.show();
		
}
		/*
		 *  extract neuriteness & eigen vecs
		 */

		t11 = System.currentTimeMillis();
		Vector<ImagePlus> nness = extractNeuritenessAndEigenVec(img, s);
		System.out.println("to get nness: "+((System.currentTimeMillis()-t11)/1000f)+"  size is "+nness.size());
		
		nness.get(0).setTitle("nness(try)");
		nness.get(0).show();
		nness.get(1).setTitle("Vx(try)");
		nness.get(1).show();
		nness.get(2).setTitle("Vy(try)");
		nness.get(2).show();
		
		inimg = Image.wrap(img); inimg.axes(Axes.X+Axes.Y);

		neuriteness = new ImagePlus("neuriteness", new FloatProcessor(W, H));
		Vx 			= new ImagePlus("neuriteness", new FloatProcessor(W, H));
		Vy 			= new ImagePlus("neuriteness", new FloatProcessor(W, H));
		
		MyHessian myhess = new MyHessian();
		
		for (int i = 0; i < s.length; i++) {
			
			System.out.println("processing scale " + t[i] + ", with sigma" + s[i]);
			
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
		
		neuriteness.setTitle("nness(orig.)");
		neuriteness.show(); 
		Vx.setTitle("Vx(orig.)");
		Vx.show();
		Vy.setTitle("Vy(orig.)");
		Vy.show();
		
		if(true) return;
		// done  with neuriteness
		/*
		 * eigenvec on original
		 */
		vector_field = new ImagePlus("vector_field", img.getProcessor());
		Overlay eigenVecs = new Overlay();
		for (int row = 0; row < H; row++) {
			for (int col = 0; col < W; col++) {
				double vx 	= Vx.getProcessor().getf(col, row);
				double vy 	= Vy.getProcessor().getf(col, row);
				//float mult 	= neuriteness.getProcessor().getf(col, row);
				if(vx*vx+vy*vy>0)
					eigenVecs.add(new Line(col+0.5, row+0.5, col+0.5+vx, row+0.5+vy));
			}
		}
		vector_field.setOverlay(eigenVecs);
		vector_field.show();

		System.out.println("click now!");
		
		
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
	
	public float[] calculateStackMinMax(ImageStack instack){
		
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
	
	public Vector<ImagePlus> extractNeuritenessAndEigenVec(ImagePlus in, double[] s){
		
		int H = in.getHeight();
		int W = in.getWidth();
		
		Image inimg = Image.wrap(in); inimg.axes(Axes.X+Axes.Y);
		
		Vector<ImagePlus> out = new Vector<ImagePlus>(3);
		
		ImagePlus neuriteness 	= new ImagePlus("neuriteness", new FloatProcessor(W, H));
		ImagePlus Vx 			= new ImagePlus("neuriteness", new FloatProcessor(W, H));
		ImagePlus Vy 			= new ImagePlus("neuriteness", new FloatProcessor(W, H));
		
		////
		
		MyHessian myhess = new MyHessian();
		
		for (int i = 0; i < s.length; i++) {
			
			System.out.println("processing scale " + (s[i]*s[i]) + ", with sigma" + s[i]);
			
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
		
		////
		
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
		
		weighted.setOverlay(ov1);
		neuriteness.setOverlay(ov2);
		img.setOverlay(ov3);
		
		
		
		double dr 		= 1;
		double darc 	= 1;
		double rratio 	= 0.2;
		
		int nrang 		= (int)Math.ceil((Math.PI*2)/(darc/radius));
		
//		// count the number of points
//		int cnt = 0;
//		for (double r = radius; r >= radius*rratio; r-=dr) {
//			for (double arc = 0; arc < 2*r*Math.PI; arc+=darc) {
//				cnt++;
//			}
//		}
		
		System.out.print("extracting " + nrang + " directional responses... ");
		
		double[] 	orts			= new double[nrang];
		
		double[] 	profile1  		= new double[nrang];
		int[]		profile1_cnt	= new int[nrang];
		double 		profile1_max 	= Double.MIN_VALUE;
		double 		profile1_min 	= Double.MAX_VALUE;
		double[] 	coeffs1  		= new double[nrang];
		
		double[] 	profile2  		= new double[nrang];
		int[]		profile2_cnt	= new int[nrang];
		double 		profile2_max 	= Double.MIN_VALUE;
		double 		profile2_min 	= Double.MAX_VALUE;
		double[] 	coeffs2  		= new double[nrang];
		
		double[] 	profile3  		= new double[nrang];
		int[]		profile3_cnt	= new int[nrang];
		double 		profile3_max 	= Double.MIN_VALUE;
		double 		profile3_min 	= Double.MAX_VALUE;
		double[] 	coeffs3  		= new double[nrang];
		
		//int cnt = 0;
		for (double r = radius; r >= radius*rratio; r-=dr) {
			for (double arc = 0; arc < 2*r*Math.PI; arc+=darc) {
				
				double ang = arc/r;
				
				double d1 = Math.sin(ang);
				double d2 = -Math.cos(ang);
				
				double x2 = mouseX + r * d1;
				double y2 = mouseY + r * d2;
				
				//ov3.add(new Line(x2+0.5, y2+0.5, x2+0.5+1*d1, y2+0.5+1*d2));
				//ov.add(new PointRoi(x2+0.5, y2+0.5));
				
				int x_loc = (int)Math.round(x2);
				int y_loc = (int)Math.round(y2);
				
				float mult1 = weighted.getProcessor().getPixelValue(x_loc, y_loc);
				float mult2 = neuriteness.getProcessor().getPixelValue(x_loc, y_loc);
				float mult3 = img.getProcessor().getPixelValue(x_loc, y_loc); 
				
				float v1 = Vx.getProcessor().getPixelValue(x_loc, y_loc);
				float v2 = Vy.getProcessor().getPixelValue(x_loc, y_loc);
				
				int idxAng = (int)Math.floor(ang/(2*Math.PI/nrang));

				// avg for the profile elements
				profile1[idxAng] 		+= Math.abs(mult1*v1*d1+mult1*v2*d2);
				profile1_cnt[idxAng]	+= 1;
				profile2[idxAng] 		+= Math.abs(mult2*v1*d1+mult2*v2*d2);
				profile2_cnt[idxAng]	+= 1;
				profile3[idxAng] 		+= Math.abs(mult3*v1*d1+mult3*v2*d2);
				profile3_cnt[idxAng]	+= 1;
				
				//orts[idxAng]			= ang*360/(2*Math.PI);
				
				coeffs1[idxAng]			+= mult1;
				coeffs2[idxAng]			+= mult2;
				coeffs3[idxAng]			+= mult3;
				
				//cnt++;
				
			}
			
		}
		
		// avg calc
		for (int k = 0; k < nrang; k++) {
			
			if(profile1_cnt[k]>1){
				profile1[k] = profile1[k] / profile1_cnt[k]; 
				coeffs1[k]  = coeffs1[k] / profile1_cnt[k];
			}
			
			if(profile2_cnt[k]>1){
				profile2[k] = profile2[k] / profile2_cnt[k]; 
				coeffs2[k]  = coeffs2[k] / profile2_cnt[k];
			}
			
			if(profile3_cnt[k]>1){
				profile3[k] = profile3[k] / profile3_cnt[k]; 
				coeffs3[k]  = coeffs3[k] / profile3_cnt[k];
			}
			
		}
		
		for (int k = 0; k < nrang; k++) {
			
			if(coeffs1[k]>profile1_max) profile1_max = coeffs1[k];
			if(coeffs1[k]<profile1_min) profile1_min = coeffs1[k];
			
			if(coeffs2[k]>profile2_max) profile2_max = coeffs2[k];
			if(coeffs2[k]<profile2_min) profile2_min = coeffs2[k];
			
			if(coeffs3[k]>profile3_max) profile3_max = coeffs3[k];
			if(coeffs3[k]<profile3_min) profile3_min = coeffs3[k];
			
		}
		
		System.out.println("done");
		
		for (int i = 0; i < nrang; i++) {
			orts[i] = i*(2*Math.PI/nrang);
		}
		
		ImageStack show_profiles = new ImageStack(800, 400);  

		Plot p1 = new Plot("", "angle", "gaborAng");
		p1.setLimits(0, 2*Math.PI, profile1_min-1, profile1_max+1);
		p1.setSize(800, 400);
		p1.setColor(Color.RED);
		p1.addPoints(orts, profile1, Plot.BOX);
		p1.setColor(Color.BLACK);
		p1.addPoints(orts, coeffs1, Plot.X);
		show_profiles.addSlice(p1.getProcessor());
		
		Plot p2 = new Plot("", "angle", "neurness");
		p2.setLimits(0, 2*Math.PI, profile2_min-1, profile2_max+1);
		p2.setSize(800, 400);
		p2.setColor(Color.RED);
		p2.addPoints(orts, profile2, Plot.BOX);
		p2.setColor(Color.BLACK);
		p2.addPoints(orts, coeffs2, Plot.X);
		show_profiles.addSlice(p2.getProcessor());
		
		Plot p3 = new Plot("", "angle", "pix.val");
		p3.setLimits(0, 2*Math.PI, 0-5, profile3_max+5);
		p3.setSize(800, 400);
		p3.setColor(Color.RED);
		p3.addPoints(orts, profile3, Plot.BOX);
		p3.setColor(Color.BLACK);
		p3.addPoints(orts, coeffs3, Plot.X);
		show_profiles.addSlice(p3.getProcessor());
		
		new ImagePlus("resp", show_profiles).show();
		
	}	
		/*
		 * radial profile (image of the profile and inside)
		 */
		
//		if(false){
//		
//	    int Q = 5;
//	    double[] r = new double[Q];
//	    for (int l = 0; l < Q; l++) {
//	    	r[l] = (l+1)*(radius/(double)Q);
//	    }
//		
//	    Overlay ov1 = new Overlay();
//	    
//	    ImageStack  	circular_profile_stack 	= new ImageStack(2*M, Q);
//	    ImageProcessor 	circular_profile1 		= new FloatProcessor(2*M, Q);
//	    ImageProcessor 	circular_profile2 		= new FloatProcessor(2*M, Q);
//	    ImageProcessor 	circular_profile3 		= new FloatProcessor(2*M, Q);
//	    
//	    for (int l = 0; l < Q; l++) {
//	    	for (int k = 0; k < 2*M; k++) {
//
//				double x2 = mouseX + r[l] * Math.sin(theta_2pi[k]);
//				double y2 = mouseY - r[l] * Math.cos(theta_2pi[k]);
//				ov1.add(new PointRoi(x2, y2));
//				
//				circular_profile1.setf(k, l, (float)gabAll.getProcessor().getInterpolatedValue(x2, y2));
//				circular_profile2.setf(k, l, (float)img.getProcessor().getInterpolatedValue(x2, y2));
//				circular_profile3.setf(k, l, (float)weighted.getProcessor().getInterpolatedValue(x2, y2));
//				
//			}
//		}
//	    
//	    circular_profile_stack.addSlice(circular_profile1);
//	    circular_profile_stack.addSlice(circular_profile2);
//	    circular_profile_stack.addSlice(circular_profile3);
//	    
//	    ImagePlus circular_profile_image = new ImagePlus("circular_profile", circular_profile_stack);
//	    circular_profile_image.show();
//	    
//	    for (int i = 0; i < 8; i++) {circular_profile_image.getCanvas().zoomIn(0, 0);}
//	    
//	    img.setOverlay(ov1);
//	    
//		}
	    
	    /*
	     * line profiles
	     */
	    
//		double[] profile1 = new double[theta_2pi.length];
//		double[] profile2 = new double[theta_2pi.length];
//		double[] profile3 = new double[theta_2pi.length];
//		// use circular_profile to obtain line profile4
//		double[] profile4 = new double[circular_profile3.getWidth()];
//		double[] profile5 = new double[circular_profile3.getWidth()];
		
//		int capture = (2*M)/4;
		
//		for (int i = 0; i < circular_profile3.getWidth(); i++) {
////			IJ.log("--->"+i+" : (capture "+capture);
//			double 	sum = 0;
//			int 	cnt = 0;
//			double 	prod = 1;
//			
//			for (int j = 0; j < circular_profile3.getHeight(); j++) {
//				
//				sum 	+= circular_profile3.getPixelValue(i, j);
//				prod 	*= circular_profile3.getPixelValue(i, j);
//				cnt++;
//				
//				/*
//				int range = (int)Math.floor(-capture*  (j/(double)circular_profile3.getHeight())  +  capture);
////				IJ.log("radial distance "+j+" : +/-  "+range);
//				for (int k = i-range; k <= i+range; k++) {
//					int s = (k<0)?(k+(2*M)):(k>=(2*M))?(k-(2*M)):k;
////					IJ.log(s+",");
//					sum += circular_profile3.getPixelValue(s, j);
//					cnt++;
//				}
//				*/
//				
//			}
//			
//			IJ.log(i+" :  sum "+sum+", prod = "+prod+" (cnt="+cnt+")");
//			
//			profile4[i] = sum;
//			profile5[i] = prod;
//			
//		}
		
//		// normalize
//		double sm = 0;
//		for (int i = 0; i < profile5.length; i++) {
//			sm += profile5[i];
//		}
//		for (int i = 0; i < profile5.length; i++) {
//			profile5[i] = profile5[i]/sm;
//		}
		
		//PointRoi center = ;
		//center.setColor(Color.GREEN);
		//ov.add(new PointRoi(mouseX+0.5, mouseY+0.5));
		

//			int k1 = (k>=theta_pi.length)? k-theta_pi.length : k ;
//			profile1[k] = gab.getStack().getProcessor(k1+1).getInterpolatedValue(x2, y2);
//			profile2[k] = gabAll.getProcessor().getInterpolatedValue(x2, y2);
//			profile3[k] = img.getProcessor().getInterpolatedValue(x2, y2);

	public void mousePressed(MouseEvent e) {}

	public void mouseReleased(MouseEvent e) {}

	public void mouseEntered(MouseEvent e) {}

	public void mouseExited(MouseEvent e) {}
	
}
