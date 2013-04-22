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

public class AlongCircleFt implements PlugInFilter, MouseListener {

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
	double 		radius; // fixed
	
	int W, H;
	
	double 		gab_min, gab_max;
	double		img_min, img_max;
	
	public void run(ImageProcessor arg0) {
		
		t1  	= Prefs.get("advantra.critpoint.start_scale", 		2.0);
		t2    	= Prefs.get("advantra.critpoint.end_scale", 		4.0);
		tn 		= (int)Prefs.get("advantra.critpoint.nr_scales", 	3);
		M		= (int)Prefs.getInt("advantra.critpoint.nr_angles", 8);
		
		GenericDialog gd = new GenericDialog("F1");
		
		gd.addNumericField("start scale", 		t1, 2);
		gd.addNumericField("end   scale", 		t2, 4);
		gd.addNumericField("nr   scales", 		tn, 0, 5, "");
		gd.addNumericField("angles(per 180 deg)",M,	0, 5, "");
		gd.addNumericField("", gamma, 2);
		gd.showDialog();
		if (gd.wasCanceled()) return;
		
		t1 	= 		gd.getNextNumber();
		t2	= 		gd.getNextNumber();
		tn	= (int)	gd.getNextNumber();
		M 	= (int)	gd.getNextNumber();
		gamma = 	gd.getNextNumber();
		
		Prefs.set("advantra.critpoint.start_scale", t1);
		Prefs.set("advantra.critpoint.end_scale", 	t2);
		Prefs.set("advantra.critpoint.nr_scales", 	tn);
		Prefs.set("advantra.critpoint.nr_angles", 	M);
		
		// scales
		t = new double[tn];
		s = new double[tn];
		for (int i = 0; i < tn; i++) {
			t[i] = (i==0)? t1 : t1+i*((t2-t1)/(tn-1));
			s[i] = Math.sqrt(t[i]);
		}
		
		radius = 3*Math.sqrt(t[t.length-1]); // xGaussianStd.
		IJ.log("surrounding radius is " + radius);
		
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
		IJ.log("converting, preparing image...");
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
		
		IJ.log("listening mouse on the original image...");
		inimg = Image.wrap(img); inimg.axes(Axes.X+Axes.Y);
		
		float[] mn_mx = calculateStackMinMax(img.getStack());
		img_min = mn_mx[0];
		img_max = mn_mx[1];
		IJ.log("image range: "+img_min+" / "+img_max);
		
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
if(true){
		IJ.log("extract directional gabor responses...");
		gab = extractDirectionalGabor(img, theta_pi, t, bandwidth, psi, gamma, isReal); 
		
		mn_mx = calculateStackMinMax(gab.getStack());
		gab_min = mn_mx[0];
		gab_max = mn_mx[1];
		IJ.log("gabor outputs: "+gab_min+" / "+gab_max);
}
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
		//gabAll.show();
		
		mn_mx = calculateStackMinMax(gabAll.getStack());
		gab_min = mn_mx[0];
		gab_max = mn_mx[1];
		IJ.log("redefining gabor extrema : "+gab_min+" / "+gab_max);
		
		weighted = imageIntensityMetric(gabAll, gab_max, 0.4, 2);
		weighted.show();
		
}
		/*
		 *  extract neuriteness
		 */

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
		
		neuriteness.show();
		
		/*
		 * eigenvec on neuriteness
		 */
		Overlay ovlyV = new Overlay();
		for (int row = 0; row < neuriteness.getHeight(); row++) {
			for (int col = 0; col < neuriteness.getWidth(); col++) {
				double vx 	= Vx.getProcessor().getf(col, row);
				double vy 	= Vy.getProcessor().getf(col, row);
				float mult 	= neuriteness.getProcessor().getf(col, row);
				if(mult>0)
				ovlyV.add(new Line(col+0.5, row+0.5, col+0.5+mult*vx, row+0.5+mult*vy));
			}
		}
		//neuriteness.setOverlay(ovlyV);

		/*
		 * eigenvec on gabAll
		 */
		Overlay ovlyV1 = new Overlay();
		for (int row = 0; row < weighted.getHeight(); row++) {
			for (int col = 0; col < weighted.getWidth(); col++) {
				double vx 	= Vx.getProcessor().getf(col, row);
				double vy 	= Vy.getProcessor().getf(col, row);
				float mult 	= weighted.getProcessor().getf(col, row);
				if(mult>0)
				ovlyV1.add(new Line(col+0.5, row+0.5, col+0.5+mult*vx, row+0.5+mult*vy));
			}
		}
		//weighted.setOverlay(ovlyV1);
		
		/*
		 * eigenvec on img
		 */
		Overlay ovlyV2 = new Overlay();
		for (int row = 0; row < img.getHeight(); row++) {
			for (int col = 0; col < img.getWidth(); col++) {
				double vx 	= Vx.getProcessor().getf(col, row);
				double vy 	= Vy.getProcessor().getf(col, row);
				float mult 	= img.getProcessor().getf(col, row);
				if(mult>0)
				ovlyV2.add(new Line(col+0.5, row+0.5, col+0.5+mult*vx, row+0.5+mult*vy));
			}
		}
//		img.setOverlay(ovlyV2);
		System.out.println("click!");
		
		
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
	
	public ImagePlus extractDirectionalGabor(
			ImagePlus input, 
			double[] angles_pi, 
			double[] t, 
			double bw, 
			double psi, 
			double gamma,
			boolean isReal){
		
		int w = input.getWidth();
		int h = input.getHeight();
		
		ImageStack 	angul = new ImageStack(w, h); 
		
		ZProjector zmax = new ZProjector();
		
		for (int i = 0; i < angles_pi.length; i++) {
			
			double current_theta = angles_pi[i];
			
			System.out.println("processing theta = "+i+" / "+ (angles_pi.length-1));
			
			ImagePlus g_theta = GaborFilt2D.run(
					input, current_theta, t, new double[t.length], bw, psi, gamma, isReal);
			
			zmax.setImage(g_theta);
			zmax.setStartSlice(1);
			zmax.setStopSlice(g_theta.getStackSize());
			zmax.setMethod(ZProjector.MAX_METHOD);
			zmax.doProjection();
			angul.addSlice("theta="+IJ.d2s(current_theta, 2), zmax.getProjection().getChannelProcessor());

		}
		
		return new ImagePlus("gabor_per_angle", angul);

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
	
	public void mouseClicked(MouseEvent e) {
		
		int mouseX = 	img.getWindow().getCanvas().offScreenX(e.getX());
		int mouseY = 	img.getWindow().getCanvas().offScreenY(e.getY());
		
		Overlay ov = new Overlay();
		
		int startX 	= mouseX-20;
		startX = startX>=0? startX : 0;
		int stopX	= mouseX+20;
		stopX = stopX>=img.getWidth()? (img.getWidth()-1) : stopX;
		int startY 	= mouseY-20;
		startY = startY>=0? startY : 0;
		int stopY	= mouseY+20;
		stopY = stopY>=img.getHeight()? (img.getHeight()-1) : stopY;
		
		for (int i = startX; i <= stopX; i++) {
			for (int j = startY; j <= stopY; j++) {
				float mult = neuriteness.getProcessor().getPixelValue(i, j);
				float v1 = Vx.getProcessor().getPixelValue(i, j);
				float v2 = Vy.getProcessor().getPixelValue(i, j);
				ov.add(new Line(i+0.5, j+0.5, i+0.5+mult*v1, j+0.5+mult*v2));
			}
		}
		
		double dr = 1;
		double darc = 1;
		double rratio = 0.2;
		
		// count the number of points
		int cnt = 0;
		for (double r = radius; r >= radius*rratio; r-=dr) {
			for (double arc = 0; arc < 2*r*Math.PI; arc+=darc) {
				cnt++;
			}
		}
		
		System.out.println("test... " + cnt);
		
		double[] orts		= new double[cnt];
		
		double[] profile1  	= new double[cnt];
		double[] coeffs1  	= new double[cnt];
		
		double[] profile2  	= new double[cnt];
		double[] coeffs2  	= new double[cnt];
		
		double[] profile3  	= new double[cnt];
		double[] coeffs3  	= new double[cnt];
		
		cnt = 0;
		for (double r = radius; r >= radius*rratio; r-=dr) {
			for (double arc = 0; arc < 2*r*Math.PI; arc+=darc) {
				
				double ang = arc/r;
//				int idxAng = (int)Math.floor(ang/(2*Math.PI/(2*M)));
				
				double d1 = Math.sin(ang);
				double d2 = -Math.cos(ang);
				
				double x2 = mouseX + r * d1;
				double y2 = mouseY + r * d2;
				
				ov.add(new Line(x2+0.5, y2+0.5, x2+0.5+1*d1, y2+0.5+1*d2));
				//ov.add(new PointRoi(x2+0.5, y2+0.5));
				
				int x_loc = (int)Math.round(x2);
				int y_loc = (int)Math.round(y2);
				
				float mult1 = weighted.getProcessor().getPixelValue(x_loc, y_loc);
				float mult2 = neuriteness.getProcessor().getPixelValue(x_loc, y_loc);
				float mult3 = img.getProcessor().getPixelValue(x_loc, y_loc); 
				
				float v1 = Vx.getProcessor().getPixelValue(x_loc, y_loc);
				float v2 = Vy.getProcessor().getPixelValue(x_loc, y_loc);
				
				profile1[cnt] = Math.abs(mult1*v1*d1+mult1*v2*d2);
				profile2[cnt] = Math.abs(mult2*v1*d1+mult2*v2*d2);
				profile3[cnt] = Math.abs(mult3*v1*d1+mult3*v2*d2);
				
				orts[cnt]		= ang*360/(2*Math.PI);
				
				coeffs1[cnt]	= mult1;
				coeffs2[cnt]	= mult2;
				coeffs3[cnt]	= mult3;
				
				cnt++;
				
			}
			
		}
		
		double profile1_max = Double.MIN_VALUE;
		double profile2_max = Double.MIN_VALUE;
		double profile3_max = Double.MIN_VALUE;
		
		double profile1_min = Double.MAX_VALUE;
		double profile2_min = Double.MAX_VALUE;
		double profile3_min = Double.MAX_VALUE;
		
		for (int k = 0; k < cnt; k++) {
			if(coeffs1[k]>profile1_max) profile1_max = coeffs1[k];
			if(coeffs1[k]<profile1_min) profile1_min = coeffs1[k];
			
			if(coeffs2[k]>profile2_max) profile2_max = coeffs2[k];
			if(coeffs2[k]<profile2_min) profile2_min = coeffs2[k];
			
			if(coeffs3[k]>profile3_max) profile3_max = coeffs3[k];
			if(coeffs3[k]<profile3_min) profile3_min = coeffs3[k];
		}
		
//		double[] angles = new double[2*M];
//		for (int i = 0; i < 2*M; i++) {
//			angles[i] = i*(2*180/(2*M));
//		}
		
		System.out.println("ready to add");
		
		ImageStack show_profiles = new ImageStack(600, 300);  

		Plot p1 = new Plot("", "angle [deg]", "gaborAng");
		p1.setLimits(0, 360, profile1_min-1, profile1_max+1);
		p1.setSize(600, 300);
		p1.setColor(Color.RED);
		p1.addPoints(orts, profile1, Plot.BOX);
		p1.setColor(Color.BLACK);
		p1.addPoints(orts, coeffs1, Plot.X);
		//p1.draw();
		//
		show_profiles.addSlice(p1.getProcessor());
		
		Plot p2 = new Plot("", "angle [deg]", "neurness");
		p2.setLimits(0, 360, profile2_min-1, profile2_max+1);
		p2.setSize(600, 300);
		p2.setColor(Color.RED);
		p2.addPoints(orts, profile2, Plot.BOX);
		p2.setColor(Color.BLACK);
		p2.addPoints(orts, coeffs2, Plot.X);
		//p2.draw();
		show_profiles.addSlice(p2.getProcessor());
		
		Plot p3 = new Plot("", "angle [deg]", "orig.");
		p3.setLimits(0, 360, profile3_min-1, profile3_max+1);
		p3.setSize(600, 300);
		p3.setColor(Color.RED);
		p3.addPoints(orts, profile3, Plot.BOX);
		p3.setColor(Color.BLACK);
		p3.addPoints(orts, coeffs3, Plot.X);
		//p3.draw();
		show_profiles.addSlice(p3.getProcessor());
		
		new ImagePlus("response", show_profiles).show();
		
		img.setOverlay(ov);
		
		System.out.println("done. ");
		
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
