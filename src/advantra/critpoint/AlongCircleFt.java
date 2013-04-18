package advantra.critpoint;

import java.awt.Color;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

import advantra.feature.GaborFilt2D;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.PointRoi;
import ij.measure.Calibration;
import ij.plugin.ZProjector;
import ij.plugin.filter.PlugInFilter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import imagescience.image.Axes;
import imagescience.image.Image;

public class AlongCircleFt implements PlugInFilter, MouseListener {

	ImagePlus 	img, gab, gabAll, gabAllWeight;
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
	
	double 		gab_min, gab_max;
	double		img_min, img_max;
	
	
	Plot p = new Plot("profile", "angle [rad]", "value");
	
	
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
		
		radius = 4*Math.sqrt(t[t.length-1]); // xGaussianStd.
		IJ.log("check radius is " + radius);
		
		// angles theta angle it makes with x axis (angle it makes with the first row)
		theta_pi 	= new double[M];
		theta_2pi	= new double[2*M];
		for (int i = 0; i < 2*M; i++) {
			theta_2pi[i] = i * (Math.PI/(double)(M));
		}
		for (int i = 0; i < M; i++) {
			theta_pi[i] = i * (Math.PI/(double)M);
		}
		
		// reset the calibration
		Calibration c = img.getCalibration();
		c.pixelWidth 	= c.pixelHeight = c.pixelDepth = 1;
		c.setUnit("pixel");
		img.setCalibration(c);
		
		// convert to float
		if(!img.getProcessor().isDefaultLut()) return;
		img.setProcessor("testImage", img.getProcessor().convertToFloat().duplicate());
	    
		img.getWindow().getCanvas().addMouseListener(this);
		inimg = Image.wrap(img); inimg.axes(Axes.X+Axes.Y);
		
		/*
		 * show gabor kernels for selected scales
		 */
		for (int i = 0; i < t.length; i++) {
			
			ImagePlus krn = GaborFilt2D.run(null, theta_pi, t[i], 0, bandwidth, psi, gamma, isReal);
			krn.setTitle("gabor_kernels_scale"+IJ.d2s(t[i],2));
			krn.show();
			for (int j = 0; j < 2; j++) {
				krn.getCanvas().zoomIn(0, 0);
			}
			
		}
		
		/*
		 * extract directional gabor responses
		 */
		gab = extractDirectionalGabor(img, theta_pi, t, bandwidth, psi, gamma, isReal); 
		
		/*
		 * extract maximal directional response in every point
		 */
		ZProjector zmax = new ZProjector();
		zmax.setImage(gab);
		zmax.setStartSlice(1);
		zmax.setStopSlice(gab.getStackSize());
		zmax.setMethod(ZProjector.MAX_METHOD);
		zmax.doProjection();
		gabAll = new ImagePlus("allThetas", zmax.getProjection().getChannelProcessor());
		gabAll.show();
		
		float[] mn_mx = calculateStackMinMax(gab.getStack());
		gab_min = mn_mx[0];
		gab_max = mn_mx[1];
		IJ.log("gabor outputs: "+gab_min+" / "+gab_max);
		
		IJ.log("redefine...");
		
		mn_mx = calculateStackMinMax(gabAll.getStack());
		gab_min = mn_mx[0];
		gab_max = mn_mx[1];
		IJ.log("gabor outputs: "+gab_min+" / "+gab_max);
		
		mn_mx = calculateStackMinMax(img.getStack());
		img_min = mn_mx[0];
		img_max = mn_mx[1];
		
		gabAllWeight = new ImagePlus("gabAllWeight", img.getProcessor().duplicate());
		gabAllWeight.duplicate().show();
		
		float[] temp = (float[])gabAllWeight.getProcessor().getPixels();
		for (int i = 0; i < temp.length; i++) {
			temp[i] = (float)Math.exp(2*Math.pow(1-(temp[i]/img_max), 2));// change wieghts
		}
		gabAllWeight.show();
		
		// extract vector field - show it on separate image
		
		
	}

	public int setup(String arg0, ImagePlus arg1) {
		img = arg1;
		return DOES_8G+NO_CHANGES;
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
		
		ov.add(new PointRoi(mouseX, mouseY));
		
		double[] profile1 = new double[theta_2pi.length];
		double[] profile2 = new double[theta_2pi.length];
		double[] profile3 = new double[theta_2pi.length];
		
		for (int k = 0; k < theta_2pi.length; k++) { //
			
			double x2 = mouseX + radius * Math.sin(theta_2pi[k]);
			double y2 = mouseY - radius * Math.cos(theta_2pi[k]);
			PointRoi pt = new PointRoi(x2, y2);
			ov.add(pt);
			
			int k1 = (k>=theta_pi.length)? k-theta_pi.length : k ;
//			IJ.log("adding "+x2+" , "+y2+" size= "+ov.size()+ "index " +k1);
			profile1[k] = gab.getStack().getProcessor(k1+1).getInterpolatedValue(x2, y2);
			profile2[k] = gabAll.getProcessor().getInterpolatedValue(x2, y2);
			profile3[k] = img.getProcessor().getInterpolatedValue(x2, y2);
		}
		
		img.setOverlay(ov);
		
		
		p.setLimits(0, 2*Math.PI, Math.min(img_min, gab_min), Math.max(img_max, gab_max));
		p.setSize(600, 300);
		
	    p.setColor(Color.RED);
	    p.addPoints(theta_2pi, profile1, Plot.BOX);
	    p.setColor(Color.GREEN);
	    p.addPoints(theta_2pi, profile2, Plot.LINE);
	    p.setColor(Color.BLUE);
	    p.addPoints(theta_2pi, profile3, Plot.CROSS);
	    p.setLineWidth(2);
	    p.show();

	    // calculate new profile
//	    double dr = radius/5;
	    int Q = 5;
	    double[] r = new double[Q];
	    for (int l = 0; l < Q; l++) {
	    	r[l] = (l+1)*(radius/(double)Q);
	    }
		
	    Overlay ov1 = new Overlay();
	    
	    
	    ImageStack  	circular_profile_stack 	= new ImageStack(2*M, Q);
	    ImageProcessor 	circular_profile1 		= new FloatProcessor(2*M, Q);
	    ImageProcessor 	circular_profile2 		= new FloatProcessor(2*M, Q);
	    
	    
	    for (int l = 0; l < Q; l++) {
	    	for (int k = 0; k < 2*M; k++) {

				double x2 = mouseX + r[l] * Math.sin(theta_2pi[k]);
				double y2 = mouseY - r[l] * Math.cos(theta_2pi[k]);
				ov1.add(new PointRoi(x2, y2));
				
				circular_profile1.setf(k, l, (float)gabAll.getProcessor().getInterpolatedValue(x2, y2));
				circular_profile2.setf(k, l, (float)img.getProcessor().getInterpolatedValue(x2, y2));
				
			}
		}
	    
	    circular_profile_stack.addSlice(circular_profile1);
	    circular_profile_stack.addSlice(circular_profile2);
	    
	    ImagePlus circular_profile_image = new ImagePlus("circular_profile", circular_profile_stack);
	    circular_profile_image.show();
	    for (int i = 0; i < 8; i++) {
	    	circular_profile_image.getCanvas().zoomIn(0, 0);
		}
	    
	    
	    
	    img.setOverlay(ov1);
	    
	    
	}

	public void mousePressed(MouseEvent e) {
		
	}

	public void mouseReleased(MouseEvent e) {
		
	}

	public void mouseEntered(MouseEvent e) {
		
	}

	public void mouseExited(MouseEvent e) {
		
	}
	
	
}
