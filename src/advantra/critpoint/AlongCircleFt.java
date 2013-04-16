package advantra.critpoint;

import java.awt.Color;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

import advantra.feature.GaborFilt2D;
import advantra.general.Sort;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.PointRoi;
import ij.measure.Calibration;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

import imagescience.image.Axes;
import imagescience.image.Image;

public class AlongCircleFt implements PlugInFilter, MouseListener {

	ImagePlus 	img;
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
	
	
	double[][] 	o;
	String 		img_path;
	double 		radius = 10.0; // fixed
	double[][] 	aDx;
	double[][] 	aDy;
	
	public void run(ImageProcessor arg0) {
		
		t1  	= Prefs.get("advantra.critpoint.start_scale", 		2.0);
		t2    	= Prefs.get("advantra.critpoint.end_scale", 		4.0);
		tn 		= (int)Prefs.get("advantra.critpoint.nr_scales", 	3);
		M		= (int)Prefs.getInt("advantra.critpoint.nr_angles", 8);
		
		GenericDialog gd = new GenericDialog("F1");
		
		gd.addNumericField("start scale", 		t1, 2);
		gd.addNumericField("end   scale", 		t2, 4);
		gd.addNumericField("nr   scales", 		tn, 3);
		gd.addNumericField("angles(per 180 deg)",M,	0, 5, "");
		gd.showDialog();
		if (gd.wasCanceled()) return;
		
		t1 	= 		gd.getNextNumber();
		t2	= 		gd.getNextNumber();
		tn	= (int)	gd.getNextNumber();
		M 	= (int)gd.getNextNumber();
		
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
		
		// angles theta angle it makes with x axis (angle it makes with the first row)
		theta_pi 	= new double[M];
		theta_2pi	= new double[2*M];
		for (int i = 0; i < 2*M; i++) {
			theta_2pi[i] = i * (Math.PI/(double)(M));
		}
		for (int i = 0; i < M; i++) {
			theta_pi[i] = i * (Math.PI/(double)M);
		}
		
//		// set of tangent direction vectors (always the same)
//		int N = 10;
//		double[][] vecs = new double[N][2];
//		for (int i = 0; i < N; i++) {
//			vecs[i][0] = Math.cos(theta_2pi[i]); // Vx
//			vecs[i][1] = Math.sin(theta_2pi[i]); // Vy
//		}
		
		// reset the calibration
		Calibration c = img.getCalibration();
		c.pixelWidth 	= c.pixelHeight = c.pixelDepth = 1;
		c.setUnit("pixel");
		img.setCalibration(c);
		
		// convert to float
		if(!img.getProcessor().isDefaultLut()) return;
		img.setProcessor("testImage", img.getProcessor().convertToFloat().duplicate());
	    
		int w = img.getWidth();
		int h = img.getHeight();
		
		img.getWindow().getCanvas().addMouseListener(this);
		inimg = Image.wrap(img); inimg.axes(Axes.X+Axes.Y);
		IJ.log("you can click!");
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
		
		
		IJ.log("you can click!");
	}

	public int setup(String arg0, ImagePlus arg1) {
		img = arg1;
		return DOES_8G+NO_CHANGES;
	}

	public void mouseClicked(MouseEvent e) {
		
		int mouseX = 	img.getWindow().getCanvas().offScreenX(e.getX());
		int mouseY = 	img.getWindow().getCanvas().offScreenY(e.getY());
		
		Overlay ov = new Overlay();
		
		ImageStack 	plot_stack 	= new ImageStack(528, 255);
		
		for (int scale_idx = 0; scale_idx < tn; scale_idx++) {
			// take locations on the circle
			
			double[] score = new double[1];
			
			for (int k = 0; k < 1; k++) {
				
				int rd = (int)Math.ceil(1 * s[scale_idx]);
				int x2 = mouseX + (int)Math.round(rd * Math.cos(theta_2pi[k]));
				int y2 = mouseY + (int)Math.round(rd * Math.sin(theta_2pi[k]));
				ov.add(new PointRoi((double)x2, (double)y2));
				//gradient value at this point 
				double Gx = aDx[y2][x2];
				//double Vx = vecs[k][0];
				double Gy = aDy[y2][x2];
				//double Vy = vecs[k][1];
				score[k] = Math.sqrt(Math.pow(Gx, 2)+Math.pow(Gy, 2));
				
			}
			
			
			double[] x_axis = new double[1];
			for (int i = 0; i < x_axis.length; i++) {
				x_axis[i] = i;
			}
			
			double min_score = Sort.findMin(score);
			double max_score = Sort.findMax(score);
			
			
			
			Plot p = new Plot("score", "ang. pos. along circle", "value");
			p.setLimits(0, 1, min_score, max_score);
		    p.setColor(Color.RED);
		    p.addPoints(x_axis, score, Plot.LINE);
		    plot_stack.addSlice(String.format("scale %.2f", s[scale_idx]), p.getProcessor());
		    
		}
		ImagePlus 	plot_score 	= new ImagePlus();
		plot_score.setStack(plot_stack);
		plot_score.show();
		
		img.getProcessor().drawOverlay(ov);
		
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
