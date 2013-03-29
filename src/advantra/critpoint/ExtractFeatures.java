package advantra.critpoint;

import java.awt.Color;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

import advantra.general.Sort;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.PointRoi;
import ij.measure.Calibration;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

import imagescience.feature.Differentiator;
import imagescience.image.Axes;
import imagescience.image.Coordinates;
import imagescience.image.Dimensions;
import imagescience.image.FloatImage;
import imagescience.image.Image;

public class ExtractFeatures implements PlugInFilter, MouseListener {

	int 		N = 50;
	double[][] 	o;
	double[] 	s;
	double[] 	theta;
	ImagePlus 	img;
	String 		img_path;
	double 		s1, s2; 
	int 		sn; 
	double 		circ = 3.0;
	double[][] 	aDx;
	double[][] 	aDy;
	
	public void run(ImageProcessor arg0) {
		
		s1  	= Prefs.get("advantra.critpoint.start_scale", 1.0);
		s2    	= Prefs.get("advantra.critpoint.end_scale", 8.0);
		sn 		= (int)Prefs.get("advantra.critpoint.nr_scales", 8);
		
		GenericDialog gd = new GenericDialog("F1");
		
		gd.addNumericField("start scale", s1, 1);
		gd.addNumericField("end   scale", s2, 1);
		gd.addNumericField("nr   scales", sn, 0);
		
		gd.showDialog();
		if (gd.wasCanceled()) return;
		
		s1 	= 		gd.getNextNumber();
		s2	= 		gd.getNextNumber();
		sn	= (int)	gd.getNextNumber();
		
		Prefs.set("advantra.critpoint.start_scale", s1);
		Prefs.set("advantra.critpoint.end_scale", s2);
		Prefs.set("advantra.critpoint.nr_scales", sn);
		
		// scales
		s = new double[sn];
		for (int i = 0; i < sn; i++) {
			s[i] = (i==0)? s1 : s1+i*((s2-s1)/(sn-1));
		}
		
		// angles theta angle it makes with x axis (angle it makes with the first row)
		theta = new double[N];
		for (int i = 0; i < theta.length; i++) {
			theta[i] = i*(2*Math.PI/(double)N);
		}
		
		// set of tangent direction vectors (always the same)
		double[][] vecs = new double[N][2];
		for (int i = 0; i < N; i++) {
			vecs[i][0] = Math.cos(theta[i]); // Vx
			vecs[i][1] = Math.sin(theta[i]); // Vy
		}
		
		// reset the calibration
		Calibration c = img.getCalibration();
		c.pixelWidth 	= c.pixelHeight = c.pixelDepth = 1;
		c.setUnit("pixel");
		img.setCalibration(c);
		img.getWindow().getCanvas().addMouseListener(this);
		
		
		Image Im, MyFeature;
		Im = Image.wrap(img); Im.axes(Axes.X+Axes.Y);
		
		Dimensions indim = Im.dimensions();
		Dimensions outdim = new Dimensions(indim.x, indim.y, indim.z, sn);
		
		MyFeature = new FloatImage(outdim); MyFeature.axes(Axes.X+Axes.Y);
		
		aDx = new double[indim.y][indim.x];
		aDy = new double[indim.y][indim.x];
		double[][] aIm = new double[indim.y][indim.x]; // rows, cols
		double[][] aMyF = new double[indim.y][indim.x];
		
		// derivatives
		Differentiator df = new Differentiator();
		Coordinates crd = new Coordinates();
		
		Im.get(crd, aIm);
		
		for (int scale_idx = 0; scale_idx < sn; scale_idx++) {
			
			Image Dx, Dy;
			
			Dx = df.run(Im.duplicate(), s[scale_idx], 1, 0, 0);
			Dy = df.run(Im.duplicate(), s[scale_idx], 0, 1, 0);
			
			Dx.axes(Axes.X+Axes.Y);
			Dy.axes(Axes.X+Axes.Y);
			
			Dx.get(crd, aDx);
			Dy.get(crd, aDy);
			
			// extract feature at this scale using gradient and orientations around the circle
			int rd = (int)Math.ceil(circ * s[scale_idx]);
			
//			int count = 0;
			
			for (int y = rd; y < aDx.length-rd; y++) {
				for (int x = rd; x < aDx[0].length-rd; x++) {
					
					double score = 0;
					// take locations on the circle
					for (int k = 0; k < N; k++) {
						
						int x2 = x + (int)Math.round(rd * Math.cos(theta[k]));
						int y2 = y + (int)Math.round(rd * Math.sin(theta[k]));
						
						//gradient valea at this point 
						double Gx = aDx[y2][x2];
						//double Vx = vecs[k][0];
						double Gy = aDy[y2][x2];
						//double Vy = vecs[k][1];
						score += Math.sqrt(Math.pow(Gx, 2)+Math.pow(Gy, 2));
						
					}
					
//					score = Math.sqrt(Math.pow(aDx[y][x], 2)+Math.pow(aDy[y][x], 2));
//					score = aDy[y][x];
					
					aMyF[y][x] = score;
				}
			}
			
//			img.getProcessor().drawOverlay(ov);

			// store it in tth t channel
			crd.t = scale_idx;
			MyFeature.set(crd, aMyF);
			crd.t = 0; // back to default
			
		}
		
		MyFeature.imageplus().show();
		
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
		
		ImagePlus plot_score 	= new ImagePlus();
		ImageStack plot_stack 	= new ImageStack(528, 255);
		
		for (int scale_idx = 0; scale_idx < sn; scale_idx++) {
			// take locations on the circle
			
			double[] score = new double[N];
			
			for (int k = 0; k < N; k++) {
				
				int rd = (int)Math.ceil(circ * s[scale_idx]);
				int x2 = mouseX + (int)Math.round(rd * Math.cos(theta[k]));
				int y2 = mouseY + (int)Math.round(rd * Math.sin(theta[k]));
				ov.add(new PointRoi((double)x2, (double)y2));
				//gradient value at this point 
				double Gx = aDx[y2][x2];
				//double Vx = vecs[k][0];
				double Gy = aDy[y2][x2];
				//double Vy = vecs[k][1];
				score[k] = Math.sqrt(Math.pow(Gx, 2)+Math.pow(Gy, 2));
				
			}
			
			
			double[] x_axis = new double[N];
			for (int i = 0; i < x_axis.length; i++) {
				x_axis[i] = i;
			}
			
			double min_score = Sort.findMin(score);
			double max_score = Sort.findMax(score);
			
			Plot p = new Plot("score", "ang. pos. along circle", "value");
			p.setLimits(0, N, min_score, max_score);
		    p.setColor(Color.RED);
		    p.addPoints(x_axis, score, Plot.LINE);
		    plot_stack.addSlice(String.format("scale %.2f", s[scale_idx]), p.getProcessor());
		    
		}
		
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
