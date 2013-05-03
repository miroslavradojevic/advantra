package advantra.feature;

import java.awt.Color;

import ij.ImageStack;
import ij.gui.Plot;
import ij.process.ImageProcessor;

public class CircularFilterConfiguration {

	public int 			angResDeg;
	public double		angResRad;
	
	public int			nrPeaks;
	public int[] 		angBtwPeakDeg;
	public double[] 	angBtwPeakRad;
	
	public int 			nrRot;
	public double 		rotStepRad;
	
	public double 		wdthRad;
	public double[][] 	peaksRad;
	
	public CircularFilterConfiguration(int angular_resolution_deg, int[] angles_between_peaks_deg){
		
		this.angResDeg = angular_resolution_deg;
		this.angResRad = (angResDeg/360f)*2*Math.PI;
		
		this.angBtwPeakDeg = new int[angles_between_peaks_deg.length];
		for (int i = 0; i < angles_between_peaks_deg.length; i++) {
			angBtwPeakDeg[i] = angles_between_peaks_deg[i];
		}
		
		this.angBtwPeakRad = new double[angles_between_peaks_deg.length];
		for (int i = 0; i < angles_between_peaks_deg.length; i++) {
			angBtwPeakRad[i] = (angles_between_peaks_deg[i]/360.0)*2*Math.PI;
		}
		this.nrPeaks = angles_between_peaks_deg.length;
		
		this.nrRot = 360/(angResDeg/2);
		this.rotStepRad = angResRad/2;
		
		this.wdthRad = angResRad/2;
		
		this.peaksRad = new double[nrRot][nrPeaks];
		
		for (int cnt_rots = 0; cnt_rots < nrRot; cnt_rots++) {
			
			double start_pos = cnt_rots*rotStepRad;
					
			for (int cnt_shfs = 0; cnt_shfs < nrPeaks; cnt_shfs++) {
				
				if(cnt_shfs==0)
					peaksRad[cnt_rots][cnt_shfs] 	= start_pos;
				else
					peaksRad[cnt_rots][cnt_shfs] 	= peaksRad[cnt_rots][cnt_shfs-1] + angBtwPeakRad[cnt_shfs-1];
				
			}
			
		}
			
	}
	
	public double calculateScore(double[] profile_2PI){
		
		/* calculate score on the configuration as 
		 * highest score on all rotated versions of  
		 * that configuration to be rotationally invariant
		 */
		
		double 	score = Double.NEGATIVE_INFINITY;//Double.MIN_VALUE;
		//int 	scoreIdx = -1;
		double[] filt = new double[profile_2PI.length];
		double[] best_filt = new double[profile_2PI.length];
		for (int i = 0; i < best_filt.length; i++) {
			best_filt[i] = -1.05;
		}
		
		int sumPos, sumNeg;
		
		double max_filt = Double.MIN_VALUE;
		
		for (int r = 0; r < nrRot; r++) {
			
			// create filt[]
			sumPos = sumNeg = 0;
			
			for (int i = 0; i < profile_2PI.length; i++) {
				
				double angle = i*((2*Math.PI)/profile_2PI.length);
				
				boolean isON = false;
				// check if it belongs to ON
				for (int p = 0; p < nrPeaks; p++) {
					if(Math.abs(wrap_PI(angle-peaksRad[r][p])) <= wdthRad/2){
						isON = true;
						break;
					}
				}
				
				if(isON){
					filt[i] = +1;
					sumPos++;
				}
				else{
					filt[i] = -1;
					sumNeg++;
				}
				
			}
			
			// normalize filt[]
			if(true){
			for (int i = 0; i < filt.length; i++) {
				if(filt[i]>0){
					filt[i] = filt[i] * ((float)sumNeg/(float)sumPos);// / (float)sumPos;// 
					if(filt[i]>max_filt) max_filt = filt[i];
				}
//				else{
//					filt[i] = filt[i];// / (float)sumNeg;
//				}
			}
			}
			
			// score for this filt[]
			double sc = 0;
			for (int i = 0; i < filt.length; i++) {
				sc += filt[i]*profile_2PI[i];
			}
			
			if(sc>score){
				score = sc;
				for (int i = 0; i < best_filt.length; i++) {
					best_filt[i] = filt[i];
				}
				
			}
				
		}
		
		double[] xValues = new double[filt.length];
		for (int i = 0; i < filt.length; i++) {
			xValues[i] = i*(360f/filt.length);
		}
		
		Plot p = new Plot("", "ang [deg]", "");
		p.setLimits(0, 360, -1.1, max_filt);
		p.addPoints(xValues, profile_2PI, Plot.X);
		p.addPoints(xValues, best_filt, Plot.LINE);
		p.setSize(400, 200);
		p.show();
		
		return score;
	}
	
	public ImageStack plotFilter(){
		
		ImageStack viz_is = new ImageStack(600, 300);
		
		int res = 128;
		
		for (int rotIdx = 0; rotIdx < nrRot; rotIdx++) {
			
			double[] filt = new double[res];
			double[] angles = new double[res];
			
			int sumPos = 0;
			int sumNeg = 0;
			
			for (int i = 0; i < res; i++) {
				
				double angle = i * (2*Math.PI / res);
				angles[i] = angle;
				
				boolean isON = false;
				// check if it belongs to ON
				for (int p = 0; p < nrPeaks; p++) {
					if(Math.abs(wrap_PI(angle-peaksRad[rotIdx][p])) <= wdthRad/2){
						isON = true;
						break;
					}
				}
				
				if(isON){
					filt[i] = +1;
					sumPos++;
				}
				else{
					filt[i] = -1;
					sumNeg++;
				}
				
			}
			
			// normalize filt[]
			if(true){
			for (int i = 0; i < filt.length; i++) {
				if(filt[i]>0){
					filt[i] = filt[i] * ((float)sumNeg/(float)sumPos);
				}
			}
			}
			
			Plot p = new Plot("filter","angle","filter value", angles, filt);
			p.setSize(600, 300);
			
			viz_is.addSlice("rot."+rotIdx+"/"+nrRot, p.getProcessor());
			
		}
		
		return viz_is;
		
	}
	
	public ImageStack plot(){
		
		ImageStack viz_is = new ImageStack(300, 300);
		int plotResolution = 128;
		for (int rotIdx = 0; rotIdx < nrRot; rotIdx++) {
			
			float[] x = new float[plotResolution];
			float[] y = new float[plotResolution];
			
			for (int k = 0; k < plotResolution; k++) {
				
				double thetaRad 	= k*((2*Math.PI)/plotResolution);
				double ro = 1.5; // OFF
				
				for (int l = 0; l < nrPeaks; l++) {
					if(Math.abs(wrap_PI(thetaRad-peaksRad[rotIdx][l])) <= wdthRad/2){
						ro=2.5;
						break;
					}
				}
				
				x[k] 	= (float) (ro * Math.sin(thetaRad));
				y[k]	= (float) (ro * Math.cos(thetaRad));
				
			}
			
			int A = (int)(360.0/angResDeg);
			double[] x1 = new double[A];
			double[] y1 = new double[A];
			for (int k = 0; k < A; k++) {
				x1[k] = 2 * Math.sin(k*angResRad);
				y1[k] = 2 * Math.cos(k*angResRad);
			}
			
			Plot p = new Plot("profile","x","y", x, y);
			p.setSize(300, 300);
			p.setLimits(-2.5, 2.5, -2.5, 2.5);
			p.setJustification(Plot.CENTER);
			p.addPoints(x1, y1, Plot.BOX);
			
			viz_is.addSlice("rot."+rotIdx+"/"+nrRot, p.getProcessor());
			
		}
		
		return viz_is;
		
	}
	
	private static double wrap_PI(double in){
		
		double out = in;
		
		while(out<=-Math.PI){
			out += 2*Math.PI;
		}
		while(out>Math.PI){
			out -= 2*Math.PI;
		}
		
		return out;
	}
	
}
