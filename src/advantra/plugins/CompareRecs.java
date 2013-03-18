package advantra.plugins;

import java.awt.Color;
import java.io.File;

import javax.swing.JFileChooser;

import flanagan.analysis.Stat;

import advantra.file.AnalyzeSWC;
import ij.IJ;
import ij.gui.Plot;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;

public class CompareRecs implements PlugIn {

	String 	path1, path2;// = "";
	File	file1, file2;// = "";
	
	public void run(String arg0) {
		
		// open file 1
		
		File dir=null;
		JFileChooser fc = null;
		try {fc = new JFileChooser();}
		catch (Throwable e) {IJ.error("This plugin requires Java 2 or Swing."); return;}
		//fc.setMultiSelectionEnabled(true);
		if (dir==null) {
			String sdir = OpenDialog.getDefaultDirectory();
			if (sdir!=null)
				dir = new File(sdir);
		}
		if (dir!=null)
			fc.setCurrentDirectory(dir);
		int returnVal = fc.showOpenDialog(IJ.getInstance());
		if (returnVal!=JFileChooser.APPROVE_OPTION)
			return;
		file1 = fc.getSelectedFile();
		path1 = fc.getSelectedFile().getAbsolutePath();
		if(!file1.exists()){
			IJ.showMessage("file "+path1+" does not exist");
			return;
		}
//		else{
//			IJ.showMessage("Chosen file path"+path1);
//		}
		
		// open file 2
		
		dir=null;
		fc = null;
		try {fc = new JFileChooser();}
		catch (Throwable e) {IJ.error("This plugin requires Java 2 or Swing."); return;}
		//fc.setMultiSelectionEnabled(true);
		if (dir==null) {
			String sdir = OpenDialog.getDefaultDirectory();
			if (sdir!=null)
				dir = new File(sdir);
		}
		if (dir!=null)
			fc.setCurrentDirectory(dir);
		returnVal = fc.showOpenDialog(IJ.getInstance());
		if (returnVal!=JFileChooser.APPROVE_OPTION)
			return;
		file2 = fc.getSelectedFile();
		path2 = fc.getSelectedFile().getAbsolutePath();
		if(!file1.exists()){
			IJ.showMessage("file "+path2+" does not exist");
			return;
		}
//		else{
//			IJ.showMessage("Chosen file path"+path2);
//		}
		
		IJ.showMessage("Compare files:\n"+path1+"\n\n with respect to \n\n"+path2);
		
		AnalyzeSWC a1 = new AnalyzeSWC(path1);
		AnalyzeSWC a2 = new AnalyzeSWC(path2);
		
		System.out.println("-----\nFILE1: "+path1);
		a1.load();
		int L1 = a1.nodes.length;
  		System.out.format("total: \t\t\t%d \n", L1);
		
		System.out.println("-----\nFILE2: "+path2);
		a2.load();
		int L2 = a2.nodes.length;
  		System.out.format("total: \t\t\t%d \n", L2);
  		
  		double[] discr_d = new double[L1];
  		double[] discr_r = new double[L1];
  		
  		double 	discr_d_max = Double.MIN_VALUE; // plot limits
  		double	discr_r_max	= Double.MIN_VALUE;
  		
  		for (int i = 0; i < L1; i++) {
			
  			double min_d2 = Double.MAX_VALUE;
  			double corr_r = Double.NaN;
  			
  			for (int j = 0; j < L2; j++) {
				
  				double dx = a1.nodes[i].getX()-a2.nodes[j].getX();
  				double dy = a1.nodes[i].getY()-a2.nodes[j].getY();
  				double dz = a1.nodes[i].getZ()-a2.nodes[j].getZ();
  				
  				double d2 = dx*dx+dy*dy+dz*dz;
  				
  				if(d2<min_d2){
  					min_d2 = d2;
  					corr_r = a2.nodes[j].getR();
  				}
  				
			}
  			
  			discr_d[i] = Math.sqrt(min_d2);
  			discr_r[i] = Math.abs(a1.nodes[i].getR()-corr_r);
  			
  			if(discr_d[i]>discr_d_max){
  				discr_d_max = discr_d[i];
  			}
  			
  			if(discr_r[i]>discr_r_max){
  				discr_r_max = discr_r[i];
  			}
  			
		}

  		// plot discrepancies in d and r
		
		double[] idx = new double[L1];
		for (int i = 0; i < idx.length; i++) {
			idx[i] = i+1;
		}
		
		Plot plot = new Plot(String.format("Euclidean distance discrepancy"), "node #", "closest gnd-tth node distance", new double[0], new double[0]);
	    plot.setLimits(0, L1+1, 0, max(discr_d));
	    plot.setLineWidth(2);
	    plot.setColor(Color.red);
	    plot.addPoints(idx, discr_d, Plot.LINE);
	    plot.show(); 
	    
	    double[][] discr_d_hist = Stat.histogramBins(discr_d, 0.5);
	    
	    plot = new Plot("Euclidean distance discrepancy histogram", "node discrepancy", "freq.", new double[0], new double[0]);
	    plot.setLimits(0, max(discr_d_hist[0]), 0, max(discr_d_hist[1]));
	    plot.setLineWidth(2);
	    plot.setColor(Color.blue);
	    plot.addPoints(discr_d_hist[0], discr_d_hist[1], Plot.LINE);
	    plot.show(); 
	    
	    
	    System.out.println("\nEuclidean distance discrepancy:");
	    System.out.println("arithmetic mean: "+Stat.mean(discr_d)+" pix.");
	    System.out.println("std: "+Stat.standardDeviation(discr_d)+" pix.");
	    
		plot = new Plot(String.format("Radius discrepancy"), "node #", "closest gnd-tth node radius", new double[0], new double[0]);
	    plot.setLimits(0, L1+1, 0, discr_r_max);
	    plot.setLineWidth(2);
	    plot.setColor(Color.red);
	    plot.addPoints(idx, discr_r, Plot.LINE);
	    plot.show(); 
	    
	    double[][] discr_r_hist = Stat.histogramBins(discr_r, 0.5);
	    
	    plot = new Plot("Radius discrepancy histogram", "radius discrepancy", "freq.", new double[0], new double[0]);
	    plot.setLimits(0, max(discr_r_hist[0]), 0, max(discr_r_hist[1]));
	    plot.setLineWidth(2);
	    plot.setColor(Color.blue);
	    plot.addPoints(discr_r_hist[0], discr_r_hist[1], Plot.LINE);
	    plot.show();
	    
	    System.out.println("\nRadius estimation discrepancy:");
	    System.out.println("arithmetic mean: "+Stat.mean(discr_r)+" pix.");
	    System.out.println("std: "+Stat.standardDeviation(discr_r)+" pix.");
		
	}
	
	private double max(double[] in) {
		double mx = Double.MIN_VALUE;
		for (int i = 0; i < in.length; i++) {
			if(in[i]>mx){
				mx = in[i];
			}
		}
		return mx;
	}
	
	

}
