package advantra.critpoint;

import java.awt.Color;
import java.io.File;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.Vector;

import advantra.feature.GaborFilt2D;
import advantra.file.AnalyzeCSV;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.PointRoi;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import ij.plugin.ZProjector;
import imagescience.feature.Differentiator;
import imagescience.feature.Hessian;
import imagescience.image.Axes;
import imagescience.image.Coordinates;
import imagescience.image.Dimensions;
import imagescience.image.FloatImage;
import imagescience.image.Image;

public class ExtractFeatures implements PlugIn {

	ImagePlus 	img, gab, gabAll, weighted, neuriteness, Vx, Vy;
	Image		inimg;
	
	String 		train_folder, test_folder;
	
	// multi-scale
	double 		t1, t2; 
	int 		tn;
	double[] 	s, t;
	
	// gabor
	int			M;
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
	
	public void run(String arg0) {
		
		
		t1  	= Prefs.get("advantra.critpoint.start_scale", 		3.0);
		t2    	= Prefs.get("advantra.critpoint.end_scale", 		11.0);
		tn 		= (int)Prefs.get("advantra.critpoint.nr_scales", 	5);
		M		= (int)Prefs.getInt("advantra.critpoint.nr_angles", 8);
		dr    	= Prefs.get("advantra.critpoint.dr", 				1.0);
		darc    = Prefs.get("advantra.critpoint.darc", 				1.0);
		dratio	= Prefs.get("advantra.critpoint.dratio", 			0.2);
		
		
		train_folder = (String)Prefs.get("advantra.critpoint.train_folder", 
				(System.getProperty("user.home")+File.separator));
		test_folder = (String)Prefs.get("advantra.critpoint.test_folder", 
				(System.getProperty("user.home")+File.separator));
		
		GenericDialog gd = new GenericDialog("ExtractFeatures");
		gd.addNumericField("start scale", t1, 1);
		gd.addNumericField("end   scale", t2, 1);
		gd.addNumericField("nr   scales", tn, 0);
		
		gd.addNumericField("angles(per 180 deg)", M,	0, 5, "");
		
		gd.addNumericField("radius step", 		dr, 1);
		gd.addNumericField("arc    step", 		darc, 1);
		gd.addNumericField("rratio step", 		dratio, 1);
		
		gd.addStringField("train folder : ", train_folder, 	50);
		gd.addStringField("test  folder : ", test_folder, 	50);
		
		gd.showDialog();
		if (gd.wasCanceled()) return;
		t1 	= 		gd.getNextNumber();
		t2	= 		gd.getNextNumber();
		tn	= (int)	gd.getNextNumber();
		
		M 	= (int)	gd.getNextNumber();
		
		dr  = 			gd.getNextNumber();
		darc  = 		gd.getNextNumber();
		dratio  = 		gd.getNextNumber();
		
		train_folder = 	gd.getNextString();
		test_folder = 	gd.getNextString();
		
		Prefs.set("advantra.critpoint.start_scale", t1);
		Prefs.set("advantra.critpoint.end_scale", 	t2);
		Prefs.set("advantra.critpoint.nr_scales", 	tn);
		
		Prefs.set("advantra.critpoint.nr_angles", 	M);
		Prefs.set("advantra.critpoint.dr", 			dr);
		Prefs.set("advantra.critpoint.darc", 		darc);
		Prefs.set("advantra.critpoint.dratio", 		dratio);
		
		Prefs.set("advantra.critpoint.train_folder", train_folder);
		Prefs.set("advantra.critpoint.test_folder", test_folder);
		
		
		// scales
		t = new double[tn];
		s = new double[tn];
		for (int i = 0; i < tn; i++) {
			t[i] = (i==0)? t1 : t1+i*((t2-t1)/(tn-1));
			s[i] = Math.sqrt(t[i]);
		}
				
		radius = 3*Math.sqrt(t[t.length-1]); // xGaussianStd.
			
		// angles theta angle it makes with x axis (angle it makes with the first row)
		theta_pi 	= new double[M];
		theta_2pi	= new double[2*M];
		for (int i = 0; i < 2*M; i++) {
			theta_2pi[i] = i * (Math.PI/(double)(M));
		}
		for (int i = 0; i < M; i++) {
			theta_pi[i] = i * (Math.PI/(double)M);
		}
		
		File dir = new File(train_folder);
		train_folder = dir.getAbsolutePath();
		if(!dir.isDirectory() ){ //|| !dir2.isDirectory()
			IJ.error("Wrong directory!");
			return;
		}
		
		File[] files_tif = listFilesEndingWith(dir, ".tif");
		File[] files_pos = new File[files_tif.length];
		File[] files_neg = new File[files_tif.length];
		Vector<double[][]> locs_pos = new Vector<double[][]>(files_tif.length);
		Vector<double[][]> locs_neg = new Vector<double[][]>(files_tif.length);
		
		int total_pos = 0;
		int total_neg = 0;
		int curr_pos = 0;
		int curr_neg = 0;
		
		System.out.println("## TRAINING ##  "+train_folder);
		
		for (int i = 0; i < files_tif.length; i++) { // for each tif file
			
			System.out.println("\n "+files_tif[i].getName()+"  ...  ");
			
			curr_pos = 0;
			curr_neg = 0;
			
			AnalyzeCSV readCSV;
			File[] check;
			String suffix;
			String file_name = files_tif[i].getName();
			file_name = file_name.substring(0, file_name.length()-4);
			
			suffix = file_name+".pos";
			check = listFilesEndingWith(dir, suffix);
			if(check!=null && check.length>0){
				files_pos[i] = check[0];
				readCSV = new AnalyzeCSV(files_pos[i].getAbsolutePath());
				double[][] bifsColRow = readCSV.readLn(2);
				locs_pos.add(bifsColRow);
				curr_neg = bifsColRow.length;
				total_pos += bifsColRow.length;
				System.out.println(bifsColRow.length+" positives ");
			}
			else{
				locs_pos.add(null);
				System.out.println("no poss");
			}
			
			suffix = file_name+".neg";
			check = listFilesEndingWith(dir, suffix);
			if(check!=null && check.length>0){
				files_neg[i] = check[0];
				readCSV = new AnalyzeCSV(files_neg[i].getAbsolutePath());
				double[][] othsColRow = readCSV.readLn(2);
				locs_neg.add(othsColRow);
				curr_pos = othsColRow.length;
				total_neg += othsColRow.length;
				System.out.println(othsColRow.length+" negatives ");
			}
			else{
				locs_neg.add(null);
				System.out.println("no negs");
			}
			
			/*
			 * actual feature extraction
			 */
			if(curr_pos>0 || curr_neg>0){
				
				System.out.println("...extract features...");
				
				img = new ImagePlus(files_tif[i].getAbsolutePath());
				resetCalibration();
				convertToFloatImage();
				IJ.freeMemory();
				
				H = img.getHeight();
				W = img.getWidth();
				
				inimg = Image.wrap(img); inimg.axes(Axes.X+Axes.Y);
				inimg.name("input_image");
				
				img.show();
				
				//gab = extractDirectionalGabor(img, theta_pi, t, bandwidth, psi, gamma, isReal); 
				
				
				img.close();
				
				
			}
			
			
		}
		
		System.out.println("////\n total (+) : "+total_pos);
		System.out.println(" total (-) : "+total_neg+"\n////\n");
		
	}
	
	private double[][] trueAdaBoost(double[][] featsP, double[][] featsN, int T){
		
		int sizeP = featsP.length;
	    int sizeN = featsN.length;
	    int fSize = featsP[0].length;
	    int m = sizeP + sizeN;
	    
	    double[] w = new double[m];
	    double[][] adaboost = new double[T][3]; // [id, alpha, threshold]
	    
	    for (int i = 0; i < sizeP; i++) {
	        w[i] = 1/(float)(2*sizeP);
	    }
	    for (int i = sizeP; i < sizeP+sizeN; i++) {
	        w[i] = 1/(float)(2*sizeN);
	    }
		
	    int tt = T;
	    
	    for (int t = 0; t < T; t++) {
	    	
	    	int steps = 500;
	        double[][] thresh = weightedWeakClassification(featsP, featsN, w, steps); 
	        // [feature_i][optimal threshold and the score /error]
	        //find minimum error 
	        double mineps = thresh[0][1];
	        int bestClassifier = 0;
	        for (int k = 1; k < fSize; k++) {
	            if (mineps > thresh[k][1]) {
	                mineps = thresh[k][1];
	                bestClassifier = k;
	            }
	        }
	        adaboost[t][0] = bestClassifier;
	        if (mineps == 0 || mineps > 0.5) {
	            adaboost[t][1] = 1;
	            adaboost[t][0] = bestClassifier;
	            adaboost[t][2] = thresh[bestClassifier][0];

	            IJ.log(t + " min eps = " + mineps + "   " +
	                    bestClassifier + "   " + thresh[bestClassifier][0]);
	            tt = t;
	            break;
	        }
	        
	        double beta = mineps / (1 - mineps);
	        adaboost[t][1] = Math.log(1 / beta);
	        adaboost[t][2] = thresh[bestClassifier][0];
	        
	        //update weights
	        for (int i = 0; i < sizeP; i++) {
	            if (applyClassifier(featsP[i][bestClassifier], thresh[bestClassifier][0])) {
	                w[i] *= beta;
	            }
	        }
	        for (int i = 0; i < sizeN; i++) {
	            if (!applyClassifier(featsN[i][bestClassifier], thresh[bestClassifier][0])) {
	                w[i + sizeP] *= beta;
	            }
	        }

	        //normalize weights
	        double sum = 0;
	        for (int i = 0; i < m; i++) {
	            sum += w[i];
	        }
	        for (int i = 0; i < m; i++) {
	            w[i] /= sum;
	        }
	        // number of runs - id of the feature - alpha - threshold - error 
	        IJ.log(
	        		"t = " + (t + 1) + " -> best classifier:" + 
	        		IJ.d2s(adaboost[t][0] + 1, 0) + "("+")"+ // nameFeatures[(int)adaboost[t][0]] +
	        		" Math.log(1 / beta):  " +
	                IJ.d2s(adaboost[t][1], 4) + "  optimal threshold: " + IJ.d2s(adaboost[t][2], 4) + "  error:   " +
	                IJ.d2s(mineps, 6));
	    }
	    
	    if (tt != T) {
	        double[][] new_adaboost = new double[tt + 1][3];
	        for (int i = 0; i < tt + 1; i++) {
	            new_adaboost[i] = adaboost[i];
	        }
	        return new_adaboost;
	    } 
	    else {
	        return adaboost;
	    }
	    
	}
		
	private double[][] weightedWeakClassification(double[][] imFeaturesP, double[][] imFeaturesN,
	        double[] w, int dt) {// featNr x [thr, score]
		
	    int sizeP = imFeaturesP.length;
	    int sizeN = imFeaturesN.length;
	    int fSize = imFeaturesN[0].length;
	    double[][] thresh = new double[fSize][2];

	    for (int k = 0; k < fSize; k++) {
	        // find max and min
	    	
	        double max = imFeaturesP[0][k];
	        double min = imFeaturesP[0][k];
	        
	        for (int i = 1; i < sizeP; i++) {
	            max = (max < imFeaturesP[i][k]) ? imFeaturesP[i][k] : max;
	            min = (min > imFeaturesP[i][k]) ? imFeaturesP[i][k] : min;
	        }
	        
	        for (int i = 0; i < sizeN; i++) {
	            max = (max < imFeaturesN[i][k]) ? imFeaturesN[i][k] : max;
	            min = (min > imFeaturesN[i][k]) ? imFeaturesN[i][k] : min;
	        }

	        double step = (max - min) / dt;
	        double thr = min;
	        double[][] count = new double[dt-1][2];
	        
	        for (int j = 0; j < dt - 1; j++) {
	            
	        	thr += step;
	            double score = 0.0;
	            
	            for (int i = 0; i < sizeP; i++) {
	                if (!applyClassifier(imFeaturesP[i][k], thr)) {
	                    score += w[i];
	                }
	            }
	            for (int i = 0; i < sizeN; i++) {
	                if (applyClassifier(imFeaturesN[i][k], thr)) {
	                    score += w[i + sizeP];
	                }
	            }
	            
	            count[j][0] = thr; 
	            count[j][1] = score;
	        
	        }
	        
	        // find optimal threshold
	        thresh[k][1] = count[0][1]; // first score
	        for (int j = 1; j < dt - 1; j++) {
	            if (thresh[k][1] > count[j][1]) {
	                thresh[k][1] = count[j][1];
	                thresh[k][0] = count[j][0];
	            }
	        }
	        
	    }

	    return thresh;

	}
	
	private boolean applyClassifier(double x, double thresh) {
	    return (x >= thresh) ? true : false;
	}
	
	private Plot plotFearutePerAllSamples(double[][] imFeaturesP, double[][] imFeaturesN, int c, double thresh) {
	    
		int sizep = imFeaturesP.length;
	    int sizen = imFeaturesN.length;
	    
	    double[] xp = new double[sizep];
	    double[] xn = new double[sizen];
	    double[] yp = new double[sizep];
	    double[] yn = new double[sizen];
	    
	    double maxyp = imFeaturesP[0][c];
	    double minyp = imFeaturesP[0][c];
	    for (int i = 0; i < sizep; i++) {
	        xp[i] = i;
	        yp[i] = imFeaturesP[i][c];
	        maxyp = (maxyp < imFeaturesP[i][c]) ? imFeaturesP[i][c] : maxyp;
	        minyp = (minyp > imFeaturesP[i][c]) ? imFeaturesP[i][c] : minyp;
	    }
	    
	    double maxyn = imFeaturesN[0][c];
	    double minyn = imFeaturesN[0][c];
	    for (int i = 0; i < sizen; i++) {
	        xn[i] = i;
	        yn[i] = imFeaturesN[i][c];
	        maxyn = (maxyn < imFeaturesN[i][c]) ? imFeaturesN[i][c] : maxyn;
	        minyn = (minyn > imFeaturesN[i][c]) ? imFeaturesN[i][c] : minyn;
	    }

	    Plot plot = new Plot("feature id = " + c, "Sample id", "Feature value", xp, yp);
	    plot.setLimits(0, Math.max(sizen, sizep), Math.min(minyn, minyp), Math.max(maxyn, maxyp));
	    plot.setColor(Color.RED);
	    
	    plot.addPoints(xn, yn, Plot.LINE);
	    plot.setColor(Color.BLACK);
	    plot.setLineWidth(5);
	    
	    double[] linex = new double[]{0, Math.max(sizen, sizep)};
	    double[] liney = new double[]{thresh, thresh};
	    plot.addPoints(linex, liney, Plot.LINE);
	    plot.setColor(Color.BLUE);
	    plot.setLineWidth(2);
	    
	    return plot;
	    
	}
	
	 int applyAdaBoost(double[][] adaboost, double[] imFeaturesT) {
		    double res = 0;
		    double object_res = 0;
		    for (int i = 0; i < adaboost.length; i++) {
		        object_res += adaboost[i][1];
		    }
		    object_res *= 0.5;

		    for (int i = 0; i < adaboost.length; i++) {
		        if (applyClassifier(imFeaturesT[i], adaboost[i][2])) {
		            res += adaboost[i][1];
		        }
		    }
		    int test = (res > object_res) ? 1 : 0;

		    return test;
		}
	
	private File[] listFilesEndingWith(File dir, String suffix){
		final String sfx = suffix;
		File[] tif_train = dir.listFiles(
			new FilenameFilter() {
				public boolean accept(File dir, String name) {
					return name.toLowerCase().endsWith(sfx);
				}
			}
		);
		return tif_train;
	}
	
	private void resetCalibration(){
		
		// reset the calibration
		Calibration c = img.getCalibration();
		c.pixelWidth 	= c.pixelHeight = c.pixelDepth = 1;
		c.setUnit("pixel");
		img.setCalibration(c);
				
	}
	
	private void convertToFloatImage(){
		
		if(!img.getProcessor().isDefaultLut()) return;
		img.setProcessor("input_image", img.getProcessor().convertToFloat().duplicate());
		
	}
	
	private ImagePlus extractDirectionalGabor(
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
		
		System.out.print("gabor filtering ");
		
		for (int i = 0; i < angles_pi.length; i++) {
			
			double current_theta = angles_pi[i];
			
			System.out.print(".");
			
			ImagePlus g_theta = GaborFilt2D.run(
					input, current_theta, t, new double[t.length], bw, psi, gamma, isReal);
			
			zmax.setImage(g_theta);
			zmax.setStartSlice(1);
			zmax.setStopSlice(g_theta.getStackSize());
			zmax.setMethod(ZProjector.MAX_METHOD);
			zmax.doProjection();
			angul.addSlice("theta="+IJ.d2s(current_theta, 2), zmax.getProjection().getChannelProcessor());

		}
		
		System.out.print(" done.");
		
		return new ImagePlus("gabor_per_angle", angul);

	}
}