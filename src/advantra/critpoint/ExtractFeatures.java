package advantra.critpoint;

import java.awt.Color;
import java.awt.image.ColorModel;
import java.io.File;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.Vector;

import advantra.feature.GaborFilt2D;
import advantra.feature.ProfileFilters;
import advantra.file.AnalyzeCSV;
import advantra.shapes.Point;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.EllipseRoi;
import ij.gui.GenericDialog;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.PointRoi;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import ij.plugin.ZProjector;
import ij.process.FloatProcessor;

public class ExtractFeatures implements PlugIn {

	ImagePlus 	img, gab, gabAll, weighted, neuriteness, Vx, Vy;
	
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
	
	// features
	double 		radius; // calculated wrt. hughest scale but can be fixed
	double 		dr, darc, rratio;
	int			surr=3;
	
	int 		W, H;
	
	int			nr_proc;
	
	// output will be examples (pos. and neg.)
	ImageStack pos_examples_profile;
	ImageStack pos_examples_patch;// not used now
	ImageStack neg_examples_profile;
	ImageStack neg_examples_patch;// not used now
	
	public void run(String arg0) {
		
		
		t1  	= Prefs.get("advantra.critpoint.start_scale", 		3.0);
		t2    	= Prefs.get("advantra.critpoint.end_scale", 		7.0);
		tn 		= (int)Prefs.get("advantra.critpoint.nr_scales", 	3);
		M		= (int)Prefs.get("advantra.critpoint.nr_angles", 	8);
		dr    	= Prefs.get("advantra.critpoint.dr", 				1.0);
		
		darc    = Prefs.get("advantra.critpoint.darc", 				1.0);
		rratio	= Prefs.get("advantra.critpoint.dratio", 			0.2);
		nr_proc = (int)Prefs.get("advantra.critpoint.nr_proc", 		4);
		
		train_folder = (String)Prefs.get("advantra.critpoint.train_folder", 
				(System.getProperty("user.home")+File.separator));
		test_folder = (String)Prefs.get("advantra.critpoint.test_folder", 
				(System.getProperty("user.home")+File.separator));
		
		GenericDialog gd = new GenericDialog("ExtractFeatures");
		gd.addNumericField("start scale", t1, 1);
		gd.addNumericField("end   scale", t2, 1);
		gd.addNumericField("nr   scales", tn, 			0, 5, "");
		gd.addNumericField("angles(per 180 deg)", M,	0, 5, "");
		
		gd.addMessage("circular extraction parameters");
		gd.addNumericField("x(lagest scale std)",	surr, 	0);
		gd.addNumericField("radius step", 		dr, 1);
		gd.addNumericField("arc    step", 		darc, 1);
		gd.addNumericField("rratio step", 		rratio, 1);
		
		gd.addStringField("train folder : ", train_folder, 	80);
		gd.addStringField("test  folder : ", test_folder, 	80);
		
		gd.addMessage("parallelization");
		gd.addNumericField("CPU #",					nr_proc, 0);
		
		gd.showDialog();
		if (gd.wasCanceled()) return;
		t1 	= 		gd.getNextNumber();
		t2	= 		gd.getNextNumber();
		tn	= (int)	gd.getNextNumber();
		M 	= (int)	gd.getNextNumber();
		
		surr = (int)gd.getNextNumber();
		dr  = 			gd.getNextNumber();
		darc  = 		gd.getNextNumber();
		rratio  = 		gd.getNextNumber();
		
		train_folder = 	gd.getNextString();
		test_folder = 	gd.getNextString();
		
		nr_proc = (int)gd.getNextNumber();
		
		Prefs.set("advantra.critpoint.start_scale", t1);
		Prefs.set("advantra.critpoint.end_scale", 	t2);
		Prefs.set("advantra.critpoint.nr_scales", 	tn);
		Prefs.set("advantra.critpoint.nr_angles", 	M);
		Prefs.set("advantra.critpoint.dr", 			dr);
		Prefs.set("advantra.critpoint.darc", 		darc);
		Prefs.set("advantra.critpoint.dratio", 		rratio);
		
		Prefs.set("advantra.critpoint.train_folder", train_folder);
		Prefs.set("advantra.critpoint.test_folder", test_folder);
		
		Prefs.set("advantra.critpoint.nr_proc", 	nr_proc);
		
		
		// scales
		t = new double[tn];
		s = new double[tn];
		for (int i = 0; i < tn; i++) {
			t[i] = (i==0)? t1 : t1+i*((t2-t1)/(tn-1));
			s[i] = Math.sqrt(t[i]);
		}
				
		radius = surr*Math.sqrt(t[t.length-1]); // xGaussianStd.
			
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
		
		
		int angular_resolution = (int)Math.ceil((Math.PI*2)/(darc/radius));// highest
//		int angular_resolution = 64;
		pos_examples_profile = new ImageStack(angular_resolution, 1);
		neg_examples_profile = new ImageStack(angular_resolution, 1);
		
		System.out.println("## TRAINING ##  "+train_folder);
		
		for (int i = 0; i < files_tif.length; i++) { // for each tif file
			
			System.out.println("\n "+files_tif[i].getName()+" "+i+"/"+(files_tif.length-1)+"  ...  ");
			
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
				double[][] A = readCSV.readLn(2);
				locs_pos.add(A);
				curr_pos = A.length;
				total_pos += A.length;
				System.out.println(A.length+" positives ");
			}
			else{
				locs_pos.add(null);
				System.out.println("no positives");
			}
			
			suffix = file_name+".neg";
			check = listFilesEndingWith(dir, suffix);
			if(check!=null && check.length>0){
				files_neg[i] = check[0];
				readCSV = new AnalyzeCSV(files_neg[i].getAbsolutePath());
				double[][] B = readCSV.readLn(2);
				locs_neg.add(B);
				curr_neg = B.length;
				total_neg += B.length;
				System.out.println(B.length+" negatives ");
			}
			else{
				locs_neg.add(null);
				System.out.println("no negatives");
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
				//img.show();
				
				H = img.getHeight();
				W = img.getWidth();
				
				System.out.print("...gabor...");
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
				for (int l = 0; l < gab_jobs.length; l++) {
					
					int start_interval 	= l*N/nr_proc;
					int end_interval	= (l+1)*N/nr_proc;
					
					gab_jobs[l] = new GaborFilt2D(start_interval,  end_interval);
					gab_jobs[l].start();
				}
				for (int l = 0; l < gab_jobs.length; l++) {
					try {
						gab_jobs[l].join();
					} catch (InterruptedException e) {
						e.printStackTrace();
					}
				}
				
				gab = new ImagePlus("", GaborFilt2D.gabor_directional_responses);
				gab.setTitle("gabor filter, directional responses");
				t1 = System.currentTimeMillis();
				System.out.println("done, "+((t1-t0)/1000f)+" seconds");
				//gab.show();
				
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
				//gabAll.show();
				
				boolean normalize = false;
				if(normalize){
					weighted = new ImagePlus("gabor filter, directional responses, max, min-max normalized", 
							VizFeatures.normalizeStackMinMax(gabAll.getStack()));
				}
				else{
					weighted = new ImagePlus("gabor filter, directional responses, max, min-max normalized", 
							gabAll.getStack());
				}
				
				/*
				 *  extract neuriteness & eigen vecs
				 */

				System.out.print("...neuriteness...");
				
				long t11 = System.currentTimeMillis();
				Vector<ImagePlus> nness = VizFeatures.extractNeuritenessAndEigenVec(img, s);
				System.out.println("done, "+((System.currentTimeMillis()-t11)/1000f)+" seconds");
				
				neuriteness = nness.get(0);
				neuriteness.setTitle("neuriteness");
				//neuriteness.show();
				
				Vx = nness.get(1);
				Vx.setTitle("Vx");
				
				Vy = nness.get(2);
				Vy.setTitle("Vy");
				
				// there are 4 options for input images and locations: neuriteness, gabAll, weighted, img
				
				// take those locations
				
				Overlay ovly = new Overlay();
				
				
				for (int k = 0; k < curr_pos; k++) {
					
					int atX = (int)locs_pos.get(i)[k][0];
					int atY = (int)locs_pos.get(i)[k][1];
					
					double[] extracted_profile = VizFeatures.extractProfile(img, Vx, Vy, atX, atY, radius, dr, darc, rratio, angular_resolution);
					
					FloatProcessor fp = new FloatProcessor(angular_resolution, 1, extracted_profile);
					
					pos_examples_profile.addSlice(fp);
					
					EllipseRoi pt = new EllipseRoi(atX-2.5, atY-2.5, atX+2.5, atY+2.5, 1);
					pt.setColor(Color.RED);
					ovly.addElement(pt);
					
				}
				
				for (int k = 0; k < curr_neg; k++) {
					
					int atX = (int)locs_neg.get(i)[k][0];
					int atY = (int)locs_neg.get(i)[k][1];
					
					double[] extracted_profile = VizFeatures.extractProfile(img, Vx, Vy, atX, atY, radius, dr, darc, rratio, angular_resolution);
					
					FloatProcessor fp = new FloatProcessor(angular_resolution, 1, extracted_profile);
					
					neg_examples_profile.addSlice(fp);
					
					PointRoi pt = new PointRoi(atX, atY);
					pt.setColor(Color.BLUE);
					ovly.addElement(pt);
					
				}
				
				ImagePlus showIt = new ImagePlus("VIZ_"+files_tif[i].getName(), img.getProcessor());
				ovly.setFillColor(Color.RED);
				showIt.setOverlay(ovly);
				showIt.show();
				
			} // if there were some
			
			
		} // loop files
		
		System.out.println("////\n total (+) : "+total_pos);
		System.out.println(" total (-) : "+total_neg+"\n////\n");
		
		ImagePlus pos_examples_image =  new ImagePlus("positive examples", pos_examples_profile);
		pos_examples_image.show();
		pos_examples_image.getCanvas().zoomIn(0, 0);
		pos_examples_image.getCanvas().zoomIn(0, 0);
		pos_examples_image.getCanvas().zoomIn(0, 0);
		pos_examples_image.getCanvas().zoomIn(0, 0);
		pos_examples_image.getCanvas().zoomIn(0, 0);
		pos_examples_image.getCanvas().zoomIn(0, 0);
		
		ImagePlus neg_examples_image =  new ImagePlus("negative examples", neg_examples_profile);
		neg_examples_image.show();
		neg_examples_image.getCanvas().zoomIn(0, 0);
		neg_examples_image.getCanvas().zoomIn(0, 0);
		neg_examples_image.getCanvas().zoomIn(0, 0);
		neg_examples_image.getCanvas().zoomIn(0, 0);
		neg_examples_image.getCanvas().zoomIn(0, 0);
		neg_examples_image.getCanvas().zoomIn(0, 0);
		
		/*
		 * generate filters to score on example profiles
		 */
		double angStep = Math.PI/8;
		ProfileFilters pft = new ProfileFilters(angular_resolution, angStep);
		pft.create();
		pft.showFilters();
		
		/*
		 * calculate features (filter scores)
		 */
		
		
		
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
	
	private int applyAdaBoost(double[][] adaboost, double[] imFeaturesT) {
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
	
}