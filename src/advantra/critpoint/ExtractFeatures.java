package advantra.critpoint;

import java.awt.Color;
import java.io.File;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.Vector;

import advantra.file.AnalyzeCSV;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.PointRoi;
import ij.plugin.PlugIn;
import imagescience.feature.Differentiator;
import imagescience.feature.Hessian;
import imagescience.image.Axes;
import imagescience.image.Coordinates;
import imagescience.image.Dimensions;
import imagescience.image.FloatImage;
import imagescience.image.Image;

public class ExtractFeatures implements PlugIn {

	String 		train_folder, test_folder;
	double 		s1, s2; 
	int 		sn;
	String[] nameFeatures;
	
	public void run(String arg0) {
		
		s1  	= Prefs.get("advantra.critpoint.start_scale", 1.0);
		s2    	= Prefs.get("advantra.critpoint.end_scale", 5.0);
		sn 		= (int)Prefs.get("advantra.critpoint.nr_scales", 5);
		train_folder = (String)Prefs.get("advantra.critpoint.train_folder", 
				(System.getProperty("user.home")+File.separator));
		test_folder = (String)Prefs.get("advantra.critpoint.test_folder", 
				(System.getProperty("user.home")+File.separator));
		
		GenericDialog gd = new GenericDialog("Extract Features");
		gd.addNumericField("start scale", s1, 1);
		gd.addNumericField("end   scale", s2, 1);
		gd.addNumericField("nr   scales", sn, 0);
		gd.addStringField("folder with images & markers: ", train_folder, 	50);
		gd.addStringField("folder with test images     : ", test_folder, 	50);
		gd.showDialog();
		if (gd.wasCanceled()) return;
		s1 	= 		gd.getNextNumber();
		s2	= 		gd.getNextNumber();
		sn	= (int)	gd.getNextNumber();
		train_folder = gd.getNextString();
		test_folder = gd.getNextString();
		
		Prefs.set("advantra.critpoint.start_scale", s1);
		Prefs.set("advantra.critpoint.end_scale", s2);
		Prefs.set("advantra.critpoint.nr_scales", sn);
		Prefs.set("advantra.critpoint.train_folder", train_folder);
		Prefs.set("advantra.critpoint.test_folder", test_folder);
		
		// scales
		double[] s = new double[sn];
		for (int i = 0; i < sn; i++) {
			s[i] = (i==0)? s1 : s1+i*((s2-s1)/(sn-1));
		}
		
		File dir = new File(train_folder);
		train_folder = dir.getAbsolutePath();
		if(!dir.isDirectory() ){ //|| !dir2.isDirectory()
			IJ.error("Wrong directory!");
			return;
		}
		
		File[] files_tif = listFilesEndingWith(dir, ".tif");
		File[] files_bif = new File[files_tif.length];
		File[] files_oth = new File[files_tif.length];
		Vector<double[][]> locs_pos = new Vector<double[][]>(files_tif.length);
		Vector<double[][]> locs_neg = new Vector<double[][]>(files_tif.length);
		
		int total_pos = 0;
		int total_neg = 0;
		
		/*
		 * 
		 */
		
		System.out.println("###TRAINING###");
		
		for (int i = 0; i < files_tif.length; i++) { // for each tif file
			
			System.out.print(files_tif[i].getName()+"  ->  ");
			
			AnalyzeCSV readCSV;
			
			File[] check;
			
			check = listFilesEndingWith(dir, files_tif[i].getName()+".bif.marker");
			if(check==null || check.length<=0){
				System.out.println("incomplete annotations"); return;
			}
			files_bif[i] = check[0];
			readCSV = new AnalyzeCSV(files_bif[i].getAbsolutePath());
			double[][] bifsColRow = readCSV.readLn(2);
			locs_pos.add(bifsColRow);
			total_pos += bifsColRow.length;
			System.out.print(bifsColRow.length+" positives ");
			
			check = listFilesEndingWith(dir, files_tif[i].getName()+".oth.marker");
			if(check==null || check.length<=0){
				System.out.println("incomplete annotations"); return;
			}
			files_oth[i] = check[0];
			readCSV = new AnalyzeCSV(files_oth[i].getAbsolutePath());
			double[][] othsColRow = readCSV.readLn(2);
			locs_neg.add(othsColRow);
			total_neg += othsColRow.length;
			System.out.print(othsColRow.length+" negatives ");
			
			System.out.println();
			
		}
		
		int nrFilters = 4;
		String[] nameFilters = new String[]{
				"Gradient", 
				"absL1", 
				"DoH", 
				"Laplacian"
				
		};
		
		nameFeatures = new String[nrFilters*sn];
		for (int scale_idx = 0; scale_idx < sn; scale_idx++) {
			for (int f = 0; f < nrFilters; f++) {
				nameFeatures[scale_idx*nrFilters+f] = 
						String.format("%s_scale_%.2f", nameFilters[f], s[scale_idx]);
			}
		}
		
		int total_feats = nrFilters*sn;
		
		double[][] featsPos = new double[total_pos][total_feats];
		double[][] featsNeg = new double[total_neg][total_feats];
		
		System.out.println("# positives "+total_pos+" from "+locs_pos.size()+" files");
		System.out.println("# negatives "+total_neg+" from "+locs_neg.size()+" files");
		System.out.println("# features  "+total_feats);
		
		int start_row_pos = 0;
		int start_row_neg = 0;
		int start_col = 0;
		
		for (int file_idx = 0; file_idx < files_tif.length; file_idx++) { // for each tif file
			
			ImagePlus img = new ImagePlus(files_tif[file_idx].getAbsolutePath());
			Overlay ovl = new Overlay();
			
			Differentiator 	diff 	= new Differentiator();
			Hessian			hess 	= new Hessian();
			
			Image Dx, Dy, Dxx, Dyy, Dxy, L1, L2;
			double[][] aDx 	= new double[img.getHeight()][img.getWidth()];
			double[][] aDy 	= new double[img.getHeight()][img.getWidth()];
			double[][] aDxx = new double[img.getHeight()][img.getWidth()];
			double[][] aDxy = new double[img.getHeight()][img.getWidth()];
			double[][] aDyy = new double[img.getHeight()][img.getWidth()];
			double[][] aL1 	= new double[img.getHeight()][img.getWidth()];
			double[][] aL2 	= new double[img.getHeight()][img.getWidth()];
			
			start_col = 0;
			
			for (int scale_idx = 0; scale_idx < sn; scale_idx++) {
				
				// process image at this scale
				Dx = diff.run(Image.wrap(img), s[scale_idx], 1, 0, 0); Dx.axes(Axes.X+Axes.Y);
				Dy = diff.run(Image.wrap(img), s[scale_idx], 0, 1, 0); Dy.axes(Axes.X+Axes.Y);
				Dxx = diff.run(Image.wrap(img), s[scale_idx], 2, 0, 0); Dxx.axes(Axes.X+Axes.Y);
				Dxy = diff.run(Image.wrap(img), s[scale_idx], 1, 1, 0); Dxy.axes(Axes.X+Axes.Y);
				Dyy = diff.run(Image.wrap(img), s[scale_idx], 0, 2, 0); Dyy.axes(Axes.X+Axes.Y);
				Vector<Image> out = hess.run(Image.wrap(img), s[scale_idx], true);
				L2 = out.get(0); L2.axes(Axes.X+Axes.Y);
				L1 = out.get(1); L1.axes(Axes.X+Axes.Y);
				Coordinates coord = new Coordinates();
				Dx.get(coord, aDx);
				Dy.get(coord, aDy);
				Dxx.get(coord, aDxx);
				Dxy.get(coord, aDxy);
				Dyy.get(coord, aDyy);
				L2.get(coord, aL2);
				L1.get(coord, aL1);
				
				/*
				 * add positives
				 */
				
				for (int k = 0; k < locs_pos.get(file_idx).length; k++) {
					
					int col = (int)Math.round(locs_pos.get(file_idx)[k][0]);
					int row = (int)Math.round(locs_pos.get(file_idx)[k][1]);
					
					ovl.add(new PointRoi(col, row));
					
					for (int l = 0; l < nrFilters; l++) {
						
						double val = Double.NaN;
						switch(l){
							case 0:
								// Gradient Intensity
								val = Math.sqrt(Math.pow(aDx[row][col], 2)+Math.pow(aDy[row][col], 2));
								break;
							case 1:
								// abs(L1)
								val = Math.abs(aL1[row][col]);
								break;
							case 2:
								// DoH
								val = aDxx[row][col]*aDyy[row][col]-aDxy[row][col]*aDxy[row][col];
								break;
							case 3:
								// Laplacian
								val = aDxx[row][col]+aDyy[row][col];
								break;
							default:
								break;
						}
						
						featsPos[start_row_pos+k][start_col+l] = val;
						
					}
					
				}
				
				/*
				 * add negatives
				 */
				
				for (int k = 0; k < locs_neg.get(file_idx).length; k++) {
					
					int col = (int)Math.round(locs_neg.get(file_idx)[k][0]);
					int row = (int)Math.round(locs_neg.get(file_idx)[k][1]);
					
					ovl.add(new PointRoi(col, row));
					
					for (int l = 0; l < nrFilters; l++) {
						
						double val = Double.NaN;
						switch(l){
							case 0:
								// Gradient Intensity
								val = Math.sqrt(Math.pow(aDx[row][col], 2)+Math.pow(aDy[row][col], 2));
								break;
							case 1:
								// abs(L1)
								val = Math.abs(aL1[row][col]);
								break;
							case 2:
								// DoH
								val = aDxx[row][col]*aDyy[row][col]-aDxy[row][col]*aDxy[row][col];
								break;
							case 3:
								// Laplacian
								val = aDxx[row][col]+aDyy[row][col];
								break;
							default:
								break;
						}
						
						featsNeg[start_row_neg+k][start_col+l] = val;

					}
					
				}
				
				start_col += nrFilters;
				
			}
			
			img.getProcessor().drawOverlay(ovl);
			img.show();
		
			start_row_pos += locs_pos.get(file_idx).length;
			start_row_neg += locs_neg.get(file_idx).length;
			
		} 
		
		System.out.println("done extracting features");
		System.out.println("positives "+featsPos.length+" x "+featsPos[0].length);
		System.out.println("negatives "+featsNeg.length+" x "+featsNeg[0].length);
		
		// train classifier
		int Tloops = 5;
		double[][] adaboost = trueAdaBoost(featsPos, featsNeg, Tloops); // T x [id, alpha, threshold]
		
	    // save the adaboost config file 
	    try {
	    	FileWriter fw;
	        fw = new FileWriter(train_folder+File.separator+"train" + "." + Tloops + ".adaboost");
	        
	        fw.write(IJ.d2s(s1, 2) + "\r\n");
	        fw.write(IJ.d2s(s2, 2) + "\r\n");
	        fw.write(IJ.d2s(sn, 0) + "\r\n");
	        fw.write(IJ.d2s(adaboost.length, 0) + "\r\n");
	        
	        for (int i = 0; i < adaboost.length; i++) {
	            fw.write(IJ.d2s(adaboost[i][0], 0) + "\r\n");
	            fw.write(IJ.d2s(adaboost[i][1], 6) + "\r\n");
	            fw.write(IJ.d2s(adaboost[i][2], 6) + "\r\n");
	        }
	        fw.close();
	    } 
	    catch (IOException e) {
	    }
	    
	    // plot selected features with thresholds 
	    ImagePlus imp1 = new ImagePlus();
	    ImageStack imstack = new ImageStack(528, 255);
	    for (int i = 0; i < adaboost.length; i++) {
	        Plot plot = plotFearutePerAllSamples(featsPos, featsNeg, (int) adaboost[i][0], adaboost[i][2]);
	        imstack.addSlice("thresh = " + IJ.d2s(adaboost[i][2], 3), plot.getProcessor());
	    }
	    imp1.setStack("adaboost train", imstack);
	    imp1.show();
	    
	    /*
	     * 
	     */
	    
	    int[] feature_idx = new int[adaboost.length];
	    int[] scale_idx = new int[adaboost.length];
	    int[] filter_idx = new int[adaboost.length];
	    
	    for (int i = 0; i < adaboost.length; i++) {
			feature_idx[i] = (int)adaboost[i][0];
			scale_idx[i] = feature_idx[i]/nrFilters;
			filter_idx[i] = feature_idx[i]%nrFilters;
			System.out.println("feature: "+feature_idx[i]+" scale: "+scale_idx[i]+" filter: "+filter_idx[i]);
		}
	    
	    /*
	     * 
	     */
	    
	    
	    System.out.println("###TESTING###");
	    dir = new File(test_folder);
	    test_folder = dir.getAbsolutePath();
		if(!dir.isDirectory() ){
			IJ.error("Wrong directory!");
			return;
		}
		
		files_tif = listFilesEndingWith(dir, ".tif");
		
		for (int i = 0; i < files_tif.length; i++) { 
			
			System.out.print(files_tif[i].getName()+"  ->  ");
			
			ImagePlus img = new ImagePlus(files_tif[i].getAbsolutePath());
			//Overlay ovl = new Overlay();
			
			Differentiator 	diff 	= new Differentiator();
			Hessian			hess 	= new Hessian();
			
			Image Dx, Dy, Dxx, Dyy, Dxy, L1, L2;
			double[][] aDx 	= new double[img.getHeight()][img.getWidth()];
			double[][] aDy 	= new double[img.getHeight()][img.getWidth()];
			double[][] aDxx = new double[img.getHeight()][img.getWidth()];
			double[][] aDxy = new double[img.getHeight()][img.getWidth()];
			double[][] aDyy = new double[img.getHeight()][img.getWidth()];
			double[][] aL1 	= new double[img.getHeight()][img.getWidth()];
			double[][] aL2 	= new double[img.getHeight()][img.getWidth()];
			
			// create image with selected features
			Dimensions dims = new Dimensions(img.getWidth(), img.getHeight(), 1, adaboost.length);
			Image featsT = new FloatImage(dims); featsT.axes(Axes.X+Axes.Y);
			double[][] aFeatsT = new double[img.getHeight()][img.getWidth()];
			
			for (int k = 0; k < adaboost.length; k++) {
				
				double choose_scale 	= s[scale_idx[k]];
				int    choose_filter 	= filter_idx[k];
				
				Dx = diff.run(Image.wrap(img), choose_scale, 1, 0, 0); Dx.axes(Axes.X+Axes.Y);
				Dy = diff.run(Image.wrap(img), choose_scale, 0, 1, 0); Dy.axes(Axes.X+Axes.Y);
				Dxx = diff.run(Image.wrap(img), choose_scale, 2, 0, 0); Dxx.axes(Axes.X+Axes.Y);
				Dxy = diff.run(Image.wrap(img), choose_scale, 1, 1, 0); Dxy.axes(Axes.X+Axes.Y);
				Dyy = diff.run(Image.wrap(img), choose_scale, 0, 2, 0); Dyy.axes(Axes.X+Axes.Y);
				Vector<Image> out = hess.run(Image.wrap(img), choose_scale, true);
				L2 = out.get(0); L2.axes(Axes.X+Axes.Y);
				L1 = out.get(1); L1.axes(Axes.X+Axes.Y);
				Coordinates coord = new Coordinates();
				Dx.get(coord, aDx);
				Dy.get(coord, aDy);
				Dxx.get(coord, aDxx);
				Dxy.get(coord, aDxy);
				Dyy.get(coord, aDyy);
				L2.get(coord, aL2);
				L1.get(coord, aL1);
				
				//
				for (int row = 0; row < img.getHeight(); row++) {
					for (int col = 0; col < img.getWidth(); col++) {
						
						double val = Double.NaN;
						switch(choose_filter){
							case 0:
								// Gradient Intensity
								val = Math.sqrt(Math.pow(aDx[row][col], 2)+Math.pow(aDy[row][col], 2));
								break;
							case 1:
								// abs(L1)
								val = Math.abs(aL1[row][col]);
								break;
							case 2:
								// DoH
								val = aDxx[row][col]*aDyy[row][col]-aDxy[row][col]*aDxy[row][col];
								break;
							case 3:
								// Laplacian
								val = aDxx[row][col]+aDyy[row][col];
								break;
							default:
								break;
						}
						
						aFeatsT[row][col] = val;
						
					}
				}
				
				featsT.set(new Coordinates(0, 0, 0, k), aFeatsT);
				
			}
			featsT.axes(Axes.T);
			
			// output image
			dims = new Dimensions(img.getWidth(), img.getHeight(), 1, 1);
			Image out = new FloatImage(dims); out.axes(Axes.X+Axes.Y);
			
			Coordinates read_coord = new Coordinates();
			double[] takeFeats = new double[adaboost.length];
			for (read_coord.x = 0; read_coord.x < dims.x; read_coord.x++) {
				for (read_coord.y = 0; read_coord.y < dims.y; read_coord.y++) {
					featsT.get(read_coord, takeFeats);
					int cls = applyAdaBoost(adaboost, takeFeats);
	                if (cls == 1) {
	                    out.set(read_coord, 255);
	                }
				}
			}
			
			out.name("afterclassification");
			out.imageplus().show();
			
		} // loop test images
	    
		
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
	        		IJ.d2s(adaboost[t][0] + 1, 0) + "(" + nameFeatures[(int)adaboost[t][0]] +")"+
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
	
}