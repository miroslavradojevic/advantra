package advantra.critpoint;

import java.awt.Button;
import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.Panel;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import javax.swing.JFileChooser;
import MLdetection.FeaturePool;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.plugin.PlugIn;
import imagescience.image.Axes;
import imagescience.image.Coordinates;
import imagescience.image.Dimensions;
import imagescience.image.Image;

public class TrainPatches implements PlugIn, ActionListener {

	Button posButton, negButton;
	int nr_runs, search_step, patch_size;
	String pos_train_file, neg_train_file, train_name;
	
	public void run(String arg0) {
		
		nr_runs 		= Prefs.getInt("advantra.critpoint.TrainPatches.nr_runs", 3);
		search_step 	= Prefs.getInt("advantra.critpoint.TrainPatches.search_step", 500);
		pos_train_file 	= Prefs.get("advantra.critpoint.TrainPatches.pos_train_file", "none");
		neg_train_file	= Prefs.get("advantra.critpoint.TrainPatches.neg_train_file", "none");
		train_name		= Prefs.getString("advantra.critpoint.TrainPatches.train_name", "train");
		
		/*
		 * generic dialog
		 */
		
		GenericDialog gd = new GenericDialog("Train Patches");
	    
		gd.addNumericField("Number of adaboost runs:", 			nr_runs, 0);
		gd.addNumericField("Search Step for weak classifier:", 	search_step, 0);
		
	    Panel buttons_panel = new Panel();
	    buttons_panel.setLayout(new FlowLayout(FlowLayout.CENTER, 5, 1));
	    posButton = new Button("Open stack with positive trainset");
	    posButton.addActionListener(this);
	    negButton = new Button("Open stack with negative trainset");
	    negButton.addActionListener(this);
	    
	    buttons_panel.add(posButton);
	    buttons_panel.add(negButton);
	    gd.addPanel(buttons_panel);
	    
	    gd.addMessage("default: "+pos_train_file);
	    gd.addMessage("default: "+neg_train_file);
	    
	    gd.showDialog();
	    if (gd.wasCanceled()) {
	        return;
	    }

	    nr_runs = (int) gd.getNextNumber();
		search_step = (int) gd.getNextNumber();
	    
	    // memorize parameters for the next execution
	    Prefs.set("advantra.critpoint.TrainPatches.nr_runs", nr_runs);
	    Prefs.set("advantra.critpoint.TrainPatches.search_step", search_step);
	    Prefs.set("advantra.critpoint.TrainPatches.pos_train_file", pos_train_file);
	    Prefs.set("advantra.critpoint.TrainPatches.neg_train_file", neg_train_file);
	    Prefs.set("advantra.critpoint.TrainPatches.train_name", train_name);
	    
	    // variables used
	    ImagePlus 		img;
	    Image 			img_in;
	    Coordinates 	coord;
	    Dimensions 		dims;
	    
	    // load positives
	    IJ.open(pos_train_file);
	    img = IJ.getImage();
	    img_in = Image.wrap(img); img_in.axes(Axes.X + Axes.Y);
	    img.close();
	    
	    // set patch size
	    if(img.getHeight()==img.getWidth()){
	    	patch_size = img.getHeight();
	    }
	    else{
	    	return;
	    }

	    // create the feature pool 
	    FeaturePool featurePool = new FeaturePool(patch_size);
	    int fSize = featurePool.getSize();
	    IJ.log("the feature pool of size " + fSize + " is created");
	    featurePool.visualizeFeaturePool();
	    
	    /*
	     *  extract positive samples
	     */
	    coord = new Coordinates();
	    dims = img_in.dimensions();
	    double[][] imageP = new double[patch_size][patch_size];
	    int[][] imFeaturesP = new int[dims.z][];
	    for (coord.z = 0; coord.z < dims.z; coord.z++) {
	    	
	    	img_in.get(coord, imageP);
	        int[][] imageIntP = getIntegralImage(imageP);
	        imFeaturesP[coord.z] = computeFeatures(imageIntP, featurePool);
	        
	    }
	    
//        // debug
//        for (int i = 0; i < dims.z; i++) {
//        	System.out.format("positive sample %d scored %d / %d at features (1147/8) \n", i, imFeaturesP[i][1146], imFeaturesP[i][1147]);
//		}
	    
	    IJ.log("features computed for "+dims.z+" (+) samples");
	    
	    // load negatives
	    IJ.open(neg_train_file);
	    img = IJ.getImage();
	    img_in = Image.wrap(img); img_in.axes(Axes.X + Axes.Y);
	    img.close();
	    // check the patch
	    if(img.getHeight()!=patch_size || img.getWidth()!=patch_size){
	    	return;
	    }
	    
	    /*
	     *  extract negative samples
	     */
	    coord = new Coordinates();
	    dims = img_in.dimensions();
	    double[][] imageN = new double[patch_size][patch_size];
	    int[][] imFeaturesN = new int[dims.z][];
	    for (coord.z = 0; coord.z < dims.z; coord.z++) {
	    	
	    	img_in.get(coord, imageN);
	        int[][] imageIntN = getIntegralImage(imageN);
	        imFeaturesN[coord.z] = computeFeatures(imageIntN, featurePool);
	        
	    }

//        // debug
//        for (int i = 0; i < dims.z; i++) {
//        	System.out.format("negative sample %d scored %d and %d at features (1147/8) \n", i, imFeaturesN[i][1146], imFeaturesN[i][1147]);
//		}
	    
	    IJ.log("features computed for "+dims.z+" (-) samples");
		
	    /*
	     * train
	     */
	    
	    double[][] adaboost = trueAdaBoost(imFeaturesP, imFeaturesN, nr_runs);
	    
	    /*
	     * save trainig results to .adaboost
	     */
	    FileWriter fw;
	    String out_file_path = System.getProperty("user.home")+File.separator+train_name + "." + nr_runs + ".adaboost";
	    try {
	        fw = new FileWriter(out_file_path);
	        fw.write(IJ.d2s(patch_size, 0) + "\r\n");
	        fw.write(IJ.d2s(adaboost.length, 0) + "\r\n");
	        for (int i = 0; i < adaboost.length; i++) {
	            fw.write(IJ.d2s(adaboost[i][0], 0) + "\r\n");
	            fw.write(IJ.d2s(adaboost[i][1], 6) + "\r\n");
	            fw.write(IJ.d2s(adaboost[i][2], 6) + "\r\n");
	        }
	        fw.close();
	    } catch (IOException e) {
	        //ij.IJ.error("Unable to save particle positions");
	    }
	    
	    IJ.log(out_file_path+"   saved.");
	    
	    // plot the features with thresholds 
	    ImagePlus imp1 = new ImagePlus();
	    ImageStack imstack = new ImageStack(528, 255);
	    for (int i = 0; i < adaboost.length; i++) {
	        Plot plot = plotFeaturePerAllSamples(imFeaturesP, imFeaturesN, (int) adaboost[i][0], adaboost[i][2]);
	        //Plot plot = plotFearutePerAllSamples(imFeaturesP, imFeaturesN, i, 0);
	        imstack.addSlice("thresh = " + IJ.d2s(adaboost[i][2], 3), plot.getProcessor());
	    }
	    imp1.setStack("training", imstack);
	    imp1.show();
	    
	}
	
	public void actionPerformed(ActionEvent e) {
		if(e.getSource()==posButton){
			JFileChooser fc = new JFileChooser();
			int returnVal = fc.showOpenDialog(null);
	        if (returnVal != JFileChooser.APPROVE_OPTION) {
	            return;
	        }
	        
	        pos_train_file = fc.getSelectedFile().getAbsolutePath();
	        train_name = fc.getSelectedFile().getName();
	        IJ.showStatus("Opened " + pos_train_file);
		}
		if(e.getSource()==negButton){
			JFileChooser fc = new JFileChooser();
			int returnVal = fc.showOpenDialog(null);
	        if (returnVal != JFileChooser.APPROVE_OPTION) {
	            return;
	        }
	        neg_train_file = fc.getSelectedFile().getAbsolutePath();
	        
	        IJ.showStatus("Opened " + neg_train_file);
		}
	}
	
	private int[][] getIntegralImage(double[][] imageZ) {
	    int dimsx = imageZ[0].length;
	    int dimsy = imageZ.length;
	    int[][] imageInt = new int[dimsy][dimsx];
	    int[][] s = new int[dimsy][dimsx];
	    for (int j = 0; j < dimsy; j++) {
	        for (int i = 0; i < dimsx; i++) {
	            s[j][i] = (j - 1 >= 0) ? s[j - 1][i] + (int) imageZ[j][i] : (int) imageZ[j][i];
	            imageInt[j][i] = (i - 1 >= 0) ? imageInt[j][i - 1] + s[j][i] : s[j][i];
	        }
	    }
	    return imageInt;
	}

	private int[] computeFeatures(int[][] imageInt, FeaturePool featurePool) {
		
	    int fSize = featurePool.getSize();
	    int[] imFeatures = new int[fSize];
	    for (int i = 0; i < fSize; i++) {
	        int[] r0 = featurePool.getFeature(i).getWhiteSquare();
	        int[] r1 = featurePool.getFeature(i).getBlackSquare();
	        double[] w = featurePool.getFeature(i).getWeights();

	        int i11 = imageInt[r1[1] + r1[3] - 1][r1[0] + r1[2] - 1];
	        int i00 = ((r1[0] - 1 < 0) || (r1[1] - 1 < 0)) ? 0 : imageInt[r1[1] - 1][r1[0] - 1];
	        int i01 = ((r1[0] - 1 < 0)) ? 0 : imageInt[r1[1] + r1[3] - 1][r1[0] - 1];
	        int i10 = ((r1[1] - 1 < 0)) ? 0 : imageInt[r1[1] - 1][r1[0] + r1[2] - 1];
	        int sum1 = i11 + i00 - i10 - i01;

	        i11 = imageInt[r0[1] + r0[3] - 1][r0[0] + r0[2] - 1];
	        i00 = ((r0[0] - 1 < 0) || (r0[1] - 1 < 0)) ? 0 : imageInt[r0[1] - 1][r0[0] - 1];
	        i01 = ((r0[0] - 1 < 0)) ? 0 : imageInt[r0[1] + r0[3] - 1][r0[0] - 1];
	        i10 = ((r0[1] - 1 < 0)) ? 0 : imageInt[r0[1] - 1][r0[0] + r0[2] - 1];
	        int sum0 = i11 + i00 - i10 - i01;

	        imFeatures[i] = (int) (w[0] * sum0 + w[1] * sum1);
	    }

	    return imFeatures;
	}
	
	private double[][] trueAdaBoost(int[][] imFeaturesP, int[][] imFeaturesN, int T) {
	    int sizeP = imFeaturesP.length;
	    int sizeN = imFeaturesN.length;
	    int fSize = imFeaturesN[0].length;
	    int m = sizeP + sizeN;
	    double[] w = new double[m];
	    double[][] adaboost = new double[T][3]; // [id, alpha, threshold]

	    //initial weights
	    for (int i = 0; i < sizeP; i++) {
	        w[i] = 1.0/sizeP;
	    }
	    for (int i = sizeP; i < m; i++) {
			w[i] = 1.0/sizeN;
		}
	    
	    int tt = T;
	    for (int t = 0; t < T; t++) {
	        double[][] thresh = weightedWeakClassification(imFeaturesP, imFeaturesN, w, search_step); // [feature_i][optimal threshold and the score /error]
	        
	        // thresh fSize x 2 (0~theta, 1~score)
	        
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

	        //update weights (multiply with beta if correct)
	        for (int i = 0; i < sizeP; i++) {
	            if (applyClassifier(imFeaturesP[i][bestClassifier], thresh[bestClassifier][0])) {
	                w[i] *= beta;
	            }
	        }
	        for (int i = 0; i < sizeN; i++) {
	            if (!applyClassifier(imFeaturesN[i][bestClassifier], thresh[bestClassifier][0])) {
	                w[i + sizeP] *= beta;
	            }
	        }
	        
	        adaboost[t][1] = Math.log(1 / beta);
	        adaboost[t][2] = thresh[bestClassifier][0];

	        //normalize weights
	        double sum = 0;
	        for (int i = 0; i < m; i++) {
	            sum += w[i];
	        }
	        for (int i = 0; i < m; i++) {
	            w[i] /= sum;
	        }
	        // number of runs - id of the feature - alpha - threshold - error 
	        IJ.log("t = " + (t + 1) + " -> best classifier:" + IJ.d2s(adaboost[t][0] + 1, 0) + " Math.log(1 / beta):  " +
	                IJ.d2s(adaboost[t][1], 4) + "  optimal threshold: " + IJ.d2s(adaboost[t][2], 4) + "  error:   " +
	                IJ.d2s(mineps, 6));
	    }
	    if (tt != T) {
	        double[][] new_adaboost = new double[tt + 1][3];
	        for (int i = 0; i < tt + 1; i++) {
	            new_adaboost[i] = adaboost[i];
	        }
	        return new_adaboost;
	    } else {
	        return adaboost;
	    }
	}
	
	private boolean applyClassifier(double x, double thresh) {
	    return (x >= thresh) ? true : false;
	}
	
	private double[][] weightedWeakClassification(int[][] imFeaturesP, int[][] imFeaturesN, double[] w, int dt) {
	    
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
	        double[][] count = new double[dt][2];
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
	            count[j][1] = score;
	            count[j][0] = thr;
	        }
	        // find optimal threshold
	        thresh[k][1] = count[0][1];
	        for (int j = 1; j < dt - 1; j++) {
	            if (thresh[k][1] > count[j][1]) {
	                thresh[k][1] = count[j][1];
	                thresh[k][0] = count[j][0];
	            }
	        }
	    }
//	    // show all features
//	    ImagePlus imp1 = new ImagePlus();
//	    ImageStack imstack = new ImageStack(528, 255);
//	    //imp1.setDimensions(1, 1, adaboost.length);
//	    for (int i = 0; i < imFeaturesN[0].length; i++) {
//	        //Plot plot = plotFearutePerAllSamples(imFeaturesP, imFeaturesN, (int)adaboost[i][0], adaboost[i][2]);
//	        Plot plot = plotFearutePerAllSamples(imFeaturesP, imFeaturesN, i, thresh[i][0]);
//	        imstack.addSlice("thresh = " + thresh[i][0], plot.getProcessor());
//	    }
//	    imp1.setStack("test", imstack);
//	    imp1.show();

	    return thresh;
	}
	
	private Plot plotFeaturePerAllSamples(int[][] imFeaturesP, int[][] imFeaturesN, int c, double thresh) {
//	    int size = imFeaturesP[0].length;
	    int sizep = imFeaturesP.length;
	    int sizen = imFeaturesN.length;
	    double[] xp = new double[sizep];
	    double[] xn = new double[sizen];
	    double[] yp = new double[sizep];
	    double[] yn = new double[sizen];
	    int maxyp = imFeaturesP[0][c];
	    int minyp = imFeaturesP[0][c];
	    for (int i = 0; i < sizep; i++) {
	        xp[i] = i;
	        yp[i] = imFeaturesP[i][c];
	        maxyp = (maxyp < imFeaturesP[i][c]) ? imFeaturesP[i][c] : maxyp;
	        minyp = (minyp > imFeaturesP[i][c]) ? imFeaturesP[i][c] : minyp;
	    }
	    int maxyn = imFeaturesN[0][c];
	    int minyn = imFeaturesN[0][c];
	    for (int i = 0; i < sizen; i++) {
	        xn[i] = i;
	        yn[i] = imFeaturesN[i][c];
	        maxyn = (maxyn < imFeaturesN[i][c]) ? imFeaturesN[i][c] : maxyn;
	        minyn = (minyn > imFeaturesN[i][c]) ? imFeaturesN[i][c] : minyn;
	    }

	    Plot plot = new Plot("Average velocity = " + c, "Sample id", "Feature", xp, yp);
	    plot.setLimits(0, Math.max(sizen, sizep), Math.min(minyn, minyp), Math.max(maxyn, maxyp));
	    plot.setColor(Color.RED);
	    plot.addPoints(xn, yn, Plot.LINE);
	    plot.setColor(Color.BLACK);
	    double[] linex = new double[]{0, Math.max(sizen, sizep)};
	    double[] liney = new double[]{thresh, thresh};
	    plot.setLineWidth(5);
	    plot.addPoints(linex, liney, Plot.LINE);
	    plot.setColor(Color.BLUE);
	    plot.setLineWidth(2);
	    //plot.show();
	    return plot;
	}
	
}