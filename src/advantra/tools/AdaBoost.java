package advantra.tools;

import ij.IJ;

public class AdaBoost {
	
	double[][] 	adaboost;	// [id, alpha, threshold]
	float[][] 	featsPos;
	float[][] 	featsNeg;
	int 		T;
	int 		search_step = 500;
	
	public AdaBoost(float[][] imFeaturesP, float[][] imFeaturesN, int T){
		
		adaboost = new double[T][3]; // [id, alpha, threshold]
		
		featsPos = new float[imFeaturesP.length][imFeaturesP[0].length];
		for (int i = 0; i < imFeaturesP.length; i++) {
			for (int j = 0; j < imFeaturesP[0].length; j++) {
				featsPos[i][j] = imFeaturesP[i][j];
			}
		}
		
		featsNeg = new float[imFeaturesN.length][imFeaturesN[0].length];
		for (int i = 0; i < imFeaturesN.length; i++) {
			for (int j = 0; j < imFeaturesN[0].length; j++) {
				featsNeg[i][j] = imFeaturesN[i][j];
			}
		}
		
		this.T= T;
		
	}

	public double[][] run() {
		
	    int sizeP = featsPos.length;
	    int sizeN = featsNeg.length;
	    int fSize = featsNeg[0].length;
	    int m = sizeP + sizeN;
	    double[] w = new double[m];
//	    double[][] adaboost = new double[T][3]; // [id, alpha, threshold]

	    //initial weights
	    for (int i = 0; i < sizeP; i++) {
	        w[i] = 1.0/sizeP;
	    }
	    for (int i = sizeP; i < m; i++) {
			w[i] = 1.0/sizeN;
		}
	    
	    int tt = T;
	    for (int t = 0; t < T; t++) {
	    	
	        double[][] thresh = weightedWeakClassification(featsPos, featsNeg, w, search_step); 
	        // [feature_i][optimal threshold and the score /error]
	        
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
	            if (applyClassifier(featsPos[i][bestClassifier], thresh[bestClassifier][0])) {
	                w[i] *= beta;
	            }
	        }
	        for (int i = 0; i < sizeN; i++) {
	            if (!applyClassifier(featsNeg[i][bestClassifier], thresh[bestClassifier][0])) {
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
	
	public int apply(float[] featTes){
		
		double res 			= 0;
	    double object_res 	= 0;
	    
	    for (int i = 0; i < adaboost.length; i++) {
	        object_res += adaboost[i][1];
	    }
	    
	    object_res *= 0.5;

	    for (int i = 0; i < adaboost.length; i++) {
	        if (applyClassifier(featTes[i], adaboost[i][2])) {
	            res += adaboost[i][1];
	        }
	    }
	    int test = (res > object_res) ? 1 : 0;

	    return test;

	}
	
	private double[][] weightedWeakClassification(
			float[][] 	imFeaturesP, 
			float[][] 	imFeaturesN, 
			double[] 	w, 
			int 		dt) {
	    
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
	
	private boolean applyClassifier(double x, double thresh) {
	    return (x >= thresh) ? true : false;
	}
	
}
