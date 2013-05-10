package advantra.critpoint;

import java.awt.Color;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.File;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.Vector;

//import advantra.feature.CircularFilterSet;
import advantra.feature.FilterSet;
import advantra.feature.GaborFilt2D;
import advantra.file.AnalyzeCSV;

import advantra.general.Sort;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.*;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import ij.plugin.ZProjector;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;

public class ExtractFeatures implements PlugIn, MouseListener {

	ImagePlus 	img, gab, gabAll, weighted, neuriteness, Vx, Vy;
	ImagePlus 	extract_from;
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
	double 		radius; // calculated wrt. highest scale but can be fixed
	double 		dr, darc, rratio;
	int			surr;
	
	int 		W, H;
	
	int			nr_proc;

    int         T;
	
	// output will be examples (pos. and neg.)
	ImageStack pos_examples_profile;
	ImageStack neg_examples_profile;

	float[][] featsP;
	float[][] featsN;

	FilterSet cft;
	float[] radiuses;

	public void run(String arg0)
	{

        /*
		 * generate filters to score on example profiles (generate features)
		 */
        int[] angleRes = new int[]{40};
        double[] radRes = new double[]{0.3, 0.6};
        cft = new FilterSet(angleRes, radRes);

        /*
        show them
        cft.showConfigs();
        */

        int nrFilters = cft.circConfs.size();

        IJ.showMessage("Formed filter bank! "+nrFilters+" filters");

		t1  	= Prefs.get("advantra.critpoint.start_scale", 		3.0);
		t2    	= Prefs.get("advantra.critpoint.end_scale", 		7.0);
		tn 		= (int)Prefs.get("advantra.critpoint.nr_scales", 	3);
		M		= (int)Prefs.get("advantra.critpoint.nr_angles", 	8);
		surr  	= (int)Prefs.get("advantra.critpoint.curr", 		3);

		dr    	= Prefs.get("advantra.critpoint.dr", 				1.0);
		darc    = Prefs.get("advantra.critpoint.darc", 				1.0);
		rratio	= Prefs.get("advantra.critpoint.dratio", 			0.2);
		nr_proc = (int)Prefs.get("advantra.critpoint.nr_proc", 		4);

        T = (int)Prefs.get("advantra.critpoint.T",                  5);

		train_folder = (String)Prefs.get("advantra.critpoint.train_folder",
				(System.getProperty("user.home")+File.separator));
		test_folder = (String)Prefs.get("advantra.critpoint.test_folder",
				(System.getProperty("user.home")+File.separator));

		String[] extract_opts = new String[3];
		extract_opts[0] = "orig.pix.";
		extract_opts[1] = "neuriteness";
		extract_opts[2] = "gabor.res.";

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

		gd.addStringField("train folder : ", train_folder, 	40);
		gd.addStringField("test  folder : ", test_folder, 	40);

		gd.addMessage("parallelization");
		gd.addNumericField("CPU #",					nr_proc, 0);

		gd.addMessage("source");
		gd.addChoice("choose: ", extract_opts, extract_opts[0]);

        gd.addMessage("Ada Boost");
        gd.addNumericField("T",					    T, 0, 6, " (total "+nrFilters+")");

		gd.showDialog();
		if (gd.wasCanceled()) return;

		t1 	= 		    gd.getNextNumber();
		t2	= 		    gd.getNextNumber();
		tn	= (int)	    gd.getNextNumber();
		M 	= (int)	    gd.getNextNumber();

		surr = (int)    gd.getNextNumber();
		dr  = 			gd.getNextNumber();
		darc  = 		gd.getNextNumber();
		rratio  = 		gd.getNextNumber();

		train_folder = 	gd.getNextString();
		test_folder = 	gd.getNextString();

		nr_proc = (int) gd.getNextNumber();

		int choose_source = gd.getNextChoiceIndex();

        T = (int)       gd.getNextNumber();

		Prefs.set("advantra.critpoint.start_scale", t1);
		Prefs.set("advantra.critpoint.end_scale", 	t2);
		Prefs.set("advantra.critpoint.nr_scales", 	tn);
		Prefs.set("advantra.critpoint.nr_angles", 	M);
		Prefs.set("advantra.critpoint.surr", 		surr);
		Prefs.set("advantra.critpoint.dr", 			dr);
		Prefs.set("advantra.critpoint.darc", 		darc);
		Prefs.set("advantra.critpoint.dratio", 		rratio);

		Prefs.set("advantra.critpoint.train_folder", train_folder);
		Prefs.set("advantra.critpoint.test_folder", test_folder);

		Prefs.set("advantra.critpoint.nr_proc", 	nr_proc);
        Prefs.set("advantra.critpoint.T", 	        T);

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

        /*
        all the params are loaded...
         */

		IJ.setTool(7);

		File dir = new File(train_folder);
		train_folder = dir.getAbsolutePath();
		if(!dir.isDirectory() ){
			IJ.error("Wrong directory!");
			return;
		}

		File[] files_tif = listFilesEndingWith(dir, ".tif");
		File[] files_pos = new File[files_tif.length];
		File[] files_neg = new File[files_tif.length];
        File[] files_tst = new File[files_tif.length];

		Vector<double[][]> locs_pos = new Vector<double[][]>(files_tif.length);
		Vector<double[][]> locs_neg = new Vector<double[][]>(files_tif.length);
        Vector<double[][]> locs_tst = new Vector<double[][]>(files_tif.length);

		int total_pos = 0;
		int total_neg = 0;
		int curr_pos = 0;
		int curr_neg = 0;

		/*
		//    set angular resolution
        //    set radiuses
		 */
        int angular_resolution = 128;

        // radiuses (experimental)
        int cnt_rads = 0;
        for (double r = radius; r >= radius*rratio; r-=dr) {
            cnt_rads++;
        }
        radiuses = new float[cnt_rads];
        cnt_rads = 0;
        for (double r = radius; r >= radius*rratio; r-=dr) {
            radiuses[cnt_rads] = (float)r;
            cnt_rads++;
        }

		System.out.println("points: "+profileLen(radiuses, darc));

        pos_examples_profile = new ImageStack(cft.circConfs.size(), 1);
		neg_examples_profile = new ImageStack(cft.circConfs.size(), 1);

		System.out.println("\n## TRAIN ##  "+train_folder);

		for (int i = 0; i < files_tif.length; i++) { // for each tif file

			System.out.print(""+files_tif[i].getName()+" "+i+"/"+(files_tif.length-1)+"  ...  ");

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
				System.out.print("\t "+A.length+" positives ");
			}
			else{
				locs_pos.add(null);
				System.out.print("\t no positives");
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
				System.out.print("\t "+B.length+" negatives ");
			}
			else{
				locs_neg.add(null);
				System.out.print("\t no negatives");
			}

			System.out.println();

			/*
			 * actual feature extraction
			 */

			if(curr_pos>0 || curr_neg>0){

				System.out.println("...extracting profiles...");

				img = new ImagePlus(files_tif[i].getAbsolutePath());
				resetCalibration();
				convertToFloatImage();
				H = img.getHeight();
				W = img.getWidth();

                if(choose_source==2){
					/*
					 * gabor
					 */
					System.out.println("...gabor...");
					gab = extractGabor(img, theta_pi, t, bandwidth, psi, gamma, isReal, nr_proc);
					//gab.show(); // gab.setTitle("gabor filter, directional responses");
					/*
					 * extract maximal directional gabor response in every point
					 */
					gabAll = maxStackZ(gab);
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
                }
                else
				{
	            	System.out.println("...skipping gabor...");
	            	gab = gabAll = weighted = img;
                }

				/*
				 *  extract neuriteness & eigen vecs
				 */

				System.out.println("...neuriteness...");
				Vector<ImagePlus> nness = VizFeatures.extractNeuritenessAndEigenVec(img, s);
				neuriteness = nness.get(0);  neuriteness.setTitle("neuriteness");
				//neuriteness.show();
				Vx = nness.get(1);// Vx.setTitle("Vx");
				Vy = nness.get(2);// Vy.setTitle("Vy");

				// there are 3 options for input images and locations: neuriteness, weighted, img
				//ImagePlus extract_from;

				switch (choose_source) {
				case 0:
					System.out.println("use: " + extract_opts[0]);
					extract_from = img;
					break;
				case 1:
					System.out.println("use: " + extract_opts[1]);
					extract_from = neuriteness;
					break;

				case 2:
					System.out.println("use: " + extract_opts[2]);
					extract_from = weighted;
					break;
				default:
					extract_from = img;
					break;
				}

				Overlay ovly = new Overlay();

				for (int k = 0; k < curr_pos; k++) {

					int atX = (int)locs_pos.get(i)[k][0];
					int atY = (int)locs_pos.get(i)[k][1];

					float[] ftPos = extractFeats(img, Vx, Vy, cft, atX, atY, radiuses, darc, false);
					FloatProcessor fp = new FloatProcessor(ftPos.length, 1, ftPos);
					pos_examples_profile.addSlice(fp);

					PointRoi pt = new PointRoi(atX-0.5, atY-0.5);
					pt.setStrokeColor(Color.RED);
					ovly.addElement(pt);

				}

				for (int k = 0; k < curr_neg; k++) {

					int atX = (int)locs_neg.get(i)[k][0];
					int atY = (int)locs_neg.get(i)[k][1];

					float[] ftNeg = extractFeats(img, Vx, Vy, cft, atX, atY, radiuses, darc, false);
					FloatProcessor fp = new FloatProcessor(ftNeg.length, 1, ftNeg);
                    neg_examples_profile.addSlice(fp);

					PointRoi pt = new PointRoi(atX-0.5, atY-0.5);
					pt.setHideLabels(false);
					pt.setStrokeColor(Color.BLUE);
					pt.setName("negative" + k);

					ovly.addElement(pt);

				}

                /*
                show loaded train image with markers
                 */
				ImagePlus showIt = new ImagePlus("VIZ_"+files_tif[i].getName(), img.getProcessor());
				showIt.setOverlay(ovly);
				showIt.show();
//
			} // if there were some

		} // loop files


		System.out.println("////\n total (+) : "+total_pos);
		System.out.println(" total (-) : "+total_neg+"\n////\n");

		/*
		 * show examples
		 */
//		ImagePlus pos_examples_image =  new ImagePlus("positive examples", pos_examples_profile);
//		pos_examples_image.show();
//		ImagePlus neg_examples_image =  new ImagePlus("negative examples", neg_examples_profile);
//		neg_examples_image.show();

		/*
		 * calculate features
		 */

		featsP = new float[pos_examples_profile.getSize()][nrFilters];
		featsN = new float[neg_examples_profile.getSize()][nrFilters];

		for (int g = 0; g < pos_examples_profile.getSize(); g++)
			featsP[g] = (float[]) pos_examples_profile.getProcessor(g + 1).getPixels();

		for (int g = 0; g < featsN.length; g++)
			featsN[g] = (float[]) neg_examples_profile.getProcessor(g + 1).getPixels();

		System.out.println(" done extracting train features (+ and -).");

		/*
		 *  train
		 */
		System.out.print("training AdaBoost..."+T);
		double[][] adaboost = trueAdaBoost(featsP, featsN, T);
		System.out.println("done.");

		for (int i = 0; i < adaboost.length; i++) {
			for (int j = 0; j < adaboost[0].length; j++) System.out.print("\t" + IJ.d2s(adaboost[i][j], 2));
			System.out.println();
		}
		System.out.println();

		// show the best features
		//for (int i = 0; i < adaboost.length; i++) cft.showConfigs((int) adaboost[i][0]);

	    /*
	     * test folder
	     */

	    File dir_test = new File(test_folder);
	    test_folder = dir_test.getAbsolutePath();
		if(!dir_test.isDirectory() ){
			IJ.error("Wrong directory!");
			return;
		}

		extract_from.show();
		extract_from.getWindow().getCanvas().addMouseListener(this);

        System.out.println("## TEST ##  "+test_folder);

        int curr_tst = 0;
        int total_tst = 0;

		boolean doTest = true;

		File[] test_files_tif = listFilesEndingWith(dir_test, ".tif");
		for (int i = 0; i < test_files_tif.length && doTest; i++) { // for each tif file

            System.out.print(""+test_files_tif[i].getName()+" "+i+"/"+(test_files_tif.length-1)+"  ...  ");

            curr_tst = 0;
            File[] check;
            String suffix;
            ImagePlus readMask;
            String file_name = test_files_tif[i].getName();
            file_name = file_name.substring(0, file_name.length()-4);
            suffix = file_name+".mask";
//            System.out.println("suffix: "+suffix);

            check = listFilesEndingWith(dir_test, suffix);
            if(check!=null && check.length>0){
                files_tst[i] = check[0];
                readMask = new ImagePlus(files_tst[i].getAbsolutePath());
                //readMask.show();
                // extract locations with logical 1
                double[][] C = extractLocations((ByteProcessor) readMask.getProcessor());
                System.out.print(
                        files_tst[i].getAbsolutePath()+
                        " -> capturing "+(100f*(float)C.length/(readMask.getWidth()*readMask.getHeight()))+
                        "% of total pix."
                );
                locs_tst.add(C);
                curr_tst = C.length;
                total_tst += C.length;

            }
            else{
                locs_tst.add(null);
                System.out.print("no mask found!");
            }

            /*
			 * actual feature extraction
			 */

            if(curr_tst>0){
//
                System.out.println("...extracting profiles...");

                img = new ImagePlus(test_files_tif[i].getAbsolutePath());
                resetCalibration();
                convertToFloatImage();
//                IJ.freeMemory();
                img.show();

                H = img.getHeight();
                W = img.getWidth();
//
                if(choose_source==2){
//
                    System.out.println("...gabor...");
					gab = extractGabor(img, theta_pi, t, bandwidth, psi, gamma, isReal, nr_proc);
                  	gabAll = maxStackZ(gab);
					boolean normalize = false;
					if(normalize){
						weighted = new ImagePlus("gabor filter, directional responses, max, min-max normalized",
														VizFeatures.normalizeStackMinMax(gabAll.getStack()));
					}
					else{
						weighted = new ImagePlus("gabor filter, directional responses, max, min-max normalized",
														gabAll.getStack());
					}
////                    weighted.show();
////                    gabAll.show();
////                    gab.show();
                }
                else{
                    System.out.println("...skipping gabor...");
                    gab = gabAll = weighted = img;
                }

                /*
				 *  extract neuriteness & eigen vecs
				 */

                System.out.println("...neuriteness...");

				Vector<ImagePlus> nness = VizFeatures.extractNeuritenessAndEigenVec(img, s);
				neuriteness = nness.get(0);  neuriteness.setTitle("neuriteness");
				//neuriteness.show();
				Vx = nness.get(1);// Vx.setTitle("Vx");
				Vy = nness.get(2);// Vy.setTitle("Vy");

				// there are 3 options for input images and locations: neuriteness, weighted, img
				switch (choose_source) {
					case 0:
						System.out.println("use test: " + extract_opts[0]);
						extract_from = img;
						break;
					case 1:
						System.out.println("use test: " + extract_opts[1]);
						extract_from = neuriteness;
						break;

					case 2:
						System.out.println("use test: " + extract_opts[2]);
						extract_from = weighted;
						break;
					default:
						extract_from = img;
						break;
				}

                Overlay ovly = new Overlay();
                System.out.print("\nclassifying...");

				for (int k = 0; k < curr_tst; k++) {
//
                    int atX = (int)locs_tst.get(i)[k][0];
                    int atY = (int)locs_tst.get(i)[k][1];

//					long t00 = System.currentTimeMillis();
//					long t01 = System.currentTimeMillis();
//					System.out.println("elapsed: "+((t01-t00)/1000f));

					float[] ftTst = extractFeats(img, Vx, Vy, cft, atX, atY, radiuses, darc, false);
                    int res = applyAdaBoost(adaboost, ftTst);

                    if(res==1){
                        PointRoi pt = new PointRoi(atX+0.5, atY+0.5);
                        pt.setStrokeColor(Color.YELLOW);
                        ovly.addElement(pt);
                    }
//                    else{
//                        PointRoi pt = new PointRoi(atX, atY);
//                        //pt.setStrokeWidth(3.0);
//                        pt.setStrokeColor(Color.GREEN);
//                        ovly.addElement(pt);
//                    }

                }

                System.out.println("done.");

                /*
                show loaded test image
                 */
				ImagePlus showIt = new ImagePlus("VIZ_"+test_files_tif[i].getName(), img.getProcessor());
				showIt.setOverlay(ovly);
				showIt.show();

            }

        System.out.println();

		}
//
////        /*
////		 * show test examples
////		 */
////
////        ImagePlus tst_examples_image =  new ImagePlus("test examples", tst_examples_profile);
////        tst_examples_image.show();
////
////        for (int i = 0; i<6; i++){
////            tst_examples_image.getCanvas().zoomIn(0, 0);
////        }
//
////		/*
////		 * generate filters to score on example profiles (calculate features)
////		 */
////
////        int[] angleRes = new int[]{60};
////        CircularFilterSet cft = new CircularFilterSet(angleRes);
////        int nrFilters = cft.filts.size();
//
//		/*
//		 * calculate features (filters are generated)
//		 */
//
////        System.out.print("calculate features...");
//
////        featsT = new double[tst_examples_profile.getSize()][nrFilters];
//
//
////        long t20 = System.currentTimeMillis();
////        for (int exampleIdx = 0; exampleIdx < tst_examples_profile.getSize(); exampleIdx++) {//
////            System.out.print(".");
////            featsT[exampleIdx] = cft.calculateScore((float[])tst_examples_profile.getProcessor(exampleIdx+1).getPixels());
////        }

		//if(true) { System.out.println("closing..."); return;}

	}

	private int profileLen(float[] radiuses, double darc){
		int len = 0;
		for (int rI = 0; rI < radiuses.length; rI++) {
			for (double arc = 0; arc < 2*radiuses[rI]*Math.PI; arc+=darc) {
				len++;
			}
		}
		return len;
	}

	public float[] extractFeats(
				ImagePlus img,
				ImagePlus Vx,
				ImagePlus Vy,
				FilterSet filterSet,
				int atX,
				int atY,
				float[] radiuses,
				double darc,
				boolean plotIt
	)
	{

		Overlay ol = new Overlay();

		int len = profileLen(radiuses, darc);

		float[] val = new float[len];
		float[] ang = new float[len];

		len = 0;
		for (int rI = 0; rI < radiuses.length; rI++) {
			for (double arc = 0; arc < 2 * radiuses[rI] * Math.PI; arc += darc) {

				ang[len] = (float) (arc / radiuses[rI]);

				double d1 = Math.sin(ang[len]);
				double d2 = -Math.cos(ang[len]);

				double x2 = atX + radiuses[rI] * d1;
				double y2 = atY + radiuses[rI] * d2;

				int x_loc = (int) Math.round(x2);
				int y_loc = (int) Math.round(y2);

				float v1 = Vx.getProcessor().getPixelValue(x_loc, y_loc);
				float v2 = Vy.getProcessor().getPixelValue(x_loc, y_loc);

				float mult = img.getProcessor().getPixelValue(x_loc, y_loc);
				val[len] = (float) Math.abs(mult * v1 * d1 + mult * v2 * d2);

				if(plotIt){
					ol.add(new PointRoi(x_loc-0.5, y_loc-0.5));
					ol.add(new Line(x_loc-0.5, y_loc-0.5, x_loc-0.5+1*v1, y_loc-0.5+1*v2));
				}

				len++;


			}
		}

		if(plotIt){

			Plot p = new Plot("", "", "");
			p.setLimits(0, 2*Math.PI, 0, 1);//Sort.findMax(val)
			p.addPoints(ang, val, Plot.CIRCLE);
			p.show();

			img.setOverlay(ol);

		}

		// allocate out
		float[] outFeat = new float[filterSet.circConfs.size()];
		for (int fI = 0; fI < filterSet.circConfs.size(); fI++){
			outFeat[fI] = filterSet.circConfs.get(fI).calculateScore(val, ang);
		}

		return outFeat;

	}

	private double[][] trueAdaBoost(float[][] featsP, float[][] featsN, int T){
		
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
	        		"t = " + (t + 1) + " -> best classifier idx:" + 
	        		IJ.d2s(adaboost[t][0], 0) + "("+")"+ // nameFeatures[(int)adaboost[t][0]] +
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
		
	private double[][] weightedWeakClassification(float[][] imFeaturesP, float[][] imFeaturesN, double[] w, int dt)
	{// featNr x [thr, score]
		
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
	
	private boolean applyClassifier(float x, double thresh)
	{
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
	
	private int applyAdaBoost(double[][] adaboost, float[] imFeaturesT) {
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

    private double[] computeSelectedFeatures(
													//double[] profile, CircularFilterSet cft,
													double[][] adaboost
	) {

//        double[] imFeatures = new double[adaboost.length];

        int[] chosen_filts = new int[adaboost.length];

        for (int i = 0; i < adaboost.length; i++) {

            chosen_filts[i] = (int) adaboost[i][0];

        }

//        return cft.calculateScore(profile, chosen_filts);
		return new double[5];
    }

    private double[][] extractLocations(ByteProcessor binary_image){

        int nrLocs = 0;
        int h = binary_image.getHeight();
        int w = binary_image.getWidth();

        for (int y = 0; y < h; y++){
            for (int x = 0; x < w; x++){

                if(binary_image.get(x, y)==255){
                    nrLocs++;
                }

            }
        }

        double[][] locations = new double[nrLocs][2];

        int cnt = 0;
        for (int y = 0; y < h; y++){
            for (int x = 0; x < w; x++){

                if(binary_image.get(x, y)==255){
                    locations[cnt][0] = x;
                    locations[cnt][1] = y;
                    cnt++;
                }

            }
        }

        return locations;

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

	private ImagePlus maxStackZ(ImagePlus inim){
				ZProjector zmax = new ZProjector();
				zmax.setImage(inim);
				zmax.setStartSlice(1);
				zmax.setStopSlice(inim.getStackSize());
				zmax.setMethod(ZProjector.MAX_METHOD);
				zmax.doProjection();
				return new ImagePlus("maxZ", zmax.getProjection().getChannelProcessor());
	}

	private ImagePlus extractGabor(ImagePlus img, double[] theta_pi, double[] t, double bandwidth, double psi, double gamma, boolean isReal, int nr_proc){
		int N 				= theta_pi.length; // parallelize
		GaborFilt2D.load(img, theta_pi, t, new double[t.length], bandwidth, psi, gamma, isReal);
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
		return new ImagePlus("extractGabor", GaborFilt2D.gabor_directional_responses);

	}


	public void mouseClicked(MouseEvent e) {

		int atX = 	extract_from.getWindow().getCanvas().offScreenX(e.getX());
		int atY = 	extract_from.getWindow().getCanvas().offScreenY(e.getY());
		float[] ftPos = extractFeats(extract_from, Vx, Vy, cft, atX, atY, radiuses, darc, true);

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