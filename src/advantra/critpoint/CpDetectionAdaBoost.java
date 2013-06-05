package advantra.critpoint;

import advantra.file.AnalyzeCSV;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.*;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.File;
import java.io.FilenameFilter;
import java.util.Random;
import java.util.Vector;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 5/28/13
 * Time: 11:59 AM
 * To change this template use File | Settings | File Templates.
 */
public class CpDetectionAdaBoost implements PlugIn, MouseListener {

    ImagePlus   img;
    ImagePlus   best_feat_img;
    String 		train_folder, test_folder;

    int 		patchRadius, patchDiameter;
    int 		W, H;
    int         T;

    ImageStack  pos_ft;
    ImageStack  neg_ft;
    float[][]   featsP;
    float[][]   featsN;

    double[][]  adaboost;

    // set of filter-kernels
    CircularConfiguration2 ccf2;
    CircularConfiguration3 ccf3;
    CircularConfiguration1 ccf1;
	SymmConfiguration scf;
	AsymmConfiguration acf;
	CircularConfiguration4 ccf4;

    int[] best_indexes;

    int nrFilters3, nrFilters2, nrFilters1, nrFiltersS, nrFiltersA, nrFilters4, nrFilters;

    static float TwoPi = (float) (2*Math.PI);

    public void run(String s) {

        patchRadius     = (int) Prefs.get("advantra.critpoint.patch_radius", 15);
        train_folder    = (String)Prefs.get("advantra.critpoint.train_folder", (System.getProperty("user.home")+ File.separator));
        test_folder     = (String)Prefs.get("advantra.critpoint.test_folder", (System.getProperty("user.home")+File.separator));

        GenericDialog gd = new GenericDialog("Critpoint Detection");
        gd.addNumericField("patch_radius", patchRadius, 0, 5, "");
        gd.addStringField("train folder : ", train_folder, 	40);
        gd.addStringField("test  folder : ", test_folder, 	40);
        gd.addCheckbox("images got same dimensions", false);
        gd.addCheckbox("equal # (+) and (-)", false);

        gd.addMessage("Feature profiles to use");

        gd.addCheckbox("1x", false);
        gd.addCheckbox("2x", false);
        gd.addCheckbox("3x", false);
        gd.addCheckbox("4x", false);
        gd.addCheckbox("Symm.", false);
        gd.addCheckbox("ASymm.", false);
        gd.showDialog();
        if (gd.wasCanceled()) return;
        patchRadius =  (int)gd.getNextNumber();
        train_folder = 	gd.getNextString();
        test_folder = 	gd.getNextString();

		System.out.println(" train folder "+train_folder);

        // check folders
        if(! new File(train_folder).exists()) {IJ.showMessage("train folder does not exist!"); return;}
        else {
            train_folder += (!train_folder.endsWith(File.separator))? File.separator : "";
        }
        if(! new File(test_folder).exists()) {IJ.showMessage("test folder does not exist!"); return;}
        else {
            test_folder += (!test_folder.endsWith(File.separator))? File.separator : "";
        }
        boolean sameSize = gd.getNextBoolean();
        boolean equal = gd.getNextBoolean();
        boolean enableF1 = gd.getNextBoolean();
        boolean enableF2 = gd.getNextBoolean();
        boolean enableF3 = gd.getNextBoolean();
        boolean enableF4 = gd.getNextBoolean();
        boolean enableSymm = gd.getNextBoolean();
        boolean enableAsymm = gd.getNextBoolean();

        Prefs.set("advantra.critpoint.train_folder",    train_folder);
        Prefs.set("advantra.critpoint.test_folder",     test_folder);
        Prefs.set("advantra.critpoint.patch_radius", 	patchRadius);

		/*
		CREATE FEATS
		 */

        patchDiameter = 2*patchRadius+1;

		ccf3 = new CircularConfiguration3(patchRadius);
        ccf2 = new CircularConfiguration2(patchRadius);
        ccf1 = new CircularConfiguration1(patchRadius);
		scf = new SymmConfiguration(patchRadius);
		acf = new AsymmConfiguration(patchRadius);
		ccf4 = new CircularConfiguration4(patchRadius);

        nrFilters3  = (enableF3)? ccf3.kernels.size() : 0;
        nrFilters2  = (enableF2)? ccf2.kernels.size() : 0;
		nrFilters1  = (enableF1)? ccf1.kernels.size() : 0;
		nrFiltersS 	= (enableSymm)? scf.kernels.size() : 0;
		nrFiltersA  = (enableAsymm)? acf.kernels.size() : 0;
		nrFilters4  = (enableF4)? ccf4.kernels.size() : 0;
        nrFilters   = nrFilters3 + nrFilters2 + nrFilters1 + nrFiltersS + nrFiltersA + nrFilters4;

		if(nrFilters==0){System.out.println("no filters set"); return;}

		/*
		SHOW FEATS
		 */

        if(nrFilters3>0) {
            ImagePlus feats3 = new ImagePlus("Y.feat."+patchRadius, ccf3.plotKernels());
            feats3.show();
            feats3.getCanvas().zoomIn(0, 0);
            feats3.getCanvas().zoomIn(0, 0);
            feats3.getCanvas().zoomIn(0, 0);
            feats3.getCanvas().zoomIn(0, 0);
        }


		if(nrFilters2>0){
        ImagePlus feats2 = new ImagePlus("V.feat."+patchRadius, ccf2.plotKernels());
        feats2.show();
        feats2.getCanvas().zoomIn(0, 0);
        feats2.getCanvas().zoomIn(0, 0);
        feats2.getCanvas().zoomIn(0, 0);
        feats2.getCanvas().zoomIn(0, 0);
		}

		if(nrFilters1>0){
		ImagePlus feats1 = new ImagePlus("I.feat."+patchRadius, ccf1.plotKernels());
		feats1.show();
		feats1.getCanvas().zoomIn(0, 0);
		feats1.getCanvas().zoomIn(0, 0);
		feats1.getCanvas().zoomIn(0, 0);
		feats1.getCanvas().zoomIn(0, 0);
		}

		if(nrFiltersS>0){
		ImagePlus featsS = new ImagePlus("S.feat."+patchRadius, scf.plotKernels());
		featsS.show();
		featsS.getCanvas().zoomIn(0, 0);
		featsS.getCanvas().zoomIn(0, 0);
		featsS.getCanvas().zoomIn(0, 0);
		featsS.getCanvas().zoomIn(0, 0);
		}

		if(nrFiltersA>0){
		ImagePlus featsA = new ImagePlus("A.feat."+patchRadius, acf.plotKernels());
		featsA.show();
		featsA.getCanvas().zoomIn(0, 0);
		featsA.getCanvas().zoomIn(0, 0);
		featsA.getCanvas().zoomIn(0, 0);
		featsA.getCanvas().zoomIn(0, 0);
		}

        if(nrFilters4>0) {
            ImagePlus feats4 = new ImagePlus("X.feat."+patchRadius, ccf4.plotKernels());
            feats4.show();
            feats4.getCanvas().zoomIn(0, 0);
            feats4.getCanvas().zoomIn(0, 0);
            feats4.getCanvas().zoomIn(0, 0);
            feats4.getCanvas().zoomIn(0, 0);
        }


        /*
        ALLOCATE STORAGE TRAIN
         */

        File dir = new File(train_folder);
        train_folder = dir.getAbsolutePath();
        if(!dir.isDirectory() ){
            IJ.error("Wrong trainset directory: " + train_folder + "   closing...");
            return;
        }

        File[] files_tif = listFilesEndingWith(dir, ".tif");
        File[] files_pos = new File[files_tif.length];
        File[] files_neg = new File[files_tif.length];

        Vector<double[][]> locs_pos = new Vector<double[][]>(files_tif.length);
        Vector<double[][]> locs_neg = new Vector<double[][]>(files_tif.length);
        Vector<double[][]> locs_tst = new Vector<double[][]>(files_tif.length);

        int curr_pos = 0;
        int curr_neg = 0;

        // to store features
        pos_ft = new ImageStack(nrFilters, 1);
        neg_ft = new ImageStack(nrFilters, 1);

        System.out.println("\n## TRAIN ##  \n"+train_folder+"\n------------------------------------\n");

        ImageStack showStk  = new ImageStack();
        boolean initialized = false;
        Overlay ovly        = new Overlay();

        for (int i = 0; i < files_tif.length; i++) { // for each tif file used to train

            System.out.print("processing "+files_tif[i].getName()+" "+i+"/"+(files_tif.length-1)+"  ...  ");

            img = convertToFloatImage(new ImagePlus(files_tif[i].getAbsolutePath()));
            H = img.getHeight();
            W = img.getWidth();

            if(sameSize && !initialized){
                showStk = new ImageStack(W, H);
                initialized = true;
            }

            curr_pos = curr_neg = 0;

            AnalyzeCSV readCSV;
            File[] check;
            String suffix;
            //ImagePlus readMask;
            String file_name = files_tif[i].getName();
            file_name = file_name.substring(0, file_name.length()-4);

            suffix = file_name+".pos";
            //suffix = file_name+".mask.pos";

            check = listFilesEndingWith(dir, suffix);
            if(check!=null && check.length>0){
                files_pos[i] = check[0];

                readCSV = new AnalyzeCSV(files_pos[i].getAbsolutePath());
                double[][] A = readCSV.readLn(2);

                /*
                readMask = new ImagePlus(files_pos[i].getAbsolutePath());
                double[][] A = extractLocations((ByteProcessor) readMask.getProcessor());
                */

                locs_pos.add(A);
                curr_pos = A.length;
//                total_pos += A.length;
                System.out.print("\t "+A.length+" positives ");
            }
            else{
                locs_pos.add(null);
                System.out.print("\t no positives");
            }

            suffix = file_name+".neg";
            //suffix = file_name+".mask.neg";

            check = listFilesEndingWith(dir, suffix);
            if(check!=null && check.length>0){
                files_neg[i] = check[0];

                readCSV = new AnalyzeCSV(files_neg[i].getAbsolutePath());
                double[][] B = readCSV.readLn(2);

                /*
                readMask = new ImagePlus(files_neg[i].getAbsolutePath());
                double[][] B = extractLocations((ByteProcessor) readMask.getProcessor());
                */

                // note: reducing B not to take those that are out of image
                Vector<double[]> Br = new Vector<double[]>();
                for (int k = 0; k < B.length; k++){
                    double take_col, take_row;
                    take_col = B[k][0];
                    take_row = B[k][1];
                    boolean isIn =
                            take_col>patchRadius &&
                                    take_col<img.getWidth()-patchRadius &&
                                    take_row>patchRadius &&
                                    take_row<img.getHeight()-patchRadius;
                    if(isIn)
                        Br.add(B[k]);
                }
                B = new double[Br.size()][2];
                for (int k = 0; k < Br.size(); k++){
                    B[k][0] = Br.get(k)[0];
                    B[k][1] = Br.get(k)[1];
                }

                locs_neg.add(B);
                curr_neg = B.length;
//                total_neg += B.length;
                System.out.print("\t "+B.length+" negatives ");
            }
            else{
                locs_neg.add(null);
                System.out.print("\t no .neg files");
            }

            System.out.println();

            // storing locations and calculating scores at locations
            int bias = 0; // notations were indexed starting from 1, and imageJ starts from 0

            if((curr_pos>0 || curr_neg>0) ){

                Overlay ovlyThisImage = new Overlay();

                for (int k = 0; k < curr_pos; k++) {

                    int atX = (int)locs_pos.get(i)[k][0] + bias;
                    int atY = (int)locs_pos.get(i)[k][1] + bias;

                    /*
                    calculate filter scores at positives
                     */
                    float[] takeScores = new float[nrFilters];
                    for (int l = 0; l < nrFilters3; l++) {
                        takeScores[l] = ccf3.score(atX, atY, l, img.getProcessor());
                    }

                    for (int l = 0; l < nrFilters2; l++) {
                        takeScores[nrFilters3+l] = ccf2.score(atX, atY, l, img.getProcessor());
                    }

					for (int l = 0; l < nrFilters1; l++) {
						takeScores[nrFilters3+nrFilters2+l] = ccf1.score(atX, atY, l, img.getProcessor());
					}

					for (int l = 0; l < nrFiltersS; l++) {
						takeScores[nrFilters3+nrFilters2+nrFilters1+l] = scf.score(atX, atY, l, img.getProcessor());
					}

					for (int l = 0; l < nrFiltersA; l++) {
						takeScores[nrFilters3+nrFilters2+nrFilters1+nrFiltersS+l] = acf.score(atX, atY, l, img.getProcessor());
					}

					for (int l = 0; l < nrFilters4; l++) {
						takeScores[nrFilters3+nrFilters2+nrFilters1+nrFiltersS+nrFiltersA+l] = ccf4.score(atX, atY, l, img.getProcessor());
					}


					pos_ft.addSlice(new FloatProcessor(nrFilters, 1, takeScores));
                    PointRoi pt = new PointRoi(atX+0.5, atY+0.5);
                    pt.setStrokeColor(Color.RED);

                    if(sameSize) {
                        pt.setPosition(showStk.getSize()+1); // indexed with 1
                        //pt.setPosition(0, showStk.getSize(), 0);
                        ovly.addElement(pt);
                    }
                    else {
                        ovlyThisImage.addElement(pt);
                    }

                }

                // fix: select same amount of negative samples
                // only in case of artificial dataset because of the disbalance pos/neg
                // each image will contribute in as many negatives as positives
                // random without sampling
                boolean[] chs;
                if(curr_neg<=0)
                    chs = new boolean[1];
                else {
                    chs = new boolean[curr_neg];
                    Random rand = new Random();
                    int cntMtches = 0;
                    while (cntMtches < curr_pos) {  // select curr_pos random ones
                        int rIdx = rand.nextInt(curr_neg);
                        if (!chs[rIdx]) {
                            chs[rIdx] = true;
                            cntMtches++;
                        }
                    }
                }
                if (!equal) {// annulate  randomizaiton
                    for (int k = 0; k < chs.length; k++){
                        chs[k] = true;
                    }
                }

                for (int k = 0; k < curr_neg; k++) {
                    if (chs[k]){

                        int atX = (int)locs_neg.get(i)[k][0] + bias;
                        int atY = (int)locs_neg.get(i)[k][1] + bias;

						/*
						 calculate filter scores at negatives input: img.getProcessor(), output: takeScores
						  */
                        float[] takeScores = new float[nrFilters];
                        for (int l = 0; l < nrFilters3; l++) {
                            takeScores[l] = ccf3.score(atX, atY, l, img.getProcessor());
                        }

                        for (int l = 0; l < nrFilters2; l++) {
                            takeScores[nrFilters3+l] = ccf2.score(atX, atY, l, img.getProcessor());
                        }

						for (int l = 0; l < nrFilters1; l++) {
							takeScores[nrFilters3+nrFilters2+l] = ccf1.score(atX, atY, l, img.getProcessor());
						}

						for (int l = 0; l < nrFiltersS; l++) {
							takeScores[nrFilters3+nrFilters2+nrFilters1+l] = scf.score(atX, atY, l, img.getProcessor());
						}

						for (int l = 0; l < nrFiltersA; l++) {
							takeScores[nrFilters3+nrFilters2+nrFilters1+nrFiltersS+l] = acf.score(atX, atY, l, img.getProcessor());
						}

						for (int l = 0; l < nrFilters4; l++) {
							takeScores[nrFilters3+nrFilters2+nrFilters1+nrFiltersS+nrFiltersA+l] = ccf4.score(atX, atY, l, img.getProcessor());
						}

                        neg_ft.addSlice(new FloatProcessor(nrFilters, 1, takeScores));
                        PointRoi pt = new PointRoi(atX+0.5, atY+0.5);
                        pt.setStrokeColor(Color.BLUE);

                        if(sameSize) {
                            pt.setPosition(showStk.getSize()+1); // indexed with 1
                            //pt.setPosition(0, showStk.getSize(), 0);
                            ovly.addElement(pt);
                        }
                        else {
                            ovlyThisImage.addElement(pt);
                        }

                    }

                }

                if(sameSize){
                    showStk.addSlice("SAMPLE"+i, img.getProcessor());
                }
                else {
                    ImagePlus showImgInd = img;
                    showImgInd.setTitle(img.getTitle());
                    System.out.println("overlay has: "+ovly.size());
                    showImgInd.setOverlay(ovlyThisImage);
                    showImgInd.show();
                }

            } // if there were some

        } // loop files

        // visualization if same size
        if(sameSize) {
            ImagePlus showImg = new ImagePlus("VIZ", showStk);
            showImg.setOverlay(ovly);
            showImg.show();
            showImg.getCanvas().zoomIn(0,0);
            showImg.getCanvas().zoomIn(0,0);
        }

        System.out.println("////\n total (+) : " + pos_ft.getSize());
        System.out.println(" total (-) : "       + neg_ft.getSize()+"\n//-------------//\n");

        /*
        STORING EXTRACTED FEATS
         */

        featsP = new float[pos_ft.getSize()][nrFilters];
        featsN = new float[neg_ft.getSize()][nrFilters];

        for (int g = 0; g < pos_ft.getSize(); g++){
            float[] getPix = (float[]) pos_ft.getProcessor(g + 1).getPixels();
            for (int g1 = 0; g1 < getPix.length; g1++){
                featsP[g][g1] = getPix[g1];
            }
        }

        for (int g = 0; g < neg_ft.getSize(); g++) {
            float[] getPix = (float[]) neg_ft.getProcessor(g + 1).getPixels();
            for (int g1 = 0; g1 < getPix.length; g1++){
                featsN[g][g1] = getPix[g1];
            }
        }

        System.out.println(" done extracting train features:\n" +
                "(+) " + featsP.length + " x " + featsP[0].length +
                "\n" +
                "(-) " + featsN.length + " x " + featsN[0].length);

        /*
        ADABOOST
         */
        T = (int)Prefs.get("advantra.critpoint.T", 5);
        GenericDialog gd1 = new GenericDialog("AdaBoost training");
        gd1.addMessage("How many feats?");
        gd1.addNumericField("T :", T, 0, 6, " (out of " + nrFilters + ")");
        gd1.showDialog();
        if (gd1.wasCanceled()) return;
        T = (int)       gd1.getNextNumber();
        Prefs.set("advantra.critpoint.T", T);

        System.out.print("training AdaBoost...");
        adaboost = trueAdaBoost(featsP, featsN, T);
        System.out.println("done.");

        best_indexes = new int[adaboost.length];
        for (int i = 0 ; i < adaboost.length; i++){
            best_indexes[i] = (int) adaboost[i][0];
        }

        // show training results
        if (adaboost.length>1){

            // show the best features
            best_feat_img = new ImagePlus();
            ImageStack best_feat = new ImageStack(patchDiameter, patchDiameter);
            for (int i = 0; i < adaboost.length; i++){
                int id = (int) adaboost[i][0];

				/*
				showing... be careful on index
				 */
				if (id < nrFilters3) {
                    best_feat.addSlice(ccf3.plotKernel(id, 0));
                }
                else if (id >= nrFilters3 && id < nrFilters2+nrFilters3) {
                    best_feat.addSlice(ccf2.plotKernel(id-nrFilters3, 0));
                }
				else if(id >= nrFilters3+nrFilters2 && id < nrFilters1+nrFilters2+nrFilters3) {
					best_feat.addSlice(ccf1.plotKernel(id-nrFilters3-nrFilters2, 0));
				}
				else if (id>= nrFilters1+nrFilters2+nrFilters3 && id < nrFiltersS+nrFilters1+nrFilters2+nrFilters3) {
					best_feat.addSlice(scf.plotKernel(id-nrFilters3-nrFilters2-nrFilters1, 0));
				}
				else if (id>= nrFiltersS+nrFilters1+nrFilters2+nrFilters3 && id < nrFiltersA+nrFiltersS+nrFilters1+nrFilters2+nrFilters3) {
					best_feat.addSlice(acf.plotKernel(id-nrFilters3-nrFilters2-nrFilters1-nrFiltersS, 0));
				}
				else if (id>= nrFiltersA+nrFiltersS+nrFilters1+nrFilters2+nrFilters3 && id < nrFilters4+nrFiltersA+nrFiltersS+nrFilters1+nrFilters2+nrFilters3) {
					best_feat.addSlice(ccf4.plotKernel(id-nrFilters3-nrFilters2-nrFilters1-nrFiltersS-nrFiltersA, 0));
				}

            }
            best_feat_img.setStack("chosen.feats", best_feat);
            best_feat_img.show();

            ImageCanvas all_feats_canvas = best_feat_img.getWindow().getCanvas();
            all_feats_canvas.setName("best");
            all_feats_canvas.addMouseListener(this);

            ImagePlus imp1 = new ImagePlus();
            ImageStack imstack = new ImageStack(528, 255);
            imp1.setDimensions(1, 1, adaboost.length);
            for (int i = 0; i < T; i++) {
                Plot plot = plotFearutePerAllSamples(featsP, featsN, (int)adaboost[i][0], adaboost[i][2]);
                imstack.addSlice("thresh = " + adaboost[i][2], plot.getProcessor());
            }
            imp1.setStack("best.scores", imstack);
            imp1.show();

        }

        /*
        ALLOCATE STORAGE TEST
         */

        File dir_test = new File(test_folder);
        test_folder = dir_test.getAbsolutePath();
        if(!dir_test.isDirectory() ){
            IJ.error("Wrong testset directory "+test_folder);
            return;
        }

        System.out.println("\n## TEST  ##  \n"+test_folder +"\n------------------------------------\n");

        int curr_tst = 0;
        boolean doTest = true;

        File[] test_files_tif = listFilesEndingWith(dir_test, ".tif");
        File[] files_tst = new File[test_files_tif.length];

        ImageStack isShow = new ImageStack();
        if (sameSize) isShow = new ImageStack(W, H);
        Overlay ovlyDetections = new Overlay();

        for (int i = 0; i < test_files_tif.length && doTest; i++) { // for each tif file

            System.out.print(""+test_files_tif[i].getName()+" "+i+"/"+(test_files_tif.length-1)+"  ...  ");

            curr_tst = 0;
            File[] check;
            String suffix;
            ImagePlus readMask;
            String file_name = test_files_tif[i].getName();
            file_name = file_name.substring(0, file_name.length()-4);
            suffix = file_name+".mask";

            check = listFilesEndingWith(dir_test, suffix);
            if(check!=null && check.length>0){
                files_tst[i] = check[0];
                readMask = new ImagePlus(files_tst[i].getAbsolutePath());

                double[][] C = extractLocations((ByteProcessor) readMask.getProcessor()); // extract locations with logical 1
                locs_tst.add(C);
                curr_tst = C.length;
            }
            else{
                locs_tst.add(null);
                System.out.print("no mask found!\n");
            }

            // check if the mask has locations then add it to the list of Overlays
            if(curr_tst>0){

                Overlay ovlyThisImage = new Overlay();

                img = convertToFloatImage(new ImagePlus(test_files_tif[i].getAbsolutePath()));
                H = img.getHeight();
                W = img.getWidth();

                System.out.println("\nclassifying " + curr_tst + " locations from " + img.getTitle() + " ... ");

                for (int k = 0; k < curr_tst; k++) {

                    int atX = (int)locs_tst.get(i)[k][0]; // these are read using imagej, no need for fixing bias
                    int atY = (int)locs_tst.get(i)[k][1];

					/*
					calculate score on selected - take care of index
					 */

                    float[] selScores = new float[best_indexes.length];
                    for (int l = 0; l < best_indexes.length; l++) {

                        int bestIdx = best_indexes[l];

                        if (bestIdx<nrFilters3){
                            selScores[l] = ccf3.score(atX, atY, bestIdx, img.getProcessor());
                        }
                        else if (bestIdx>=nrFilters3 && bestIdx<nrFilters2+nrFilters3) {
                            selScores[l] = ccf2.score(atX, atY, bestIdx-nrFilters3, img.getProcessor());
                        }
						else if (bestIdx>=nrFilters2+nrFilters3 && bestIdx<nrFilters1+nrFilters2+nrFilters3) {
							selScores[l] = ccf1.score(atX, atY, bestIdx-nrFilters3-nrFilters2, img.getProcessor());
						}
						else if (bestIdx>=nrFilters1+nrFilters2+nrFilters3 && bestIdx<nrFiltersS+nrFilters1+nrFilters2+nrFilters3) {
							selScores[l] = scf.score(atX, atY, bestIdx-nrFilters3-nrFilters2-nrFilters1, img.getProcessor());
						}
						else if (bestIdx>=nrFiltersS+nrFilters1+nrFilters2+nrFilters3 && bestIdx<nrFiltersA+nrFiltersS+nrFilters1+nrFilters2+nrFilters3) {
							selScores[l] = acf.score(atX, atY, bestIdx-nrFilters3-nrFilters2-nrFilters1-nrFiltersS, img.getProcessor());
						}
						else if (bestIdx>=nrFiltersA+nrFiltersS+nrFilters1+nrFilters2+nrFilters3 && bestIdx<nrFilters4+nrFiltersA+nrFiltersS+nrFilters1+nrFilters2+nrFilters3) {
							selScores[l] = ccf4.score(atX, atY, bestIdx-nrFilters3-nrFilters2-nrFilters1-nrFiltersS-nrFiltersA, img.getProcessor());
						}

                    }

                    int res = applyAdaBoost(adaboost, selScores);

                    if(res==1){

                        PointRoi pt = new PointRoi(atX+0.5, atY+0.5);
                        pt.setStrokeColor(Color.YELLOW);

                        if(sameSize) {
                            pt.setPosition(isShow.getSize()+1); // indexed with 1
                            //pt.setPosition(0, showStk.getSize(), 0);
                            ovlyDetections.addElement(pt);
                        }
                        else {
                            ovlyThisImage.addElement(pt);
                        }

                    }

                }

                if(sameSize)
                    isShow.addSlice(img.getTitle(), img.getProcessor());
                else {
                    ImagePlus showImgInd = img;
                    showImgInd.setTitle(img.getTitle());
//                    System.out.println(ovlyDetections.size()+ " detections ");
                    showImgInd.setOverlay(ovlyThisImage);
                    showImgInd.show();
                }

            }

        }
        if (sameSize){
            ImagePlus a23 = new ImagePlus("det", isShow);
            //System.out.println("finally "+ovlyDetections.size()+ " detections ");
            a23.setOverlay(ovlyDetections);
            a23.show();
            a23.getCanvas().zoomIn(0,0);
            a23.getCanvas().zoomIn(0,0);
        }

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

                IJ.log(t + " min eps = " + mineps + "   best.cl.idx" +
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

    private boolean applyClassifier(float x, double thresh)
    {
        return (x >= thresh) ? true : false;
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

    private Plot plotFearutePerAllSamples(float[][] imFeaturesP, float[][] imFeaturesN, int c, double thresh) {

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
        plot.setColor(Color.BLUE);

        plot.addPoints(xn, yn, Plot.LINE);
        plot.setColor(Color.BLACK);
        plot.setLineWidth(4);

        double[] linex = new double[]{0, Math.max(sizen, sizep)};
        double[] liney = new double[]{thresh, thresh};
        plot.addPoints(linex, liney, Plot.LINE);
        plot.setColor(Color.RED);
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

    private File[] listFilesEndingWith(
            File dir,
            String suffix
    )
    {
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

    private ImagePlus convertToFloatImage(
            ImagePlus inim
    )
    {

        int W = inim.getWidth();
        int H = inim.getHeight();

        ImageProcessor ip = new FloatProcessor(W, H);
        for (int i = 0; i < H*W; i++ ) {
            ip.setf(i, inim.getProcessor().getPixelValue(i%W, i/W));
        }

        ImagePlus outim = new ImagePlus(inim.getTitle()+"ToFloat", ip);
        return outim;

    }

    public void mouseClicked(MouseEvent e)
    {

        ImageCanvas srcCanv = (ImageCanvas) e.getSource();
        String source  =  srcCanv.getName();
        if (source=="best") {

            int mouseZ = best_feat_img.getCurrentSlice()-1;

            ImagePlus imTest = convertToFloatImage(IJ.openImage());

            System.out.println("caluclate score on "+imTest.getTitle()+"\nwith extracted feature of idx. " + best_indexes[mouseZ]);

            if (best_indexes[mouseZ]<nrFilters3){
                new ImagePlus("max", ccf3.score(best_indexes[mouseZ], imTest.getProcessor())).show();
                new ImagePlus("per.rot", ccf3.scoreAllRot(best_indexes[mouseZ], imTest.getProcessor())).show();
            }
            else if (best_indexes[mouseZ]>=nrFilters3 && best_indexes[mouseZ]<nrFilters3+nrFilters2) {
                new ImagePlus("max", ccf2.score(best_indexes[mouseZ]-nrFilters3, imTest.getProcessor())).show();
                new ImagePlus("per.rot", ccf2.scoreAllRot(best_indexes[mouseZ]-nrFilters3, imTest.getProcessor())).show();
            }
			else if (best_indexes[mouseZ]>=nrFilters3+nrFilters2 && best_indexes[mouseZ]<nrFilters3+nrFilters2+nrFilters1) {
				new ImagePlus("max", ccf1.score(best_indexes[mouseZ]-nrFilters3-nrFilters2, imTest.getProcessor())).show();
				new ImagePlus("per.rot", ccf1.scoreAllRot(best_indexes[mouseZ]-nrFilters3-nrFilters2, imTest.getProcessor())).show();
			}
			else if (best_indexes[mouseZ]>=nrFilters3+nrFilters2+nrFilters1 && best_indexes[mouseZ]<nrFilters3+nrFilters2+nrFilters1+nrFiltersS) {
				new ImagePlus("max", scf.score(best_indexes[mouseZ]-nrFilters3-nrFilters2-nrFilters1, imTest.getProcessor())).show();
				new ImagePlus("per.rot", scf.scoreAllRot(best_indexes[mouseZ]-nrFilters3-nrFilters2-nrFilters1, imTest.getProcessor())).show();
			}
			else if (best_indexes[mouseZ]>=nrFilters3+nrFilters2+nrFilters1+nrFiltersS && best_indexes[mouseZ]<nrFilters3+nrFilters2+nrFilters1+nrFiltersS+nrFiltersA) {
				new ImagePlus("max", acf.score(best_indexes[mouseZ]-nrFilters3-nrFilters2-nrFilters1-nrFiltersS, imTest.getProcessor())).show();
				new ImagePlus("per.rot", acf.scoreAllRot(best_indexes[mouseZ]-nrFilters3-nrFilters2-nrFilters1-nrFiltersS, imTest.getProcessor())).show();
			}

			else if (best_indexes[mouseZ]>=nrFilters3+nrFilters2+nrFilters1+nrFiltersS+nrFiltersA && best_indexes[mouseZ]<nrFilters3+nrFilters2+nrFilters1+nrFiltersS+nrFiltersA+nrFilters4) {
				new ImagePlus("max", ccf4.score(best_indexes[mouseZ]-nrFilters3-nrFilters2-nrFilters1-nrFiltersS-nrFiltersA, imTest.getProcessor())).show();
				new ImagePlus("per.rot", ccf4.scoreAllRot(best_indexes[mouseZ]-nrFilters3-nrFilters2-nrFilters1-nrFiltersS-nrFiltersA, imTest.getProcessor())).show();
			}

        }

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
