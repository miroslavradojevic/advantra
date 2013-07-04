import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.PointRoi;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;

/*
import weka.classifiers.trees.RandomForest;
import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
*/

import java.awt.*;
import java.io.File;
import java.util.ArrayList;
import java.util.Random;
import java.util.Vector;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 6/23/13
 * Time: 5:38 PM
 */
public class BifDetectML implements PlugIn {

	int neuronDiamMin, neuronDiamMax;
	double scaleMin, scaleMax, scaleStep;
	String train_folder, test_folder;

    int         H, W;

	public void run(String s) {

		neuronDiamMin     	= (int) Prefs.get("advantra.critpoint.neuron_diam_min", 2);
		neuronDiamMax     	= (int) Prefs.get("advantra.critpoint.neuron_diam_max", 4);

		scaleMin     	= Prefs.get("advantra.critpoint.scale_min", 2);
		scaleMax     	= Prefs.get("advantra.critpoint.scale_max", 4);
        scaleStep       = Prefs.get("advantra.critpoint.scale_step", 0.5);

		train_folder    = (String)Prefs.get("advantra.critpoint.train_folder", (System.getProperty("user.home")+ File.separator));
		test_folder     = (String)Prefs.get("advantra.critpoint.test_folder", (System.getProperty("user.home")+File.separator));

		GenericDialog gd = new GenericDialog("Bif. Detection");
		gd.addNumericField("n. diam. (min)", neuronDiamMin, 0, 5, "");
		gd.addNumericField("n. diam. (max)", neuronDiamMax, 0, 5, "");

        gd.addNumericField("scale (min)", scaleMin, 0, 5, "");
        gd.addNumericField("scale (max)", scaleMax, 0, 5, "");

		gd.addStringField("train folder : ", train_folder, 	40);
		gd.addStringField("test  folder : ", test_folder, 	40);
		gd.addCheckbox("images got same dimensions", false);
		gd.addCheckbox("equal # (+) and (-)", false);

		gd.showDialog();
		if (gd.wasCanceled()) return;

		neuronDiamMin =  (int)gd.getNextNumber();
        neuronDiamMax =  (int)gd.getNextNumber();
		Prefs.set("advantra.critpoint.neuron_diam_min", 	neuronDiamMin);
        Prefs.set("advantra.critpoint.neuron_diam_max", 	neuronDiamMax);

        scaleMin = gd.getNextNumber();
        scaleMax = gd.getNextNumber();
        Prefs.set("advantra.critpoint.scale_min", 	scaleMin);
        Prefs.set("advantra.critpoint.scale_max", 	scaleMax);

		train_folder = 	gd.getNextString();
		test_folder = 	gd.getNextString();
		if(! new File(train_folder).exists())	{IJ.showMessage("train folder does not exist!"); 	return;}
		if(! new File(test_folder).exists()) 	{IJ.showMessage("test folder does not exist!"); 	return;}
		Prefs.set("advantra.critpoint.train_folder",    train_folder);
		Prefs.set("advantra.critpoint.test_folder",     test_folder);

		boolean sameSize = gd.getNextBoolean();
		boolean equal = gd.getNextBoolean();

        ArrayList<Feat> feats = new ArrayList<Feat>();

        int patchRadiusMax = Integer.MIN_VALUE;
        for (int nd = neuronDiamMin; nd<=neuronDiamMax; nd++) {
            for (double scl = scaleMin; scl<=scaleMax; scl+=1) {

                Feat ft = new Feat(nd, scl);

                if (ft.r>patchRadiusMax) patchRadiusMax=ft.r;

                feats.add(ft);

            }
        }

        double[] angles = new double[]{0, 2*Math.PI/3, 4*Math.PI/3};
        for (int i=0; i<feats.size(); i++) {
            new ImagePlus("n.diam="+feats.get(i).diam+",r="+feats.get(i).r, feats.get(i).exportTemplate(angles)).show();
        }

        File dir = new File(train_folder);
        train_folder = dir.getAbsolutePath();
        if(!dir.isDirectory() ){
            IJ.error("Wrong trainset directory: " + train_folder + "   closing...");
            return;
        }

        File[] files_tif = Tools.listFilesEndingWith(dir, ".tif");
        File[] files_pos = new File[files_tif.length];
        File[] files_neg = new File[files_tif.length];

        Vector<double[][]> locs_pos = new Vector<double[][]>(files_tif.length);
        Vector<double[][]> locs_neg = new Vector<double[][]>(files_tif.length);
        Vector<double[][]> locs_tst = new Vector<double[][]>(files_tif.length);

        int curr_pos = 0;
        int curr_neg = 0;

        // to store features
        int nrFilters = feats.size()*9;
        ImageStack  pos_ft = new ImageStack(nrFilters, 1);
        ImageStack  neg_ft = new ImageStack(nrFilters, 1);

        System.out.println("\n## TRAIN ##  \n"+train_folder+"\n------------------------------------\n");

        ImageStack showStk  = new ImageStack();
        boolean initialized = false;
        Overlay ovly        = new Overlay();

        ImagePlus   img;


        for (int i = 0; i < files_tif.length; i++) { // for each tif file used to train

            System.out.print("processing "+files_tif[i].getName()+" "+i+"/"+(files_tif.length-1)+"  ...  ");

            img = Tools.convertToFloatImage(new ImagePlus(files_tif[i].getAbsolutePath()));
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

            check = Tools.listFilesEndingWith(dir, suffix);
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

            check = Tools.listFilesEndingWith(dir, suffix);
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
                                    take_col>patchRadiusMax &&
                                    take_col<img.getWidth()-patchRadiusMax &&
                                    take_row>patchRadiusMax &&
                                    take_row<img.getHeight()-patchRadiusMax;
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
            else {
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

                    for (int l = 0; l < feats.size(); l++) {

                        // fill in takeScores

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

                        for (int l = 0; l < feats.size(); l++) {

                            // fill in takeScores score

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
        float[][]   featsP;
        float[][]   featsN;
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


//        double[][]  adaboost;
//        int T = 5;
//        T = (int)Prefs.get("advantra.critpoint.T", 5);
//        GenericDialog gd1 = new GenericDialog("AdaBoost training");
//        gd1.addMessage("How many feats?");
//        gd1.addNumericField("T :", T, 0, 6, " (out of " + nrFilters + ")");
//        gd1.showDialog();
//        if (gd1.wasCanceled()) return;
//        T = (int)       gd1.getNextNumber();
//        Prefs.set("advantra.critpoint.T", T);
//
//        System.out.print("training AdaBoost...");
//        adaboost = trueAdaBoost(featsP, featsN, T);
//        System.out.println("done.");
//
//        int[] best_indexes = new int[adaboost.length];
//        for (int i = 0 ; i < adaboost.length; i++){
//            best_indexes[i] = (int) adaboost[i][0];
//        }
//
//        // show training results
//        if (adaboost.length>1){
//
//            ImagePlus imp1 = new ImagePlus();
//            ImageStack imstack = new ImageStack(528, 255);
//            imp1.setDimensions(1, 1, adaboost.length);
//            for (int i = 0; i < T; i++) {
//                Plot plot = plotFearutePerAllSamples(featsP, featsN, (int)adaboost[i][0], adaboost[i][2]);
//                imstack.addSlice("thresh = " + adaboost[i][2], plot.getProcessor());
//            }
//            imp1.setStack("best.scores", imstack);
//            imp1.show();
//
//        }

        /*
        RF
         */

/*
        System.out.print("training RF...");
        // form weka instances object
        int totalFeats = (featsP[0].length == featsN[0].length)? featsP[0].length : 0;
        FastVector attList = new FastVector();
        for (int i = 0; i < totalFeats; i++) {
            attList.addElement(new Attribute("feat."+i));
        }
        FastVector categ = new FastVector();
        categ.addElement("yes");
        categ.addElement("no");
        Attribute class_attribute = new Attribute("class", categ);
        attList.addElement(class_attribute);

        Instances trn = new Instances("Trainset", attList, (featsP.length + featsN.length));

        double[] values;
        for (int i = 0; i < featsP.length; i++){
            values = new double[featsP[0].length+1];
            for (int j =0; j < featsP[0].length; j++){
                values[j] = featsP[i][j];
            }
            values[featsP[0].length] = class_attribute.indexOfValue("yes");
            trn.add(new Instance(1.0, values));
        }

        for (int i = 0; i < featsN.length; i++){
            values = new double[featsN[0].length+1];
            for (int j =0; j < featsN[0].length; j++){
                values[j] = featsN[i][j];
            }
            values[featsN[0].length] = class_attribute.indexOfValue("no");
            trn.add(new Instance(1.0, values));
        }

        trn.setClassIndex(trn.numAttributes()-1);
        System.out.println("weka trainset formed " + trn.numInstances() + ", atts " + trn.numAttributes());

        int nrTrees = 10;
        int nrFeatsToConsider = (int)Math.sqrt(trn.numAttributes());
        GenericDialog gd1 = new GenericDialog("RF training");
        gd1.addNumericField("# of trees", nrTrees, 0, 6, "");
        gd1.addNumericField("# of features to consider", nrFeatsToConsider, 0, 6, " (out of " + trn.numAttributes() + ")");
        gd1.showDialog();
        if (gd1.wasCanceled()) return;
        nrTrees = (int)       gd1.getNextNumber();
        nrFeatsToConsider = (int)       gd1.getNextNumber();

        RandomForest rf = new RandomForest();
        rf.setNumTrees(nrTrees);
        rf.setNumFeatures(nrFeatsToConsider);
        try {

            System.out.println("building RF... " + nrTrees+" trees and "+nrFeatsToConsider+" features considered");
            rf.buildClassifier(trn);

        } catch (Exception e) {

        }

        System.out.println("done.");
        */


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

        File[] test_files_tif = Tools.listFilesEndingWith(dir_test, ".tif");
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

            check = Tools.listFilesEndingWith(dir_test, suffix);
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

                img = Tools.convertToFloatImage(new ImagePlus(test_files_tif[i].getAbsolutePath()));
                H = img.getHeight();
                W = img.getWidth();

                System.out.println("\nclassifying " + curr_tst + " locations from " + img.getTitle() + " ... ");

                float[][] full_profile = new float[curr_tst][];


                /*
                Instances tst = new Instances("Test", attList, curr_tst);
*/


                for (int k = 0; k < curr_tst; k++) {

                    int atX = (int)locs_tst.get(i)[k][0]; // these are read using imagej, no need for fixing bias
                    int atY = (int)locs_tst.get(i)[k][1];

					/*
					calculate score on selected - take care of the index
					 */

                    float[] takeScores = new float[nrFilters];
                    for (int l = 0; l < feats.size(); l++) {
                        // fill in takeScores
                    }

                    // taking values for RF
                    full_profile[k] = new float[nrFilters];
                    for (int w = 0; w < nrFilters; w++)
                        full_profile[k][w] = takeScores[w];

                }

                // store image features to weka test dataset
                for (int k1 = 0; k1 < curr_tst; k1++){
                    double[] val = new double[nrFilters+1];
                    for (int k2 = 0; k2 < nrFilters; k2++) {
                        val[k2] = full_profile[k1][k2];
                    }

                /*
                    val[nrFilters] = trn.classAttribute().indexOfValue("no"); //negative class for all


                    tst.add(new Instance(1.0, val));
                */

                }

                /*
                tst.setClassIndex(tst.numAttributes()-1);
*/

                // classify each
                try {
                    for (int ti = 0; ti < curr_tst; ti++){ // tst.numInstances()

                       /*
                        double clsLabel = rf.classifyInstance(tst.instance(ti));


                        if(clsLabel==trn.classAttribute().indexOfValue("yes")){

                            int atX = (int)locs_tst.get(i)[ti][0];
                            int atY = (int)locs_tst.get(i)[ti][1];
                            PointRoi pt = new PointRoi(atX+0.5, atY+0.5);
                            pt.setStrokeColor(Color.BLUE);

                            if(sameSize) {
                                pt.setPosition(isShow.getSize()+1); // indexed with 1
                                //pt.setPosition(isShow.getSize());  // if there was only one image
                                //pt.setPosition(0, showStk.getSize(), 0);
                                ovlyDetections.addElement(pt);
                            }
                            else {
                                ovlyThisImage.addElement(pt);
                            }

                        }
                        */
                    }
                } catch (Exception e) { }

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
        if (sameSize && isShow!=null){
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

}


