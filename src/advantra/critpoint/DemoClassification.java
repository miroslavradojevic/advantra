package advantra.critpoint;

import java.awt.Color;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.File;
import java.io.FilenameFilter;
import java.util.Random;
import java.util.Vector;

import advantra.feature.GaborFilt2D;
import advantra.file.AnalyzeCSV;

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
import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instances;

public class DemoClassification implements PlugIn, MouseListener {

	ImagePlus 	img;
    //, gab, gabAll, weighted, neuriteness, Vx, Vy;
	ImagePlus 	imgSrc;
	String 		train_folder, test_folder;

//    ImagePlus   img_click; // this one just to see how it performs on manually chosen locations (visualisation)

    // feature parameters (filter parameters)
    int         angScale1, angScale2;
    int[]       angScale;

    double      ring1, ring2;
    int         nr_ring;
    double[]    rings;

    double      radial1, radial2;
    int         nr_radial;
    double[]    radials;

	double 		innerRadius;

    int 		patchRadius, patchDiameter;
    int 		W, H;
    int         T;

	ImageStack posPatches;
	ImageStack negPatches;

	ImageStack pos_ft;
	ImageStack neg_ft;

	float[][] featsP;
	float[][] featsN;

	FilterSet fs;

	float[] vals;
	float[] angs;
	float[] rads;

	double[][] adaboost;

    static float TwoPi = (float) (2*Math.PI);

	public void run(String arg0)
	{

        patchRadius     = (int) Prefs.get("advantra.critpoint.patch_radius", 15);

        String[] angleChoices       =  new String[]{"20", "40", "60", "80"};
        String angleStartChoice  	=  (String) Prefs.get("advantra.critpoint.start_ang_scale", 	"40");
		String angleEndChoice       =  (String) Prefs.get("advantra.critpoint.end_ang_scale", 	    "60");

        ring1  	    = Prefs.get("advantra.critpoint.start_ring", 	0.3);
		ring2    	= Prefs.get("advantra.critpoint.end_ring", 		0.7);
		nr_ring		= (int)Prefs.get("advantra.critpoint.nr_ring",  2);

        radial1  	    = Prefs.get("advantra.critpoint.start_radial",      0.4);
        radial2    	    = Prefs.get("advantra.critpoint.end_radial", 	    0.8);
        nr_radial		= (int)Prefs.get("advantra.critpoint.nr_radial",    3);

		innerRadius     = Prefs.get("advantra.critpoint.inner_rad", 0.2);

		train_folder    = (String)Prefs.get("advantra.critpoint.train_folder",
				(System.getProperty("user.home")+File.separator));
		test_folder     = (String)Prefs.get("advantra.critpoint.test_folder",
				(System.getProperty("user.home")+File.separator));

        T               = (int)Prefs.get("advantra.critpoint.T", 5);

		GenericDialog gd = new GenericDialog("Critpoint Detection");

        gd.addNumericField("patch_radius", patchRadius, 0, 5, "");

        gd.addChoice("alfa_1",      angleChoices, angleStartChoice);
        gd.addChoice("alfa_2",      angleChoices, angleEndChoice);

        gd.addNumericField("ring_1:",    ring1,       1);
        gd.addNumericField("ring_2:",    ring2,       1);
        gd.addNumericField("nr_rings:",  nr_ring, 	  0,  5, "");

        gd.addNumericField("radial_1:",    radial1,       1); //radial1
        gd.addNumericField("radial_2:",    radial2,       1); //radial2
        gd.addNumericField("nr_radials:",  nr_radial, 	  0,  5, "");//nr_radials

		gd.addNumericField("inner radius:",    innerRadius,       1);

		gd.addStringField("train folder : ", train_folder, 	40);
		gd.addStringField("test  folder : ", test_folder, 	40);

		gd.addNumericField("T: ",					T, 0, 5, "");

		gd.showDialog();
		if (gd.wasCanceled()) return;

        patchRadius =  (int)gd.getNextNumber();
        System.out.println("patch radius = " + patchRadius);

        switch ((int) gd.getNextChoiceIndex()) {
            case 0: angScale1 = 20; angleStartChoice = "20"; break;
            case 1: angScale1 = 40; angleStartChoice = "40"; break;
            case 2: angScale1 = 60; angleStartChoice = "60"; break;
            case 3: angScale1 = 80; angleStartChoice = "80"; break;
            default: angScale1 = 40;angleStartChoice = "40"; break;
        }

        switch ((int) gd.getNextChoiceIndex()) {
            case 0: angScale2 = 20; angleEndChoice = "20"; break;
            case 1: angScale2 = 40; angleEndChoice = "40"; break;
            case 2: angScale2 = 60; angleEndChoice = "60"; break;
            case 3: angScale2 = 80; angleEndChoice = "80"; break;
            default: angScale2 = 40;angleEndChoice = "40"; break;
        }

        angScale = new int[(angScale2-angScale1)/20+1];
        int cnt = 0;
        for (int i = angScale1; i <= angScale2; i+=20) angScale[cnt++] = i;

        ring1 = gd.getNextNumber();
        ring2 = gd.getNextNumber();
        nr_ring = (int) gd.getNextNumber();

        rings = new double[nr_ring];
        for (int i = 0; i <nr_ring; i++) rings[i] = (i == 0) ? ring1 : ring1 + i * ((ring2 - ring1) / (nr_ring - 1));
        // TODO: block inputs higher than 1.0 for rings

        radial1 = gd.getNextNumber();
        radial2 = gd.getNextNumber();
        nr_radial = (int) gd.getNextNumber();

        radials = new double[nr_radial];
        for (int i = 0; i <nr_radial; i++) radials[i] = (i == 0) ? radial1 : radial1 + i * ((radial2 - radial1) / (nr_radial - 1));
        // TODO: block inputs higher than 1.0 for raidals

		innerRadius = gd.getNextNumber();

		train_folder = 	gd.getNextString();
		test_folder = 	gd.getNextString();

		T = (int) gd.getNextNumber();

		Prefs.set("advantra.critpoint.start_ang_scale", angleStartChoice);
		Prefs.set("advantra.critpoint.end_ang_scale", 	angleEndChoice);

		Prefs.set("advantra.critpoint.start_ring", 	    ring1);
		Prefs.set("advantra.critpoint.end_ring", 	    ring2);
		Prefs.set("advantra.critpoint.nr_ring",         nr_ring);

        Prefs.set("advantra.critpoint.start_radial",    radial1);
        Prefs.set("advantra.critpoint.end_radial", 	    radial2);
        Prefs.set("advantra.critpoint.nr_radial",       nr_radial);

        Prefs.set("advantra.critpoint.inner_rad",       innerRadius);

		Prefs.set("advantra.critpoint.train_folder",    train_folder);
		Prefs.set("advantra.critpoint.test_folder",     test_folder);
		Prefs.set("advantra.critpoint.patch_radius", 	patchRadius);

        Prefs.set("advantra.critpoint.T", 	            T);

        /*
		 * generate filters to score on example profiles (generate features)
		 */
        fs = new FilterSet(angScale, rings, radials, innerRadius);
        int nrFilters = fs.circConfs.size()+fs.radlConfs.size();
        //fs.print();
        System.out.println(nrFilters+" filters formed");
        /*
        show them
        */
        patchDiameter = 2*patchRadius+1;
        ImagePlus all_feats = new ImagePlus("ALL", fs.plot(patchDiameter));
        all_feats.setTitle("All_Features");
        all_feats.show();

		int toAlloc = Calc.circularProfileSize(patchRadius);
		vals = new float[toAlloc];
		angs = new float[toAlloc];
		rads = new float[toAlloc];

        /*
        	all the params are loaded...
         */

		File dir = new File(train_folder);
		train_folder = dir.getAbsolutePath();
		if(!dir.isDirectory() ){
			IJ.error("Wrong train-set directory: "+train_folder+"   closing...");
			return;
		}

		File[] files_tif = listFilesEndingWith(dir, ".tif");
		File[] files_pos = new File[files_tif.length];
		File[] files_neg = new File[files_tif.length];

		Vector<double[][]> locs_pos = new Vector<double[][]>(files_tif.length);
		Vector<double[][]> locs_neg = new Vector<double[][]>(files_tif.length);
        Vector<double[][]> locs_tst = new Vector<double[][]>(files_tif.length);

		int total_pos = 0;
		int total_neg = 0;
		int curr_pos = 0;
		int curr_neg = 0;

        // to store features
        pos_ft = new ImageStack(nrFilters, 1);
		neg_ft = new ImageStack(nrFilters, 1);

        // to store examples
        posPatches = new ImageStack(patchDiameter, patchDiameter);
		negPatches = new ImageStack(patchDiameter, patchDiameter);

		System.out.println("\n## TRAIN ##  \n"+train_folder+"\n------------------------------------\n");

        // img1, W1, H1, showStk, showImg, ovly are just for visualization
        ImagePlus img1 = new ImagePlus(files_tif[0].getAbsolutePath());
        int W1 = img1.getWidth();
        int H1 = img1.getHeight();
        ImageStack showStk = new ImageStack(W1, H1);
        Overlay ovly = new Overlay();

		for (int i = 0; i < files_tif.length; i++) { // for each tif file

			System.out.print(""+files_tif[i].getName()+" "+i+"/"+(files_tif.length-1)+"  ...  ");

			img = new ImagePlus(files_tif[i].getAbsolutePath());
			resetCalibration();
			convertToFloatImage();
			H = img.getHeight();
			W = img.getWidth();

			curr_pos = 0;
			curr_neg = 0;

			AnalyzeCSV readCSV;
			File[] check;
			String suffix;
            ImagePlus readMask;
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
				total_pos += A.length;
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

                // fix: reducing B take those that are out of image
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
				total_neg += B.length;
				System.out.print("\t "+B.length+" negatives ");
			}
			else{
				locs_neg.add(null);
				System.out.print("\t no negatives");
			}

			System.out.println();

            if((curr_pos>0 || curr_neg>0) && true){

                showStk.addSlice("SAMPLE"+i, img.getProcessor());

				for (int k = 0; k < curr_pos; k++) {

					int atX = (int)locs_pos.get(i)[k][0];
					int atY = (int)locs_pos.get(i)[k][1];

                    Calc.getProfile(img, atX, atY, patchRadius, vals, angs, rads);
					pos_ft.addSlice(new FloatProcessor(nrFilters, 1,    Calc.getProfileResponse(fs, vals, angs, rads)));
					posPatches.addSlice("pos_train_sample_"+i,          Calc.getProfilePatch(img, atX, atY, patchRadius));
					PointRoi pt = new PointRoi(atX+0.5, atY+0.5);
					pt.setStrokeColor(Color.RED);
//                    pt.setPosition(showStk.getSize());
                    pt.setPosition(0, showStk.getSize(), 0);
					ovly.addElement(pt);

				}

                // fix: select same amount of negative samples
                // only in case of artificial dataset because of the disbalance pos/neg
                // each image will contribute in as many negatives as positives
                // random without sampling
                boolean[] chs = new boolean[curr_neg];
                Random rand = new Random();
                int cntMtches = 0;
                while (cntMtches < curr_pos) {  // select curr_pos random ones
                    int rIdx = rand.nextInt(curr_neg);
                    if (!chs[rIdx]) {
                        chs[rIdx] = true;
                        cntMtches++;
                    }
                }

                System.out.println("feature extraction for negatives "+img.getTitle());
				for (int k = 0; k < curr_neg; k++) {

                    if (chs[k]){

                        int atX = (int)locs_neg.get(i)[k][0];
                        int atY = (int)locs_neg.get(i)[k][1];

                        Calc.getProfile(img, atX, atY, patchRadius, vals, angs, rads);
                        neg_ft.addSlice(new FloatProcessor(nrFilters, 1,    Calc.getProfileResponse(fs, vals, angs, rads)));
                        negPatches.addSlice("neg_train_sample_"+i,          Calc.getProfilePatch(img, atX, atY, patchRadius));
                        PointRoi pt = new PointRoi(atX+0.5, atY+0.5);
                        pt.setStrokeColor(Color.BLUE);
                        pt.setName("negative" + k);
//                        pt.setPosition(showStk.getSize());
                        pt.setPosition(0, showStk.getSize(), 0);
                        ovly.addElement(pt);

                    }

				}

			} // if there were some

		} // loop files

        // visualization
        ImagePlus showImg = new ImagePlus("VIZ", showStk);
        showImg.setOverlay(ovly);
        showImg.show();
        showImg.getCanvas().zoomIn(0,0);
        showImg.getCanvas().zoomIn(0,0);

		System.out.println("////\n total (+) : " + posPatches.getSize());
		System.out.println(" total (-) : "       + negPatches.getSize()+"\n//-------------//\n");

        new ImagePlus("POS", posPatches).show();
		new ImagePlus("NEG", negPatches).show();

        // equalize number of positives and negatives
		featsP = new float[pos_ft.getSize()][nrFilters];
		featsN = new float[neg_ft.getSize()][nrFilters];

		for (int g = 0; g < pos_ft.getSize(); g++){
			float[] getPix = (float[]) pos_ft.getProcessor(g + 1).getPixels();
			for (int g1 = 0; g1 < getPix.length; g1++){
				featsP[g][g1] = getPix[g1];
			}
		}

		for (int g = 0; g < featsN.length; g++) {
            // randomize when choosing negative
            //int chooseIdx =   rd.nextInt(neg_ft.getSize());
			float[] getPix = (float[]) neg_ft.getProcessor(g + 1).getPixels();
			for (int g1 = 0; g1 < getPix.length; g1++){
				featsN[g][g1] = getPix[g1];
			}
		}

		System.out.println(" done extracting train features:\n" +
                "(+) " + featsP.length + " x " + featsP[0].length +
                "(-) " + featsN.length + " x " + featsN[0].length);

        // form weka instances object
        int totalFeats = (featsP[0].length == featsN[0].length)? featsP[0].length : 0;
        FastVector attList = new FastVector();
        for (int i = 0; i < totalFeats; i++) {
            attList.addElement(new Attribute("feat."+i));
        }
        FastVector categ = new FastVector(); categ.addElement("yes"); categ.addElement("no");
        attList.addElement(new Attribute("class", categ));

        Instances trn = new Instances("Trainset", attList, (featsP.length + featsN.length));

        IJ.log(trn.numInstances()+", atts "+trn.numAttributes());



        if(true) { System.out.println("closing... "); return; }


        GenericDialog gd1 = new GenericDialog("AdaBoost training");
        gd1.addMessage("How many feats to keep? stored: " + T);
        gd1.addNumericField("T", T, 0, 6, " (total " + nrFilters + ")");
        gd1.showDialog();
        if (gd1.wasCanceled()) return;
        T = (int)       gd1.getNextNumber();

		System.out.print("training AdaBoost..."+T);
		adaboost = trueAdaBoost(featsP, featsN, T);
		System.out.println("done.");

		for (int i = 0; i < adaboost.length; i++) {
			for (int j = 0; j < adaboost[0].length; j++) System.out.print("\t" + IJ.d2s(adaboost[i][j], 2));
			System.out.println();
		}
		System.out.println();

		int[] best_indexes = new int[adaboost.length];
		for (int i = 0 ; i < adaboost.length; i++){
			best_indexes[i] = (int) adaboost[i][0];
		}

		// show training results
		if (adaboost.length>1){

			// show the best features
			ImagePlus imp2 = new ImagePlus();
            int imp2_size = patchDiameter;
			ImageStack best_feat = new ImageStack(imp2_size, imp2_size);
			for (int i = 0; i < adaboost.length; i++){
				best_feat.addSlice(fs.plotOne((int) adaboost[i][0], imp2_size));
			}
			imp2.setStack("best_feats", best_feat);
			imp2.show();

			ImagePlus imp1 = new ImagePlus();
			ImageStack imstack = new ImageStack(528, 255);
			imp1.setDimensions(1, 1, adaboost.length);
			for (int i = 0; i < T; i++) {
				Plot plot = plotFearutePerAllSamples(featsP, featsN, (int)adaboost[i][0], adaboost[i][2]);
				imstack.addSlice("thresh = " + adaboost[i][2], plot.getProcessor());
			}
			imp1.setStack("best_feats_scores", imstack);
			imp1.show();

		}

	    File dir_test = new File(test_folder);
	    test_folder = dir_test.getAbsolutePath();
		if(!dir_test.isDirectory() ){
			IJ.error("Wrong testset directory "+test_folder);
			return;
		}

        System.out.println("\n## TEST  ##  \n"+test_folder +"\n------------------------------------\n");

        int curr_tst = 0;
        int total_tst = 0;

		boolean doTest = true;

		File[] test_files_tif = listFilesEndingWith(dir_test, ".tif");
        File[] files_tst = new File[test_files_tif.length];

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
                System.out.println("i: "+i);
                files_tst[i] = check[0];
                readMask = new ImagePlus(files_tst[i].getAbsolutePath());
                //readMask.show();
                // extract locations with logical 1
                double[][] C = extractLocations((ByteProcessor) readMask.getProcessor());
                locs_tst.add(C);
                curr_tst = C.length;
                total_tst += C.length;

            }
            else{
                locs_tst.add(null);
                System.out.print("no mask found!");
            }

            // check if the mask has anything and add it to list ov Overlays and add to ImageStack

            if(curr_tst>0){

                img = new ImagePlus(test_files_tif[i].getAbsolutePath());
                resetCalibration();
                convertToFloatImage();

                H = img.getHeight();
                W = img.getWidth();

                ovly = new Overlay();
                System.out.println("\nclassifying " + curr_tst + " locations from " + img.getTitle());

				for (int k = 0; k < curr_tst; k++) {

                    int atX = (int)locs_tst.get(i)[k][0];
                    int atY = (int)locs_tst.get(i)[k][1];

                    Calc.getProfile(img, atX, atY, patchRadius, vals, angs, rads);
                    int res = applyAdaBoost(adaboost, Calc.getProfileResponseSelection(fs, best_indexes, vals, angs, rads));

                    if(res==1){
                        PointRoi pt = new PointRoi(atX+0.5, atY+0.5);
                        pt.setStrokeColor(Color.YELLOW);
                        ovly.addElement(pt);
                    }

                }

                System.out.println("done.");

				ImagePlus showIt = new ImagePlus("RES_"+test_files_tif[i].getName(), img.getProcessor());
				showIt.setOverlay(ovly);

				showIt.show();
                showIt.getCanvas().zoomIn(0,0);
                showIt.getCanvas().zoomIn(0,0);
                showIt.getCanvas().zoomIn(0,0);
                showIt.getCanvas().zoomIn(0,0);

            }

        System.out.println();

		}

        //show ImageStack with test images and each has overlay (if possible)

	}

	private int profileLength(int rin)
	{

		int cnt = 0;

		for (int x = -rin; x <= rin; x++){
			for (int y = -rin; y <= rin; y++){
                if (x*x+y*y<=rin*rin)
				    cnt++;
			}
		}

		return cnt;

	}

	public void profile(
								ImagePlus img,
								int xin,
								int yin,
								int rin
	)
	{

		// will extract out vals[], rads[], and angs[]
		// those will be used to score() on filterSet

		int xs = xin-rin;
		int xe = xin+rin;

		int ys = yin-rin;
		int ye = yin+rin;

		int cnt = 0;
		for (int x = xs; x <= xe; x++){
			for (int y = ys; y <= ye; y++){

				float d = (float) (Math.pow(x-xin, 2) + Math.pow(y-yin, 2));

				if(d <= rin*rin){

					vals[cnt] = img.getProcessor().getPixelValue(x, y);
					rads[cnt] = (float) (Math.sqrt(d) / rin);
					angs[cnt] = (float) (Math.atan2(y-yin, x-xin) + Math.PI);
					angs[cnt] = (angs[cnt]>=(float)(2*Math.PI))? 0 : angs[cnt];
					angs[cnt] = (angs[cnt]<0)? 0 : angs[cnt];
					cnt++;

				}

			}
		}
    }

	public void profile(
							   ImagePlus img,
							   ImagePlus vx,
							   ImagePlus vy,
							   int xin,
							   int yin,
							   int rin
	)
	{

		// will extract out vals[], rads[], and angs[], those will be used to score() on filterSet

		int xs = xin-rin;
		int xe = xin+rin;

		int ys = yin-rin;
		int ye = yin+rin;

		int cnt = 0;
		for (int x = xs; x <= xe; x++){
			for (int y = ys; y <= ye; y++){

				float d = (float) (Math.pow(x-xin, 2) + Math.pow(y-yin, 2));

				if(d <= rin*rin){

					float dx = x-xin;	//Math.sin(ang[len]);
					float dy = y-yin;	//-Math.cos(ang[len]);

					// unit vector pointing towards the center
					dx = (float) (dx / Math.sqrt((x - xin) * (x - xin) + (y - yin) * (y - yin)));
					dy = (float) (dy / Math.sqrt((x - xin) * (x - xin) + (y - yin) * (y - yin)));

					float I 	= img.getProcessor().getPixelValue(x, y);
					float VX 	= vx.getProcessor().getPixelValue(x, y);
					float VY 	= vy.getProcessor().getPixelValue(x, y);

					vals[cnt] = (float) Math.abs(I * VX * dx + I * VY * dy);
					//vals[cnt] = img.getProcessor().getPixelValue(x, y);
					rads[cnt] = (float) (Math.sqrt(d) / rin);
					angs[cnt] = (float) (Math.atan2(y-yin, x-xin) + Math.PI);

					cnt++;

				}

			}
		}

	}


	public ImageStack plotProfile()
	{
		ImageStack viz = new ImageStack(400, 200);

		// find max for plotting
		float max_val = vals[0];
		for (int i = 1; i < vals.length; i++){
			if (vals[i] > max_val) max_val = vals[i];
		}

		Plot p = new Plot("circular_profile", "angle[rad]", "value");
		p.setLimits(0, 2*Math.PI, 0, max_val);
		p.setSize(400, 200);
		p.addPoints(angs, vals, Plot.BOX);
		viz.addSlice("circular_profile", p.getProcessor());

		Plot p1 = new Plot("radial_profile", "radius", "value");
		p1.setLimits(0, 1, 0, max_val);
		p1.setSize(400, 200);
		p1.addPoints(rads, vals, Plot.CIRCLE);
		viz.addSlice("radial_profile", p1.getProcessor());

		// check circular score
		float[] c_idx = new float[fs.circConfs.size()];
		float[] c_sco = new float[fs.circConfs.size()];
		for (int i = 0; i < fs.circConfs.size(); i++){
			c_idx[i] = i;
			c_sco[i] = fs.score[i];
		}
		Plot p3 = new Plot("filtering_score", "circular_configuration", "score", c_idx, c_sco);
		p3.setSize(400, 200);
		viz.addSlice("scores", p3.getProcessor());

		// check radial score
		float[] r_idx = new float[fs.radlConfs.size()];
		float[] r_sco = new float[fs.radlConfs.size()];
		for (int i = 0; i < fs.radlConfs.size(); i++){
			r_idx[i] = i;
			r_sco[i] = fs.score[fs.circConfs.size()+i];
		}
		Plot p4 = new Plot("filtering_score", "radial_configuration", "score", r_idx, r_sco);
		p4.setSize(400, 200);
		viz.addSlice("scores", p4.getProcessor());

		return viz;
	}

	private void addProfilePoints(
		ImagePlus img,
		int xin,
		int yin,
		int rin
	)
	{
		Overlay ol = new Overlay();

		int xs = xin-rin;
		int xe = xin+rin;

		int ys = yin-rin;
		int ye = yin+rin;

		for (int x = xs; x <= xe; x++){
			for (int y = ys; y <= ye; y++){

				if(Math.pow(x-xin, 2) + Math.pow(y-yin, 2) <= rin*rin){

					ol.add(new PointRoi(x-0.5, y-0.5));

				}

			}
		}

		img.setOverlay(ol);

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
	    plot.setColor(Color.RED);
	    
	    plot.addPoints(xn, yn, Plot.LINE);
	    plot.setColor(Color.BLACK);
	    plot.setLineWidth(4);
	    
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

		int atX = 	img.getWindow().getCanvas().offScreenX(e.getX());
		int atY = 	img.getWindow().getCanvas().offScreenY(e.getY());

		addProfilePoints(img, atX, atY, patchRadius);
        img.getCanvas().unzoom();
        img.getCanvas().zoomIn(atX, atY);

//        if(useDotProd) profile(img_click, Vx, Vy, atX, atY, radius);
//        else profile(img_click, atX, atY, radius);

//		fs.score(vals, angs, rads);
//        float[] ft_id = new float[fs.circConfs.size()+fs.radlConfs.size()];
//        for (int i = 0; i < ft_id.length; i++) ft_id[i] = i+1;

        int choose_ft;
        GenericDialog gd = new GenericDialog("Choose Feature");
        gd.addNumericField("choose feature:" ,	1, 0);
        gd.showDialog();
        if (gd.wasCanceled()) return;
        choose_ft = (int)       gd.getNextNumber();

        System.out.print("chosen ft: "+choose_ft);
        float[] ft_id = new float[fs.circConfs.get(choose_ft).nrRot];
        for (int i = 0; i < ft_id.length; i++) ft_id[i] = i;

        new ImagePlus("chosen_feature_", plotProfile()).show();

//        for (int i = 0; i < vals.length; i++){
//            System.out.println(""+i+" -> "+vals[i]);
//        }

        new ImagePlus("chosen_ft"+choose_ft+"/"+(fs.circConfs.size()-1), fs.circConfs.get(choose_ft).plotAllRotations(65)).show();
        Plot p = new Plot("ft"+choose_ft+"/"+(fs.circConfs.size()-1), "rot#", "scores per rot.", ft_id, fs.circConfs.get(choose_ft).scoreAllRot(vals, angs, rads));
        p.show();


	}

	public void mousePressed(MouseEvent e) {}

	public void mouseReleased(MouseEvent e) {}

	public void mouseEntered(MouseEvent e) {}

	public void mouseExited(MouseEvent e) {}

}

/*	public float[] extract(
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
		// allocate out
		float[] outFeat = new float[filterSet.circConfs.size()];
		for (int fI = 0; fI < filterSet.circConfs.size(); fI++){
			outFeat[fI] = filterSet.circConfs.get(fI).calculateScore(val, ang);
		}

		return outFeat;

	}*/
