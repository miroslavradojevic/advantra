package advantra.critpoint;

import java.awt.*;
import java.io.File;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Random;
import java.util.Vector;

import advantra.file.AnalyzeCSV;
import advantra.general.CreateDirectory;
import advantra.tools.BranchModel2D;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.io.FileSaver;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

public class GenerateBranch implements PlugIn {
	
	/*
	 * generates the branch, shows it, and saves it in a output .tif file
	 */

    String outDirPath;
    String outDirTrain;
    String outDirTest;

    int bgBias, bgRange, fgBias, fgRange;

	boolean addPoisson;

	int minStd, rngStd;

	public void run(
						   String arg0
	)
	{

		GenericDialog gd = new GenericDialog("Generate Random Bifurcations");

		gd.addMessage("generate bifurcations \n export locations of critical points");

		gd.addNumericField("Image Width:",      129,  0,  5, "");
		gd.addNumericField("Image Height:", 	129,  0,  5, "");
		gd.addNumericField("train images:", 	200,  0,  5, "");
        gd.addNumericField("test  images:", 	5,    0,  5, "");

        String moment = (new SimpleDateFormat("dd-MM-yyyy-HH-mm-ss")).format(Calendar.getInstance().getTime());
        String def_dir = System.getProperty("user.home")+File.separator+
                "cp_"+moment+File.separator;

        gd.addStringField("destination folder: ", def_dir, 40);

        gd.addNumericField("background (poisson):",  50, 0, 5, "");
        gd.addNumericField("foreground (poisson):", 70, 0, 5, "");

		gd.addNumericField("min.std:", 1, 0, 5, "");
		gd.addNumericField("max.std:", 2, 0, 5, "(random int 0 -- max.std-1)");

		gd.addCheckbox("add poisson", true);

		gd.showDialog();
		if (gd.wasCanceled()) return;

        int H = (int) gd.getNextNumber();
		int W = (int) gd.getNextNumber();
		int N = (int) gd.getNextNumber();
        int L = (int) gd.getNextNumber();
        outDirPath = gd.getNextString();
        if (!outDirPath.endsWith(File.separator)) {
            outDirPath += File.separator;
        }
        bgBias = (int) gd.getNextNumber();
        fgBias = (int) gd.getNextNumber();

		minStd = (int) gd.getNextNumber();
		rngStd = (int) gd.getNextNumber();

		gd.getNextBoolean();

//        bgBias  = 20;
        bgRange = 5;
//        fgBias  = 150;
        fgRange = 5;

        FileWriter 			fw;
        String              file_name;

        outDirTrain = outDirPath+"train"+File.separator;
        outDirTest = outDirPath+"test"+File.separator;

        System.out.println("will be storing train in: "+outDirTrain);

        CreateDirectory.createOneDir(outDirPath);
        CreateDirectory.createOneDir(outDirTrain);
        CreateDirectory.createOneDir(outDirTest);

        BranchModel2D bm2d = new BranchModel2D(W, H);

        for (int i = 0; i < N; i ++){

			bm2d.generateRandomBranch(bgBias, bgRange, fgBias, fgRange, minStd, rngStd);
			ImageProcessor ip_to_add = new ByteProcessor(W, H);

			for (int j = 0; j < W*H; j++) ip_to_add.setf(j, bm2d.values[j]);

            file_name = "bch_"+i;
            try{
                /*
                 */
                fw = new FileWriter(outDirTrain+file_name+".pos");

				for (int a = 0; a < 1; a++){ // take it several times
					fw.write(IJ.d2s(bm2d.pc[1], 2)+", "+IJ.d2s(bm2d.pc[0], 2)+"\n");
				}

//                fw.write(IJ.d2s(bm2d.pc[1]-1, 2)+", "+IJ.d2s(bm2d.pc[0], 2)+"\n");
//                fw.write(IJ.d2s(bm2d.pc[1]+1, 2)+", "+IJ.d2s(bm2d.pc[0], 2)+"\n");
//                fw.write(IJ.d2s(bm2d.pc[1], 2)+", "+IJ.d2s(bm2d.pc[0]-1, 2)+"\n");
//                fw.write(IJ.d2s(bm2d.pc[1], 2)+", "+IJ.d2s(bm2d.pc[0]+1, 2)+"\n");

                fw.close();
                /*
                 */
                fw = new FileWriter(outDirTrain+file_name+".neg");
                fw.write(IJ.d2s(bm2d.p1[1], 2)+", "+IJ.d2s(bm2d.p1[0], 2)+"\n");
                fw.write(IJ.d2s(bm2d.p2[1], 2)+", "+IJ.d2s(bm2d.p2[0], 2)+"\n");
                fw.write(IJ.d2s(bm2d.p3[1], 2)+", "+IJ.d2s(bm2d.p3[0], 2)+"\n");

                // add endpoints
                fw.write(IJ.d2s(bm2d.p1[1], 2)+", "+IJ.d2s(bm2d.p1[0], 2)+"\n");
                fw.write(IJ.d2s(bm2d.p2[1], 2)+", "+IJ.d2s(bm2d.p2[0], 2)+"\n");
                fw.write(IJ.d2s(bm2d.p3[1], 2)+", "+IJ.d2s(bm2d.p3[0], 2)+"\n");
                // add midpoints
                int[] mid1 = new int[2];
                mid1[0] = (bm2d.pc[0]+bm2d.p1[0])/2;
                mid1[1] = (bm2d.pc[1]+bm2d.p1[1])/2;
                int[] mid2 = new int[2];
                mid2[0] = (bm2d.pc[0]+bm2d.p2[0])/2;
                mid2[1] = (bm2d.pc[1]+bm2d.p2[1])/2;
                int[] mid3 = new int[2];
                mid3[0] = (bm2d.pc[0]+bm2d.p3[0])/2;
                mid3[1] = (bm2d.pc[1]+bm2d.p3[1])/2;
                fw.write(IJ.d2s(mid1[1], 2)+", "+IJ.d2s(mid1[0], 2)+"\n");
                fw.write(IJ.d2s(mid2[1], 2)+", "+IJ.d2s(mid2[0], 2)+"\n");
                fw.write(IJ.d2s(mid3[1], 2)+", "+IJ.d2s(mid3[0], 2)+"\n");



                int bor_row, bor_col;

				for (float k = 0.1f; k <=0.5; k+=0.1f) {      // till +/-0.5 std

					// cross-line

					// bor1
					bor_row = (int) ((bm2d.p1[0]+bm2d.pc[0])/2 + k*bm2d.rstd*(-bm2d.v1[1]));
					bor_col = (int) ((bm2d.p1[1]+bm2d.pc[1])/2 + k*bm2d.rstd*bm2d.v1[0]);
					fw.write(IJ.d2s(bor_col, 2)+", "+IJ.d2s(bor_row, 2)+"\n");
					bor_row = (int) ((bm2d.p1[0]+bm2d.pc[0])/2 - k*bm2d.rstd*(-bm2d.v1[1]));
					bor_col = (int) ((bm2d.p1[1]+bm2d.pc[1])/2 - k*bm2d.rstd*bm2d.v1[0]);

					fw.write(IJ.d2s(bor_col, 2)+", "+IJ.d2s(bor_row, 2)+"\n");
					// bor2
					bor_row = (int) ((bm2d.p2[0]+bm2d.pc[0])/2 + k*bm2d.rstd*(-bm2d.v2[1]));
					bor_col = (int) ((bm2d.p2[1]+bm2d.pc[1])/2 + k*bm2d.rstd*bm2d.v2[0]);
					fw.write(IJ.d2s(bor_col, 2)+", "+IJ.d2s(bor_row, 2)+"\n");
					bor_row = (int) ((bm2d.p2[0]+bm2d.pc[0])/2 - k*bm2d.rstd*(-bm2d.v2[1]));
					bor_col = (int) ((bm2d.p2[1]+bm2d.pc[1])/2 - k*bm2d.rstd*bm2d.v2[0]);
					fw.write(IJ.d2s(bor_col, 2)+", "+IJ.d2s(bor_row, 2)+"\n");

					// bor3
					bor_row = (int) ((bm2d.p3[0]+bm2d.pc[0])/2 + k*bm2d.rstd*(-bm2d.v3[1]));
					bor_col = (int) ((bm2d.p3[1]+bm2d.pc[1])/2 + k*bm2d.rstd*bm2d.v3[0]);
					fw.write(IJ.d2s(bor_col, 2)+", "+IJ.d2s(bor_row, 2)+"\n");
					bor_row = (int) ((bm2d.p3[0]+bm2d.pc[0])/2 - k*bm2d.rstd*(-bm2d.v3[1]));
					bor_col = (int) ((bm2d.p3[1]+bm2d.pc[1])/2 - k*bm2d.rstd*bm2d.v3[0]);
					fw.write(IJ.d2s(bor_col, 2)+", "+IJ.d2s(bor_row, 2)+"\n");

				}

                fw.close();

            }
            catch(IOException exIO){}

			// imagej
			ImagePlus im_to_add = new ImagePlus(file_name, ip_to_add);

			if (addPoisson) {

				// imagescience
				//IJ.run(new ImagePlus(file_name, ip_to_add), "RandomJ Poisson", "mean=10 insertion=additive");

				IJ.run(im_to_add, "Poisson Noise", "");
				//the other one was more convenient giving 8-bit output
				// and storing it directly into opened image
				// but you have to be careful setting bg and fg values

			}

			// save image
			FileSaver fs = new FileSaver(im_to_add);
			fs.saveAsTiff(outDirTrain+file_name+".tif");

        }

        // test set with masks
        for (int i = 0; i < L; i ++){

            bm2d.generateRandomBranch(bgBias, bgRange, fgBias, fgRange, minStd, rngStd);

            ImageProcessor ip_to_add = new ByteProcessor(W, H);
            ImageProcessor ipmask_to_add = new ByteProcessor(W, H);

            for (int j = 0; j < W*H; j++) {
                ip_to_add.setf(j, bm2d.values[j]);
                ipmask_to_add.setf(j, bm2d.mask[j]);
            }
            file_name = "bch_"+i;
            ImagePlus im_to_add = new ImagePlus(file_name, ip_to_add);

			if (addPoisson) {
				IJ.run(im_to_add, "Poisson Noise", "");
			}

            //save
            FileSaver fs = new FileSaver(im_to_add);
            fs.saveAsTiff(outDirTest+file_name+".tif");

            ImagePlus immask_to_add = new ImagePlus(file_name, ipmask_to_add);

            //save
            fs = new FileSaver(immask_to_add);
            fs.saveAsTiff(outDirTest+file_name+".mask");

        }

        IJ.log("trainset dir. path:\n\n" + outDirTrain + "\n\ntestset dir. path:\n\n" + outDirTest);
		//GenericDialog gdOut = new GenericDialog("Directories");
		gd.removeAll();
		gd.addStringField("Train: ", outDirTrain, 40);
		gd.addStringField("Test: ", outDirTest, 40);
		gd.showDialog();
		//if (gd.wasCanceled()) return;

		/*
		show them
		 */

		File dir = new File(outDirTrain);
		String train_folder = dir.getAbsolutePath();
		if(!dir.isDirectory() ){
			IJ.error("Wrong trainset directory: " + train_folder + "   closing...");
			return;
		}

		int curr_pos = 0;
		int curr_neg = 0;

		File[] files_tif = listFilesEndingWith(dir, ".tif");
		File[] files_pos = new File[files_tif.length];
		File[] files_neg = new File[files_tif.length];

		Vector<double[][]> locs_pos = new Vector<double[][]>(files_tif.length);
		Vector<double[][]> locs_neg = new Vector<double[][]>(files_tif.length);

		ImagePlus img;
		ImageStack showStk  = new ImageStack();
		boolean initialized = false;
		Overlay ovly        = new Overlay();
		boolean equal = false;
		boolean sameSize = true;

		for (int i = 0; i < files_tif.length; i++){ // for each tif file used to train

			System.out.print("reading "+files_tif[i].getName()+" "+i+"/"+(files_tif.length-1)+"  ...  ");

			img = convertToFloatImage(new ImagePlus(files_tif[i].getAbsolutePath()));
			H = img.getHeight();
			W = img.getWidth();

			if(true && !initialized){
				showStk = new ImageStack(W, H);
				initialized = true;
			}

			curr_pos = curr_neg = 0;

			AnalyzeCSV readCSV;
			File[] check;
			String suffix;
			//ImagePlus readMask;
			file_name = files_tif[i].getName();
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
//					boolean isIn =
//							take_col>patchRadius &&
//									take_col<img.getWidth()-patchRadius &&
//									take_row>patchRadius &&
//									take_row<img.getHeight()-patchRadius;
//					if(isIn)
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

					PointRoi pt = new PointRoi(atX+0.5, atY+0.5);
					pt.setStrokeColor(Color.RED);

					if(true) {
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

						PointRoi pt = new PointRoi(atX+0.5, atY+0.5);
						pt.setStrokeColor(Color.BLUE);

						if(sameSize) { // they are same size
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


	// terminal command option (empty)
    public static void main(
            String[] args
    )
    {
    }

}
