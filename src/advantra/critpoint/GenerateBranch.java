package advantra.critpoint;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Calendar;

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
import ij.process.ImageProcessor;

public class GenerateBranch implements PlugIn {
	
	/*
	 * generates the branch, shows it, and saves it in a output .tif file
	 */

    String outDirPath;
    String outDirTrain;
    String outDirTest;

    int bgBias, bgRange, fgBias, fgRange;

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
        gd.addNumericField("foreground (poisson):", 100, 0, 5, "");

		gd.showDialog();
		if (gd.wasCanceled()) return;

        int H = (int) gd.getNextNumber();
		int W = (int) gd.getNextNumber();
		int N = (int) gd.getNextNumber();
        int L = (int) gd.getNextNumber();
        outDirPath = gd.getNextString();
        bgBias = (int) gd.getNextNumber();
        fgBias = (int) gd.getNextNumber();

//        bgBias  = 20;
        bgRange = 1;
//        fgBias  = 150;
        fgRange = 1;

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

			bm2d.generateRandomBranch(bgBias, bgRange, fgBias, fgRange);
			ImageProcessor ip_to_add = new ByteProcessor(W, H);

			for (int j = 0; j < W*H; j++) ip_to_add.setf(j, bm2d.values[j]);

            file_name = "bch_"+i;
            try{
                /*
                 */
                fw = new FileWriter(outDirTrain+file_name+".pos");
                fw.write(IJ.d2s(bm2d.pc[1], 2)+", "+IJ.d2s(bm2d.pc[0], 2)+"\n");
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
                // add points at the border

                int bor_row, bor_col;

                // bor1
                bor_row = (int) ((bm2d.p1[0]+bm2d.pc[0])/2 + 2*bm2d.rstd*(-bm2d.v1[1]));
                bor_col = (int) ((bm2d.p1[1]+bm2d.pc[1])/2 + 2*bm2d.rstd*bm2d.v1[0]);
                fw.write(IJ.d2s(bor_col, 2)+", "+IJ.d2s(bor_row, 2)+"\n");
                bor_row = (int) ((bm2d.p1[0]+bm2d.pc[0])/2 - 2*bm2d.rstd*(-bm2d.v1[1]));
                bor_col = (int) ((bm2d.p1[1]+bm2d.pc[1])/2 - 2*bm2d.rstd*bm2d.v1[0]);
                fw.write(IJ.d2s(bor_col, 2)+", "+IJ.d2s(bor_row, 2)+"\n");

                // bor2
                bor_row = (int) ((bm2d.p2[0]+bm2d.pc[0])/2 + 2*bm2d.rstd*(-bm2d.v2[1]));
                bor_col = (int) ((bm2d.p2[1]+bm2d.pc[1])/2 + 2*bm2d.rstd*bm2d.v2[0]);
                fw.write(IJ.d2s(bor_col, 2)+", "+IJ.d2s(bor_row, 2)+"\n");
                bor_row = (int) ((bm2d.p2[0]+bm2d.pc[0])/2 - 2*bm2d.rstd*(-bm2d.v2[1]));
                bor_col = (int) ((bm2d.p2[1]+bm2d.pc[1])/2 - 2*bm2d.rstd*bm2d.v2[0]);
                fw.write(IJ.d2s(bor_col, 2)+", "+IJ.d2s(bor_row, 2)+"\n");

                // bor3
                bor_row = (int) ((bm2d.p3[0]+bm2d.pc[0])/2 + 2*bm2d.rstd*(-bm2d.v3[1]));
                bor_col = (int) ((bm2d.p3[1]+bm2d.pc[1])/2 + 2*bm2d.rstd*bm2d.v3[0]);
                fw.write(IJ.d2s(bor_col, 2)+", "+IJ.d2s(bor_row, 2)+"\n");
                bor_row = (int) ((bm2d.p3[0]+bm2d.pc[0])/2 - 2*bm2d.rstd*(-bm2d.v3[1]));
                bor_col = (int) ((bm2d.p3[1]+bm2d.pc[1])/2 - 2*bm2d.rstd*bm2d.v3[0]);
                fw.write(IJ.d2s(bor_col, 2)+", "+IJ.d2s(bor_row, 2)+"\n");

                fw.close();

            }
            catch(IOException exIO){}

            // imagescience
            //IJ.run(new ImagePlus(file_name, ip_to_add), "RandomJ Poisson", "mean=10 insertion=additive");
            // imagej
            ImagePlus im_to_add = new ImagePlus(file_name, ip_to_add);
            IJ.run(im_to_add, "Poisson Noise", "");
            //the other one was more convenient giving 8-bit output
            // and storing it directly into opened image
            // but you have to be careful setting bg and fg values
            // save image
            FileSaver fs = new FileSaver(im_to_add);
            fs.saveAsTiff(outDirTrain+file_name+".tif");

        }

        // test set with masks
        for (int i = 0; i < L; i ++){

            bm2d.generateRandomBranch(bgBias, bgRange, fgBias, fgRange);

            ImageProcessor ip_to_add = new ByteProcessor(W, H);
            ImageProcessor ipmask_to_add = new ByteProcessor(W, H);

            for (int j = 0; j < W*H; j++) {
                ip_to_add.setf(j, bm2d.values[j]);
                ipmask_to_add.setf(j, bm2d.mask[j]);
            }
            file_name = "bch_"+i;
            ImagePlus im_to_add = new ImagePlus(file_name, ip_to_add);
            IJ.run(im_to_add, "Poisson Noise", "");

            //save
            FileSaver fs = new FileSaver(im_to_add);
            fs.saveAsTiff(outDirTest+file_name+".tif");

            ImagePlus immask_to_add = new ImagePlus(file_name, ipmask_to_add);

            //save
            fs = new FileSaver(immask_to_add);
            fs.saveAsTiff(outDirTest+file_name+".mask");

        }

        IJ.log("trainset dir. path:\n\n" + outDirTrain + "\n\ntestset dir. path:\n\n" + outDirTest);

	}

    // terminal command option (empty)
    public static void main(
            String[] args
    )
    {
    }

}
