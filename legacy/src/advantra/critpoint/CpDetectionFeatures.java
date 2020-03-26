package advantra.critpoint;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 5/29/13
 * Time: 5:54 PM
 * To change this template use File | Settings | File Templates.
 */
public class CpDetectionFeatures implements PlugIn, MouseListener {

	int 		patchRadius, patchDiameter;

	CircularConfiguration2 ccf2;
	ImagePlus feats2;

	CircularConfiguration3 ccf3;
	ImagePlus feats3;

	CircularConfiguration1 ccf1;
	ImagePlus feats1;

	SymmConfiguration scf;
	ImagePlus featsS;

	AsymmConfiguration acf;
	ImagePlus featsA;

	CircularConfiguration4 ccf4;
	ImagePlus feats4;

	public void run(String s) {

		patchRadius     = (int) Prefs.get("advantra.critpoint.patch_radius", 15);

		GenericDialog gd = new GenericDialog("Critpoint Detection");
		gd.addNumericField("patch_radius", patchRadius, 0, 5, "");
		gd.showDialog();
		if (gd.wasCanceled()) return;
		patchRadius =  (int)gd.getNextNumber();
		patchDiameter = 2*patchRadius+1;

		Prefs.set("advantra.critpoint.patch_radius", 	patchRadius);

		ccf3 = new CircularConfiguration3(patchRadius);
		feats3 = new ImagePlus("Y.feat."+patchRadius, ccf3.plotKernels());
		feats3.show();
		ImageCanvas feats3_cnv = feats3.getCanvas();
		feats3_cnv.setName("feats3");
		feats3_cnv.zoomIn(0, 0);
		feats3_cnv.addMouseListener(this);


//		ccf2 = new CircularConfiguration2(patchRadius);
//		feats2 = new ImagePlus("V.feat."+patchRadius, ccf2.plotKernels());
//		feats2.show();
//		ImageCanvas feats2_cnv = feats2.getCanvas();
//		feats2_cnv.setName("feats2");
//		feats2_cnv.zoomIn(0, 0);
//		feats2_cnv.addMouseListener(this);
//
//		ccf1 = new CircularConfiguration1(patchRadius);
//		feats1 = new ImagePlus("I.feat."+patchRadius, ccf1.plotKernels());
//		feats1.show();
//		ImageCanvas feats1_cnv = feats1.getCanvas();
//		feats1_cnv.setName("feats1");
//		feats1_cnv.zoomIn(0, 0);
//		feats1_cnv.addMouseListener(this);
//
//		scf = new SymmConfiguration(patchRadius);
//		featsS = new ImagePlus("S.feat."+patchRadius, scf.plotKernels());
//		featsS.show();
//		ImageCanvas featsS_cnv = featsS.getCanvas();
//		featsS_cnv.setName("featsS");
//		featsS_cnv.zoomIn(0, 0);
//		featsS_cnv.addMouseListener(this);
//
//		acf = new AsymmConfiguration(patchRadius);
//		featsA = new ImagePlus("A.feat."+patchRadius, acf.plotKernels());
//		featsA.show();
//		ImageCanvas featsA_cnv = featsA.getCanvas();
//		featsA_cnv.setName("featsA");
//		featsA_cnv.zoomIn(0, 0);
//		featsA_cnv.addMouseListener(this);

		ccf4 = new CircularConfiguration4(patchRadius);
		feats4 = new ImagePlus("X.feat."+patchRadius, ccf4.plotKernels());
		feats4.show();
		ImageCanvas feats4_cnv = feats4.getCanvas();
		feats4_cnv.setName("feats4");
		feats4_cnv.zoomIn(0, 0);
		feats4_cnv.addMouseListener(this);

	}

	public void mouseClicked(MouseEvent e)
	{

		ImageCanvas srcCanv = (ImageCanvas) e.getSource();
		String source  =  srcCanv.getName();

		if (source=="feats3") {
			int mouseZ = feats3.getCurrentSlice()-1;
			ImagePlus imTest = convertToFloatImage(IJ.openImage());
			System.out.println("caluclate score on "+imTest.getTitle()+"\nwith extracted feature of idx. " + mouseZ);
			imTest.show();

            new ImagePlus("max", ccf3.score(mouseZ, imTest.getProcessor())).show();
			new ImagePlus("per.rot", ccf3.scoreAllRot(mouseZ, imTest.getProcessor())).show();

            //ccf3.score_Experimental(mouseZ, 0, imTest.getProcessor());

//            for (int i = 0; i < CircularConfiguration3.nRot; i++) {
//                new ImagePlus("ft", new FloatProcessor(ccf3.d, ccf3.d, ccf3.kernels.get(mouseZ)[i]));
//            }
		}

		if (source=="feats2") {
			int mouseZ = feats2.getCurrentSlice()-1;
			ImagePlus imTest = convertToFloatImage(IJ.openImage());
			System.out.println("caluclate score on "+imTest.getTitle()+"\nwith extracted feature of idx. " + mouseZ);
			imTest.show();
			new ImagePlus("max", ccf2.score(mouseZ, imTest.getProcessor())).show();
			new ImagePlus("per.rot", ccf2.scoreAllRot(mouseZ, imTest.getProcessor())).show();
		}

		if (source=="feats1") {
			int mouseZ = feats1.getCurrentSlice()-1;
			ImagePlus imTest = convertToFloatImage(IJ.openImage());
			System.out.println("caluclate score on "+imTest.getTitle()+"\nwith extracted feature of idx. " + mouseZ);
			imTest.show();
			new ImagePlus("max", ccf1.score(mouseZ, imTest.getProcessor())).show();
			new ImagePlus("per.rot", ccf1.scoreAllRot(mouseZ, imTest.getProcessor())).show();
		}

		if (source=="featsS") {
			int mouseZ = featsS.getCurrentSlice()-1;
			ImagePlus imTest = convertToFloatImage(IJ.openImage());
			System.out.println("caluclate score on "+imTest.getTitle()+"\nwith extracted feature of idx. " + mouseZ);
			imTest.show();
			new ImagePlus("max", scf.score(mouseZ, imTest.getProcessor())).show();
			new ImagePlus("per.rot", scf.scoreAllRot(mouseZ, imTest.getProcessor())).show();
		}

		if (source=="featsA") {
			int mouseZ = featsA.getCurrentSlice()-1;
			ImagePlus imTest = convertToFloatImage(IJ.openImage());
			System.out.println("caluclate score on "+imTest.getTitle()+"\nwith extracted feature of idx. " + mouseZ);
			imTest.show();
			new ImagePlus("max", acf.score(mouseZ, imTest.getProcessor())).show();
			new ImagePlus("per.rot", acf.scoreAllRot(mouseZ, imTest.getProcessor())).show();
		}

		if (source=="feats4") {
			int mouseZ = feats4.getCurrentSlice()-1;
			ImagePlus imTest = convertToFloatImage(IJ.openImage());
			System.out.println("caluclate score on "+imTest.getTitle()+"\nwith extracted feature of idx. " + mouseZ);
			imTest.show();
			new ImagePlus("max", ccf4.score(mouseZ, imTest.getProcessor())).show();
			new ImagePlus("per.rot", ccf4.scoreAllRot(mouseZ, imTest.getProcessor())).show();
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


}
