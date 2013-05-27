package advantra.critpoint;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

public class DemoCircularConfigurations implements PlugIn {

	int         angScale1, angScale2;
	int[]       angScale;

	double      ring1, ring2;
	int         nr_ring;
	double[]    rings;

	double 		innerRadius;

	int 		patchRadius;

	// profile values
	float[] vals, angs, rads;

	FilterSet fs;

	public void run(String arg0) {

		CircularConfiguration3 ccf = new CircularConfiguration3(15);

		new ImagePlus("all_feat", ccf.plotKernels()).show();

		int chIdx = 3;
		new ImagePlus("feat.idx."+chIdx, ccf.plotKernel(chIdx)).show();

		// open image to convolve
		ImagePlus im = convertToFloatImage(new ImagePlus("/home/miroslav/tst.tif"));
		im.show();

		new ImagePlus("0,all", ccf.scoreAllRot(chIdx, im.getProcessor())).show();
		new ImagePlus("0,all,max", ccf.score(chIdx, im.getProcessor())).show();

		System.out.println("calculating all features...");
		ImageStack all_scores = new ImageStack(im.getWidth(), im.getHeight());
		for (int featIdx = 0; featIdx < ccf.kernels.size(); featIdx++ ) {
			all_scores.addSlice(ccf.score(featIdx, im.getProcessor()));
		}
		System.out.println("done.");
		new ImagePlus("all_scores", all_scores).show();


	}

	private ImagePlus convertToFloatImage(ImagePlus inim){

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


//GenericDialog gd = new GenericDialog("Demo filter");
//		gd.addMessage("WEAK LEARNERS (feature set)");
//		gd.addChoice("alfa_1", new String[]{"20", "40", "60", "80"}, "40");
//		gd.addChoice("alfa_2",      new String[]{"20", "40", "60", "80"}, "40");
//		gd.addNumericField("ring_1:", 0.4, 1); //ring1
//		gd.addNumericField("ring_2:",    0.7,       1); //ring2
//		gd.addNumericField("nr_rings:", 2, 0, 5, "");//nr_rings
//		gd.addNumericField("inner ring:",   0.2,    1);
//		gd.addNumericField("patch_size:", 20, 0, 5, "");
//
//		gd.showDialog();
//		if (gd.wasCanceled()) return;
//
//		switch ((int) gd.getNextChoiceIndex()) {
//			case 0: angScale1 = 20; 	break;
//			case 1: angScale1 = 40;		break;
//			case 2: angScale1 = 60;		break;
//			case 3: angScale1 = 80;		break;
//			default: angScale1 = 40; 	break;
//		}
//
//		switch ((int) gd.getNextChoiceIndex()) {
//			case 0: angScale2 = 20; 	break;
//			case 1: angScale2 = 40;		break;
//			case 2: angScale2 = 60;		break;
//			case 3: angScale2 = 80;		break;
//			default: angScale2 = 40;	break;
//		}
//
//		angScale = new int[(angScale2-angScale1)/20+1];
//		int cnt = 0;
//		for (int i = angScale1; i <= angScale2; i+=20) angScale[cnt++] = i;
//
//		ring1 = gd.getNextNumber();
//		ring2 = gd.getNextNumber();
//		nr_ring = (int) gd.getNextNumber();
//
//        System.out.println("nr_ring = "+ nr_ring);
//
//		rings = new double[nr_ring];
//		for (int i = 0; i <nr_ring; i++) rings[i] = (i == 0) ? ring1 : ring1 + i * ((ring2 - ring1) / (nr_ring - 1));
//
//        for (int i = 0; i <rings.length; i++) System.out.println("ring["+i+"]="+rings[i]);
//
//		innerRadius = gd.getNextNumber();
//
//		patchRadius = (int) gd.getNextNumber();
//
//		int toAlloc = Calc.circularProfileSize(patchRadius);
////		IJ.showMessage("profile will have "+toAlloc+" elements.");
//		vals = new float[toAlloc];
//		angs = new float[toAlloc];
//		rads = new float[toAlloc];
//
//		System.out.println("allocate "+toAlloc);
//
//		/*
//		 * generate filters to score on example profiles (generate features)
//		 */
//		fs = new FilterSet(angScale, rings, new double[]{0.4, 0.6, 0.8}, innerRadius);
//		int nrFilters = fs.circConfs.size()+fs.radlConfs.size();
//		System.out.println(nrFilters + " filters (weak classifiers) formed!");
//        /*
//        show them
//        */
//		int dispSize = 2*patchRadius+1;
//		all_feats = new ImagePlus("FEATURES", fs.plot(dispSize));
//		all_feats.setTitle("All_Features");
//		all_feats.show();
//
//		ImageCanvas all_feats_canvas = all_feats.getWindow().getCanvas();
//		all_feats_canvas.setName("feats");
//		all_feats_canvas.addMouseListener(this);
//
//		img.show();
//		ImageCanvas img_canvas = img.getWindow().getCanvas();
//		img_canvas.setName("image");
//		img_canvas.addMouseListener(this);
