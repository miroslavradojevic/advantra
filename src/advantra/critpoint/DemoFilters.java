package advantra.critpoint;

import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

import advantra.feature.FilterSet;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.measure.Calibration;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

public class DemoFilters implements PlugInFilter, MouseListener {

	/*
		creates filter set (set of weak classifiers, weak learners),
		visualises them and shows their scores
		on selected locations (mouse)
		it is possible to choose a feature
		and plot its score on the whole image

	 */

	ImagePlus img, all_feats;

	// parameters here
	int         angScale1, angScale2;
	int[]       angScale;

	double      ring1, ring2;
	int         nr_ring;
	double[]    rings;

	int 		patchRadius;

	// profile values
	float[] vals, angs, rads;

	FilterSet fs;

	public void run(ImageProcessor arg0) {
		
		GenericDialog gd = new GenericDialog("Demo filter");
		gd.addMessage("WEAK LEARNERS (feature set)");
		gd.addChoice("alfa_1",      new String[]{"20", "40", "60", "80"}, "20");
		gd.addChoice("alfa_2",      new String[]{"20", "40", "60", "80"}, "40");
		gd.addNumericField("ring_1:",    0.4,       1); //ring1
		gd.addNumericField("ring_2:",    0.7,       1); //ring2
		gd.addNumericField("nr_rings:", 2, 0, 5, "");//nr_rings
		gd.addNumericField("patch_size:", 15, 0, 5, "");
		gd.showDialog();
		if (gd.wasCanceled()) return;

		switch ((int) gd.getNextChoiceIndex()) {
			case 0: angScale1 = 20; 	break;
			case 1: angScale1 = 40;		break;
			case 2: angScale1 = 60;		break;
			case 3: angScale1 = 80;		break;
			default: angScale1 = 40; 	break;
		}

		switch ((int) gd.getNextChoiceIndex()) {
			case 0: angScale2 = 20; 	break;
			case 1: angScale2 = 40;		break;
			case 2: angScale2 = 60;		break;
			case 3: angScale2 = 80;		break;
			default: angScale2 = 40;	break;
		}

		angScale = new int[(angScale2-angScale1)/20+1];
		int cnt = 0;
		for (int i = angScale1; i <= angScale2; i+=20) angScale[cnt++] = i;

		ring1 = gd.getNextNumber();
		ring2 = gd.getNextNumber();
		nr_ring = (int) gd.getNextNumber();

		rings = new double[nr_ring];
		for (int i = 0; i <nr_ring; i++) rings[i] = (i == 0) ? ring1 : ring1 + i * ((ring2 - ring1) / (nr_ring - 1));

		patchRadius = (int) gd.getNextNumber();

		int toAlloc = Calc.circularProfileSize(patchRadius);
		vals = new float[toAlloc];
		angs = new float[toAlloc];
		rads = new float[toAlloc];

		System.out.println("allocate "+toAlloc);

		/*
		 * generate filters to score on example profiles (generate features)
		 */
		fs = new FilterSet(angScale, rings, new double[]{0.3, 0.5, 0.7});
		int nrFilters = fs.circConfs.size()+fs.radlConfs.size();
		System.out.println(nrFilters + " filters (weak classifiers) formed!");
        /*
        show them
        */
		int dispSize = 2*patchRadius+1;
		all_feats = new ImagePlus("FEATURES", fs.plot(dispSize));
		all_feats.setTitle("All_Features");
		all_feats.show();

		ImageCanvas all_feats_canvas = all_feats.getWindow().getCanvas();
		all_feats_canvas.setName("feats");
		all_feats_canvas.addMouseListener(this);

		img.show();
		ImageCanvas img_canvas = img.getWindow().getCanvas();
		img_canvas.setName("image");
		img_canvas.addMouseListener(this);


	}

	public int setup(String arg0, ImagePlus imp) {
		if(imp!=null) {
			img = imp;
					// reset calibration before going further 
					Calibration cal = new Calibration(img);
					cal.pixelWidth = cal.pixelHeight = cal.pixelDepth = 1;
					cal.setUnit("pixel");
					img.setCalibration(cal);
					 //new FloatImage(Image.wrap(imp));
		}
		return DOES_8G+NO_CHANGES;
	}

	public void mouseClicked(MouseEvent e) {

		ImageCanvas srcCanv = (ImageCanvas) e.getSource();

		String source  =  srcCanv.getName();

		if (source=="image"){

			int mouseX = img.getWindow().getCanvas().offScreenX(e.getX());
			int mouseY = img.getWindow().getCanvas().offScreenX(e.getY());
			System.out.println("X: " + mouseX + "," + "Y: " + mouseY);

			if(mouseX<patchRadius || mouseX>=img.getWidth()-patchRadius || mouseY<patchRadius || mouseY>=img.getHeight()-patchRadius )
				System.out.println("it's out, click again...");
			else {

				Overlay o = new Overlay();
				int N = Calc.circularProfileSize(patchRadius);
				int[] xloc = new int[N];
				int[] yloc = new int[N];
				Calc.getProfileLocations(mouseX, mouseY, patchRadius, xloc, yloc);
				for (int i = 0; i < N; i++) o.addElement(new PointRoi(xloc[i], yloc[i]));
				img.setOverlay(o);

				Calc.getProfile(img, mouseX, mouseY, patchRadius, vals, rads, angs);
				new ImagePlus("extracted_values", Calc.plotProfile(vals, rads, angs)).show();
				new ImagePlus("filter_responses", Calc.plotAllResponses(fs, vals, angs, rads)).show();
			}

		}

		if (source=="feats") {

			int mouseZ = all_feats.getCurrentSlice()-1;
			System.out.println("caluclate score on feature with idx. " + mouseZ);
			new ImagePlus("", Calc.filterResponse(img, fs, mouseZ, patchRadius)).show();

		}



	}

	public void mousePressed(MouseEvent e) {}

	public void mouseReleased(MouseEvent e) {}

	public void mouseEntered(MouseEvent e) {}

	public void mouseExited(MouseEvent e) {}

}
