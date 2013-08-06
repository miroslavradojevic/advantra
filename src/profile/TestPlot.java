package profile;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.plugin.PlugIn;

import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;

public class TestPlot implements PlugIn {

//	double[] start = new double[10];
//	double[] finish = new double[10];

    public void run(String s) {

		GenericDialog gd = new GenericDialog("POSITIVE?");
		gd.enableYesNoCancel();
		gd.showDialog();
		if (gd.wasCanceled()) return;
		if (gd.wasOKed()) {
			IJ.log("YES");
		}
		else {
			IJ.log("NO");
		}

/*        int a = 21;
        IJ.log("before wrap: "+a);
        a = Tools.wrap(a, 20);
        IJ.log("after wrap: "+a);


        float[] x = new float[]{Float.NaN, Float.NaN, Float.NaN, 1, 2, 3};
        float[] y = new float[]{Float.NaN, Float.NaN, Float.NaN, 1, 2, 3};

        Plot p = new Plot("", "", "", x, y);

        PlotWindow pw = p.show();*/


/*

		int N = 4;
		int k = 3;

		ArrayList<int[]> c = Tools.comb(N, k);

		IJ.log("done, "+c.size()+" combinations");

		for (int i =0; i< c.size(); i++) {
			IJ.log("comb. " + i + " :" + Arrays.toString(c.get(i)) );
		}
*/



//        start[0] = 1;

//		fn(start, finish);

//        IJ.log("start after: "+Arrays.toString(start));
//        IJ.log("finish after: "+Arrays.toString(finish));

//        double[] x = new double[]{1, 2, 3, 4, 5};
//        double[] y1 = new double[]{4, 6, 9, 2, 5};
//        double[] y2 = new double[]{6, 8, 1, 1, 1};
//
//        Plot p = new Plot("name", "x", "y", x, y1);
//        p.setColor(Color.BLACK);
//        //p.draw();
//        //p.setColor(Color.BLUE);
//        p.addPoints(x, y2, Plot.BOX);
//
//        ImageStack plotStack = new ImageStack(p.getProcessor().getWidth(), p.getProcessor().getHeight());
//        plotStack.addSlice(p.getProcessor());
//
//        p = new Plot("name", "x", "y", x, y1);
//        p.setColor(Color.BLACK);
//        p.draw();
//        p.setColor(Color.BLUE);
//        p.addPoints(x, y2, Plot.BOX);
//
//        plotStack.addSlice(p.getProcessor());
//
//        new ImagePlus("example", plotStack).show();

    }

	public static void fn (double[] start, double[] finish) {
		start[0]++;
        for (int i=0; i<start.length; i++) {
			finish[i] = start[i] + 4;
		}
	}

}
