import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Plot;
import ij.plugin.PlugIn;

import java.awt.*;
import java.util.Arrays;

public class TestPlot implements PlugIn {

	double[] start = new double[10];
	double[] finish = new double[10];

    public void run(String s) {

		IJ.log("start before: "+Arrays.toString(start));
		IJ.log("finish before: "+Arrays.toString(finish));


        start[0] = 1;

		fn(start, finish);

        IJ.log("start after: "+Arrays.toString(start));
        IJ.log("finish after: "+Arrays.toString(finish));

        System.out.println("value: "+Double.MIN_VALUE);

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
