import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Plot;
import ij.plugin.PlugIn;

import java.awt.*;

public class TestPlot implements PlugIn {

    public void run(String s) {


        double[] x = new double[]{1, 2, 3, 4, 5};
        double[] y1 = new double[]{4, 6, 9, 2, 5};
        double[] y2 = new double[]{6, 8, 1, 1, 1};

        Plot p = new Plot("name", "x", "y", x, y1);
        p.setColor(Color.BLACK);
        //p.draw();
        //p.setColor(Color.BLUE);
        p.addPoints(x, y2, Plot.BOX);

        ImageStack plotStack = new ImageStack(p.getProcessor().getWidth(), p.getProcessor().getHeight());
        plotStack.addSlice(p.getProcessor());

        p = new Plot("name", "x", "y", x, y1);
        p.setColor(Color.BLACK);
        p.draw();
        p.setColor(Color.BLUE);
        p.addPoints(x, y2, Plot.BOX);

        plotStack.addSlice(p.getProcessor());

        new ImagePlus("example", plotStack).show();

    }
}
