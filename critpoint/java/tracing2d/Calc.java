package tracing2d;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Line;
import ij.gui.Overlay;
import ij.plugin.ZProjector;
import ij.process.ByteProcessor;
import imagescience.image.*;
import imagescience.utility.ImageScience;

import java.util.Arrays;
import java.util.Vector;

/**
 * Created by miroslav on 31-8-14.
 * some handy image processing methods
 */
public class Calc {

    public static ImagePlus neuriteness(ImagePlus inimg, float[] scales)
    {

        Image im = new FloatImage(Image.wrap(inimg));
        Dimensions dims 		= im.dimensions();
        Dimensions new_dims 	= new Dimensions(inimg.getWidth(), inimg.getHeight(), 1);

        MyHessian my_hess = new MyHessian();

        ImageStack isout = new ImageStack(inimg.getWidth(), inimg.getHeight());

        for (int i = 0; i < scales.length; i++) {

            float scale = scales[i];

            Image L1scales = new FloatImage(new_dims); L1scales.axes(Axes.X);
            Image L2scales = new FloatImage(new_dims); L2scales.axes(Axes.X);
            Image v1scales = new FloatImage(new_dims); v1scales.axes(Axes.X);
            Image v2scales = new FloatImage(new_dims); v2scales.axes(Axes.X);

            Image nness = new FloatImage(new_dims);     nness.axes(Axes.X);

            double[] aL1 	= new double[dims.x];
            double[] aL2 	= new double[dims.x];
            double[] aV11 	= new double[dims.x];
            double[] aV12 	= new double[dims.x];

            double Lmin = Double.MAX_VALUE;

            Vector<Image> hess = my_hess.eigs(im.duplicate(), scale, false);

            // assign values to layers of L1, L2, v1, v2
            Image L2 	= hess.get(0); L2.axes(Axes.X); // higher
            Image L1 	= hess.get(1); L1.axes(Axes.X); // lower

            Coordinates coords 	= new Coordinates();

            coords.z = 0;
            for (coords.y=0; coords.y<dims.y; ++coords.y) {
                for (coords.x = 0; coords.x < dims.x; ++coords.x) {

                    if (Math.abs(L1.get(coords))>Math.abs(L2.get(coords))) {

                        // L1 to add
                        if (L1.get(coords)>=0) {
                            nness.set(coords, 0);
                        }
                        else {
                            nness.set(coords, L1.get(coords));
                            if (L1.get(coords)<Lmin) Lmin = L1.get(coords);
                        }

                    }
                    else {

                        // L2 to add
                        if (L2.get(coords)>=0) {
                            nness.set(coords, 0);
                        }
                        else {
                            nness.set(coords, L2.get(coords));
                            if (L2.get(coords)<Lmin) Lmin = L2.get(coords);
                        }

                    }
                }
            }

            for (coords.y=0; coords.y<dims.y; ++coords.y) {
                for (coords.x = 0; coords.x < dims.x; ++coords.x) {
                    if (nness.get(coords)!=0) {
                        double value = nness.get(coords) / Lmin;
//                        ov.add(new Line(coords.x+0.5, coords.y+0.5, coords.x+0.5+value*V11.get(coords), coords.y+0.5+value*V12.get(coords)));
                        nness.set(coords, value);
                    }
                }
            }



            // loop once more to set nness & create vector overlay
//            hess.clear();
//            hess = my_hess.eigs(im.duplicate(), scale, true);

//            nness.name();
            ImagePlus neuriteness = nness.imageplus();
//            neuriteness.show();

            isout.addSlice("nness,s="+ IJ.d2s(scale, 1), neuriteness.getProcessor());

        }

        ZProjector zproj = new ZProjector(new ImagePlus("", isout));
        zproj.setMethod(ZProjector.MAX_METHOD);
        zproj.doProjection();
        ImagePlus imout = zproj.getProjection();

//        IJ.run(new ImagePlus("", isout), "Z Project...", "projection=[Max Intensity]");
//        ImagePlus imout = IJ.getImage();
        imout.setTitle("Neuriteness"+ Arrays.toString(scales));
        return imout;

    }

}
