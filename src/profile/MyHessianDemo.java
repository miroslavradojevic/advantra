package profile;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.Line;
import ij.gui.Overlay;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import imagescience.image.*;
import profile.MyHessian;

import java.util.Vector;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 6/29/13
 * Time: 10:29 AM
 */
public class MyHessianDemo implements PlugInFilter {

    ImagePlus   img;
    Image       im;

    public int setup(String s, ImagePlus imagePlus) {

        if (imagePlus==null) return DONE;
        img   = imagePlus;
        im = new FloatImage(Image.wrap(img));
        return DOES_8G+DOES_32+NO_CHANGES;

    }

    public void run(ImageProcessor imageProcessor) {

		double scale = 2;

        // extract hessian eigen values and eigen vectors for different scales

		scale               = Prefs.get("advantra.critpoint.scale", 1.5);

		GenericDialog gd 	= new GenericDialog("NEURITENESS");
		gd.addNumericField("scale ",  scale, 1, 10, " pix");

		gd.showDialog();
		if (gd.wasCanceled()) return;

		scale               =  gd.getNextNumber();
		Prefs.set("advantra.critpoint.scale",   scale);

        Dimensions dims 		= im.dimensions();
        Dimensions new_dims 	= new Dimensions(img.getWidth(), img.getHeight(), 1);

        Image L1scales = new FloatImage(new_dims); L1scales.axes(Axes.X);
        Image L2scales = new FloatImage(new_dims); L2scales.axes(Axes.X);
        Image v1scales = new FloatImage(new_dims); v1scales.axes(Axes.X);
        Image v2scales = new FloatImage(new_dims); v2scales.axes(Axes.X);

        Image nness = new FloatImage(new_dims);     nness.axes(Axes.X);

        MyHessian my_hess = new MyHessian();

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

        // loop once more to set nness & create vector overlay
        hess.clear();
        hess = my_hess.eigs(im.duplicate(), scale, true);

        Image V11 	= hess.get(4); V11.axes(Axes.X);
        Image V12 	= hess.get(5); V12.axes(Axes.X);

        Overlay ov = new Overlay();
        coords.z = 0;
        for (coords.y=0; coords.y<dims.y; ++coords.y) {
            for (coords.x = 0; coords.x < dims.x; ++coords.x) {
                if (nness.get(coords)!=0) {
                    double value = nness.get(coords) / Lmin;
                    ov.add(new Line(coords.x+0.5, coords.y+0.5, coords.x+0.5+value*V11.get(coords), coords.y+0.5+value*V12.get(coords)));
                    nness.set(coords, value);
                }
            }
        }
//            L2.get(coords,aL2);
//            L1.get(coords,aL1);
//            for (coords.x = 0; coords.x < dims.x; ++coords.x) {
//                if (aL1[coords.x]<0) {
//                    nness.set(coords, aL1[coords.x]/Lmin);
//                }
//                else {
//                    nness.set(coords, 0);
//                }
//            }
        //}

        nness.name("nness,s="+IJ.d2s(scale,1));
        ImagePlus neuriteness = nness.imageplus();
        neuriteness.setOverlay(ov);
        neuriteness.show();

//      L1scales.name("L1");
//		L1scales.imageplus().show();
//      L2scales.name("L2");
//		L2scales.imageplus().show();
//		v1scales.imageplus().show();
//		v2scales.imageplus().show();
        // neuriteness calculation

    }

}
