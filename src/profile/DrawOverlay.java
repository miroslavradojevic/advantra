package profile;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Overlay;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 7/17/13
 * Time: 8:14 PM
 */
public class DrawOverlay implements PlugInFilter {

    ImagePlus inimg;

    public int setup(String s, ImagePlus imagePlus) {
        if(imagePlus==null) return DONE;
        inimg = imagePlus;//Tools.convertToFloatImage(imagePlus);
        return DOES_8G+DOES_32+NO_CHANGES;
    }

    public void run(ImageProcessor imageProcessor) {

        //inimg.show();
        if (inimg.getOverlay()!=null) {
            Overlay ov = inimg.getOverlay();
            IJ.log("has overlay -> "+ov.size()+" elements");
            for (int i=0; i<ov.size(); i++) {
                // draw
                ov.get(i).drawPixels(inimg.getProcessor());
            }
        }
        else {
            IJ.log("does not have overlay!");
        }

        //inimg.setHideOverlay(true);
        inimg.show();

    }
}
