package detection3d;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.ImageCanvas;
import ij.gui.ImageWindow;
import ij.gui.Overlay;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;

import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 11/22/13
 * Time: 6:03 AM
 */
public class JunctionDetection3D implements PlugInFilter, MouseListener, MouseMotionListener {

	ImagePlus       imp;
	ImageCanvas     cnv;

    ImagePlus       pfl_im;
    ImageProcessor  pfl_ip;
    ImageWindow     pfl_iw;

    Detector3D det3D;

	// parameters necessary for the detection
	float 	D;
	float 	iDiff;

//	float 	minCosAngle;
//	float 	minFuzzyScore;
//	int 	MIN_SIZE;
//	float 	scatterDistSquared = 5f;
//	static float 	 	LOCATION_TOLERANCE_SCALE 	= 1.8f;

	private static float 	Deg2Rad = (float) (Math.PI/180f);
	private static float 	Rad2Deg = (float) (180f/Math.PI);

	// in case it is necessary to downsample
	boolean downsample = false;
	float   neuriteDiameter = Float.NaN;

    public int setup(String s, ImagePlus imagePlus) {
		if(imagePlus==null) return DONE;
		imp = imagePlus;
        cnv = imagePlus.getCanvas();
        return DOES_8G+NO_CHANGES;
    }

    public void run(ImageProcessor imageProcessor) {

		float zDist     = 3.0f;
		float D         = 4;
		float iDiff     = 10;

		// uses Detector3D
		det3D = new Detector3D();
		det3D.run(imp, zDist, D, iDiff);

		System.out.println("DONE!");

        cnv.addMouseListener(this);
        cnv.addMouseMotionListener(this);

    }

    public static void main(String[] args){

        ImagePlus inimg  = new ImagePlus("/home/miroslav/test11.tif");
        float zDist     = 3.4f;
        float D         = 3;
        float iDiff     = 10;

        Detector3D det3D = new Detector3D();
        det3D.run(inimg, zDist, D, iDiff);

		IJ.log("DONE!");

    }

    public void mouseClicked(MouseEvent e) {

        int offscreenX = cnv.offScreenX(e.getX());  // zero indexed location x y z
        int offscreenY = cnv.offScreenY(e.getY());
        int offscreenZ = cnv.getImage().getZ()-1;

        // overlay the skeleton
        Overlay skeleton = det3D.getLocalSkeleton(offscreenX, offscreenY, offscreenZ);
        cnv.getImage().setOverlay(skeleton);

        // update the profile image
        if (pfl_iw==null) {
            pfl_ip = new ShortProcessor(1, 1);
            pfl_im = new ImagePlus("", pfl_ip);
            pfl_iw = new ImageWindow(pfl_im);
            pfl_iw.setSize(600, 300);
        }
        else { // it exists, just update it
            pfl_ip = det3D.getLocalProfile(offscreenX, offscreenY, offscreenZ);

            if (pfl_ip==null) {
                pfl_ip = new ShortProcessor(1, 1);
            }

            pfl_im.setProcessor(pfl_ip);
            pfl_iw.updateImage(pfl_im);
            pfl_iw.setSize(600, 300);
            pfl_iw.getCanvas().fitToWindow();//setMagnification(8);

        }

        // plot detected & selected peaks on the profile
        Overlay local_profile_peaks = det3D.getLocalProfilePeaks(offscreenX, offscreenY, offscreenZ);
        Overlay selected_local_profile_peaks = det3D.getSelectedLocalProfilePeaks(offscreenX, offscreenY, offscreenZ);

        Overlay joint = new Overlay();
        for (int a=0; a<local_profile_peaks.size(); a++) {
            joint.add(local_profile_peaks.get(a));
        }
        for (int a=0; a<selected_local_profile_peaks.size(); a++) {
            joint.add(selected_local_profile_peaks.get(a));
        }

        if (joint.size()>0) {
            pfl_iw.getCanvas().getImage().setOverlay(joint);
        }

        det3D.debug(offscreenX, offscreenY, offscreenZ);

    }

    public void mouseMoved(MouseEvent e) {
        //mouseClicked(e);
    }

    @Override
    public void mousePressed(MouseEvent e) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void mouseReleased(MouseEvent e) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void mouseEntered(MouseEvent e) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void mouseExited(MouseEvent e) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void mouseDragged(MouseEvent e) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

}
