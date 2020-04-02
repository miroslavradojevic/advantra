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
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;

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

		float zDist     = 3.4f;
		float D         = 4;
		float iDiff     = 10;

		// uses Detector3D
		det3D = new Detector3D();
//		ArrayList<float[]>
		float[][] det_list_xyzr =  det3D.run(imp, zDist, D, iDiff);
		if (det_list_xyzr!=null) {
//			System.out.println(det_list_xyzr.length + " els.");
			exportSwcWithLocs("/home/miroslav/test_output.swc", det_list_xyzr);
			for (int w=0; w<det_list_xyzr.length; w++) {
//				System.out.println(Arrays.toString(det_list_xyzr[w]));
			}
		}

//		cnv.getImage().setOverlay(det3D.exportOverlayForegroundLocations());
		//cnv.getImage().setOverlay(det3D.exportOverlayPeakLocations()); // debug only

        cnv.addMouseListener(this);
        cnv.addMouseMotionListener(this);

    }

    public static void main(String[] args){

        ImagePlus inimg  = new ImagePlus("/home/miroslav/test11.tif");
        float zDist     = 3.4f;
        float D         = 4;
        float iDiff     = 10;

        Detector3D det3D = new Detector3D();
		float[][] det_list_xyzr =  det3D.run(inimg, zDist, D, iDiff);

    }

	public void exportSwcWithLocs(String out_path, float[][] _xyzr_list){ // ArrayList<float[]>

		PrintWriter     logWriter = null;

		try {
			logWriter = new PrintWriter(out_path);
			logWriter.print("");
			logWriter.close();
		}
		catch (FileNotFoundException ex) {}

		try {
			logWriter = new PrintWriter(new BufferedWriter(new FileWriter(out_path, true)));
			logWriter.println("# SWC WITH LOCATIONS ");
		} catch (IOException e) {}

		for (int l=0; l<_xyzr_list.length; l++) {
			logWriter.println((l+1)+" "+0+" "+_xyzr_list[l][0]+" "+_xyzr_list[l][1]+" "+_xyzr_list[l][2]+" "+_xyzr_list[l][3]+" "+(-1));

		}


		logWriter.close();

	}

    public void mouseClicked(MouseEvent e) {

        int offscreenX = cnv.offScreenX(e.getX());  // zero indexed location x y z
        int offscreenY = cnv.offScreenY(e.getY());
        int offscreenZ = cnv.getImage().getZ()-1;

        // overlay the skeleton
        //Overlay skeleton = det3D.getLocalSkeleton(offscreenX, offscreenY, offscreenZ);
        //cnv.getImage().setOverlay(skeleton);

        pfl_ip = det3D.getLocalProfile(offscreenX, offscreenY, offscreenZ);

        if (pfl_ip==null) {
			IJ.log("point in backgr");
            pfl_ip = new ShortProcessor(1, 1);
			pfl_im = new ImagePlus("BACK", pfl_ip);

        }
		else {
			pfl_im = new ImagePlus("FGR", pfl_ip);
//				pfl_im.setProcessor(pfl_ip);
		}

		if (pfl_iw==null) {
			pfl_iw = new ImageWindow(pfl_im);
		}
		else {
//			pfl_im.show();
			pfl_iw.setImage(pfl_im);
			pfl_iw.setSize(600, 300);
			pfl_iw.getCanvas().fitToWindow();
//			pfl_iw.updateImage(pfl_im);
		}

		/*
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
		*/


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
