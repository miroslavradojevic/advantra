package detection2d;

import aux.ReadSWC;
import conn.Find_Connected_Regions;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.*;
import ij.io.FileSaver;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.io.*;
import java.util.ArrayList;

/**
 * Created by miroslav on 1/6/14.
 *
 */

// hint on zooming: imageplus.show(); imageplus.getWindow().setSize(200, 200); imageplus.getCanvas().fitToWindow();
public class Critpoint2D implements PlugIn, MouseListener, MouseMotionListener {

    String      image_path;

    /*
    interface elements - all the windows that pop up as you click/move with mouse
     */
    ImagePlus       pfl_im  = new ImagePlus();
    ImagePlus       pfl_im1 = new ImagePlus();
    ImagePlus       pfl_im2 = new ImagePlus();
    ImagePlus       pfl_im3 = new ImagePlus();
    ImagePlus       pfl_im4 = new ImagePlus();

    ImageStack      pfl_is  = null;
    ImageStack      pfl_is1 = null;
    ImageStack      pfl_is2 = null;
    ImageStack      pfl_is3 = null;
    ImageStack      pfl_is4 = null;

    ImageCanvas cnv;

	public void run(String sss){

        /*
        load the image
         */
        String in_folder = Prefs.get("id.folder", System.getProperty("user.home"));
        OpenDialog.setDefaultDirectory(in_folder);
        OpenDialog dc = new OpenDialog("Select file");
        in_folder = dc.getDirectory();
        image_path = dc.getPath();
        if (image_path==null) return;
        Prefs.set("id.folder", in_folder);

        ImagePlus ip_load = new ImagePlus(image_path);
        if(ip_load==null) return;


        /*
        load the detection parameters - generic dialog
         */
        boolean 	_show_bifpoints 		= Prefs.get("critpoint.detection2d.show_bifpoints", true);
        boolean 	_show_endpoints 		= Prefs.get("critpoint.detection2d.show_endpoints", true);
        boolean 	_enable_interactive     = Prefs.get("critpoint.detection2d.enable_interactive", true);
		float 		_s						= (float)	Prefs.get("critpoint.detection2d.s", 1.2f);
        String 		_Dlist 					= Prefs.get("critpoint.detection2d.d", "4");
        int 		_L						= (int) 	Prefs.get("critpoint.detection2d.l", 10);
        float 		_ncc_high               = (float)   Prefs.get("critpoint.detection2d.ncc_high", 1f);
        float 		_ncc_low				= (float) 	Prefs.get("critpoint.detection2d.ncc_low", 0.7f);
        float 		_likelihood_high        = (float)   Prefs.get("critpoint.detection2d.likelihood_high", 0.9f);
        float 		_likelihood_low			= (float) 	Prefs.get("critpoint.detection2d.likelihood_low", 0.0f);
        float 		_output_sigma			= (float) 	Prefs.get("critpoint.detection2d.output_sigma", 0.4f);

        GenericDialog gd = new GenericDialog("DETECTOR2D");

        gd.addCheckbox("BIFURCATIONS ", 		_show_bifpoints);
        gd.addCheckbox("ENDPOINTS    ", 		_show_endpoints);
        gd.addCheckbox("INTERACTIVE  ",         _enable_interactive);
		gd.addNumericField("s:", 				_s,					1,	10,	"(scale)");
        gd.addStringField("Dlist", 				_Dlist);
        gd.addNumericField("L",                 _L, 		        0,  10, "");
        gd.addNumericField("NCC_HIGH", 	        _ncc_high, 			2,  10, "");
        gd.addNumericField("NCC_LOW",           _ncc_low, 		    2,  10, "");
        gd.addNumericField("LIKELIHOOD_HIGH", 	_likelihood_high, 	2,  10, "");
        gd.addNumericField("LIKELIHOOD_LOW",    _likelihood_low, 	2,  10, "");
        gd.addNumericField("OUT",    			_output_sigma, 		2,  10, "SIGMA");

        gd.showDialog();
        if (gd.wasCanceled()) return;

        _show_bifpoints = gd.getNextBoolean();      	    	Prefs.set("critpoint.detection2d.show_bifpoints", 		_show_bifpoints);
        _show_endpoints = gd.getNextBoolean();      	    	Prefs.set("critpoint.detection2d.show_endpoints", 		_show_endpoints);
        _enable_interactive = gd.getNextBoolean();				Prefs.set("critpoint.detection2d.enable_interactive", 	_enable_interactive);
		_s          	= (float) gd.getNextNumber(); 			Prefs.set("critpoint.detection2d.s", _s);
        _Dlist       	= gd.getNextString(); 				    Prefs.set("critpoint.detection2d.d", _Dlist);
        _L		    	= (int) gd.getNextNumber();			    Prefs.set("critpoint.detection2d.l", _L);
        _ncc_high   	= (float) gd.getNextNumber();   	    Prefs.set("critpoint.detection2d.ncc_high", _ncc_high);
        _ncc_low    	= (float) gd.getNextNumber();  		    Prefs.set("critpoint.detection2d.ncc_low", _ncc_low);
        _likelihood_high= (float) gd.getNextNumber();   		Prefs.set("critpoint.detection2d.likelihood_high", 	_likelihood_high);
        _likelihood_low = (float) gd.getNextNumber();   		Prefs.set("critpoint.detection2d.likelihood_low", 	_likelihood_low);
        _output_sigma 	= (float) gd.getNextNumber();           Prefs.set("critpoint.detection2d.output_sigma", 	_output_sigma);

        String[] dd = _Dlist.split(",");
        float[] _D = new float[dd.length];
		for (int i=0; i<dd.length; i++) _D[i] = Float.valueOf(dd[i]);

		// detection
		Detector2D det2d = new Detector2D(
										ip_load,
										_s,
										_D,
										_L,
										_ncc_high,
										_ncc_low,
										_likelihood_high,
										_likelihood_low,
										_output_sigma
		);

		det2d.run();
		det2d.doEvaluation(); // append results if there was a ground truth
		Overlay ov_bifurcations = det2d.getBifurcationOverlay();
		Overlay ov_endpoints 	= det2d.getEndpointOverlay();

		Overlay ov_all = new Overlay();
		if (_show_endpoints)
			for (int i_ov=0; i_ov<ov_endpoints.size(); i_ov++) ov_all.add(ov_endpoints.get(i_ov));
		if (_show_bifpoints)
			for (int i_ov=0; i_ov<ov_bifurcations.size(); i_ov++) ov_all.add(ov_bifurcations.get(i_ov));

		ip_load.setOverlay(ov_all);

		// save as tif with detection overlay
		FileSaver fs = new FileSaver(ip_load);
		String save_path = det2d.output_dir_name + File.separator + det2d.image_name+".tif";
		fs.saveAsTiff(save_path);

		if (_enable_interactive) {
			ip_load.show();
			cnv = ip_load.getCanvas();
			cnv.addMouseListener(this);
			cnv.addMouseMotionListener(this);
			IJ.setTool("hand");
		}


	}

//        /*
//        if (false) {
//
//            IJ.log("export features ");
//            t1 = System.currentTimeMillis();
//
//            String export_path = image_dir+image_name+".feat";
//            Delineator2D.exportFeats(export_path);
//            IJ.log("done exporting to: \t" + export_path + "\n");
//
//            export_path = image_dir+image_name+".desc";
//            Delineator2D.exportDescripts(export_path);
//            IJ.log("done exporting to: \t" + export_path + "\n");
//
//            export_path = image_dir + image_name + ".frame";
//            Delineator2D.exportFrames(export_path);
//            IJ.log("done exporting to: \t" + export_path + "\n");
//
//            export_path = image_dir + image_name + ".i2xy";
//            Masker2D.exportI2xyCsv(export_path);
//            IJ.log("done exporting to: \t" + export_path + "\n");
//
//            t2 = System.currentTimeMillis();
//            IJ.log("done. " + ((t2 - t1) / 1000f) + "sec.");
//
//        }
//        */
//    }

    public void mouseClicked(MouseEvent e)
    {
//        mouseMoved(e);
        int clickX = cnv.offScreenX(e.getX());
        int clickY = cnv.offScreenY(e.getY());

		/*
        pfl_is1 = PeakAnalyzer2D.plotDelineationProfiles(clickX, clickY);
        pfl_im1.setStack("cross-profiles", pfl_is1);
        pfl_im1.show();
        */

        IJ.setTool("hand");

    }

    public void mouseMoved(MouseEvent e)
    {
        //mouseClicked(e);

        int clickX = cnv.offScreenX(e.getX());
        int clickY = cnv.offScreenY(e.getY());

        /*
            output Overlay & update canvas with the original
         */
        Overlay ov_to_add = FeatureExtractor2D.getDelineationOverlay(clickX, clickY);

        // add the whole image on top as overlay each time mouse is moved
        //ImageRoi fgroi = new ImageRoi(0, 0, Masker2D.getMask());    // !!! not very efficient to be done each time
        //ImageRoi simple_det_roi = new ImageRoi(0, 0, new ByteProcessor(inimg_xy.length, inimg_xy[0].length, SimpleDetsector2D.score2));

        //simple_det_roi.setOpacity(0.1);     // add foreground to know what was removed always
        //ov_to_add.add(simple_det_roi);
        //simple_det_roi.setOpacity(0.1);     // add foreground to know what was removed always
		//ov.add(simple_det_roi);

        cnv.setOverlay(ov_to_add);

		// print extracted features
        FeatureExtractor2D.print(clickX, clickY);

        /*
         output stack with local profile
        */
//        pfl_is3 = Profiler2D.getProfile(clickX, clickY);
//        pfl_im3.setStack("local_profile", pfl_is3);
//        pfl_im3.show();
		/*
			output stack: local profile with peaks detected
		 */

        pfl_is = FeatureExtractor2D.getFuzzyAgg(clickX, clickY);
        pfl_im.setStack(pfl_is);
        pfl_im.show();


        pfl_is4 = PeakExtractor2D.getProfileWithPeaks(clickX, clickY);
        pfl_im4.setStack("local_profile_with_peaks", pfl_is4);
        pfl_im4.show();

		pfl_is2 = FeatureExtractor2D.plotDelineationFeatures(clickX, clickY);
		pfl_im2.setStack("features", pfl_is2);
		pfl_im2.show();

		IJ.setTool("hand");

	}

    public void mousePressed(MouseEvent e)  {}
    public void mouseReleased(MouseEvent e) {}
    public void mouseEntered(MouseEvent e)  {}
    public void mouseExited(MouseEvent e)   {}
    public void mouseDragged(MouseEvent e)  {}

}