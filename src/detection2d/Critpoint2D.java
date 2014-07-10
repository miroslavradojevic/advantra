package detection2d;

import aux.ReadSWC;
import conn.Find_Connected_Regions;
import ij.*;
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

	ImageWindow   	pfl_iw4 = null;
	ImageWindow   	pfl_iw3 = null;
	ImageWindow   	pfl_iw2 = null;

	boolean first4 = true;
	boolean first3 = true;
	boolean first2 = true;

    ImageCanvas cnv;

	int plot_w = 528;
	int plot_h = 255;
	int upper_left_x = 70;
	int upper_left_y = 50;

    Detector2D det2d = null;

	public void run(String sss){

        // load the image through the menu
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
        	load the detection parameters
         */
        boolean _show_junctions, _show_endpoints, _enable_interactive, _save_midresults, _auto_smoothness;
        float _s;
		float _sigma_ratio;
        String _Dlist="";
        float _ncc_high, _ncc_low, _likelihood_high, _likelihood_low, _smoothness_low, _smoothness_high;//, _output_sigma;

        // check if the parameters were submitted through the macro before rising up the Generic Dialog
        // enables calling plugin from the macro without opening the graphical window
        // useful to call plugins with parameters submitted by ij macro in fiji headless mode
        if (Macro.getOptions()==null) {

            // generic dialog (graphic)
            _show_junctions 		= 			Prefs.get("critpoint.detection2d.show_junctions", true);
            _show_endpoints 		= 			Prefs.get("critpoint.detection2d.show_endpoints", true);
            _enable_interactive     = 			Prefs.get("critpoint.detection2d.enable_interactive", true);
            _auto_smoothness        =           Prefs.get("critpoint.detection2d.auto_smoothness", false);
            _save_midresults        =           Prefs.get("critpoint.detection2d.save_midresults", false);
            _s						= (float)	Prefs.get("critpoint.detection2d.s", 1.2f);
            _Dlist 					= 			Prefs.get("critpoint.detection2d.d", "6");
			_sigma_ratio			= (float) 	Prefs.get("critpoint.detection2d.sigma_ratio", 		0.25f);
            _ncc_high               = (float)   Prefs.get("critpoint.detection2d.ncc_high", 		0.90f);
            _ncc_low				= (float) 	Prefs.get("critpoint.detection2d.ncc_low",	 		0.60f);
            _likelihood_high        = (float)   Prefs.get("critpoint.detection2d.likelihood_high", 	0.90f);
            _likelihood_low			= (float) 	Prefs.get("critpoint.detection2d.likelihood_low", 	0.30f);
            _smoothness_high        = (float)   Prefs.get("critpoint.detection2d.smoothness_high", 	10f);
            _smoothness_low         = (float)   Prefs.get("critpoint.detection2d.smoothness_low", 	20f);

            GenericDialog gd = new GenericDialog("DETECTOR2D");

            gd.addCheckbox("JUNCTIONS", 		_show_junctions);
            gd.addCheckbox("ENDPOINTS", 		_show_endpoints);
            gd.addCheckbox("INTERACTIVE",       _enable_interactive);
            gd.addCheckbox("AUTO_SMOOTHNESS",       _auto_smoothness);
            gd.addCheckbox("SAVE_MIDRESULTS",       _save_midresults);

            gd.addStringField("Dlist", 				_Dlist);
			gd.addNumericField("sigma_ratio",       _sigma_ratio, 	    2,  10, "");
			gd.addNumericField("s", 				_s,					1,	10,	"");
			gd.addNumericField("NCC_HIGH", 	        _ncc_high, 			2,  10, "");
            gd.addNumericField("NCC_LOW",           _ncc_low, 		    2,  10, "");
            gd.addNumericField("LIKELIHOOD_HIGH", 	_likelihood_high, 	2,  10, "");
            gd.addNumericField("LIKELIHOOD_LOW",    _likelihood_low, 	2,  10, "");
            gd.addNumericField("SMOOTHNESS_HIGH", 	_smoothness_high, 	2,  10, "");
            gd.addNumericField("SMOOTHNESS_LOW",    _smoothness_low, 	2,  10, "");

            gd.showDialog();
            if (gd.wasCanceled()) return;

            _show_junctions = gd.getNextBoolean();      	    	Prefs.set("critpoint.detection2d.show_junctions", 		_show_junctions);
            _show_endpoints = gd.getNextBoolean();      	    	Prefs.set("critpoint.detection2d.show_endpoints", 		_show_endpoints);
            _enable_interactive = gd.getNextBoolean();				Prefs.set("critpoint.detection2d.enable_interactive", 	_enable_interactive);
			_auto_smoothness =  gd.getNextBoolean();				Prefs.set("critpoint.detection2d.auto_smoothness", 		_auto_smoothness);
			_save_midresults = gd.getNextBoolean();					Prefs.set("critpoint.detection2d.save_midresults", 		_save_midresults);

			_Dlist       	= gd.getNextString(); 				    Prefs.set("critpoint.detection2d.d", _Dlist);
			_sigma_ratio	= (float) gd.getNextNumber();			Prefs.set("critpoint.detection2d.sigma_ratio", _sigma_ratio);
			_s          	= (float) gd.getNextNumber(); 			Prefs.set("critpoint.detection2d.s", _s);
            _ncc_high   	= (float) gd.getNextNumber();   	    Prefs.set("critpoint.detection2d.ncc_high", _ncc_high);
            _ncc_low    	= (float) gd.getNextNumber();  		    Prefs.set("critpoint.detection2d.ncc_low", _ncc_low);
            _likelihood_high= (float) gd.getNextNumber();   		Prefs.set("critpoint.detection2d.likelihood_high", 	_likelihood_high);
            _likelihood_low = (float) gd.getNextNumber();   		Prefs.set("critpoint.detection2d.likelihood_low", 	_likelihood_low);
            _smoothness_high= (float) gd.getNextNumber();   		Prefs.set("critpoint.detection2d.smoothness_high", 	_smoothness_high);
            _smoothness_low = (float) gd.getNextNumber();   		Prefs.set("critpoint.detection2d.smoothness_low", 	_smoothness_low);

        }
        else { // continue with macro arguments without rising graphic window

            _show_junctions = Boolean.valueOf(Macro.getValue(Macro.getOptions(), "junctions", String.valueOf(false)));
            _show_endpoints = Boolean.valueOf(Macro.getValue(Macro.getOptions(), "endpoints", String.valueOf(false)));
            _enable_interactive = Boolean.valueOf(Macro.getValue(Macro.getOptions(), "interactive", String.valueOf(false)));
            _auto_smoothness = Boolean.valueOf(Macro.getValue(Macro.getOptions(), "auto_smoothness", String.valueOf(false)));
            _save_midresults = Boolean.valueOf(Macro.getValue(Macro.getOptions(), "save_midresults", String.valueOf(false)));
			_Dlist = Macro.getValue(Macro.getOptions(), "Dlist", String.valueOf(4));
			_sigma_ratio = Float.valueOf(Macro.getValue(Macro.getOptions(), "sigma_ratio", String.valueOf(0.25)));
			_s = Float.valueOf(Macro.getValue(Macro.getOptions(), "s", String.valueOf(1.1)));

			_ncc_high           = Float.valueOf(Macro.getValue(Macro.getOptions(), "NCC_HIGH", String.valueOf(0.95)));
            _ncc_low            = Float.valueOf(Macro.getValue(Macro.getOptions(), "NCC_LOW", String.valueOf(0.2)));
            _likelihood_high    = Float.valueOf(Macro.getValue(Macro.getOptions(), "LIKELIHOOD_HIGH", String.valueOf(0.5)));
            _likelihood_low     = Float.valueOf(Macro.getValue(Macro.getOptions(), "LIKELIHOOD_LOW", String.valueOf(0.0)));
            _smoothness_high    = Float.valueOf(Macro.getValue(Macro.getOptions(), "SMOOTHNESS_HIGH", String.valueOf(10)));
            _smoothness_low     = Float.valueOf(Macro.getValue(Macro.getOptions(), "SMOOTHNESS_LOW", String.valueOf(20)));

        }

        String[] dd = _Dlist.split(",");
        float[] _D = new float[dd.length];
		for (int i=0; i<dd.length; i++) _D[i] = Float.valueOf(dd[i]);

		// detection
		System.out.print("Detector2D... ");
		det2d = new Detector2D(
										ip_load,
										_s,
										_D,
										_sigma_ratio,
										_ncc_high,
										_ncc_low,
										_likelihood_high,
										_likelihood_low,
                                        _smoothness_high,
                                        _smoothness_low
		);

		// they are initialized by default already , this just overwrites the
		det2d.saveMidresults(_save_midresults); 	 // t
        det2d.auto_smoothness = _auto_smoothness;    // f
		det2d.do_endpoints = _show_endpoints;        // t
		det2d.do_junctions = _show_junctions;        // t

		String save_path = "";
		if (true) {
			det2d.run();
			Overlay ov_critpoints = det2d.getDetectionOverlay();
			if (ov_critpoints.size()>0) ip_load.setOverlay(ov_critpoints);
			else System.out.print("Empty overlay with detections.");

			// save as tif with detection overlay (this will be input for evaluation, comparison with gndtth swc file)
			save_path = det2d.output_dir_name + File.separator + det2d.image_name+".tif";
			IJ.saveAs(ip_load, "Tiff", save_path);
			//System.out.println("exported: " + save_path);

			// save swc with dets
			save_path = det2d.output_dir_name + File.separator + det2d.image_name+".det";
			det2d.saveDetection(save_path);
			//System.out.println("exported: " + save_path);
		}
		System.out.println("DONE.");

        if (_enable_interactive) {
			ip_load.show();
			cnv = ip_load.getCanvas();
			cnv.addMouseListener(this);
			cnv.addMouseMotionListener(this);
			IJ.setTool("hand");
		}

	}

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

    public void mouseClicked(MouseEvent e)
    {
//        mouseMoved(e);
        int clickX = cnv.offScreenX(e.getX());
        int clickY = cnv.offScreenY(e.getY());


		pfl_is3 = Delineator2D.getClusterMask(clickX, clickY);
		if (pfl_is3.getSize()>0) {
            pfl_im3.setStack("cluster_mask", pfl_is3);
            if (first3) {
                pfl_im3.show();
                pfl_iw3 = pfl_im3.getWindow();
                pfl_iw3.setLocation(upper_left_x+plot_w+20, upper_left_y);
                first3 = false;
            }
            pfl_im3.setStack("fitting_scores", pfl_is3);
            //ImageWindow iw = pfl_im3.getWindow();
            pfl_im3.updateAndDraw();
//            pfl_im3.show();
        }



		pfl_is2 = (det2d!=null)? det2d.getFuzzyLogicResponse(clickX, clickY) : null;
        if (pfl_is2!=null && pfl_is2.getSize()>0) {
            pfl_im2.setStack("what happened in fuzzy", pfl_is2);
            if (first2) {
                pfl_im2.show();
                pfl_iw2 = pfl_im2.getWindow();
                pfl_iw2.setLocation(upper_left_x, upper_left_y+plot_h+50);
                first2 = false;
            }
            pfl_im2.setStack("what happened in fuzzy", pfl_is2);
//            ImageWindow iw = pfl_im2.getWindow();
            pfl_im2.updateAndDraw();
        }



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
        Overlay ov_to_add = Delineator2D.getDelineationOverlay(clickX, clickY);

        // add the whole image on top as overlay each time mouse is moved
        //ImageRoi fgroi = new ImageRoi(0, 0, Masker2D.getMask());    // !!! not very efficient to be done each time
        //ImageRoi simple_det_roi = new ImageRoi(0, 0, new ByteProcessor(inimg_xy.length, inimg_xy[0].length, SimpleDetsector2D.score2));

        //simple_det_roi.setOpacity(0.1);     // add foreground to know what was removed always
        //ov_to_add.add(simple_det_roi);
        //simple_det_roi.setOpacity(0.1);     // add foreground to know what was removed always
		//ov.add(simple_det_roi);

        cnv.setOverlay(ov_to_add);

        Detector2D.print(clickX, clickY); // print extracted features

/*
        pfl_is = FeatureExtractor2D.getFuzzyAgg(clickX, clickY);
		if (pfl_is.getSize()>0) {
        pfl_im.setStack("fuzzy_agg", pfl_is);
        pfl_im.show();
		}
*/

        pfl_is4 = PeakExtractor2D.getProfileWithPeaks(clickX, clickY);
		if (pfl_is4.getSize()>0) pfl_im4.setStack("local_profile_with_peaks", pfl_is4);
		if (first4) {
			pfl_im4.show();
			pfl_iw4 = pfl_im4.getWindow();
			pfl_iw4.setLocation(upper_left_x, upper_left_y);
			first4 = false;
		}

		IJ.setTool("hand");

	}

    public void mousePressed(MouseEvent e)  {}
    public void mouseReleased(MouseEvent e) {}
    public void mouseEntered(MouseEvent e)  {}
    public void mouseExited(MouseEvent e)   {}
    public void mouseDragged(MouseEvent e)  {}

}