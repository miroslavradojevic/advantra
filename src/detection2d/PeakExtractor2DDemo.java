package detection2d;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.*;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 12/17/13
 * Time: 4:10 PM
 * Demo usage of PeakExtractor2D: uses also Masker2D, Profiler2D, and shows peaks extracted on all foreground locations as a joint Overlay on the image
 */
public class PeakExtractor2DDemo implements PlugInFilter, MouseListener, MouseMotionListener {

    float[][] 	inimg_xy; // store input image as an array

    /*
    parameters
     */
    float       iDiff, D, s=1.5f; // maskNhoodRadius
    int         CPU_NR;

    /*
    interface things
     */
	ImagePlus       pfl_im = new ImagePlus();
	ImageProcessor  pfl_ip = null;
	ImageWindow 	pfl_iw;
	ImageCanvas 	cnv;

	public int setup(String s, ImagePlus imagePlus) {

        if(imagePlus==null) return DONE;

        inimg_xy = new float[imagePlus.getWidth()][imagePlus.getHeight()]; // x~column, y~row

        if (imagePlus.getType()== ImagePlus.GRAY8) {
            byte[] read = (byte[]) imagePlus.getProcessor().getPixels();
            for (int idx=0; idx<read.length; idx++) {
                inimg_xy[idx%imagePlus.getWidth()][idx/imagePlus.getWidth()] = (float) (read[idx] & 0xff);
            }

        }
        else if (imagePlus.getType()==ImagePlus.GRAY32) {
            float[] read = (float[]) imagePlus.getProcessor().getPixels();
            for (int idx=0; idx<read.length; idx++) {
                inimg_xy[idx%imagePlus.getWidth()][idx/imagePlus.getWidth()] = read[idx];
            }
        }
        else {
            IJ.log("image type not recognized");
            return DONE;
        }

        /******************************
         Generic Dialog
         *****************************/
        //this.maskNhoodRadius             	= (float)   Prefs.get("advantra.critpoint.mask.nhoodRadius", 5);
        this.iDiff 					    	= (float)   Prefs.get("advantra.critpoint.mask.iDiff", 5);
        this.D 					        	= (float)   Prefs.get("advantra.critpoint.profile.d", 4);

        GenericDialog gd = new GenericDialog("PROFILER2DDEMO");
        //gd.addNumericField("radius ", 	maskNhoodRadius, 	0, 10, "spatial neighbourhood");
        gd.addNumericField("iDiff ", 	iDiff, 			0, 10, "intensity margin");
        gd.addNumericField("D ", 	D, 			0, 10, "neuron diameter[pix]");

        gd.showDialog();
        if (gd.wasCanceled()) return DONE;

        //maskNhoodRadius      = (float) gd.getNextNumber();
        //Prefs.set("advantra.critpoint.mask.nhoodRadius",	maskNhoodRadius);

        iDiff       	= (float) gd.getNextNumber();
        Prefs.set("advantra.critpoint.mask.iDiff", 	    	iDiff);

        CPU_NR = Runtime.getRuntime().availableProcessors();

        cnv = imagePlus.getCanvas();

        return DOES_8G+DOES_32+NO_CHANGES;

    }

    public void run(ImageProcessor imageProcessor) {

//		System.out.println(String.format("range: %4d", (int)(6.5f/0.7f)));

        Sphere2D sph2d = new Sphere2D(D, s);
        ImagePlus samplingScheme =  sph2d.showSampling();
        samplingScheme.show();
        samplingScheme.getWindow().setSize(600, 600);
        samplingScheme.getCanvas().fitToWindow();
        sph2d.showWeights().show();

        /*
        main
         */
        long t1, t2;
        t1 = System.currentTimeMillis();

        Masker2D.loadTemplate(inimg_xy, 0, sph2d.getOuterRadius(), iDiff);
        int totalLocs = inimg_xy.length * inimg_xy[0].length;

        Masker2D ms_jobs[] = new Masker2D[CPU_NR];
        for (int i = 0; i < ms_jobs.length; i++) {
            ms_jobs[i] = new Masker2D(i*totalLocs/CPU_NR,  (i+1)*totalLocs/CPU_NR);
            ms_jobs[i].start();
        }
        for (int i = 0; i < ms_jobs.length; i++) {
            try {
                ms_jobs[i].join();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }

        Masker2D.formRemainingOutputs();

        Profiler2D.loadTemplate(sph2d, Masker2D.i2xy, inimg_xy);
        int totalProfileComponents = sph2d.getProfileLength();

        Profiler2D pf_jobs[] = new Profiler2D[CPU_NR];
        for (int i = 0; i < pf_jobs.length; i++) {
            pf_jobs[i] = new Profiler2D(i*totalProfileComponents/CPU_NR,  (i+1)*totalProfileComponents/CPU_NR);
            pf_jobs[i].start();
        }
        for (int i = 0; i < pf_jobs.length; i++) {
            try {
                pf_jobs[i].join();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }

        PeakExtractor2D.loadTemplate(sph2d, Masker2D.i2xy, Profiler2D.prof2, inimg_xy, Masker2D.xy2i);
        //System.out.println("managed to initialize");
        int totalPeakExtrComponents = Profiler2D.prof2.length; // number of profiles == number of locations

        PeakExtractor2D pe_jobs[] = new PeakExtractor2D[CPU_NR];
        for (int i = 0; i < pe_jobs.length; i++) {
            pe_jobs[i] = new PeakExtractor2D(i*totalPeakExtrComponents/CPU_NR, (i+1)*totalPeakExtrComponents/CPU_NR);
            pe_jobs[i].start();
        }
        for (int i = 0; i < pe_jobs.length; i++) {
            try {
                pe_jobs[i].join();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }


        t2 = System.currentTimeMillis();
		        System.out.println("DEBUG: finished peak extraction");
        IJ.log("done. "+((t2-t1)/1000f)+"sec.");

		/*
			mouse listeners after processing
		 */
		cnv.addMouseListener(this);
		cnv.addMouseMotionListener(this);

        /* DEBUG: show the peaks (all peaks in one overlay at the current image)
        Overlay ov = new Overlay();
        float R_MAX = 1; // ~ corresponds highest weight
        float R = R_MAX;
        for (int l = 0; l<PeakExtractor2D.peaks_xy.length; l++) { // loop locations

            for (int t = 0; t<PeakExtractor2D.peaks_xy[l].length; t++) {  // loop threads

                int pk_x = PeakExtractor2D.peaks_xy[l][t][0];
                int pk_y = PeakExtractor2D.peaks_xy[l][t][1];

                if (pk_x != -1 && pk_y != -1) {
                    // R ~ size - weight of the peak
                    OvalRoi ovroi = new OvalRoi(pk_x+.5 - (R/2), pk_y+.5 - (R/2), R, R);
					ovroi.setFillColor(Color.BLUE);
                    ov.add(ovroi);

                }

            }

        }

        cnv.getImage().setOverlay(ov);
        */

    }

	public void mouseClicked(MouseEvent e) {

		int clickX = cnv.offScreenX(e.getX());
		int clickY = cnv.offScreenY(e.getY());

//		IJ.log(String.format("clicked:\t %4d, %4d", clickX, clickY));

		pfl_ip = Profiler2D.getProfile(clickX, clickY, Masker2D.xy2i);

		if (pfl_ip == null) {
			pfl_ip = new ShortProcessor(1, 1);
			pfl_im.setTitle("background");
		}
		else {
			pfl_im.setTitle("foreground");
		}

		pfl_im.setProcessor(pfl_ip);

		if (pfl_iw==null) {
			pfl_iw = new ImageWindow(pfl_im);
		}
		else {
			pfl_iw.setImage(pfl_im);
		}

		//pfl_iw.setSize(600, 300);
		//pfl_iw.getCanvas().fitToWindow();

		// add an overlay with peaks
		Overlay ov = new Overlay();
		float R = 2;
		OvalRoi ovalroi = new OvalRoi(clickX-(R/2)+.5f, clickY-(R/2)+.5f, R, R);
		//ovalroi.setFillColor(Color.BLUE);
		ov.add(ovalroi);

		// read extracted peaks at this location
		int idx = Masker2D.xy2i[clickX][clickY];
		if (idx!=-1) {
			int[][] pk_locs_xy = PeakExtractor2D.peaks_xy[idx];
			for (int i = 0; i<pk_locs_xy.length; i++) {

				int pk_x = pk_locs_xy[i][0];
				int pk_y = pk_locs_xy[i][1];

				if (pk_x!=-1) {
					// peak exists
					ovalroi = new OvalRoi(pk_x-(R/2)+.5f, pk_y-(R/2)+.5f, R, R);
					Color c = Color.BLACK;
					if(i==0) {
						c = Color.RED;
					}
					else if (i==1) {
						c = Color.YELLOW;
					}
					else if (i==2) {
						c = Color.GREEN;
					}
					else if (i==3) {
						c = Color.BLUE;
					}
					ovalroi.setFillColor(c);
				}

				ov.add(ovalroi);
			}

		}

		cnv.setOverlay(ov);

	}

	public void mousePressed(MouseEvent e) {}
	public void mouseReleased(MouseEvent e) {}
	public void mouseEntered(MouseEvent e) {}
	public void mouseExited(MouseEvent e) {}
	public void mouseDragged(MouseEvent e) {}
	public void mouseMoved(MouseEvent e) {}
}
