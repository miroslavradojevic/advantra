package detection2d;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 * Created by miroslav on 1/10/14.
 * simplified version of PeakAnalyzer2DDemo stripped down to extract features, no debug things and user interface
 * input is an image, the outputs are exported into output files:
 *
 * image.feat   (csv) with extracted features: nr. locations x 7 features,
 * image.feat.description (txt) file with the description of the features used
 * image.feat.params (txt) extraction log, parameters used for feature extraction
 * image.i2xy   (csv) with index (i) to zero-indexed coordinate (x,y) mapping, lookup table
 * image.xy2i   (csv) with zero-indexed coordinate (x,y) to index (i) mapping, lookup table
 * image.frame  (csv) with frame map for each location (4xM index map describing the streamlines going from every location)
 * image.back.tif (tif image) estimated background
 * image.mask.tif (tif image) selection of foreground pixels
 * image.sampling.tif (tif image) imagej visualization of the sampling scheme for filtering
 * image.weights.tif (tif image) weights of the filter profile
 *
 * this is designed to be ij plugin, can be automatically called using macros, outputs have
 * the same name as the input image and are stored in the same folder
 */

public class FeatureExtractor2D implements PlugInFilter {

    float[][] 	inimg_xy;               // store input image as an array
    String      image_name;
    String      image_dir;
    String      output_dir;             // parameter signature will be in folder name

    /*
    feature extraction parameters
     */
    float       iDiff, D, minCos, scatterDist;
    int         M = 2;
    float       s = 1.5f;
	float 		threshold; 				// used for feature score normalization

    int         CPU_NR;

    public int setup(String str, ImagePlus imagePlus) {

        if(imagePlus==null) return DONE;

        image_name = imagePlus.getShortTitle();
        image_dir  = imagePlus.getOriginalFileInfo().directory + image_name;

        File image_dir_f = new File(image_dir);
        if (!image_dir_f.exists()) {
            boolean img_dir_created = new File(image_dir).mkdirs();
            if (!img_dir_created) return DONE;
        }

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
        this.iDiff 					    	= (float)   Prefs.get("advantra.critpoint.mask.iDiff", 10);
        this.D 					        	= (float)   Prefs.get("advantra.critpoint.profile.d", 4);           // scale
        this.M 					        	= (int)     Prefs.get("advantra.critpoint.analyze.m", 2);
        this.minCos 					    = (float)   Prefs.get("advantra.critpoint.analyze.min_cos", 0.2f);
        this.scatterDist                    = (float)   Prefs.get("advantra.critpoint.analyze.scatter_d", 2*this.D);
        this.s                              = (float)   Prefs.get("advantra.critpoint.profile.s", 1.5f);
		this.threshold						= (float) 	Prefs.get("advantra.critpoint.analyze.threshold", 10f);

        GenericDialog gd = new GenericDialog("FEATURE_EXTRACTION");
        gd.addNumericField("idiff", 	iDiff, 			1, 20, "intensity margin");
        gd.addNumericField("d", 	    D, 			    1, 20, "neuron diameter[pix]");
        gd.addNumericField("m", 	    M, 			    0, 20, "branch len.");
        gd.addNumericField("mincos", 	minCos, 	    1, 20, "alignment parameter.");
        gd.addNumericField("scatter_d",scatterDist, 	1, 20, "max allowed scatter (robustness test)");
        gd.addNumericField("scale",    s, 	            1, 20, "profiler scale");
		gd.addMessage("feat. calculation");
		gd.addNumericField("threshold ",threshold, 		1, 10, "intensity margin");

        gd.showDialog();
        if (gd.wasCanceled()) return DONE;

        iDiff       	= (float) gd.getNextNumber();
        Prefs.set("advantra.critpoint.mask.iDiff", 	    	iDiff);

        D       	    = (float) gd.getNextNumber();
        Prefs.set("advantra.critpoint.profile.d", 	    	D);

        M       	    = (int) gd.getNextNumber();
        Prefs.set("advantra.critpoint.analyze.m", 	    	M);

        minCos       	= (float) gd.getNextNumber();
        Prefs.set("advantra.critpoint.analyze.min_cos", 	minCos);

        scatterDist     = (float) gd.getNextNumber();
        Prefs.set("advantra.critpoint.analyze.scatter_d", 	scatterDist);

        s               = (float) gd.getNextNumber();
        Prefs.set("advantra.critpoint.profile.s", 	        s);

		threshold		= (float) gd.getNextNumber();
		Prefs.set("advantra.critpoint.analyze.threshold", 	threshold);

        output_dir =
                image_dir + File.separator +
                "idiff_"+Float.toString(iDiff)+"_"+
                "D_"+Float.toString(D)+"_"+
                "M_"+Float.toString(M)+"_"+
                "minCos_"+Float.toString(minCos)+"_"+
                "scatterDist_"+Float.toString(scatterDist)+"_"+
                "s_"+Float.toString(s)+File.separator
        ;

        boolean delete_success = deleteDir(new File(output_dir));
        if (!delete_success) { IJ.log("delete " + output_dir + " failed."); return DONE; }

        boolean out_dir_created = new File(output_dir).mkdirs();
        if (!out_dir_created) { IJ.log("creating " + output_dir + " failed."); return DONE; }

        CPU_NR = Runtime.getRuntime().availableProcessors();

        return DOES_8G+DOES_32+NO_CHANGES;

    }

    public void run(ImageProcessor imageProcessor) {

        Sphere2D sph2d = new Sphere2D(D, s);

		IJ.log("start extraction...");

		long t1, t2;
        t1 = System.currentTimeMillis();

        /********************************/
        Masker2D.loadTemplate(inimg_xy, 0, sph2d.getOuterRadius());
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
		new ImagePlus("MASK", Masker2D.getMask()).show();
        /********************************/
        Profiler2D.loadTemplate(sph2d, Masker2D.i2xy, Masker2D.xy2i, inimg_xy);
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
        /********************************/
        PeakExtractor2D.loadTemplate(sph2d, Masker2D.i2xy, Profiler2D.prof2, inimg_xy, Masker2D.xy2i);
        int totalPeakExtrComponents = Profiler2D.prof2.length; // number of profiles == number of locations i2xy.length

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
        /********************************/
        PeakAnalyzer2D.loadTemplate(Masker2D.i2xy, Masker2D.xy2i, PeakExtractor2D.peaks_i, PeakExtractor2D.peaks_w, inimg_xy, Masker2D.back_xy, M, minCos, scatterDist, threshold, D); // initialize peak analyzer parameters
        int totalPeakAnalyzeComponents = Masker2D.i2xy.length; // number of locations

        PeakAnalyzer2D pa_jobs[] = new PeakAnalyzer2D[CPU_NR];
        for (int i = 0; i < pa_jobs.length; i++) {
            pa_jobs[i] = new PeakAnalyzer2D(i*totalPeakAnalyzeComponents/CPU_NR, (i+1)*totalPeakAnalyzeComponents/CPU_NR);
            pa_jobs[i].start();
        }
        for (int i = 0; i < pa_jobs.length; i++) {
            try {
                pa_jobs[i].join();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }
        /********************************/

        t2 = System.currentTimeMillis();
        IJ.log(((t2-t1)/1000f)+"sec.");

		IJ.log("exporting...");
        PeakAnalyzer2D.exportFeats(      output_dir+image_name+".feat");
		PeakAnalyzer2D.exportDescripts(	 output_dir+image_name+".desc");
        PeakAnalyzer2D.exportFeatsLegend(output_dir+image_name+".feat.legend");
        exportExtractionLegend(          output_dir+image_name+".feat.params");
        Masker2D.exportI2xyCsv(          output_dir+image_name+".i2xy"); // export csv lookup table  from Masker2D.i2xy to image_name.i2xy
        Masker2D.exportXy2iCsv(			 output_dir+image_name+".xy2i"); // export csv lookup from Masker2D.xy2i to image_name.xy2i
		PeakAnalyzer2D.exportFrames(     output_dir+image_name+".frame");
        Masker2D.exportEstBackground(    output_dir+image_name+".back.tif");
		Masker2D.exportForegroundMask(   output_dir+image_name+".mask.tif");
		sph2d.exportSampling(            output_dir+image_name+".sampling.tif");
		sph2d.exportWeights(             output_dir+image_name+".weights.tif");
		IJ.log("done.");

    }

    private void exportExtractionLegend(String file_path) {

        PrintWriter logWriter = null; //initialize writer

        try {
            logWriter = new PrintWriter(file_path);
            logWriter.print("");
            logWriter.close();
        } catch (FileNotFoundException ex) {}   // empty the file before logging...

        try {                                   // initialize detection log file
            logWriter = new PrintWriter(new BufferedWriter(new FileWriter(file_path, true)));
        } catch (IOException e) {}


        logWriter.println("extraction parameters:");

        logWriter.println("advantra.critpoint.mask.iDiff:       \t"+        this.iDiff);
        logWriter.println("advantra.critpoint.profile.d:        \t"+        this.D);
        logWriter.println("advantra.critpoint.profile.s:        \t"+        this.s);
        logWriter.println("advantra.critpoint.analyze.m:        \t"+        PeakAnalyzer2D.M            + this.M);
        logWriter.println("advantra.critpoint.analyze.min_cos:  \t"+        PeakAnalyzer2D.minCos       + this.minCos);
        logWriter.println("advantra.critpoint.analyze.scatter_d:\t"+        PeakAnalyzer2D.scatterDist  + this.scatterDist);
        logWriter.close(); // close log

    }

    private static boolean deleteDir(File dir)
    {
        if (dir.isDirectory())
        {
            String[] children = dir.list();
            for (int i=0; i<children.length; i++)
            {
				File f = new File(dir, children[i]);
                boolean success = deleteDir(f);
                if (!success)
                {
					IJ.log("could not delete " + f.getAbsolutePath());
                    return false;
                }
            }
            // The directory is now empty so delete it
            return dir.delete();
        }
		else {
			// if it was a file - empty the dir anyway
			dir.delete();
		}

        return true;

    }

}
