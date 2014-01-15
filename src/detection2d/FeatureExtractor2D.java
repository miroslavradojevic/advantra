package detection2d;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

import java.io.*;

/**
 * Created by miroslav on 1/10/14.
 * simplified version of PeakAnalyzer2DDemo stripped down to extract features, no debug things and user interface
 * input is an image, the output is
 *
 * image.feat (csv file) with extracted features,
 * image.i2xy ( also .csv) with index to (x,y) mapping,
 * image.param (regular txt output) with parameters used for f.e. and f. description
 *
 * this is designed to be ij plugin, can be automatically called using macros, outputs have
 * the same name as the input image and are stored int the same folder
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

    int         CPU_NR;

    public int setup(String str, ImagePlus imagePlus) {

        if(imagePlus==null) return DONE;

        image_name = imagePlus.getShortTitle();
        image_dir  = imagePlus.getOriginalFileInfo().directory;

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


        GenericDialog gd = new GenericDialog("FEATURE_EXTRACTION");
        gd.addNumericField("idiff", 	iDiff, 			1, 20, "intensity margin");
        gd.addNumericField("d", 	    D, 			    1, 20, "neuron diameter[pix]");
        gd.addNumericField("m", 	    M, 			    0, 20, "branch len.");
        gd.addNumericField("mincos", 	minCos, 	    1, 20, "alignment parameter.");
        gd.addNumericField("scatter_d",scatterDist, 	1, 20, "max allowed scatter (robustness test)");
        gd.addNumericField("scale",    s, 	            1, 20, "profiler scale");

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

        output_dir =
                image_dir + image_name + "__" +
                "idiff_"+Float.toString(iDiff)+"_"+
                "D_"+Float.toString(D)+"_"+
                "M_"+Float.toString(M)+"_"+
                "minCos_"+Float.toString(minCos)+"_"+
                "scatterDist_"+Float.toString(scatterDist)+"_"+
                "s_"+Float.toString(s)+File.separator
        ;

        boolean out_dir_created = new File(output_dir).mkdirs();
        if (!out_dir_created) {

        }

        CPU_NR = Runtime.getRuntime().availableProcessors();

        return DOES_8G+DOES_32+NO_CHANGES;

    }

    public void run(ImageProcessor imageProcessor) {

        Sphere2D sph2d = new Sphere2D(D, s);

        long t1, t2;
        t1 = System.currentTimeMillis();

        /********************************/
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
        /********************************/
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
        PeakAnalyzer2D.loadTemplate(Masker2D.i2xy, Masker2D.xy2i, PeakExtractor2D.peaks_xy, inimg_xy, Masker2D.back_xy, M, minCos, scatterDist); // initialize peak analyzer parameters
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
        IJ.log("done. "+((t2-t1)/1000f)+"sec.");

        PeakAnalyzer2D.exportFeatsCsv(   output_dir+image_name+".feat");    // export csv features      from PeakAnalyzer2D.feat2 to image_name.feat
        PeakAnalyzer2D.exportFeatsLegend(output_dir+image_name+".feat.description");
        Masker2D.exportI2xyCsv(          output_dir+image_name+".i2xy"); // export csv lookup table  from Masker2D.i2xy to image_name.i2xy
        exportExtractionLegend(          output_dir+image_name+".feat.params");

    }

    private void exportExtractionLegend(String file_path) {

        IJ.log("exporting feature extraction parameters...");

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

        logWriter.println("advantra.critpoint.mask.iDiff: \t"+          this.iDiff);
        logWriter.println("advantra.critpoint.profile.d:  \t"+          this.D);
        logWriter.println("advantra.critpoint.analyze.m: \t"+           this.M);
        logWriter.println("advantra.critpoint.analyze.min_cos: \t"+     this.minCos);
        logWriter.println("advantra.critpoint.analyze.scatter_d: \t"+   this.scatterDist);
        logWriter.println("advantra.critpoint.profile.s: \t"+           this.s);

        logWriter.close(); // close log
        IJ.log("Saved in "+file_path);

    }


}
