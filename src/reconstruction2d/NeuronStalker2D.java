package reconstruction2d;

import aux.Tools;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.io.OpenDialog;
import ij.measure.Measurements;
import ij.measure.ResultsTable;
import ij.plugin.ImageCalculator;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import tracing2d.BayesianTracer;

import java.awt.*;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * main plugin for neuron reconstruction in 2d
 * Created by miroslav on 10-2-15.
 */

public class NeuronStalker2D implements PlugIn {

    ImagePlus   inimg;
    String      inimg_path;
    String      inimg_dir_path;
    String      inimg_name;
    float[][]   inimg_xy;

    int W, H;

    // soma extraction
    float   sigma_min   = 10;
    float   sigma_max   = 20;
    int     sigma_step  = 4;
    float   tolerance   = 5;
    int     minSz       = 5; // min size of soma area

    // params
    float new_masker_percentile = 80;
    float relative_correlation_threshold = 0.6f;

    float neuron_radius = 4; // neurite diameter
    float sigma_degrees = 80;
    float step = 2f;

    int Ni = 100; // hardcoded
    int Ns  = 50; // hardcoded

    int counter;
    // tag map for book-keeping the traces
    int[] tag_map = null; // will be used to book-keep the tags - value at some location will correspond to the

    // trace list
    ArrayList<Trace2D> tlist = new ArrayList<Trace2D>(); // some traces will be soma, some branches, etc.

    // soma centroid list
    ArrayList<float[]> slist = new ArrayList<float[]>(); // soma list (used when calcualting the distances)

    // computation module
    Zncc2D zncc;

    // bayesian tracer
    BayesianTracer2D bt;
    //

    int CPU_NR;
    boolean save_midresults = true;

    String output_dir_name;
    String output_dir_midresults;

    public void run(String s) {

//        float[] test = new float[5];
//        test[0] = 2;
//        test[1] = 5;
//        test[2] = 1;
//        test[3] = 4;
//        System.out.println(Arrays.toString(test));
//        System.out.println(Arrays.toString(Tools.descending(test)));
//        if (true) return;

        CPU_NR = Runtime.getRuntime().availableProcessors() + 1; //5;

        // load the image
        String in_folder = Prefs.get("id.folder", System.getProperty("user.home"));
        OpenDialog.setDefaultDirectory(in_folder);
        OpenDialog dc = new OpenDialog("Select image file");
        in_folder = dc.getDirectory();
        inimg_path = dc.getPath();
        if (inimg_path==null) return;
        Prefs.set("id.folder", in_folder);

        inimg = new ImagePlus(inimg_path);//.duplicate();
        inimg.setTitle("inimg");

        inimg_dir_path      = inimg.getOriginalFileInfo().directory; //  + File.separator  + image_name
        inimg_name          = inimg.getShortTitle();
        inimg_name          = Tools.removeExtension(inimg.getOriginalFileInfo().fileName);

        if (inimg.getNSlices()>1) {
            System.out.println("\nInput image needs to be 2d...\n");
            return;
        }

        if (inimg.getType()!=ImagePlus.GRAY8) {
            System.out.println("Input image needs to be BYTE8");
            return;
        }

        // image ok, go on...
        output_dir_name = inimg_dir_path + "NS2D";
        output_dir_midresults = output_dir_name + File.separator + "mids_" + inimg_name + File.separator;

        // create output directories
        File f = new File(output_dir_name);
        if (!f.exists()) f.mkdirs();

        f = new File(output_dir_midresults);
        if (!f.exists()) f.mkdirs();

        W = inimg.getWidth();
        H = inimg.getHeight();

        inimg_xy = new float[W][H]; 	// x~column, y~row
        if (inimg.getType()== ImagePlus.GRAY8) {
            byte[] read = (byte[]) inimg.getProcessor().getPixels();
            for (int idx=0; idx<read.length; idx++) {
                inimg_xy[idx%W][idx/W] = (float) (read[idx] & 0xff);
            }
        }
        else {
            System.out.println("image type needs to be GRAY8");
            return;
        }


        // allocate tag map with trace tags
        tag_map = new int[W*H]; // zeros at the beginnning
        counter = 1;

        System.out.println("_____________________\nNeuronStalker2D...");

        // load params
        sigma_min = (float) Prefs.get("neuronstalker2d.sigmaMin", 20);
        sigma_max = (float) Prefs.get("neuronstalker2d.sigmaMax", 30);
        tolerance = (float) Prefs.get("neuronstalker2d.tolerance", 5);
        neuron_radius = (float) Prefs.get("neuronstalker2d.neuron_radius", 6);
        new_masker_percentile = (float) Prefs.get("neuronstalker2d.percentile", 90);
        sigma_degrees = (float) Prefs.get("neuronstalker2d.sigma_degrees", 90);

        GenericDialog gd = new GenericDialog("NeuronStalker2D");

        gd.addNumericField("sigmaMin",      sigma_min, 0);
        gd.addNumericField("sigmaMax",      sigma_max, 0);
        gd.addNumericField("tolerance", tolerance, 0);
        gd.addNumericField("neuron_radius", neuron_radius, 0);
        gd.addNumericField("percentile", new_masker_percentile, 0);
        gd.addNumericField("angularSigmaDeg", sigma_degrees, 0);

        gd.showDialog();
        if (gd.wasCanceled()) return;

        // take params
        sigma_min = (float) gd.getNextNumber();         Prefs.set("neuronstalker2d.sigmaMin", sigma_min);
        sigma_max = (float) gd.getNextNumber();         Prefs.set("neuronstalker2d.sigmaMax", sigma_max);
        tolerance = (float) gd.getNextNumber();         Prefs.set("neuronstalker2d.tolerance", tolerance);
        neuron_radius = (float) gd.getNextNumber();     Prefs.set("neuronstalker2d.neuron_radius", neuron_radius);
        new_masker_percentile = (float) gd.getNextNumber();Prefs.set("neuronstalker2d.percentile", new_masker_percentile);
        sigma_degrees = (float) gd.getNextNumber(); Prefs.set("neuronstalker2d.sigma_degrees", sigma_degrees);

        //
        //
        // MASKER2D
        //
        //

        System.out.print("_____________________\nforeground extraction...\n");
        float new_masker_radius = 2 * neuron_radius;// 1.5f*sph2d.getOuterRadius();   	// important that it is outer radius of the sphere
                           		// used to have these two as argument but not necessary
        Masker2D.loadTemplate(
                inimg_xy,
                (int)Math.ceil(new_masker_radius),
                new_masker_radius,
                new_masker_percentile); //image, margin, check, percentile
        int totalLocs = W*H;
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
        Masker2D.defineThreshold();
        Masker2D.formRemainingOutputs();
        if (save_midresults) {
            ImagePlus mask = new ImagePlus("mask", Masker2D.getMask());
            IJ.saveAs(mask, "Tiff", output_dir_midresults+"fg.tif");
        }
        System.out.println("done. " + (Masker2D.i2xy.length/1000f) + "k locations (" +
                IJ.d2s((Masker2D.i2xy.length / (float) (W * H)) * 100f, 2) + "%%)");

        System.out.print("_____________________\nsoma extraction...\n");

        // create duplicate image that will be blurred
        ImagePlus blur_temp         = new ImagePlus("", inimg.getProcessor().duplicate()); // temp variable for the
        ImagePlus blur_temp_prev    = new ImagePlus("", inimg.getProcessor().duplicate()); // temp variable

        for (int i = 0; i < sigma_step; i++) {

            float sigma = sigma_min + i * ((sigma_max-sigma_min)/(float)(sigma_step-1));

            System.out.print("GaussianBlur(" + IJ.d2s(sigma,2) + ")...");
            blur_temp.setProcessor(inimg.getProcessor().duplicate());
            IJ.run(blur_temp, "Gaussian Blur...", "sigma="+sigma);
            System.out.print("done.");

            System.out.print("MaximumFinder...");
            MaximumFinder mf = new MaximumFinder();
            mf.setup("", blur_temp);
            blur_temp = mf.extractMaxWithTolerance(blur_temp.getProcessor(), tolerance);
            System.out.println("done.");

            if (i>0) logic_and((ByteProcessor) blur_temp.getProcessor(), (ByteProcessor) blur_temp_prev.getProcessor());

            blur_temp_prev.setProcessor(blur_temp.getProcessor());

        }

        int opts = ParticleAnalyzer.SHOW_NONE;
        int meas = Measurements.CENTROID|Measurements.CENTER_OF_MASS|Measurements.AREA|Measurements.PERIMETER;
        int maxSz = Math.round(W*H*0.25f);
        ResultsTable rt = new ResultsTable();
        ParticleAnalyzer pa = new ParticleAnalyzer(opts, meas, rt, minSz, maxSz);
        pa.setup("", blur_temp);
        pa.run(blur_temp.getProcessor());
        System.out.print("->" + rt.getCounter() + " soma(s)... ");



        // adding soma traces to the list and to the tag_map
        for (int i = 0; i < rt.getCounter(); i++) {

            float xx = rt.getColumn(ResultsTable.X_CENTER_OF_MASS)[i];
            float yy = rt.getColumn(ResultsTable.Y_CENTER_OF_MASS)[i];
            float rr = (float) Math.sqrt(rt.getColumn(ResultsTable.AREA)[i] / Math.PI);

            Trace2D soma_trace = new Trace2D(counter);
            soma_trace.add(xx, yy, rr);
            soma_trace.set_prev(-1);
            counter += soma_trace.locs_x.size();

            tlist.add(soma_trace);
            slist.add(new float[]{xx, yy});

            addToMap(soma_trace, tag_map, W, H);

        }

        System.out.print("done.\n");

        if (save_midresults) { // same soma overlayed on the original image

            Overlay ovl = new Overlay();

            for (int i = 0; i < rt.getCounter(); i++) {

                float xx = rt.getColumn(ResultsTable.X_CENTER_OF_MASS)[i];
                float yy = rt.getColumn(ResultsTable.Y_CENTER_OF_MASS)[i];
                float rr = (float) Math.sqrt(rt.getColumn(ResultsTable.AREA)[i] / Math.PI);

                PointRoi p = new PointRoi(xx+.5f, yy+.5f);
                OvalRoi  o = new OvalRoi(xx-rr, yy-rr, 2*rr, 2*rr);
                o.setFillColor(new Color(1,1,0,0.3f));
                ovl.add(p);
                ovl.add(o);


                ImagePlus outip = inimg.duplicate();
                outip.setOverlay(ovl);
                IJ.saveAs(outip, "Tiff", output_dir_midresults+"somas.tif");

            }

        }

        if (save_midresults) {
            IJ.saveAs(getTagMap(), "Tiff", output_dir_midresults+"tagMap_00.tif");
        }

        //
        //
        //  define ZNCC2D
        //
        //

        zncc = new Zncc2D(neuron_radius);
        if (save_midresults) {
            IJ.saveAs(zncc.getSampling(),   "Tiff", output_dir_midresults+"sampling.tif");
            IJ.saveAs(zncc.getTemplates(),  "Tiff", output_dir_midresults+"wgts.tif");
        }

        //
        //
        //
        //     PROFILER2D
        //
        //
        //

        System.out.print("_____________________\nlocal correlations (Profiler2D)...");
        long t1 = System.currentTimeMillis();
        Profiler2D.loadTemplate(
                zncc,
                Masker2D.i2xy,
                Masker2D.xy2i,
                slist, // current soma list to calcualate the distances towards somas
                inimg_xy);
        int totalProfileComponents = Masker2D.i2xy.length;//sph2d.getProfileLength();
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

        if (save_midresults) {
            IJ.saveAs(Profiler2D.getScores(), "Tiff", output_dir_midresults + "zncc_map.tif");
        }

        Profiler2D.createQueue(relative_correlation_threshold);

        if (save_midresults){
            IJ.saveAs(Profiler2D.getQueueMap(), "Tiff", output_dir_midresults + "queue_map.tif");
        }


        long t2 = System.currentTimeMillis();
        System.out.println("done ["+ ((t2-t1)/1000f) +" sec.].");


        //
        //
        //
        // BAYESIANTRACER2D
        //
        //

        System.out.print("_____________________\ninitializing BayesianTracer2D...");
        bt = new BayesianTracer2D(Ni, Ns, step, neuron_radius, sigma_degrees);

        if (save_midresults) {
            IJ.saveAs(bt.getTemplates(), "Tiff", output_dir_midresults + "trace_templates.tif");
        }
        System.out.println("done.");

        //
        //
        // STALKING...
        //
        //

        System.out.print("_____________________\nstalking...");

        for (int i = 0; i < Profiler2D.queue.length; i++) {

            System.out.print(".");

            if (i==1) break; // test case

            float startx = Masker2D.i2xy[Profiler2D.queue[i]][0];


        }

        System.out.println("done.");

        System.out.print("_____________________\nswc export...");
        exportSwc();
        System.out.println("done.");

        System.out.println("FINISHED.");
        IJ.setTool("hand");

    }

    private float distToSoma(float x, float y){

        float min_dst = Float.MAX_VALUE; // will return min dist to soma

        for (int i = 0; i < slist.size(); i++) {
            float dst = (float) (Math.pow(x-slist.get(i)[0],2) + Math.pow(y-slist.get(i)[1],2));
            if (dst<min_dst)
                min_dst = dst;
        }

        return (float) Math.sqrt(min_dst);
    }

    public void logic_and(ByteProcessor p1, ByteProcessor p2) { // the result will be stored in p1

        byte[] p1array = (byte[]) p1.getPixels();
        byte[] p2array = (byte[]) p2.getPixels();

        for (int i = 0; i < p1array.length; i++) {
            if ( (p1array[i]&0xff) > (p2array[i]&0xff) ) p1array[i] = p2array[i];
        }

    }

    private void showDuplicate(ImagePlus im, String tag){
        ImagePlus ii = im.duplicate();
        ii.setTitle(tag);
        ii.show();
    }

    private void printTraces(){

    }

    private void exportSwc(){ // will export all traces

//        String.format(
//                "REC.D.scalesList.priorSigmaDeg.Nt.maxIter.minIter.useOrig.nrSampMargin_"+"%.2f_%s_%.2f_%d_%d_%d_%b_%d",
//                _D,
//                _scales_list,
//                _prior_sigma_deg,
//                _Nt,
//                _MAX_ITER,
//                _MIN_ITER,
//                _USE_ORIGINAL,
//                _NR_SAMPLES_MARGIN,
//                _show_it
//        )

        String det_path = output_dir_name + File.separator + inimg_name + ".swc";

        PrintWriter logWriter = null;
        try {
            logWriter = new PrintWriter(det_path);
            logWriter.print("");
            logWriter.close();
        } catch (FileNotFoundException ex) {}

        try {
            logWriter = new PrintWriter(new BufferedWriter(new FileWriter(det_path, true)));
            logWriter.println("#NeuronStalker2D\n#ADVATNTRA project\n#author: miroslavR");
        } catch (IOException e) {}

        for (int i = 0; i < tlist.size(); i++) {

            Trace2D curr_trace = tlist.get(i);

            int cnt = tlist.get(i).tag_init;

            for (int j = curr_trace.locs_x.size()-1; j >= 0; j--) {

                if (j==curr_trace.locs_x.size()-1) { // last in the trace list, first in tagging
                    logWriter.println(String.format(
                            "%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f %-4d",
                            cnt, 1, curr_trace.locs_x.get(j), curr_trace.locs_y.get(j), 0f, curr_trace.rads.get(j), tlist.get(i).tag_prev));
                }
                else {
                    logWriter.println(String.format( // somewhere along
                            "%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f %-4d",
                            cnt, 1, curr_trace.locs_x.get(j), curr_trace.locs_y.get(j), 0f, curr_trace.rads.get(j), cnt-1));

                }

                cnt++;

            }

        }

        logWriter.close();
        System.out.println("done.\nExported .SWC: " + det_path);

    }

    private void addToMap(Trace2D _trc_to_add, int[] _map, int _mapW, int _mapH) {

        // will update tag map over the selected trace _trc_to_add
        // update with tag indexes incremented starting from the start tag
        // that's how it also works in swc export...

        int curr_tag = _trc_to_add.tag_init;
//        System.out.println("curr_tag = " + curr_tag);

        for (int tidx = _trc_to_add.locs_x.size()-1; tidx >= 0; tidx--) {

            int xx_min = (int) Math.floor(_trc_to_add.locs_x.get(tidx)-_trc_to_add.rads.get(tidx));
            int xx_max = (int) Math.ceil(_trc_to_add.locs_x.get(tidx) + _trc_to_add.rads.get(tidx));

            int yy_min = (int) Math.floor(_trc_to_add.locs_y.get(tidx)-_trc_to_add.rads.get(tidx));
            int yy_max = (int) Math.ceil(_trc_to_add.locs_y.get(tidx) + _trc_to_add.rads.get(tidx));

//            System.out.println(xx_min + " -- " + xx_max + " , " + yy_min + " -- " + yy_max);

            for (int xx = xx_min; xx <= xx_max; xx++) {
                for (int yy = yy_min; yy <= yy_max; yy++) {
                    if (xx>0 && xx<_mapW && yy>0 && yy<_mapH) {
//                        System.out.println(curr_tag+" ---> checking... x="+xx+", y="+yy+", i(x,y)="+(yy * _mapW + xx)+" map value = "+_map[yy * _mapW + xx]);
                        if (_map[yy * _mapW + xx] == 0) {// fill new value if it has not been filled already

                            _map[yy * _mapW + xx] = curr_tag;
//                            System.out.println("x="+xx+", y="+yy+", i(x,y)="+(yy * _mapW + xx)+_map[yy * _mapW + xx]);

                        }
                    }
                }
            }

            curr_tag++;

        }
    }

    private ImagePlus getTagMap() {
        return new ImagePlus("tag_map", new FloatProcessor(W, H, tag_map));
    }

}
