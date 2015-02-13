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
import tracing2d.BayesianTracer;

import java.awt.*;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * main plugin for neuron reconstruction in 2d
 * Created by miroslav on 10-2-15.
 */

public class NeuronReconstructor2D implements PlugIn {

    ImagePlus inimg;
    String inimg_path;
    String inimg_dir_path;
    String inimg_name;
    float[][] inimg_xy;

    int W, H;

    // soma extraction
    float   sigma_min   = 10;
    float   sigma_max   = 20;
    int     sigma_step  = 4;
    float   tolerance   = 5;
    int     minSz       = 5; // min size of soma area

    float neuron_radius = 4; // neurite diameter
    float sigma_degrees = 80;

    int Ni = 100; // hardcoded
    int Ns  = 50; // hardcoded

    int counter;
    // tag map for book-keeping the traces
    int[] tag_map = null; // will be used to book-keep the tags - value at some location will correspond to the

    // trace list
    ArrayList<Trace2D> tlist = new ArrayList<Trace2D>();

    // bayesian tracer
    BayesianTracer2D bt;
    //


    public void run(String s) {

        // some test things
//        int L = 15;
//        System.out.println("L2 = " + (L/2));
//        if (true )return;

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

        W = inimg.getWidth();
        H = inimg.getHeight();

        tag_map = new int[W*H]; // zeros at the beginnning
        counter = 1;

        System.out.println("\n_____________________\nNeuronStalker2D...");

        // load params
        GenericDialog gd = new GenericDialog("NeuronStalker2D");
        gd.addMessage("SOMA");
        gd.addNumericField("sigmaMin",  sigma_min, 0);
        gd.addNumericField("sigmaMax",  sigma_max, 0);
        gd.addNumericField("tolerance", tolerance, 0);

        gd.addMessage("TRACE");
        gd.addNumericField("neuronRadius", neuron_radius, 0);
        gd.addNumericField("angularSigmaDeg", sigma_degrees, 0);

        gd.showDialog();
        if (gd.wasCanceled()) {
            return;
        }

        sigma_min = (float) gd.getNextNumber();
        sigma_max = (float) gd.getNextNumber();
        tolerance = (float) gd.getNextNumber();
        neuron_radius = (float) gd.getNextNumber();
        sigma_degrees = (float) gd.getNextNumber();

        System.out.print("\n_____________________\nsoma extraction...\n");

        // create duplicate image that will be blurred
        ImagePlus blur_temp         = new ImagePlus("", inimg.getProcessor().duplicate()); // temp variable for the
        ImagePlus blur_temp_prev    = new ImagePlus("", inimg.getProcessor().duplicate()); // temp variable

        for (int i = 0; i < sigma_step; i++) {

            float sigma = sigma_min + i * ((sigma_max-sigma_min)/(float)(sigma_step-1));

            System.out.print("GaussianBlur(" + IJ.d2s(sigma,2) + ")...");
            blur_temp.setProcessor(inimg.getProcessor().duplicate());
            IJ.run(blur_temp, "Gaussian Blur...", "sigma="+sigma);
            System.out.println("done");

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

        Overlay ovl = new Overlay();

        // adding soma traces to the list and to the tag_map
        for (int i = 0; i < rt.getCounter(); i++) {

            float xx = rt.getColumn(ResultsTable.X_CENTER_OF_MASS)[i];
            float yy = rt.getColumn(ResultsTable.Y_CENTER_OF_MASS)[i];
            float rr = (float) Math.sqrt(rt.getColumn(ResultsTable.AREA)[i] / Math.PI);

            Trace2D soma_trace = new Trace2D(counter++);
            soma_trace.add(xx, yy, rr);
            soma_trace.terminate(-1);

            tlist.add(soma_trace);

            soma_trace.add_to_map(tag_map);

            PointRoi p = new PointRoi(xx+.5f, yy+.5f);
            OvalRoi  o = new OvalRoi(xx-rr, yy-rr, 2*rr, 2*rr);
            o.setFillColor(new Color(1,1,0,0.3f));
            ovl.add(p);
            ovl.add(o);

        }

        System.out.println("\n->"+tlist.size() + " trace(s).");
        System.out.println("done.");

        System.out.print("\n_____________________\nextracting foreground Masker2D...");
        System.out.println("done.");

        System.out.print("\n_____________________\nextracting correlation Zncc2D...");
        System.out.println("done.");

        System.out.print("\n_____________________\ninitializing BayesianTracer2D...");
        bt = new BayesianTracer2D(Ni, Ns, neuron_radius, sigma_degrees);
        bt.showTemplates();
        System.out.println("done.");

        System.out.print("\n_____________________\nswc export...");
        exportSwc();
        System.out.println("done.");

        // show overlay with side-results
        inimg.show();
        inimg.setOverlay(ovl);


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

    private void exportSwc(){

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

        String output_dir_name = inimg_dir_path + "NS2D";

        // create output directory
        File f = new File(output_dir_name);
        if (!f.exists()) f.mkdirs();

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
            int start_tag = tlist.get(i).tag_start;
            int end_tag = tlist.get(i).tag_end;

            int cnt = start_tag;

            System.out.println("DBG: " + curr_trace.locs_x.size());


            for (int j = 0; j < curr_trace.locs_x.size(); j++) {

                if (j==curr_trace.locs_x.size()-1) { // last
                    logWriter.println(String.format(
                            "%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f %-4d",
                            cnt, 1, curr_trace.locs_x.get(j), curr_trace.locs_y.get(j), 0f, curr_trace.rads.get(j), end_tag));
                }
                else if (j==0) {
                    logWriter.println(String.format( // first
                            "%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f %-4d",
                            cnt, 1, curr_trace.locs_x.get(j), curr_trace.locs_y.get(j), 0f, curr_trace.rads.get(j), -1));
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

}
