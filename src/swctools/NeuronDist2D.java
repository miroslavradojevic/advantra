package swctools;


import aux.ReadSWC;
import aux.Tools;
import ij.Macro;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;

import java.io.File;

/**
 * Created by miroslav on 27-10-14.
 * compares the distance between two .swc reconstructions
 * it is recommendable that they are resampled before being compared
 * that is possible with vaa3d, using resample_swc plugin which
 * can be called as a terminal command ./vaa3d -x neuron_distance -f neuron_distance
 *
 * distance measure is explained in, there are three outputs SD, and SSD with Ddiv(A,B) and Ddiv(B,A) as partial results
 * Peng, Hanchuan, et al. "V3D enables real-time 3D visualization and quantitative analysis of large-scale biological image data sets." Nature biotechnology 28.4 (2010): 348-353.
 */
public class NeuronDist2D implements PlugIn {

    String _swc_path_A;
    String _swc_path_B;
    float  _dst;                // used for SSD distance

    int CPU_NR;

    public void run(String s) {

        CPU_NR = Runtime.getRuntime().availableProcessors() + 1;

        if (Macro.getOptions()==null) {

            _swc_path_A 		= 			Prefs.get("critpoint.tracing2d.swc_path_a", "");
            _swc_path_B 		= 			Prefs.get("critpoint.tracing2d.swc_path_b", "");
            _dst                = (float)   Prefs.get("critpoint.tracing2d.dst", 2);

            GenericDialog gd = new GenericDialog("NeuronDist2D");

            gd.addStringField( "Swc_A", 	_swc_path_A,     60);
            gd.addStringField( "Swc_B", 	_swc_path_B,     60);
            gd.addNumericField("dst",       _dst,            0, 10, "");

            gd.showDialog();
            if (gd.wasCanceled()) return;

            _swc_path_A       	= gd.getNextString(); 				Prefs.set("critpoint.tracing2d.swc_path_a",     _swc_path_A);
            _swc_path_B       	= gd.getNextString(); 				Prefs.set("critpoint.tracing2d.swc_path_b",     _swc_path_B);
            _dst                = (float) gd.getNextNumber();       Prefs.set("critpoint.tracing2d.dst",            _dst);

        }
        else {

            _swc_path_A  =                          Macro.getValue(Macro.getOptions(),  "swc_a", "");
            _swc_path_B  =                          Macro.getValue(Macro.getOptions(),  "swc_b", "");
            _dst         =         Float.valueOf(Macro.getValue(Macro.getOptions(), "dst", ""));

        }

        _swc_path_A = new File(_swc_path_A).getAbsolutePath();
        _swc_path_B = new File(_swc_path_B).getAbsolutePath();

        if (! new File(_swc_path_A).exists()) return;
        if (! new File(_swc_path_B).exists()) return;

        if (
                !Tools.getFileExtension(_swc_path_A).equals("swc") ||
                        !Tools.getFileExtension(_swc_path_B).equals("swc")
                ) {
            return;
        }

        ReadSWC rswc_A = new ReadSWC(_swc_path_A);
        ReadSWC rswc_B = new ReadSWC(_swc_path_B);

        // threaded implementation of neuron distance
        SwcDistanceComputer2D.load(rswc_A.nodes, rswc_B.nodes, _dst);
        int total = rswc_A.nodes.size();

        SwcDistanceComputer2D jobs[] = new SwcDistanceComputer2D[CPU_NR];

        for (int i = 0; i < jobs.length; i++) {
            jobs[i] = new SwcDistanceComputer2D(i*total/CPU_NR,  (i+1)*total/CPU_NR);
            jobs[i].start();
        }

        for (int i = 0; i < jobs.length; i++) {
            try {
                jobs[i].join();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }

//        SwcDistanceComputer2D.remainder(); // expelled!!

        // done calculating, results are stored in static class, can be extracted from there

        //  output the result - printout
        System.out.println("Ddiv(A,B)=\t"+  SwcDistanceComputer2D.dAB());
        System.out.println("Ddiv(B,A)=\t"+  SwcDistanceComputer2D.dBA());
        System.out.println("SD       =\t"+  SwcDistanceComputer2D.SD());
        System.out.println("SSD      =\t"+  SwcDistanceComputer2D.SSD(_dst));
        System.out.println("%SSD     =\t"+  SwcDistanceComputer2D.percSSD(_dst));

        // output the result in text file eval.csv in the same folder as swc A with appended scores


    }

}
