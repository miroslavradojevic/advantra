import ij.Macro;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;

import java.io.*;

/**
 * measures the distance between two swc reconstructions
 * first choose swc file reconstruction then submit the path to the other one that will be used for the comparison
 * based on measuring the distance measures between two swc-s (SD,SSD,%SSD), we evaluate the reconstruction
 * appends the result (line with the: filename annot_tag sd ssd %ssd) into eval.csv file that's located in the same directory as
 * chosen swc file (the one that is submitted is normally used as a )
 * Created by miroslav on 16-3-15.
 */
public class Evaluator implements PlugIn {

    public void run(String s) {

        String rec_swc_path;

        /**
         *  load the image with detections through the menu
          */
        String in_folder = Prefs.get("id.folder", System.getProperty("user.home"));
        OpenDialog.setDefaultDirectory(in_folder);
        OpenDialog dc = new OpenDialog("Select reconstruction file");
        in_folder = dc.getDirectory();
        rec_swc_path = dc.getPath();
        if (rec_swc_path==null) return;
        Prefs.set("id.folder", in_folder);

        if (!Toolbox.getFileExtension(rec_swc_path).equals("swc")) {
            System.out.println("file needs to be .swc");
            return;
        }
        System.out.println(in_folder);
        String  output_dir_name = in_folder;
        String  output_log_name = output_dir_name + "eval.csv"; // output (file append) will be stored in the same folder as the chosen swc file
        String 	gndtth_path;			                        // ground truth swc file
        String	gndtth_tag;
        float   dst;
        System.out.println(output_log_name);
        if (Macro.getOptions()==null) {

            GenericDialog gd = new GenericDialog("GROUND TRUTH?");
            gd.addStringField("gndtth_path", 	new File(output_dir_name).getAbsolutePath(), 70);
            gd.addStringField("gndtth_tag", 	"LABEL", 50);
            gd.addNumericField("dst", 2f, 0, 10, "pix");
//            gd.addCheckbox("show",  false);
            gd.showDialog();
            if (gd.wasCanceled()) return;
            gndtth_path	= gd.getNextString();
            gndtth_tag	= gd.getNextString();
            dst = (float) gd.getNextNumber();
//            show = gd.getNextBoolean();

        }
        else {
            gndtth_path = Macro.getValue(Macro.getOptions(), "gndtth_path", "ground_truth_path");
            gndtth_tag 	= Macro.getValue(Macro.getOptions(), "gndtth_tag", 	"ground_truth_tag");
            dst = Float.valueOf(Macro.getValue(Macro.getOptions(), "dst", 	Float.toString(2)));
        }

        /**
         * read ground truth
         */
        File f_gndtth = new File(gndtth_path);

        if (!f_gndtth.exists()) {
            System.out.println("file does not exist");
            return;
        }

        if (!Toolbox.getFileExtension(gndtth_path).equals("swc")) {
            System.out.println("file needs to be .swc");
            return;
        }

        /* calculate scores */
        ReadSWC rswc_A = new ReadSWC(rec_swc_path);
        ReadSWC rswc_B = new ReadSWC(gndtth_path);

        // it can happen that they are empty - then it is not necessary to go further
        if (rswc_A.nodes.size()>0 && rswc_B.nodes.size()>0) {

            // threaded implementation of neuron distance
            SwcDistanceComputer.load(rswc_A.nodes, rswc_B.nodes, dst);
            int total = rswc_A.nodes.size();

            int CPU_NR = Runtime.getRuntime().availableProcessors() + 1;

            SwcDistanceComputer jobs[] = new SwcDistanceComputer[CPU_NR];

            for (int i = 0; i < jobs.length; i++) {
                jobs[i] = new SwcDistanceComputer(i * total / CPU_NR, (i + 1) * total / CPU_NR);
                jobs[i].start();
            }

            for (int i = 0; i < jobs.length; i++) {
                try {
                    jobs[i].join();
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }

//            System.out.println(" : "+ Arrays.toString(SwcDistanceComputer.dBA));

            PrintWriter logWriter = null;
            String legend = String.format("%10s,%10s,%6s,%6s,%6s",
                    "NAME", "ANNOT_TAG",
                    "SD", "SSD", "percSSD"
            );

            File f = new File(output_log_name);
            if (!f.exists()) {
                try {
                    logWriter = new PrintWriter(output_log_name);
                    logWriter.println(legend);
                    logWriter.close();
                } catch (FileNotFoundException ex) {
                }
            }
            // if it exists already in the folder, just prepare to append on the existing file
            try {
                logWriter = new PrintWriter(new BufferedWriter(new FileWriter(output_log_name, true)));
            } catch (IOException e) {
            }

            String eval = String.format("%10s,%10s,%6.2f,%6.2f,%6.2f",
                    Toolbox.getFileName(rec_swc_path), gndtth_tag,
                    SwcDistanceComputer.SD(),
                    SwcDistanceComputer.SSD(dst),
                    SwcDistanceComputer.percSSD(dst)
            );


            logWriter.println(eval);
            logWriter.close();

            System.out.println(legend);
            System.out.println(eval);
        }
        else {
            if (rswc_A.nodes.size()<=0) System.out.println("Empty Swc:\t" + rec_swc_path);
            if (rswc_B.nodes.size()<=0) System.out.println("Empty Swc:\t" + gndtth_path);
        }

    }

}


