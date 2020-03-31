import ij.IJ;
import ij.ImagePlus;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.io.OpenDialog;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

import java.awt.*;
import java.io.*;
import java.util.ArrayList;

public class ViewLog implements PlugInFilter {

    ImagePlus imp;

    public int setup(String s, ImagePlus imagePlus) {
        this.imp = imagePlus;
        return DOES_8G;
    }

    public void run(ImageProcessor imageProcessor) {

        IJ.log(imp.getTitle());

        OpenDialog dc = new OpenDialog("Select file");
        String path = dc.getPath();
        if (path==null) return;
        IJ.log(path);

        if (getFileExtension(path).equalsIgnoreCase("log")) {
            ArrayList<Rectangle> r = new ArrayList<Rectangle>();

            read_log(path, r); // x, y, w, h

            Overlay ov  = new Overlay();

            for (int i = 0; i < r.size(); i++) {
                ov.add(new Roi(r.get(i)));
            }

            imp.setOverlay(ov);
            imp.updateAndDraw();

        }

    }

    public static String getFileExtension(String file_path) {
        String extension = "";

        int i = file_path.lastIndexOf('.');
        if (i >= 0) {
            extension = file_path.substring(i+1);
        }

        return extension;
    }

    public static void read_log(String log_file_path, ArrayList<Rectangle> r){ // ArrayList<Float> x, ArrayList<Float> y, ArrayList<Float> w, ArrayList<Float> h) {

        String log_path = new File(log_file_path).getAbsolutePath(); // path to log file

        if (!(new File(log_path).exists())) {
            IJ.log(log_path + " does not exist!");
            return;
        }

        r.clear();

        try { // scan the file

            FileInputStream fstream 	= new FileInputStream(log_path);
            BufferedReader br 			= new BufferedReader(new InputStreamReader(new DataInputStream(fstream)));
            String read_line;

            while ( (read_line = br.readLine()) != null ) {

                if (read_line.isEmpty()) continue;

                if(!read_line.trim().startsWith("#")) { // # are comments

                    String[] 	readLn = 	read_line.trim().replaceAll("," , ".").split("\\s+");

                    if (readLn.length>=5) {

                        int x = Integer.valueOf(readLn[1].trim()).intValue();//floatValue();
                        int y = Integer.valueOf(readLn[2].trim()).intValue();//floatValue();
                        int w = Integer.valueOf(readLn[3].trim()).intValue();//floatValue();
                        int h = Integer.valueOf(readLn[4].trim()).intValue();//floatValue();

                        r.add(new Rectangle(x, y, w, h));

                    }

                }

            }

            br.close();
            fstream.close();

        }
        catch (Exception e) {
            System.err.println("Error: " + e.getMessage());
        }

    }

}
