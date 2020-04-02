package swctools;

import aux.ReadSWC;
import aux.Tools;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;

/**
 * Created by miroslav on 11-11-14.
 */
public class Shifter implements PlugIn {

    static String                           path_swc;

    public void run(String s) {

        // load the .swc through the menu
        String in_folder = Prefs.get("swc.folder", System.getProperty("user.home"));
        OpenDialog.setDefaultDirectory(in_folder);
        OpenDialog dc = new OpenDialog("Select SWC reconstruction");
        in_folder = dc.getDirectory();
        path_swc = dc.getPath();
        if (path_swc==null) return;
        Prefs.set("swc.folder", in_folder);

        if (!Tools.getFileExtension(path_swc).equals("swc")) {
            System.out.println("file needs to be .swc");
            return;
        }

        GenericDialog gd = new GenericDialog("Shifter");
//        gd.addChoice("CritpointTYPE", new String[]{"ALL", "JUN", "END"}, "ALL");
//        gd.addChoice("ExportFORMAT", new String[]{"DET_2D", "SWC"}, "DET_2D");
        gd.addNumericField("dx", 0, 1, 10, "pix");
        gd.addNumericField("dy", 0, 1, 10, "pix");
        gd.addNumericField("dz", 0, 1, 10, "pix");
        gd.addNumericField("dr", 0, 1, 10, "pix");

        gd.showDialog();
        if (gd.wasCanceled()) return;

        float dx = (float) gd.getNextNumber();
        float dy = (float) gd.getNextNumber();
        float dz = (float) gd.getNextNumber();
        float dr = (float) gd.getNextNumber();

        ReadSWC reader = new ReadSWC(path_swc);

        reader.modifySwc(Tools.removeExtension(path_swc) + "_shifted.swc", dx, dy, dz, dr);

    }

}
