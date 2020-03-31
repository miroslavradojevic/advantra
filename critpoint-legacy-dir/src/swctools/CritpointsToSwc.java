package swctools;

import aux.ReadDET;
import aux.Tools;
import ij.Prefs;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;

import java.io.*;

/**
 * Created by miroslav on 11-11-14.
 */
public class CritpointsToSwc implements PlugIn {

    String path_det;

    public void run(String s) {

        // choose .det file
        String in_folder = Prefs.get("det.folder", System.getProperty("user.home"));
        OpenDialog.setDefaultDirectory(in_folder);
        OpenDialog dc = new OpenDialog("Select DET detection");
        in_folder = dc.getDirectory();
        path_det = dc.getPath();
        if (path_det==null) return;
        Prefs.set("swc.folder", in_folder);

        if (!Tools.getFileExtension(path_det).equals("det")) {
            System.out.println("file needs to be .det (2 dimensional)");
            return;
        }

        ReadDET reader = new ReadDET(path_det);

        System.out.println(reader.x.size());
        reader.exportSwc(Tools.removeExtension(path_det) + ".swc");


    }
}
