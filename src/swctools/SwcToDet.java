package swctools;

import aux.ReadSWC;
import aux.Tools;
import ij.ImagePlus;
import ij.Prefs;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;

import java.io.File;

/**
 * Created by miroslav on 4-10-14.
 * terminal call
 * java -cp critpoint_.jar:/home/miroslav/ImageJ/ij.jar swctools.SwcToDet AAA
 */
public class SwcToDet implements PlugIn {

    static String path_swc;
//    static String path_det;

    public static void main(String[] args) {

        File f;

        if (args.length == 1) {
            f = new File(args[0]);
            if (f.exists()) path_swc = f.getAbsolutePath();
            else System.out.println(args[0] + " file does not exist.");
        }
        else System.out.println("one argument necessary");

        swc_to_det(path_swc);

    }

    public void run(String s) {

        // load the .swc through the menu
        String in_folder = Prefs.get("swc.folder", System.getProperty("user.home"));
        OpenDialog.setDefaultDirectory(in_folder);
        OpenDialog dc = new OpenDialog("Select swc reconstruction");
        in_folder = dc.getDirectory();
        path_swc = dc.getPath();
        if (path_swc==null) return;
        Prefs.set("swc.folder", in_folder);

        swc_to_det(path_swc);

    }

    // convert reconstruction swc to the .det format used by Critpoint2D for detections in 2d
    private static void swc_to_det(String _path_swc){

        System.out.println("converting:\n" + _path_swc);
        // save it to the same location with .det extension
//        File 	fileSWC = new File(_path_swc);
        ReadSWC reader = new ReadSWC(_path_swc);

        reader.exportDetCritpoint(Tools.removeExtension(_path_swc) + ".det");
        reader.exportSwcCritpoint(Tools.removeExtension(_path_swc) + ".det.swc");

    }

}
