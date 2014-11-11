package swctools;

import aux.ReadSWC;
import aux.Tools;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;

import java.io.File;

/**
 * Created by miroslav on 4-10-14.
 */
public class TreeToCritpoints implements PlugIn {

    static String                           path_swc;
    static ReadSWC.ExportCritpointType      what_to_extract;
    static ReadSWC.ExportCritpointFormat    which_format;

    public static void main(String[] args)
    { // terminal input

        File f;

        if (args.length == 3) {

            f = new File(args[0]);
            if (f.exists()) path_swc = f.getAbsolutePath();
            else System.out.println(args[0] + " file does not exist.");

            if (args[1].equals("ALL"))      what_to_extract = ReadSWC.ExportCritpointType.ALL;
            else if (args[1].equals("JUN")) what_to_extract = ReadSWC.ExportCritpointType.JUN;
            else if (args[1].equals("END")) what_to_extract = ReadSWC.ExportCritpointType.END;
            else                            what_to_extract = null;

            if (args[2].equals("DET_2D"))       which_format = ReadSWC.ExportCritpointFormat.DET_2D;
            else if (args[2].equals("SWC"))     which_format = ReadSWC.ExportCritpointFormat.SWC;
            else                                which_format = null;

        }
        else System.out.println(
                "Usage: \n" +
                "java -cp critpoint_.jar:ij.jar swctools.TreeToCritpoints  swc_path  what_to_extract  which_format\n" +
                        "(what_to_extract: ALL, END, JUN)\n" +
                        "(which_format: DET_2D, SWC)\n");

        swc_to_critpoint(path_swc, which_format, what_to_extract);

    }

    public void run(String s)
    {

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

        GenericDialog gd = new GenericDialog("Swc->Critpoint(.det,.swc)");
        gd.addChoice("CritpointTYPE", new String[]{"ALL", "JUN", "END"}, "ALL");
        gd.addChoice("ExportFORMAT", new String[]{"DET_2D", "SWC"}, "DET_2D");

        gd.showDialog();
        if (gd.wasCanceled()) return;

        String take = gd.getNextChoice();
        if (take.equals("ALL"))         what_to_extract = ReadSWC.ExportCritpointType.ALL;
        else if (take.equals("JUN"))    what_to_extract = ReadSWC.ExportCritpointType.JUN;
        else if (take.equals("END"))    what_to_extract = ReadSWC.ExportCritpointType.END;
        else                            what_to_extract = null;

        take = gd.getNextChoice();

        if (take.equals("DET_2D"))          which_format = ReadSWC.ExportCritpointFormat.DET_2D;
        else if (take.equals("SWC"))        which_format = ReadSWC.ExportCritpointFormat.SWC;
        else                                which_format = null;

        swc_to_critpoint(path_swc, which_format, what_to_extract);

    }

    private static void swc_to_critpoint(
            String                              _path_swc,
            ReadSWC.ExportCritpointFormat       _export_format,
            ReadSWC.ExportCritpointType         _export_type
    )
    {

        System.out.println("converting:\n" + _path_swc);

        ReadSWC reader = new ReadSWC(_path_swc);

        System.out.println("# nodes = " + reader.nodes.size());

        String export_path = "";
        if (_export_format == ReadSWC.ExportCritpointFormat.DET_2D)
            export_path = Tools.removeExtension(_path_swc) + ".det";
        else if (_export_format == ReadSWC.ExportCritpointFormat.SWC)
            export_path = Tools.removeExtension(_path_swc) + "_critpoints.swc";

        reader.exportCritpoint(export_path, _export_format, _export_type);

    }

}
