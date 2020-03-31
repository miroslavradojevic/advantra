import ij.IJ;
import ij.Prefs;
import ij.WindowManager;
import ij.io.DirectoryChooser;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import ij.text.TextPanel;
import ij.text.TextWindow;

import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by miroslav on 5/9/14.
 */
public class TableTool implements PlugIn {

    public void run(String s) {

        TextWindow tw = (TextWindow) WindowManager.getFrame("NeuronJ: Tracings");
        if (tw==null) {IJ.showMessage("'NeuronJ: Tracings' was not found."); return;}

        TextPanel  tp = tw.getTextPanel();
		if (tp==null) {IJ.showMessage("'NeuronJ: Tracings' does not have TextPanel."); return;}

		String[] legend = tp.getColumnHeadings().split("\t");

        ArrayList<ResultsTable> res_tabs = new ArrayList<ResultsTable>();

		String first_row = tp.getLine(0);
        String[] first_row_parts = first_row.split("\t");
        String first_cls = first_row_parts[2];
        String first_len = first_row_parts[5];

		ResultsTable rt_to_add = new ResultsTable();
		rt_to_add.getFreeColumn(legend[2]);
		rt_to_add.getFreeColumn(legend[5]);
		rt_to_add.setValue(0, 0, first_cls);
		rt_to_add.setValue(1, 0, Float.valueOf(first_len));

		res_tabs.add(rt_to_add);

        for (int row=1; row<tp.getLineCount(); row++) {

            String curr_row = tp.getLine(row);
            String[] curr_row_parts = curr_row.split("\t");
            String curr_row_cls = curr_row_parts[2];
            String curr_row_len = curr_row_parts[5];

            // check if it is in the list and append if not
            boolean isAvailable = false;
            for (int ll = 0; ll<res_tabs.size(); ll++) {

                if (curr_row_cls.equals(res_tabs.get(ll).getStringValue(0,0))) {
                    isAvailable = true;

					// append
					int curr_size = res_tabs.get(ll).getColumn(0).length;
					res_tabs.get(ll).setValue(0, curr_size, curr_row_cls);
					res_tabs.get(ll).setValue(1, curr_size, Float.valueOf(curr_row_len));

					break; // stop looping further
                }
            }
            if (!isAvailable) {

				rt_to_add = new ResultsTable();
				rt_to_add.getFreeColumn(legend[2]);
				rt_to_add.getFreeColumn(legend[5]);
				rt_to_add.setValue(0, 0, curr_row_cls);
				rt_to_add.setValue(1, 0, Float.valueOf(curr_row_len));

				res_tabs.add(rt_to_add);

			}

        }

		// destination folder
		String out_folder = Prefs.get("neuronj_add.out_folder", System.getProperty("user.home"));
		DirectoryChooser.setDefaultDirectory(out_folder);
		DirectoryChooser dc = new DirectoryChooser("Select destination folder");
		out_folder = dc.getDirectory();
		if (out_folder==null) return;

		File dir = new File(out_folder);
		dir.mkdir();

		Prefs.set("neuronj_add.out_folder", out_folder);

		IJ.log("saving:");
		for (int ii=0; ii<res_tabs.size(); ii++) {
			String out_name = String.format("root%02d.xls", ii+1);
			res_tabs.get(ii).show(out_name);
			try {
				IJ.log(out_folder+out_name);
				res_tabs.get(ii).saveAs(out_folder+out_name);
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

	}

}
