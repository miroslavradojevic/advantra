package tracing2d;

import aux.ReadDET;
import conn.Find_Connected_Regions;
import ij.IJ;
import ij.ImagePlus;
import ij.Macro;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by miroslav on 6-10-14.
 * takes critical point regions and applies the bayesian tracking between the points
 * ISBI15 submission code
 * inputs are .soma and .det to define critical regions and Bayesian Tracking (copied from Bayesian Tracking class)
 * is used to trace between the points
 */
public class TreeReconstructor2D implements PlugIn {

    BayesianTracer2D btracer;
    ReadDET rdet;

    public void run(String s) {

        // load the image through the menu
        String in_folder = Prefs.get("id.folder", System.getProperty("user.home"));
        OpenDialog.setDefaultDirectory(in_folder);
        OpenDialog dc = new OpenDialog("Select image file");
        in_folder = dc.getDirectory();
        String image_path = dc.getPath();
        if (image_path==null) return;
        Prefs.set("id.folder", in_folder);

        ImagePlus curr_img = new ImagePlus(image_path);
        int W = curr_img.getWidth();
        int H = curr_img.getHeight();

        // variables to read through the menu
        String _det_path;    // path to the detection image
        String _soma_path;   // path to the soma binary image file

        if (Macro.getOptions()==null) {

            _det_path 		= 			Prefs.get("critpoint.tracing2d.det_path", "");
            _soma_path 		= 			Prefs.get("critpoint.tracing2d.soma_path", "");

            GenericDialog gd = new GenericDialog("TreeReconstructor2D");
            gd.addStringField("DetectionFilePath", 			_det_path, 50);
            gd.addStringField("SomaFilePath", 				_soma_path, 50);
            gd.showDialog();
            if (gd.wasCanceled()) return;
            _det_path       	= gd.getNextString(); 				Prefs.set("critpoint.tracing2d.det_path", _det_path);
            _soma_path       	= gd.getNextString(); 				Prefs.set("critpoint.tracing2d.soma_path", _soma_path);

        }
        else {

//            System.out.println(Macro.getOptions());
            _det_path = Macro.getValue(Macro.getOptions(),  "DetectionFilePath", "");
            _soma_path = Macro.getValue(Macro.getOptions(),  "SomaFilePath", "");

        }

        //













//        curr_img.show();

        // load the critical point regions: .det
        rdet = new ReadDET(_det_path);
        rdet.print();

        // create critical point region map, fill it with zeros
        int[][] region_map = new int[W][H];

        int label = 1;

        // loop locations from the det file and fill the labels in the region_map
        int nr_reg = rdet.x.size(); // numner of x coords

        for (int i = 0; i < nr_reg; i++) {


            int xc = Math.round(rdet.x.get(i));
            int yc = Math.round(rdet.y.get(i));
            int rc = (int) Math.ceil(rdet.r.get(i));



            // loop throughout the sphere to get the locations that are within
            for (int xreg = xc-rc; xreg <=xc+rc; xreg++) {
                for (int yreg = yc-rc; yreg <=yc+rc; yreg++) {
                    //
                    if (xreg>=0 && xreg<W && yreg>=0 && yreg < H) { // the point is within the dims

                        // here if you wanna make circles
                        region_map[xreg][yreg] = label;
                        label++;
                    }
                    //
                }
            }

        }

        // load soma regions: binary .tif file
        ImageProcessor somap = new ImagePlus(_soma_path).getProcessor();

        //check
        if ((somap.getWidth()!=W) || (somap.getHeight()!=H)) return;

        ByteProcessor score = new ByteProcessor(somap.getWidth(), somap.getHeight());
        for (int ii=0; ii<somap.getWidth()*somap.getHeight(); ii++) if (somap.getf(ii) > 0) score.set(ii, 255);
        Find_Connected_Regions conn_reg = new Find_Connected_Regions(new ImagePlus("", score), true);  // true means save locations
        conn_reg.run("");
        System.out.println(" done.");
        ArrayList<ArrayList<int[]>> soma_conn =  conn_reg.getConnectedRegions();

        // add somas top the map - each connected region will be assigned a negative label: -1, -2, -3...
        label = -1;

        for (int soma_idx = 0; soma_idx < soma_conn.size(); soma_idx++) {

            for (int loc_idx = 0; loc_idx < soma_conn.get(soma_idx).size(); loc_idx++) {
                int xcoord = soma_conn.get(soma_idx).get(loc_idx)[1];
                int ycoord = soma_conn.get(soma_idx).get(loc_idx)[0];
                region_map[xcoord][ycoord] = label;
            }

            label = label - 1;
        }

        new ImagePlus("region_map", new FloatProcessor(region_map)).show();

        // form trace starting points




        // tracker parameters
        String          scales_list         = "3,5,7,9";
        float           D                   = 5f;
        float           prior_sigma_deg     = 35f;
        int             Nt                  = 100;

        btracer = new BayesianTracer2D(
                scales_list,
                curr_img,
                D,
                prior_sigma_deg,
                Nt
        );
        System.out.println("initialized bayesian tracer");

    }



}
