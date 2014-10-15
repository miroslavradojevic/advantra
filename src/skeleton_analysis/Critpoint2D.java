package skeleton_analysis;

import aux.Stat;
import conn.Find_Connected_Regions;
import detection2d.CritpointRegion;
import ij.IJ;
import ij.ImagePlus;
import ij.Macro;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.Line;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import javafx.scene.paint.*;
import jdk.nashorn.internal.objects.NativeUint16Array;

import java.awt.*;
import java.awt.Color;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by miroslav on 13-10-14.
 */
public class Critpoint2D implements PlugIn {

    String      image_path;
    ArrayList<CritpointRegion> detected_regions;

    String		image_name;
    String 		image_dir;
    String      output_dir_name;

    public void run(String s) {

        // load the image through the menu
        String in_folder = Prefs.get("id.folder", System.getProperty("user.home"));
        OpenDialog.setDefaultDirectory(in_folder);
        OpenDialog dc = new OpenDialog("Select file");
        in_folder = dc.getDirectory();
        image_path = dc.getPath();
        if (image_path==null) return;
        Prefs.set("id.folder", in_folder);

        ImagePlus ip_load = new ImagePlus(image_path);

        if(ip_load==null) return;
        if (ip_load.getType()!=ImagePlus.GRAY8) return;

        ImagePlus origIP = ip_load.duplicate();
        image_name = ip_load.getShortTitle();
        image_dir = ip_load.getOriginalFileInfo().directory;


        /*
        	load the detection parameters
         */

        boolean _show_det = false;
        int     _sigma_smooth = 0;
        int     _threshold = 125; // 0 - 255 8bit
        int     _pruneIndex = 0;// AnalyzeSkeleton_.LOWEST_INTENSITY_VOXEL;// NONE    SHORTEST_BRANCH  LOWEST_INTENSITY_VOXEL LOWEST_INTENSITY_BRANCH
        boolean _pruneEnds = false;//, _show_endpoints, _enable_interactive, _save_midresults, _calculate_directions;


        if (Macro.getOptions()==null) {

            GenericDialog gd = new GenericDialog("DETECTOR2D");

            gd.addCheckbox("SHOW_DET", 		    _show_det);
            gd.addNumericField("SIGMA_SMOOTH",  _sigma_smooth, 0, 10,   "pix.");
            gd.addNumericField("THRESHOLD",     _threshold, 0, 10,      "grey.");
            gd.addNumericField("PRUNE_INDEX", 	_pruneIndex,   0, 10,   "tag.");
            gd.addCheckbox("PRUNE_ENDS", 		_pruneEnds);

            gd.showDialog();
            if (gd.wasCanceled()) return;

            _show_det       = gd.getNextBoolean();
            _sigma_smooth   = (int) Math.round(gd.getNextNumber());
            _threshold      = (int) Math.round(gd.getNextNumber());
            _pruneIndex     = (int) Math.round(gd.getNextNumber());
            _pruneEnds      = gd.getNextBoolean();

        }
        else {

//            System.out.println(Macro.getOptions());
            _show_det = Boolean.valueOf(Macro.getValue(Macro.getOptions(),          "show_det", String.valueOf(_show_det)));
            _sigma_smooth = Integer.valueOf(Macro.getValue(Macro.getOptions(),      "sigma_smooth", String.valueOf(_sigma_smooth)));
            _threshold = Integer.valueOf(Macro.getValue(Macro.getOptions(),         "threshold", String.valueOf(_threshold)));
            _pruneIndex = Integer.valueOf(Macro.getValue(Macro.getOptions(),        "prune_index", String.valueOf(_pruneIndex)));
            _pruneEnds = Boolean.valueOf(Macro.getValue(Macro.getOptions(),         "prune_ends", String.valueOf(_pruneEnds)));

        }

//        System.out.println("show_det\t" + _show_det);
//        System.out.println("sigma_smooth\t" + _sigma_smooth);
//        System.out.println("THRESHOLD\t" + _threshold);
//        System.out.println("PRUNE_INDEX\t" + _pruneIndex);
//        System.out.println("PRUNE_END\t" + _pruneEnds);

//        if (true) return;

        output_dir_name = image_dir+String.format(
                "DET.sigma.th.pruneIndex.pruneEnds_"+
                        "%d_%d_%d_%s",
                _sigma_smooth,
                _threshold,
                _pruneIndex,
                Boolean.toString(_pruneEnds)
        );

        if (_sigma_smooth>0)
            IJ.run(ip_load, "Gaussian Blur...", "sigma="+Integer.toString(_sigma_smooth));
//        else{
////            System.out.println("skipping the blurring...");
//        }

        IJ.setThreshold(ip_load, _threshold, 255);
//        IJ.setAutoThreshold(ip_load, threshold_method+" dark");

        Prefs.blackBackground = true;

        IJ.run(ip_load, "Convert to Mask", ""); // mask
        if (_show_det) ip_load.duplicate().show();

        IJ.run(ip_load, "Skeletonize (2D/3D)", ""); // threshold
        if (_show_det) ip_load.duplicate().show();

        // analyse skeleton
        AnalyzeSkeleton_ skel = new AnalyzeSkeleton_(); // Initialize AnalyzeSkeleton_
        skel.calculateShortestPath = true;
        skel.setup("", ip_load);

        // (work on a copy of the ImagePlus if you don't want it displayed)
        // run(int pruneIndex, boolean pruneEnds, boolean shortPath, ImagePlus origIP, boolean silent, boolean verbose)
        SkeletonResult skelResult = skel.run(_pruneIndex, _pruneEnds, true, origIP, true, false); // Perform analysis in silent mode

        ArrayList<Point> slab_locs = skelResult.getListOfSlabVoxels();
        boolean[][] slab_map = new boolean[ip_load.getWidth()][ip_load.getHeight()];
        for (int i = 0; i < slab_locs.size(); i++) {
            int xloc = slab_locs.get(i).x;
            int yloc = slab_locs.get(i).y;
            slab_map[xloc][yloc] = true;
        }

        ByteProcessor region_map = new ByteProcessor(ip_load.getWidth(), ip_load.getHeight());

        detected_regions = new ArrayList<CritpointRegion>(); // initialize

        int[] dx = new int[]{-1, -1, -1,  0,     0, +1, +1, +1};
        int[] dy = new int[]{-1,  0, +1, -1,    +1, -1,  0, +1};

//        Overlay dbgov = new Overlay();

        ArrayList<Point> locs = skelResult.getListOfJunctionVoxels();

        // loop junction locations and add them to the junciton region in case they got at least one neighbouring slab point
        // having one slab in 8-neighbourhood confirms the critical point
        for (int i = 0; i < locs.size(); i++) {

            int xloc = locs.get(i).x;
            int yloc = locs.get(i).y;

            // check if it is inside within the 8-neighbourhood margin
            if (xloc>0 && xloc<ip_load.getWidth()-1 && yloc>0 && yloc<ip_load.getHeight()-1) {

                for (int nb = 0; nb < 8; nb++) {// loop neighbours

                    int xloc8 = xloc+dx[nb];
                    int yloc8 = yloc+dy[nb];

                    if (slab_map[xloc8][yloc8]) {

                        region_map.set(xloc8, yloc8, 255);
                        region_map.set(xloc, yloc, 255);

//                        OvalRoi ovroi = new OvalRoi(xloc, yloc, 1, 1);
//                        ovroi.setFillColor(Color.RED);
//                        ovroi.setStrokeColor(Color.RED);
//                        ovroi.setStrokeWidth(0.1f);
//                        dbgov.add(ovroi);
//
//                        ovroi = new OvalRoi(xloc8, yloc8, 1, 1);
//                        ovroi.setFillColor(Color.RED);
//                        ovroi.setStrokeColor(Color.RED);
//                        ovroi.setStrokeWidth(0.1f);
//                        dbgov.add(ovroi);

                    }
                }

            }

        }

        appendDetectedRegions(region_map, CritpointRegion.RegionType.BIF_CROSS, detected_regions);

        for (int i = 0; i < region_map.getWidth()*region_map.getHeight(); i++) region_map.set(i, 0); // reset

        locs = skelResult.getListOfEndPoints();

        for (int i = 0; i < locs.size(); i++) {

            int xloc = locs.get(i).x;
            int yloc = locs.get(i).y;

            // check if it is inside within the 8-neighbourhood margin
            if (xloc>0 && xloc<ip_load.getWidth()-1 && yloc>0 && yloc<ip_load.getHeight()-1) {

                for (int nb = 0; nb < 8; nb++) {// loop neighbours

                    int xloc8 = xloc+dx[nb];
                    int yloc8 = yloc+dy[nb];

                    if (slab_map[xloc8][yloc8]) {

                        region_map.set(xloc8, yloc8, 255);
                        region_map.set(xloc, yloc, 255);

//                        OvalRoi ovroi = new OvalRoi(xloc, yloc, 1, 1);
//                        ovroi.setFillColor(Color.YELLOW);
//                        ovroi.setStrokeColor(Color.YELLOW);
//                        ovroi.setStrokeWidth(0.1f);
//                        dbgov.add(ovroi);
//
//                        ovroi = new OvalRoi(xloc8, yloc8, 1, 1);
//                        ovroi.setFillColor(Color.YELLOW);
//                        ovroi.setStrokeColor(Color.YELLOW);
//                        ovroi.setStrokeWidth(0.1f);
//                        dbgov.add(ovroi);

                    }
                }

            }

        }

        appendDetectedRegions(region_map, CritpointRegion.RegionType.END, detected_regions);

        Overlay ov = saveDetection(); // export centroids and radiuses to .det

        if (_show_det) {
            origIP.setOverlay(ov);
            origIP.show();
//            ip_load.setOverlay(dbgov);
//            ip_load.show();
        }

        System.out.println("DONE.");

    }




    private void appendDetectedRegions( // appends to the list of CritpointRegion, without directions
                                        ByteProcessor               _region_map,
                                        CritpointRegion.RegionType  _choose_type,
                                        ArrayList<CritpointRegion>  region_list
    ) {

        ///// connected components
        Find_Connected_Regions conn_reg = new Find_Connected_Regions(new ImagePlus("", _region_map), true);
        conn_reg.run("");
        ArrayList<ArrayList<int[]>> regs = conn_reg.getConnectedRegions();

        // regs
        for (int i=0; i<regs.size(); i++) {

//            if (regs.get(i).size() > Amax) continue; // go to the next region
//            if (regs.get(i).size() < 2) continue; // go to the next region

//            nr_regs++;

            float Cx=0, Cy=0; 	// centroid
//            float C=0;        	// score
            float Cr=3f; 			// radius

            CritpointRegion.RegionType Ctype;   	// type
            float[][] Cdirections;           		// directions

            for (int aa=0; aa<regs.get(i).size(); aa++) {  // loop elements of the region

                int xcoord = regs.get(i).get(aa)[1];
                int ycoord = regs.get(i).get(aa)[0];

                Cx += regs.get(i).get(aa)[1];
                Cy += regs.get(i).get(aa)[0];
//                C += _score_map.getf(xcoord, ycoord); //_critpoint_det[xcoord][ycoord];

            }

            Cx /= regs.get(i).size();
            Cy /= regs.get(i).size();   // centroid

//            C /= regs.get(i).size();    // score (calculate it here regardless of the type, says how much average fuzzy score was confident on the output)
//            Cr = 1f;//(float) Math.ceil(             1f*Math.sqrt(regs.get(i).size()/3.14f)        ); // radius is wrt to the area
//            Cr = (Cr<1)? 1 : Cr;

            if (_choose_type == CritpointRegion.RegionType.END) {

                Ctype = CritpointRegion.RegionType.END;
                Cdirections = new float[1][2];
                Arrays.fill(Cdirections[0], Float.NaN);

                region_list.add(new CritpointRegion(Ctype, Cx, Cy, Cr, 1f, Cdirections, 1)); // take one from Cdirections

            }
            else if (_choose_type == CritpointRegion.RegionType.BIF_CROSS) {

                Ctype = CritpointRegion.RegionType.BIF;
                Cdirections = new float[3][2];
                Arrays.fill(Cdirections[0], Float.NaN);
                Arrays.fill(Cdirections[1], Float.NaN);
                Arrays.fill(Cdirections[2], Float.NaN);

                region_list.add(new CritpointRegion(Ctype, Cx, Cy, Cr, 1f, Cdirections, 3)); // take three from Cdirections

            }
        }
    }

    public Overlay saveDetection() // will save the detection: ArrayList<CritpointRegion> -> file.det
    {
        // saves detections in predefined output folders
        // 2 outputs:
        // 1. image_name.det textual file with the description of critpoint regions
        //      format: x, y, radius, score, type{BIF,END,CROSS}, dir{vx,vy; vx,vy; ...}
        // 2. return Overlay array that was extracted

        Overlay ov = new Overlay();

        //// initialize file
        String det_path = output_dir_name + File.separator + image_name+".det";

        // create output dir
            File f = new File(output_dir_name);
            if (!f.exists()) f.mkdirs();

            try {
                PrintWriter logWriter = new PrintWriter(det_path);
                logWriter.print("");
                logWriter.close();
            } catch (FileNotFoundException ex) {}

        try {

            PrintWriter logWriter;

            logWriter = new PrintWriter(new BufferedWriter(new FileWriter(det_path, true)));
            logWriter.println("# DETECTOR2D...\n# X, Y, RADIUS, SCORE, TYPE, DIRECTIONS(vx,vy, vx,vy...) ");

            // write the regions
            for (int ii=0; ii<detected_regions.size(); ii++) {

                if (detected_regions.get(ii)!=null) { // just in case if the CritpointRegion was null

                    // write one component of the Overlay
                    float cx = detected_regions.get(ii).centroid[0];
                    float cy = detected_regions.get(ii).centroid[1];
                    float cr = detected_regions.get(ii).radius;
                    float sc = detected_regions.get(ii).score;
                    CritpointRegion.RegionType ctype = detected_regions.get(ii).type;

                    Color region_color = null;
                    switch (ctype) {
                        case BIF:
                            region_color = Color.RED;//new Color(1, 0, 0, sc);
                            break;

                        case END:
                            region_color = Color.YELLOW;//new Color(1, 1, 0, sc);
                            break;

                        case CROSS:
                            region_color = Color.GREEN;//new Color(0, 1, 0, sc);
                            break;

                        default:
                            IJ.log("non valid critical point");
                            break;
                    }

                    if (ctype==null) continue; // skip adding the overlay and line in .det

                    // write line to .det for each region
                    String curr_detection = "";

                    //add region as OvalRoi to output array of Overlays[]
                    OvalRoi ovroi = new OvalRoi(cx-cr, cy-cr, 2*cr, 2*cr);
                    ovroi.setStrokeWidth(1);
                    ovroi.setStrokeColor(region_color);
                    ovroi.setFillColor(region_color);
                    ov.add(ovroi);

                    // add line to .det  (what's loaded so far) // HERE IS HOW IT IS WRITTTEN!!!
                    curr_detection +=
                            IJ.d2s(cx,2)+", "+IJ.d2s(cy,2)+", "+
                                    IJ.d2s(cr,2)+", "+
                                    IJ.d2s(sc,2)+", "+detected_regions.get(ii).type+", ";

                    // add directions (outward_directions) os Line  Overlay component
                    float scale_direction = 1.3f; // just for visualization
                    int nr_directions = detected_regions.get(ii).outward_directions.length;

                    for (int j = 0; j < nr_directions; j++) {

                        float dx = detected_regions.get(ii).outward_directions[j][0];
                        float dy = detected_regions.get(ii).outward_directions[j][1];

                        curr_detection += IJ.d2s(dx,2)+", "+IJ.d2s(dy,2); // line for .det file
                        if (j<nr_directions-1) curr_detection += ", ";

                        if (!Float.isNaN(dx) && !Float.isNaN(dy)) {

                            Line l = new Line(cx, cy, cx+scale_direction*cr*dx, cy+scale_direction*cr*dy);
                            l.setStrokeWidth(2);
                            l.setStrokeColor(region_color);
                            l.setFillColor(region_color);
                            ov.add(l);

                        }

                    }

                    logWriter.println(curr_detection);

                }

            }

            logWriter.close();
            System.out.println("exported .DET: " + det_path);

        } catch (IOException e) {
            System.out.println("crash!!! \n");
            e.printStackTrace();
        }

        return ov;

    }

}
