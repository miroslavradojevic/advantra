import ij.*;
import ij.gui.GenericDialog;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;

import java.io.*;
import java.util.ArrayList;

/**
 * Created by miroslav on 8-3-15.
 */
public class NeuronStalker implements PlugIn {

    String      inimg_path;
    ImagePlus   inimg;
    String      inimg_dir_path;
    String      inimg_name;
    String      output_dir_midresults;

    float       neuron_diameter;
    float       percentile;
    float       correlation_boundary;
    float       std_angle_deg;
    float       std_gcsstd_pix;
    int         Ni;
    int         Ns;
    float       step = 2f; // fix this one
    float       zDist;

    int         search_mode;
    boolean     save_midresults = false;

    ArrayList<Node>                 stalker     = new ArrayList<Node>();
    ArrayList<Trace>                trace_list  = new ArrayList<Trace>();

    public void run(String s) {

        /**
         * load the image
         */
        String in_folder = Prefs.get("id.folder", System.getProperty("user.home"));
        OpenDialog.setDefaultDirectory(in_folder);
        OpenDialog dc = new OpenDialog("Select image file");
        in_folder = dc.getDirectory();
        inimg_path = dc.getPath();
        if (inimg_path==null) return;
        Prefs.set("id.folder", in_folder);

        inimg = new ImagePlus(inimg_path);
        if (inimg.getType()!=ImagePlus.GRAY8) return;
        inimg.setTitle("inimg");

        int W = inimg.getWidth();
        int H = inimg.getHeight();
        int L = inimg.getNSlices();

        inimg_dir_path      = inimg.getOriginalFileInfo().directory;
        inimg_name          = inimg.getShortTitle();
        inimg_name          = Toolbox.removeExtension(inimg.getOriginalFileInfo().fileName);


        /**
         * load the parameters
         */

        String[] mode = new String[]{"SOMA_CLOSE", "SOMA_AWAY", "HIGH_CORR"};

        if (Macro.getOptions()==null) {

            neuron_diameter         = (float)   Prefs.get("neuronstalker.neuron_diameter", 8f);
            percentile              = (float)   Prefs.get("neuronstalker.percentile", 85f);
            correlation_boundary    = (float)   Prefs.get("neuronstalker.correlation_boundary", 0.75f);
            std_angle_deg           = (float)   Prefs.get("neuronstalker.std_angle_deg", 70f);
            std_gcsstd_pix          = (float)   Prefs.get("neuronstalker.std_gcsstd_pix", 3f);
            Ni                      = (int)     Prefs.get("neuronstalker.Ni", 100);
            Ns                      = (int)     Prefs.get("neuronstalker.Ns", 50);
            zDist                   = (float)   Prefs.get("neuronstalker.zDist", 2f);
            search_mode             = (int)     Prefs.get("neuronstalker.search_mode", 2);
            save_midresults         = (boolean) Prefs.get("neuronstalker.save_midresults", false);

            GenericDialog gd = new GenericDialog("NeuronStalker");

        gd.addNumericField("neuron_diameter",       neuron_diameter, 0);
        gd.addNumericField("percentile",            percentile, 0);
        gd.addNumericField("correlation_boundary",  correlation_boundary, 2);
        gd.addNumericField("std_angle_deg",         std_angle_deg, 0);
        gd.addNumericField("std_gcsstd_pix",        std_gcsstd_pix, 0);
        gd.addNumericField("Ni", Ni, 0);
        gd.addNumericField("Ns", Ns, 0);
        gd.addNumericField("zDist", zDist, 1);
        gd.addChoice("search_mode", mode, mode[search_mode]);
        gd.addCheckbox("save_midresults", save_midresults);

        gd.showDialog();
        if (gd.wasCanceled()) return;

        neuron_diameter = (float) gd.getNextNumber();       Prefs.set("neuronstalker.neuron_diameter", neuron_diameter);            //System.out.println(neuron_diameter);
        percentile = (float) gd.getNextNumber();            Prefs.set("neuronstalker.percentile", percentile);                      //System.out.println(percentile);
        correlation_boundary = (float) gd.getNextNumber();  Prefs.set("neuronstalker.correlation_boundary", correlation_boundary);  //System.out.println(correlation_boundary);
        std_angle_deg = (float) gd.getNextNumber();         Prefs.set("neuronstalker.std_angle_deg", std_angle_deg);                //System.out.println(std_angle_deg);
        std_gcsstd_pix = (float) gd.getNextNumber();        Prefs.set("neuronstalker.std_gcsstd_pix", std_gcsstd_pix);              //System.out.println(std_gcsstd_pix);
        Ni = (int) gd.getNextNumber();                      Prefs.set("neuronstalker.Ni", Ni);                                      //System.out.println(Ni);
        Ns = (int) gd.getNextNumber();                      Prefs.set("neuronstalker.Ns", Ns);                                      //System.out.println(Ns);
        zDist = (float) gd.getNextNumber();                 Prefs.set("neuronstalker.zDist", zDist);                                //System.out.println(zDist);
        search_mode = gd.getNextChoiceIndex();              Prefs.set("neuronstalker.search_mode", search_mode);                    //System.out.println(search_mode);
        save_midresults = gd.getNextBoolean();              Prefs.set("neuronstalker.save_midresults", save_midresults);            //System.out.println(save_midresults);

        }
        else {

            neuron_diameter         = Float.valueOf(Macro.getValue(Macro.getOptions(),  "neuron_diameter", String.valueOf(8)));
            percentile              = Float.valueOf(Macro.getValue(Macro.getOptions(),  "percentile", String.valueOf(85)));
            correlation_boundary    = Float.valueOf(Macro.getValue(Macro.getOptions(),  "correlation_boundary", String.valueOf(0.75)));
            std_angle_deg           = Float.valueOf(Macro.getValue(Macro.getOptions(),  "std_angle_deg", String.valueOf(70)));
            std_gcsstd_pix          = Float.valueOf(Macro.getValue(Macro.getOptions(),  "std_gcsstd_pix", String.valueOf(3)));
            Ni                      = Integer.valueOf(Macro.getValue(Macro.getOptions(),"Ni", String.valueOf(100)));
            Ns                      = Integer.valueOf(Macro.getValue(Macro.getOptions(),"Ns", String.valueOf(50)));
            zDist                   = Float.valueOf(Macro.getValue(Macro.getOptions(),  "zDist", String.valueOf(2)));
            String read_mode        = Macro.getValue(Macro.getOptions(),"search_mode", mode[2]);
            if (read_mode.equalsIgnoreCase(mode[0])) {
                search_mode = 0;
            }
            else if (read_mode.equalsIgnoreCase(mode[1])) {
                search_mode = 1;
            }
            else if (read_mode.equalsIgnoreCase(mode[2])){
                search_mode = 3;
            }
            else return;

            System.out.println("mode idx. = " + search_mode);
            save_midresults         = false;

        }


        output_dir_midresults = inimg_dir_path + inimg_name + "_midresults" + File.separator;
        File f = new File(output_dir_midresults);
        if (!f.exists() && save_midresults) f.mkdirs();

        long t1 = System.currentTimeMillis();

        /**
         * initialize classes used for reconstruction stages
         */
        SomaExtractor se = new SomaExtractor();
        ForegroundExtractor fe = new ForegroundExtractor();
        GuidepointExtractor gpe = new GuidepointExtractor();


        /**
         * extract soma(s)
         */
        float open_radius = .4f*neuron_diameter;
        System.out.print("_____________________\nextracting soma(s) (R=" + IJ.d2s(open_radius,1) + ")... ");
        se.work(inimg, open_radius, (int) neuron_diameter, (int) (0.25*inimg.getWidth()*inimg.getHeight()), save_midresults?output_dir_midresults:"");
        System.out.println(se.soma_list.size() + " soma(s) found");

        if (L==1) { // 2d

            float[][] inimg_xy = new float[W][H]; 	// x~column, y~row
            byte[] read = (byte[]) inimg.getProcessor().getPixels();
            for (int idx=0; idx<read.length; idx++) inimg_xy[idx % W][idx / W] = (float) (read[idx] & 0xff);

            /**
             * extract foreground
             */
            System.out.print("_____________________\nextracting foreground (R=" + IJ.d2s(neuron_diameter,0) + ")... ");
            fe.work2D(inimg_xy, (int) Math.ceil(Math.sqrt(2)*neuron_diameter), neuron_diameter, percentile, save_midresults?output_dir_midresults:"");
            System.out.println(IJ.d2s(fe.i2xy.get(0).length/1000f, 0) + "k foreground locs, " + IJ.d2s((100f * fe.i2xy.get(0).length) / (float) (W * H * L), 1) + "%");


            /**
             * extract guidepoints
             */
            System.out.println("_____________________\nextracting guidepoints");
            gpe.work2D(inimg_xy, .5f*neuron_diameter, fe.i2xy.get(0), se.soma_list, correlation_boundary, search_mode, save_midresults?output_dir_midresults:"");
            System.out.println(IJ.d2s(gpe.queue.get(0).length/1000f,0) + "k guidepoints");

            /**
             * stalking in 2d...
             */
            System.out.println("_____________________\nstalking... ");
            stalker.clear();
            stalker.add(new Node());    // stalker elements will be referenced by tags maintained through the tag map, tags start from 1
                                        // corresponding index in the list is initialized with a dummy node due to such indexing

            int[] queue = gpe.queue.get(0);
            boolean[][] queue_map = gpe.isqueue_loc.get(0);
            boolean[] queue_success = new boolean[queue.length];
            int[] tag_map = new int[W*H];

            // soma addition
            for (int i = 0; i < se.soma_list.size(); i++) {

                float xx = (float) se.soma_list.get(i)[0];
                float yy = (float) se.soma_list.get(i)[1];
                float rr = (float) se.soma_list.get(i)[2];

                Node soma_node = new Node(xx, yy, rr);
                stalker.add(soma_node);

                addToMap(stalker, stalker.size()-1, tag_map, W, H);
                if (save_midresults) IJ.saveAs(getTagMap(tag_map, W, H), "Tiff", output_dir_midresults + "tagMap" + String.format("%05d", stalker.size()) + ".tif");

            }

            BayesianTracer2D bt = new BayesianTracer2D(Ni, Ns, step, .5f*neuron_diameter);
            int new_traces_found = Integer.MAX_VALUE;
            int cnt_evol = 0;

            while (new_traces_found>0) {

                new_traces_found = 0;
                int LOG_EVERY = queue.length / 10;

                for (int i = 0; i < queue.length; i++) {

                    if (i % LOG_EVERY == 0) System.out.print(IJ.d2s((i/LOG_EVERY)*10,0) + "%\t");

                    if (queue_success[i]) continue;

                    // trace in 1 direction
                    int startpx         =  fe.i2xy.get(0)[    queue[i]][0];
                    int startpy         =  fe.i2xy.get(0)[    queue[i]][1];
                    float startvx       = gpe.i2dir.get(0)[   queue[i]][0];
                    float startvy       = gpe.i2dir.get(0)[   queue[i]][1];
                    float startgcsstd   = gpe.i2sigma.get(0)[ queue[i]];

                    if (tag_map[startpy * W + startpx] != 0) continue; // no retracing

                    // trace 1
                    int outcome_1 = bt.trace(startpx, startpy, startvx, startvy, startgcsstd, inimg_xy, std_angle_deg, std_gcsstd_pix, tag_map, fe.mask_xy.get(0), queue_map);

                    int prev_tag, next_tag, curr_tag, lim; // aux. variables for storing the neighboring indexes
                    boolean found_t1 = false;
                    boolean found_tag = outcome_1>0 && bt.iter_counter>=1;

                    if (found_tag || bt.last_queue_element_checked>0) {

                        found_t1 = true;

                        lim = (found_tag)? bt.iter_counter-1 : bt.last_queue_element_checked-1;

                        stalker.add(new Node(bt.xc[lim][0], bt.xc[lim][1], bt.rc[lim]));

                        addToMap(stalker, stalker.size()-1, tag_map, W, H);
                        if (save_midresults) IJ.saveAs(getTagMap(tag_map, W, H), "Tiff", output_dir_midresults + "tagMap" + String.format("%05d", stalker.size()) + ".tif");

                        if (found_tag) { // additional linking with another trace
                            prev_tag = outcome_1;
                            curr_tag = stalker.size()-1;
                            stalker.get(prev_tag).neighbors.add(curr_tag);
                            stalker.get(curr_tag).neighbors.add(prev_tag);
                        }

                        for (int j = lim-1; j >= 0; j--) { // trace 1 backwards

                            stalker.add(new Node(bt.xc[j][0], bt.xc[j][1], bt.rc[j]));

                            addToMap(stalker, stalker.size()-1, tag_map, W, H);
                            if (save_midresults) IJ.saveAs(getTagMap(tag_map, W, H), "Tiff", output_dir_midresults + "tagMap" + String.format("%05d", stalker.size()) + ".tif");

                            prev_tag = stalker.size()-2;
                            curr_tag = stalker.size()-1;
                            stalker.get(prev_tag).neighbors.add(curr_tag);
                            stalker.get(curr_tag).neighbors.add(prev_tag);
                        }
                    }

                    // trace 2
                    int outcome_2 = bt.trace(startpx, startpy, -startvx, -startvy, startgcsstd, inimg_xy, std_angle_deg, std_gcsstd_pix, tag_map, fe.mask_xy.get(0), queue_map);

                    found_tag = outcome_2>0 && bt.iter_counter>=1;

                    if (found_tag || bt.last_queue_element_checked>0) {

                        lim = (found_tag)? bt.iter_counter-1 : bt.last_queue_element_checked-1;

                        stalker.add(new Node(bt.xc[0][0], bt.xc[0][1], bt.rc[0]));

                        addToMap(stalker, stalker.size()-1, tag_map, W, H);
                        if (save_midresults) IJ.saveAs(getTagMap(tag_map, W, H), "Tiff", output_dir_midresults + "tagMap" + String.format("%05d", stalker.size()) + ".tif");

                        if (found_t1) { // additional linking with trace 1
                            prev_tag = stalker.size()-2; // the last one from trace 1
                            curr_tag = stalker.size()-1;
                            stalker.get(prev_tag).neighbors.add(curr_tag);
                            stalker.get(curr_tag).neighbors.add(prev_tag);
                        }

                        for (int j = 1; j <= lim; j++) { // trace 2 forward

                            stalker.add(new Node(bt.xc[j][0], bt.xc[j][1], bt.rc[j]));

                            addToMap(stalker, stalker.size()-1, tag_map, W, H);
                            if (save_midresults) IJ.saveAs(getTagMap(tag_map, W, H), "Tiff", output_dir_midresults + "tagMap" + String.format("%05d", stalker.size()) + ".tif");

                            prev_tag = stalker.size()-2;
                            curr_tag = stalker.size()-1;
                            stalker.get(prev_tag).neighbors.add(curr_tag);
                            stalker.get(curr_tag).neighbors.add(prev_tag);
                        }

                        if (found_tag) { // additional linking with another trace
                            next_tag = outcome_2;
                            curr_tag = stalker.size()-1;
                            stalker.get(next_tag).neighbors.add(curr_tag);
                            stalker.get(curr_tag).neighbors.add(next_tag);
                        }

                    }

                }

                cnt_evol++;
                System.out.println(" " + cnt_evol + " -> " + new_traces_found + " new traces" );

            }

            System.out.println(stalker.size() + " nodes");

        }
        else { // 3d

//            System.out.println("3d in development! Exiting...");
//            System.exit(-1);

            ArrayList<float[][]> inimg_xyz = new ArrayList<float[][]>();    // x~column, y~row, z~slice
            for (int i = 0; i < L; i++) {                                   // loop through the number of slices
                byte[] read = (byte[]) inimg.getStack().getProcessor(i+1).getPixels();
                float[][] tt = new float[W][H];
                for (int idx=0; idx<read.length; idx++) tt[idx % W][idx / W] = (float) (read[idx] & 0xff);
                inimg_xyz.add(tt);
            }

            fe.work3D(inimg_xyz, (int) Math.ceil(Math.sqrt(2)*neuron_diameter), neuron_diameter, percentile, output_dir_midresults);
            int nr_locs = 0;
            for (int i = 0; i < fe.i2xy.size(); i++) nr_locs += fe.i2xy.get(i).length;
            System.out.println(IJ.d2s(nr_locs/1000f, 0) + "k locations of foreground, " + IJ.d2s((100f*nr_locs) / (float) (W * H * L), 1) + "%");

            gpe.work3D(inimg_xyz, .5f*neuron_diameter, fe.i2xy, se.soma_list, correlation_boundary, search_mode, output_dir_midresults);


            /**
             * stalking in 3d...
             */
            stalker.clear();
            stalker.add(new Node());    // elements of the list will be referenced by tags maintained through the tag map
            // the tags start from 1, tags zero are background
            // corresponding index in the list is initialized with a dummy node due to such indexing

            ArrayList<int[]> queue = gpe.queue;

            ArrayList<boolean[][]> queue_map = gpe.isqueue_loc;

            ArrayList<boolean[]> queue_success = new ArrayList<boolean[]>();
            for (int i = 0; i < queue.size(); i++) queue_success.add(new boolean[queue.get(i).length]);

//            int counter = 1; // counting nodes in swc starts from 1

            ArrayList<int[]> tag_map = new ArrayList<int[]>();
            for (int i = 0; i < queue.size(); i++) tag_map.add(new int[W * H]);

            // soma addition
            for (int i = 0; i < se.soma_list.size(); i++) {

                float xx = (float) se.soma_list.get(i)[0];
                float yy = (float) se.soma_list.get(i)[1];
                float zz = (float) se.soma_list.get(i)[2];
                float rr = (float) se.soma_list.get(i)[3];

                Node soma_node = new Node(xx, yy, zz, rr);
                stalker.add(soma_node);

                addToMap(stalker, stalker.size()-1, tag_map, W, H, zDist);
                IJ.saveAs(getTagMap(tag_map, W, H), "Tiff", output_dir_midresults + "tagMap" + String.format("%05d", stalker.size()) + ".tif");

            }

            BayesianTracer3D bt = new BayesianTracer3D(Ni, Ns, step, .5f*neuron_diameter);
            bt.getTemplates().show();
            int new_traces_found = Integer.MAX_VALUE;
            int cnt_evol = 0;


        }

        /**
         *  build tree using breath-first search
         */
        System.out.print("_____________________\ngenerating tree... ");

        // first elements of the stalker list are soma's and they are taken as seeds, stalker.get(1),...,stalker.get(nr_soma)
        int[] seed_indexes = new int[se.soma_list.size()];
        for (int i = 0; i <seed_indexes.length; i++) seed_indexes[i] = i + 1;

        ArrayList[] discovered = init_discovered_list(stalker);
        trace_list = build_tree(stalker, discovered, seed_indexes, Trace.SOMA, Trace.AXON, Trace.LOOSE); // start with soma as seeds

        System.out.println(trace_list.size() + " traces.");

        /**
         * exporting .swc reconstruction
         **/
        System.out.print("_____________________\nswc export... ");

        File swc_dir = new File(inimg_dir_path +
                String.format("NS.D.perc.corr.ang.gcs.Ni.zDist.mode_%.0f_%.0f_%.2f_%.0f_%.0f_%d_%.1f_%d",
                        neuron_diameter, percentile, correlation_boundary, std_angle_deg, std_gcsstd_pix, Ni, zDist, search_mode));
        if (!swc_dir.exists()) swc_dir.mkdirs();
        String swc_path = swc_dir.getAbsolutePath() + File.separator + inimg_name + ".swc";

        exportSwc(trace_list, swc_path);
        System.out.println(swc_path); // "done.\nExport .SWC: " +

        System.out.println("FINISHED.");

        long t2 = System.currentTimeMillis();
        System.out.println("elapsed " + ((t2-t1)/1000f) + " s.");

        if (Macro.getOptions()==null) IJ.setTool("hand");

    }

    private int get_undiscovered(ArrayList[] discovered){

        for (int i = 1; i < discovered.length; i++) { // first element in stalker list of nodes is nothing
            for (int j = 0; j < discovered[i].size(); j++) {
                if (!(Boolean) discovered[i].get(j)) {
                    return i;
                }
            }
        }

        return -1; // all are discovered

    }

    private ArrayList[] init_discovered_list(ArrayList<Node> _stalker) {
        // will be used to book-keep discovered traces (graph vertices), all are not discovered at the beginning
        ArrayList[] discovered = new ArrayList[_stalker.size()];
        for (int i = 1; i < discovered.length; i++) { // skip the first index as the first element of stalker is dummy
            discovered[i] = new ArrayList<Boolean>(_stalker.get(i).neighbors.size());
            for (int j = 0; j < _stalker.get(i).neighbors.size(); j++) {
                discovered[i].add(false);
            }
        }
        return discovered;
    }

    private ArrayList<Trace> build_tree(
            ArrayList<Node> _stalker,
            ArrayList[] discovered,
            int[] _seed_indexes,
            int seed_trace_type,
            int tree_trace_type,
            int loose_trace_type) {

        /**
         *  breadth-first search (BFS) to traverse the tree from extracted node list
         */
        boolean is3D = _stalker.get(1).loc.length==3;

        ArrayList<int[]> Q = new ArrayList<int[]>(); // queue, element will have [[prev][next]] node index
        ArrayList<Trace> tlist = new ArrayList<Trace>(); // output trace list

        for (int i = 0; i < _seed_indexes.length; i++) { // add seed indexes to queue, mark the as discovered

            int seed = _seed_indexes[i];

            Trace t = new Trace(seed_trace_type); // init seed trace

            float x_node = _stalker.get(seed).loc[0];
            float y_node = _stalker.get(seed).loc[1];
            float z_node = (is3D)?_stalker.get(seed).loc[2]:0;
            float r_node = _stalker.get(seed).r;

            t.add(-1, seed, x_node, y_node, r_node);

            tlist.add(t);

            // add the neighbors to the queue and label them as discovered
            for (int j = 0; j <_stalker.get(seed).neighbors.size(); j++) {

                int neighbor_seed = _stalker.get(seed).neighbors.get(j);

                Q.add(new int[]{seed, neighbor_seed}); // [0, 1] index
                discovered[seed].set(j, true);
                discovered[neighbor_seed].set(_stalker.get(neighbor_seed).neighbors.indexOf(seed), true);

            }

        }

        while (Q.size()>0) {

            // Q.pop()
            int prev = Q.get(Q.size()-1)[0];
            int curr = Q.get(Q.size()-1)[1];
            Q.remove(Q.size()-1);

            Trace t = new Trace(tree_trace_type); // init branch trace

            // always add the first node (it exists since this one was stored in the queue)
            t.add(
                    prev, curr,
                    _stalker.get(curr).loc[0], _stalker.get(curr).loc[1], (is3D)?_stalker.get(curr).loc[2]:0,
                    _stalker.get(curr).r
            );

            while(_stalker.get(curr).neighbors.size()==2) {

                int next = (_stalker.get(curr).neighbors.get(0)!=prev)?
                        _stalker.get(curr).neighbors.get(0):
                        _stalker.get(curr).neighbors.get(1);

                prev = curr;

                curr = next;

                t.add(
                        prev, curr,
                        _stalker.get(curr).loc[0], _stalker.get(curr).loc[1], (is3D)?_stalker.get(curr).loc[2]:0,
                        _stalker.get(curr).r
                );

                // mark as discovered the connections
                discovered[curr].set(_stalker.get(curr).neighbors.indexOf(prev), true);
                discovered[prev].set(_stalker.get(prev).neighbors.indexOf(curr), true);

            }

            tlist.add(t);

            // means that it was not on the neurite anymore, check adjacent traces
            for (int i = 0; i < discovered[curr].size(); i++) {
                boolean isDiscovered =  (Boolean) discovered[curr].get(i);
                if (!isDiscovered) { // if it was not discovered

                    int next = _stalker.get(curr).neighbors.get(i);

                    // add it to queue
                    Q.add(new int[]{curr, next});

                    // label as discovered
                    discovered[curr].set(i, true);
                    discovered[next].set(_stalker.get(next).neighbors.indexOf(curr), true);

                }
            }

        }

        // before returning do bfs from those  that were not discovered, seed each time from an undiscovered one until all are discovered
        int next_seed = get_undiscovered(discovered);

        while (next_seed!=-1) {

            // do BFS from first undiscovered seed
            Trace t = new Trace(loose_trace_type); // init seed trace

            float x_node = _stalker.get(next_seed).loc[0];
            float y_node = _stalker.get(next_seed).loc[1];
            float z_node = (is3D)?_stalker.get(next_seed).loc[2]:0;
            float r_node = _stalker.get(next_seed).r;

            t.add(-1, next_seed, x_node, y_node, r_node);

            tlist.add(t);

            // add the neighbors to the queue and label them as discovered
            for (int j = 0; j <_stalker.get(next_seed).neighbors.size(); j++) {

                int neighbor_seed = _stalker.get(next_seed).neighbors.get(j);

                Q.add(new int[]{next_seed, neighbor_seed}); // [0, 1] index
                discovered[next_seed].set(j, true);
                discovered[neighbor_seed].set(_stalker.get(neighbor_seed).neighbors.indexOf(next_seed), true);

            }

            while (Q.size()>0) {

                // Q.pop()
                int prev = Q.get(Q.size()-1)[0];
                int curr = Q.get(Q.size()-1)[1];
                Q.remove(Q.size()-1);

                t = new Trace(loose_trace_type); // init branch trace

                // always add the first node (it exists since this one was stored in the queue)
                t.add(
                        prev, curr,
                        _stalker.get(curr).loc[0], _stalker.get(curr).loc[1], (is3D)?_stalker.get(curr).loc[2]:0,
                        _stalker.get(curr).r
                );

                while(_stalker.get(curr).neighbors.size()==2) {

                    int next = (_stalker.get(curr).neighbors.get(0)!=prev)?
                            _stalker.get(curr).neighbors.get(0):
                            _stalker.get(curr).neighbors.get(1);

                    prev = curr;

                    curr = next;

                    t.add(
                            prev, curr,
                            _stalker.get(curr).loc[0], _stalker.get(curr).loc[1], (is3D)?_stalker.get(curr).loc[2]:0,
                            _stalker.get(curr).r
                    );

                    // mark as discovered the connections
                    discovered[curr].set(_stalker.get(curr).neighbors.indexOf(prev), true);
                    discovered[prev].set(_stalker.get(prev).neighbors.indexOf(curr), true);

                }

                tlist.add(t);

                // means that it was not on the neurite anymore, check adjacent traces
                for (int i = 0; i < discovered[curr].size(); i++) {
                    boolean isDiscovered =  (Boolean) discovered[curr].get(i);
                    if (!isDiscovered) { // if it was not discovered

                        int next = _stalker.get(curr).neighbors.get(i);

                        // add it to queue
                        Q.add(new int[]{curr, next});

                        // label as discovered
                        discovered[curr].set(i, true);
                        discovered[next].set(_stalker.get(next).neighbors.indexOf(curr), true);

                    }
                }

            }

            next_seed = get_undiscovered(discovered);

        }

        return tlist;

    }

    private ImagePlus getTagMap(int[] tag_map, int w, int h) {

        return new ImagePlus("tag_map", new FloatProcessor(w, h, tag_map.clone()));

    }

    private ImagePlus getTagMap(ArrayList<int[]> tag_map, int w, int h) {

        ImageStack isout = new ImageStack(w, h);

        for (int i = 0; i < tag_map.size(); i++) {
            isout.addSlice(new FloatProcessor(w, h, tag_map.get(i)));
        }

        return new ImagePlus("tag_map", isout);

    }

    private void addToMap(ArrayList<Node> _stalker, int _tag, int[] _tag_map, int _mapW, int _mapH) {

        Node _node_to_add = stalker.get(_tag);

        float curr_x = _node_to_add.loc[0];
        float curr_y = _node_to_add.loc[1]; //s.get(tidx)[1];
        float curr_r = _node_to_add.r; //rads.get(tidx);

        int xx_min = (int) Math.floor(curr_x - curr_r);
        int xx_max = (int) Math.ceil( curr_x + curr_r);

        int yy_min = (int) Math.floor(curr_y - curr_r);
        int yy_max = (int) Math.ceil( curr_y + curr_r);

        for (int xx = xx_min; xx <= xx_max; xx++) {
                for (int yy = yy_min; yy <= yy_max; yy++) {
                    if (xx>0 && xx<_mapW && yy>0 && yy<_mapH) {
                        if (  Math.pow(xx-curr_x,2) + Math.pow(yy-curr_y,2) <= Math.pow(curr_r,2) ) {
                            if (_tag_map[yy * _mapW + xx] == 0) {// fill new value if it has not been filled already
                                _tag_map[yy * _mapW + xx] = _tag;
                            }
                        }
                    }
                }
        }

//            curr_tag++;
//        }
    }

    private void addToMap(ArrayList<Node> _stalker, int _tag, ArrayList<int[]> _tag_map, int _mapW, int _mapH, float _zDist) {

            float curr_x = _stalker.get(_tag).loc[0];   // _trc_to_add.locs.get(tidx)[0];
            float curr_y = _stalker.get(_tag).loc[1];   // _trc_to_add.locs.get(tidx)[1];
            float curr_z = _stalker.get(_tag).loc[2];   // _trc_to_add.locs.get(tidx)[2];
            float curr_r = _stalker.get(_tag).r;        // _trc_to_add.rads.get(tidx);

            int xx_min = (int) Math.floor(curr_x - curr_r);
            int xx_max = (int) Math.ceil( curr_x + curr_r);

            int yy_min = (int) Math.floor(curr_y - curr_r);
            int yy_max = (int) Math.ceil( curr_y + curr_r);

            int zz_min = (int) Math.floor(curr_z - curr_r/_zDist);
            int zz_max = (int) Math.ceil( curr_z + curr_r/_zDist);

            for (int xx = xx_min; xx <= xx_max; xx++) {
                for (int yy = yy_min; yy <= yy_max; yy++) {
                    for (int zz = zz_min; zz <= zz_max; zz++) {

                        if (xx>0 && xx<_mapW && yy>0 && yy<_mapH && zz>0 && zz<_tag_map.size()) {

                            if ( Math.pow(xx-curr_x,2) + Math.pow(yy-curr_y,2) + Math.pow((zz-curr_z)*zDist,2) <= Math.pow(curr_r,2) ) {
                                if (_tag_map.get(zz)[yy*_mapW+xx] == 0) { // fill new value if it has not been filled already
                                    _tag_map.get(zz)[yy*_mapW+xx] = _tag; // _map[zz*(_mapW*_mapH) + yy * _mapW + xx]
                                }
                            }

                        }

                    }
                }
            }

    }

    private void exportSwc(ArrayList<Trace> trace_list, String swc_path){

        if (trace_list.size()==0) {
            System.out.println("trace list was empty");
            return;
        }

        boolean is3dtrace = trace_list.get(0).locs.get(0).length==3;

        PrintWriter logWriter = null;
        try {
            logWriter = new PrintWriter(swc_path);
            logWriter.print("");
            logWriter.close();
        } catch (FileNotFoundException ex) {}

        try {
            logWriter = new PrintWriter(new BufferedWriter(new FileWriter(swc_path, true)));
            logWriter.println("#NeuronStalker\n#ADVANTRA\n#author: Miro\n# " + trace_list.size() + " traces");
        } catch (IOException e) {}

        exportTrace(trace_list, logWriter, is3dtrace);

        logWriter.close();

    }

    private void exportTrace(ArrayList<Trace> trace_list, PrintWriter logWriter, boolean is3dtrace){

        for (int i = 0; i < trace_list.size(); i++) {

            Trace tt = trace_list.get(i);

            for (int j = 0; j < tt.locs.size(); j++) {

                logWriter.println(String.format(
                        "%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f %-4d",
                        tt.curr.get(j), tt.type, tt.locs.get(j)[0], tt.locs.get(j)[1], (is3dtrace)?tt.locs.get(j)[2]:0f, tt.rads.get(j), tt.prev.get(j)));

            }

        }

    }

}



/////////////////////////////////////// experimental
//    public static int[] clustering(boolean[][] dists)
//    {
//        // indxs represent indexes of values that need to be clustered,
//        // dists are the distances,
//        // threshold_dists is the distance limit
//        // output is list of unique labels
//        int[] labels = new int[dists.length];
//        for (int i = 0; i < labels.length; i++) labels[i] = i;
//
//        for (int i = 0; i < dists.length; i++) {
//
//            // one versus the rest
//            for (int j = 0; j < dists[0].length; j++) {
//
//                if (i != j) {
//
//                    if (dists[i][j]) {
//
//                        if (labels[j] != labels[i]) {
//
//                            int currLabel = labels[j];
//                            int newLabel  = labels[i];
//
//                            labels[j] = newLabel;
//
//                            //set all that also were currLabel to newLabel
//                            for (int k = 0; k < labels.length; k++)
//                                if (labels[k]==currLabel)
//                                    labels[k] = newLabel;
//
//                        }
//
//                    }
//
//                }
//
//            }
//
//        }
//
//        return labels;
//
//    }
//    /*
//        use output of clustering to give out the final cluster centroids
//     */
//    public static ArrayList<float[]> extracting(int[] labels, int[] vals)
//    { // int[] idxs,
//
//        boolean[] checked = new boolean[labels.length];
//        ArrayList<float[]> out = new ArrayList<float[]>();
//
//        for (int i = 0; i < labels.length; i++) {
//            if (!checked[i]) {
//
//                float centroid = vals[ i ]; // idxs[i]
//                int count = 1;
//                checked[i] = true;
//
//                // check the rest
//                for (int j = i+1; j < labels.length; j++) {
//                    if (!checked[j]) {
//                        if (labels[j]==labels[i]) {
//
//                            centroid += vals[ j ]; // idxs[j]
//                            count++;
//                            checked[j] = true;
//
//                        }
//                    }
//                }
//
//                out.add(new float[]{centroid/count, count});
//
//            }
//        }
//
//        return out;
//
//    }


//                    if (outcome_1!=0) {
//                        for (int j = bt.iter_counter-1; j >= 0; j--) {
//                            if (j==bt.iter_counter-1) {
//                                stalker.add(new Node(bt.xc[j][0], bt.xc[j][1], bt.rc[j], outcome_1));
//                            }
//                            else {
//                                int prev_idx = stalker.size()-1;
//                                stalker.add(new Node(bt.xc[j][0], bt.xc[j][1], bt.rc[j], prev_idx));
//                                int curr_idx = stalker.size()-1;
//                                stalker.get(prev_idx).neighbors.add(curr_idx);
//                            }
//                        }
//                    }
//                    else if (bt.last_queue_element_checked > 0){
//                            for (int j = bt.last_queue_element_checked-1; j >= 0; j--) {
//                                if (j==bt.last_queue_element_checked-1) {
//                                    stalker.add(new Node(bt.xc[j][0], bt.xc[j][1], bt.rc[j]));
//                                }
//                                else {
//                                    int prev_idx = stalker.size()-1;
//                                    stalker.add(new Node(bt.xc[j][0], bt.xc[j][1], bt.rc[j], prev_idx));
//                                    int curr_idx = stalker.size()-1;
//                                    stalker.get(prev_idx).neighbors.add(curr_idx);
//                                }
//                            }
//                            outcome_1=-1; // loose end
//                    }
//                    if (outcome_2!=0) {
//                        if (outcome_1!=0) {
//                            // there was a trace and the other side reached TAG
//                            for (int j = 0; j < bt.iter_counter; j++) {
//                                int prev_idx = stalker.size()-1;
//                                stalker.add(new Node(bt.xc[j][0], bt.xc[j][1], bt.rc[j], prev_idx));
//                                int curr_idx = stalker.size()-1;
//                                stalker.get(prev_idx).neighbors.add(curr_idx);
//                                if (j==bt.iter_counter-1) {
//                                    stalker.get(curr_idx).neighbors.add(outcome_2);
//                                    stalker.get(outcome_1).neighbors.add(curr_idx);
//                                }
//                            }
//                        }
//                        else {
//                            // there was NO trace and the other side reached TAG
//                            for (int j = 0; j < bt.iter_counter; j++) {
//                                if (j==bt.iter_counter-1) {
//                                }
//                                else {
//                                }
//                            }
//                        }
//                    }
//                    ArrayList<Integer> trace_neighbour = new ArrayList<Integer>();
//                    Trace trace = null;
//                    // ---------------------------
//                    // different cases
//                    // ---------------------------
//                    if (outcome_1==-1   &&  outcome_2==-1) {
//                        trace = new Trace(counter, -1, Trace.AXON);
//                        for (int j = trace_1.locs.size()-1; j >= 0; j--)    {trace.add(trace_1.locs.get(j)[0], trace_1.locs.get(j)[1], trace_1.rads.get(j)); trace_idx.add(trace_list.size());}
//                        for (int j = 0; j < trace_2.locs.size(); j++)       {trace.add(trace_2.locs.get(j)[0], trace_2.locs.get(j)[1], trace_2.rads.get(j)); trace_idx.add(trace_list.size());}
//                    }
//
//                    if (outcome_1==-1   && outcome_2==0) {
//                        trace = new Trace(counter, -1, Trace.AXON);
//                        for (int j = trace_1.locs.size()-1; j >= 0; j--)    {trace.add(trace_1.locs.get(j)[0], trace_1.locs.get(j)[1], trace_1.rads.get(j)); trace_idx.add(trace_list.size());}
//                    }
//
//                    if (outcome_1==-1  &&  outcome_2>0) {
//                        trace_neighbour.add(trace_idx.get(outcome_2));
//                        trace = new Trace(counter, outcome_2, Trace.AXON);
//                        for (int j = trace_2.locs.size()-1; j >= 0; j--)    {trace.add(trace_2.locs.get(j)[0], trace_2.locs.get(j)[1], trace_2.rads.get(j)); trace_idx.add(trace_list.size());}
//                        for (int j = 0; j < trace_1.locs.size(); j++)       {trace.add(trace_1.locs.get(j)[0], trace_1.locs.get(j)[1], trace_1.rads.get(j)); trace_idx.add(trace_list.size());}
//                    }
//
//                    // ---------------------------
//
//                    if (outcome_1==0  && outcome_2==-1) {
//                        trace = new Trace(counter, -1, Trace.AXON);
//                        for (int j = 0; j < trace_2.locs.size(); j++)       {trace.add(trace_2.locs.get(j)[0], trace_2.locs.get(j)[1], trace_2.rads.get(j)); trace_idx.add(trace_list.size());}
//                    }
//
//                    if (outcome_1==0  && outcome_2>0) {
//                        trace_neighbour.add(trace_idx.get(outcome_2));
//                        trace = new Trace(counter, outcome_2, Trace.AXON);
//                        for (int j = trace_2.locs.size()-1; j >= 0; j--)    {trace.add(trace_2.locs.get(j)[0], trace_2.locs.get(j)[1], trace_2.rads.get(j)); trace_idx.add(trace_list.size());}
//                    }
//
//                    // ---------------------------
//
//                    if (outcome_1>0  && outcome_2==-1) {
//                        trace_neighbour.add(trace_idx.get(outcome_1));
//                        trace = new Trace(counter, outcome_1, Trace.AXON);
//                        for (int j = trace_1.locs.size()-1; j >= 0; j--)    {trace.add(trace_1.locs.get(j)[0], trace_1.locs.get(j)[1], trace_1.rads.get(j)); trace_idx.add(trace_list.size());}
//                        for (int j = 0; j < trace_2.locs.size(); j++)       {trace.add(trace_2.locs.get(j)[0], trace_2.locs.get(j)[1], trace_2.rads.get(j)); trace_idx.add(trace_list.size());}
//                    }
//
//                    if (outcome_1>0  && outcome_2==0) {
//                        trace_neighbour.add(trace_idx.get(outcome_1));
//                        trace = new Trace(counter, outcome_1, Trace.AXON);
//                        for (int j = trace_1.locs.size()-1; j >= 0; j--)    {trace.add(trace_1.locs.get(j)[0], trace_1.locs.get(j)[1], trace_1.rads.get(j)); trace_idx.add(trace_list.size());}
//                    }
//
//                    if (outcome_1>0  && outcome_2>0) { // (this is actually a closed loop here!)
//                        trace_neighbour.add(trace_idx.get(outcome_1));
//                        trace_neighbour.add(trace_idx.get(outcome_2));
//                        trace = new Trace(counter, outcome_1, Trace.AXON);
//                        for (int j = trace_1.locs.size()-1; j >= 0; j--)    {trace.add(trace_1.locs.get(j)[0], trace_1.locs.get(j)[1], trace_1.rads.get(j)); trace_idx.add(trace_list.size());}
//                        for (int j = 0; j < trace_2.locs.size(); j++)       {trace.add(trace_2.locs.get(j)[0], trace_2.locs.get(j)[1], trace_2.rads.get(j)); trace_idx.add(trace_list.size());}
//                    }
//
//                    if (outcome_1!=0 || outcome_2!=0) { // means one of the previous cases worked out
//                        if (trace.locs.size()>0) {
//                            trace_list.add(trace);
//                            trace_nbr_idx.add(trace_neighbour);
//                            counter += trace.locs.size();
//                            addToMap(trace, tag_map, W, H);
//                            IJ.saveAs(getTagMap(tag_map, W, H), "Tiff", output_dir_midresults + "tagMap" + String.format("%05d", trace_list.size()) + ".tif");
//                            new_traces_found++;
//                            queue_success[i]=true; // so that this location is not examined again
//                        }
//                    }

//            // denote loose segments
//            boolean[][] trace_nbr_map = new boolean[trace_list.size()][trace_list.size()];
//            for (int i = 0; i < trace_nbr_idx.size(); i++) {
//                for (int j = 0; j < trace_nbr_idx.get(i).size(); j++) {
//                    trace_nbr_map[i][trace_nbr_idx.get(i).get(j)] = true;
//                    trace_nbr_map[trace_nbr_idx.get(i).get(j)][i] = true;
//                }
//            }
//            int[] labs = clustering(trace_nbr_map);
//            // we'll denote all labels that are connected with soma as AXON and the rest will be LOOSE, change will be visible in the type
//            int nr_somas = se.soma_list.size(); // first soma will be trace(0) etc...
//            for (int i = nr_somas; i < labs.length; i++) {
//                int curr_label = labs[i];
//                boolean is_conn = false;
//                for (int j = 0; j < nr_somas; j++) {
//                    if (curr_label==labs[j]) {
//                        is_conn = true;
//                        break;
//                    }
//                }
//                if (!is_conn) trace_list.get(i).type=Trace.LOOSE;           // this wil be taken into account when exporting
//            }
//                Trace soma_trace = new Trace(counter, -1, Trace.SOMA);
//                soma_trace.add(xx, yy, rr); // here trace has only one element, without known neighbors
//                trace_idx.add(trace_list.size());
//                counter += soma_trace.locs.size(); // update count
//                trace_list.add(soma_trace);
//                trace_nbr_idx.add(new ArrayList<Integer>()); // there are no neighbours
//                addToMap(soma_trace, tag_map, W, H);


//    private void addToMap(Trace _trc_to_add, int[] _map, int _mapW, int _mapH) {
//
//        int curr_tag = _trc_to_add.tag_init;
//
//        for (int tidx = 0; tidx < _trc_to_add.locs.size(); tidx++) {
//
//            float curr_x = _trc_to_add.locs.get(tidx)[0];
//            float curr_y = _trc_to_add.locs.get(tidx)[1];
//            float curr_r = _trc_to_add.rads.get(tidx);
//
//            int xx_min = (int) Math.floor(curr_x - curr_r);
//            int xx_max = (int) Math.ceil( curr_x + curr_r);
//
//            int yy_min = (int) Math.floor(curr_y - curr_r);
//            int yy_max = (int) Math.ceil( curr_y + curr_r);
//
//            for (int xx = xx_min; xx <= xx_max; xx++) {
//                for (int yy = yy_min; yy <= yy_max; yy++) {
//                    if (xx>0 && xx<_mapW && yy>0 && yy<_mapH) {
//                        if (  Math.pow(xx-curr_x,2) + Math.pow(yy-curr_y,2) <= Math.pow(curr_r,2) ) {
//                            if (_map[yy * _mapW + xx] == 0) {// fill new value if it has not been filled already
//                                _map[yy * _mapW + xx] = curr_tag;
//                            }
//                        }
//                    }
//                }
//            }
//
//            curr_tag++;
//
//        }
//    }

//    private void addToMap(Trace _trc_to_add, ArrayList<int[]> _map, int _mapW, int _mapH, float _zDist) {
//
//        int curr_tag = _trc_to_add.tag_init;
//
//        for (int tidx = 0; tidx < _trc_to_add.locs.size(); tidx++) {
//
//            float curr_x = _trc_to_add.locs.get(tidx)[0];
//            float curr_y = _trc_to_add.locs.get(tidx)[1];
//            float curr_z = _trc_to_add.locs.get(tidx)[2];
//            float curr_r = _trc_to_add.rads.get(tidx);
//
//
//            int xx_min = (int) Math.floor(_trc_to_add.locs.get(tidx)[0] - _trc_to_add.rads.get(tidx));
//            int xx_max = (int) Math.ceil( _trc_to_add.locs.get(tidx)[0] + _trc_to_add.rads.get(tidx));
//
//            int yy_min = (int) Math.floor(_trc_to_add.locs.get(tidx)[1] - _trc_to_add.rads.get(tidx));
//            int yy_max = (int) Math.ceil( _trc_to_add.locs.get(tidx)[1] + _trc_to_add.rads.get(tidx));
//
//            int zz_min = (int) Math.floor(_trc_to_add.locs.get(tidx)[2] - _trc_to_add.rads.get(tidx)/_zDist);
//            int zz_max = (int) Math.ceil( _trc_to_add.locs.get(tidx)[2] + _trc_to_add.rads.get(tidx)/_zDist);
//
//            for (int xx = xx_min; xx <= xx_max; xx++) {
//                for (int yy = yy_min; yy <= yy_max; yy++) {
//                    for (int zz = zz_min; zz <= zz_max; zz++) {
//
//                        if (xx>0 && xx<_mapW && yy>0 && yy<_mapH && zz>0 && zz<_map.size()) {
//
//                            if ( Math.pow(xx-curr_x,2) + Math.pow(yy-curr_y,2) + Math.pow((zz-curr_z)*zDist,2) <= Math.pow(curr_r,2) ) {
//                                if (_map.get(zz)[yy*_mapW+xx] == 0) { // fill new value if it has not been filled already
//                                    _map.get(zz)[yy*_mapW+xx] = curr_tag; // _map[zz*(_mapW*_mapH) + yy * _mapW + xx]
//                                }
//                            }
//
//                        }
//
//                    }
//                }
//            }
//
//            curr_tag++;
//
//        }
//    }

