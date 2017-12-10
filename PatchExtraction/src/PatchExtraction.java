import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.io.FileSaver;
import ij.plugin.PlugIn;
import ij.plugin.frame.RoiManager;
import ij.process.ByteProcessor;

import java.awt.*;
import java.io.*;
import java.util.ArrayList;
import java.util.Random;

public class PatchExtraction implements PlugIn {

    String mosaic_dir_path="";
    int D_out=128;
    float step=1.0f, O_pos=0.75f, O_neg=0.00f;

//    byte[] patch_readout;
    int patch_x, patch_y;
//    ByteProcessor patch_in;     // before rescaling
    ByteProcessor patch_out;    // after rescaling
    ImagePlus     patch_out_ip; // wrapper for the re-scaled image

    public void run(String s) {
        IJ.log("test..");

        GenericDialog gd = new GenericDialog("com.braincadet.ndet.PatchExtraction");
        gd.addStringField("dir",                Prefs.get("com.braincadet.ndet.dir",                    mosaic_dir_path), 20);
//        gd.addNumericField("Din",               Prefs.get("com.braincadet.ndet.din",                    D_in),  0, 5, "pix");
        gd.addNumericField("Dout",              Prefs.get("com.braincadet.ndet.dout",                   D_out), 0, 5, "pix");
        gd.addNumericField("step",              Prefs.get("com.braincadet.ndet.step",                   step),  2, 5, "");
        gd.addNumericField("Opos",              Prefs.get("com.braincadet.ndet.opos",                   O_pos), 2, 5, "");
        gd.addNumericField("Oneg",              Prefs.get("com.braincadet.ndet.oneg",                   O_neg), 2, 5, "");

        gd.showDialog();
        if (gd.wasCanceled()) return;

        mosaic_dir_path = gd.getNextString();   Prefs.set("com.braincadet.ndet.dir",                mosaic_dir_path);
//        D_in =  (int)gd.getNextNumber();        Prefs.set("com.braincadet.ndet.din",                D_in);
        D_out = (int)gd.getNextNumber();        Prefs.set("com.braincadet.ndet.dout",               D_out);
        step = (float)gd.getNextNumber();       Prefs.set("com.braincadet.ndet.step",               step);
        O_pos = (float)gd.getNextNumber();      Prefs.set("com.braincadet.ndet.opos",               O_pos);
        O_neg = (float)gd.getNextNumber();      Prefs.set("com.braincadet.ndet.oneg",               O_neg);


        String outdir = "NDET.step.Dout.Opos.Oneg_"+IJ.d2s(step,2)+"_"+D_out+"_"+O_pos+"_"+O_neg;

        patch_out = new ByteProcessor(D_out, D_out, new byte[D_out * D_out]);
        patch_out_ip = new ImagePlus("rescaled_patch", new ByteProcessor(D_out, D_out, new byte[D_out*D_out]));

        // read mosaic tif files and corresponding annotation rectangles
        // nomenclature: m01.tif, m01.log for the mosaic 01
        File mosaic_dir = new File(mosaic_dir_path);
        if (!mosaic_dir.isDirectory()) return;

        try {
            mosaic_dir_path = mosaic_dir.getCanonicalPath();
        } catch (IOException e) {
            e.printStackTrace();
        }

        mosaic_dir = new File(mosaic_dir_path);
        File[] list = mosaic_dir.listFiles();

        createAndCleanDir(mosaic_dir.getParent()+File.separator+outdir+File.separator+"annotat_maps");
        createAndCleanDir(mosaic_dir.getParent()+File.separator+outdir+File.separator+"overlap_maps");
        createAndCleanDir(mosaic_dir.getParent()+File.separator+outdir+File.separator+"sampling_background");
        createAndCleanDir(mosaic_dir.getParent()+File.separator+outdir+File.separator+"sampling_neuron");
        createAndCleanDir(mosaic_dir.getParent()+File.separator+outdir+File.separator+"NEURONS");
        createAndCleanDir(mosaic_dir.getParent()+File.separator+outdir+File.separator+"BACKGROUND");

        IJ.log("dir=\t"+mosaic_dir_path);
//        IJ.log("Din=\t"+D_in);
        IJ.log("step=\t"+step);
        IJ.log("Dout=\t"+D_out);
        IJ.log("Opos=\t"+O_pos);
        IJ.log("Oneg=\t"+O_neg);

        ArrayList<String>               mosaic_paths        = new ArrayList<String>();                          // mosaic paths used to sample the rectangles
        ArrayList<Integer>              D                   = new ArrayList<Integer>();
        ArrayList<ByteProcessor>        mosaic_overlap_maps = new ArrayList<ByteProcessor>();          //
        ArrayList<Integer>              neuron_cnt          = new ArrayList<Integer>();          // count number of annotated neurons per mosaic
        ArrayList<ArrayList<Integer>>   neuron_idx          = new ArrayList<ArrayList<Integer>>();   // add locs with overlap > O_pos
        ArrayList<ArrayList<Float>>     neuron_wgt          = new ArrayList<ArrayList<Float>>();     // corresponding weights used for importance sampling
        int count_background = 0;                                                           // total # background rectangles
        int count_neuron_annots = 0;                                                        // total # annotated neurons

        // go through the files that fit the nomenclature and calculate overlap maps
        for (int i = 0; i < list.length; i++) {
            if (list[i].isFile()) {
                if (list[i].getName().matches("m[0-9][0-9]\\.tif")) {

                    mosaic_paths.add(list[i].getAbsolutePath());

                    ImagePlus mXY = new ImagePlus(list[i].getAbsolutePath());

                    IJ.log(list[i].getAbsolutePath()+"\t w="+
                            IJ.d2s(mXY.getWidth() /1000f,2)+"k \t h="+
                            IJ.d2s(mXY.getHeight()/1000f,2)+"k \t l="+
                            IJ.d2s(mXY.getStackSize(),0));

                    byte[] mXYarray = (byte[])mXY.getProcessor().getPixels();

                    String mXYannot = list[i].getParent()+File.separator+getFileName(list[i].getAbsolutePath())+".log";
                    IJ.log(mXYannot);

                    ArrayList<Rectangle> r = new ArrayList<Rectangle>();

                    read_log(mXYannot, r); // x (pixels), y (pixels), w (pixels), h (pixels)
                    IJ.log(r.size() + " ");
                    int[] sizes = new int[2*r.size()];
                    for (int j = 0; j < r.size(); j++) {
                        sizes[2*j]   = r.get(j).width;
                        sizes[2*j+1] = r.get(j).height;
                    }

                    int D_in = median(sizes);

                    D.add(D_in); // D annotation will be read from the annotation file and set for each mosaic
                    byte[] patch_readout = new byte[D_in*D_in];
                    ByteProcessor patch_in  = new ByteProcessor(D_in,  D_in,  new byte[D_in * D_in]);

                    neuron_cnt.add(r.size());
                    count_neuron_annots += r.size();

                    // generate the overlap map
                    ByteProcessor fp = new ByteProcessor(mXY.getWidth(), mXY.getHeight()); // initiate with zeros
                    byte[] ovlp_map  = (byte[]) fp.getPixels();
                    int val;

                    Overlay ov_annot = new Overlay();
                    RoiManager rm_annot = new RoiManager();

                    for (int j = 0; j < r.size(); j++) {

                        // annotation overlay elements
                        Roi rr = new Roi(r.get(j));
                        rr.setStrokeColor(Color.RED);
                        ov_annot.add(rr);
                        rm_annot.addRoi(rr);

                        // calculate overlap map
                        for (int xp = r.get(j).x - r.get(j).width; xp <= r.get(j).x + r.get(j).width; xp++) {
                            for (int yp = r.get(j).y - r.get(j).height; yp <= r.get(j).y + r.get(j).height; yp++) {
                                // calculate the overlap and store value in the overlap map
                                if (xp>=0 && xp<fp.getWidth() && yp>=0 && yp<fp.getHeight()) {
                                    val = Math.round((
                                            overlap(
                                                    xp, yp,                 r.get(j).width, r.get(j).height,
                                                    r.get(j).x, r.get(j).y, r.get(j).width, r.get(j).height
                                            )/(r.get(j).width * r.get(j).height)
                                    )*255f);
                                    val = (val>255)? 255 : (val<0)? 0 : val ;

                                    if (val > (ovlp_map[yp*fp.getWidth()+xp] & 0xff)) {
                                        ovlp_map[yp*fp.getWidth()+xp] = (byte) val;
                                    }

                                }
                            }
                        }
                    } // now the overlap map is complete!

//                    IJ.log("done.");

                    // save it to a list
                    mosaic_overlap_maps.add(fp);

//                    IJ.log("saving annotation roi...");
                    rm_annot.runCommand("Save", mosaic_dir.getParent()+File.separator+outdir+File.separator+"annotat_maps"+File.separator+list[i].getName()+".zip");
                    rm_annot.reset();
                    rm_annot.close();
//                    IJ.log("done");

//                    IJ.log("save overlap map...");
                    ImagePlus ip = new ImagePlus(list[i].getName(), fp); // add annots to the overlap map and save
                    ip.setOverlay(ov_annot);
                    FileSaver fs = new FileSaver(ip);
                    fs.saveAsZip(mosaic_dir.getParent()+File.separator+outdir+File.separator+"overlap_maps"+File.separator+list[i].getName()+".zip");
//                    IJ.log("done.");

//                    mXY.setOverlay(ov_annot);   // add annots to the originals and save
//                    fs = new FileSaver(mXY);
//                    fs.saveAsZip();

                    // go through the mosaic to sample the negatives using step
                    // skip if the overlap is above O_neg at sampled location - take only those that overlap below threshold
                    // needs overlap map to be filly computed first
                    Overlay ov_background = new Overlay();
                    RoiManager rm_background = new RoiManager();

                    for (int xneg = 0; xneg < fp.getWidth(); xneg+=step*D_in) {

                        for (int yneg = 0; yneg < fp.getHeight(); yneg+=step*D_in) {

                            if (xneg+D_in>=fp.getWidth()) continue;
                            if (yneg+D_in>=fp.getHeight()) continue;

                            if ((ovlp_map[yneg*fp.getWidth()+xneg] & 0xff) <= (O_neg*255f)) {


                                // extract the background patch and save to the destination folder
                                Roi background_patch = new Roi(new Rectangle(xneg, yneg, D_in, D_in));
                                background_patch.setStrokeColor(Color.BLUE);
                                ov_background.add(background_patch);
                                rm_background.addRoi(background_patch);

                                count_background++; // global background patch count

                                // patch_readout store values
                                for (int j = 0; j < D_in * D_in; j++) {
                                    patch_x = j % D_in;
                                    patch_y = j / D_in;
                                    patch_readout[patch_y * D_in + patch_x] = mXYarray[(yneg+patch_y) * mXY.getWidth() + (xneg+patch_x)];
                                }

                                patch_in.setPixels(patch_readout);
                                patch_out_ip.setProcessor(patch_in.resize(D_out, D_out, true));
                                IJ.saveAs(patch_out_ip, "Tiff",
                                        mosaic_dir.getParent()+File.separator+outdir+File.separator+"BACKGROUND"+File.separator+
                                                mosaic_paths.size()+"_"+String.format("%010d",count_background)+".tif");

                            }
                        }
                    }

                    // save background sampling visualization
                    rm_background.runCommand("Save", mosaic_dir.getParent()+File.separator+outdir+File.separator+"sampling_background"+File.separator+list[i].getName()+".zip");
                    rm_background.reset();
                    rm_background.close();

                    // go through the overlap map (all pixels) - and add to the pool of positive sampling locations and corresponding importance weights
                    ArrayList<Integer> neuron_idx_I = new ArrayList<Integer>();
                    ArrayList<Float>   neuron_wgt_I = new ArrayList<Float>();

                    for (int idx = 0; idx < ovlp_map.length; idx++) {

                        int x = idx % mXY.getWidth();
                        int y = idx / mXY.getWidth();

                        if ((ovlp_map[idx] & 0xff) >= (O_pos*255f) &&  // add to sampling pool
                                x+D_in < mXY.getWidth() &&             // if overlap higher than specified threshold for the positives
                                y+D_in < mXY.getHeight()               // if the whole patch/kernel would be within the mosaic boundaries
                                ) {  // if overlap higher than specified threshold and if the rectangle would fit within the frame

                            // add the sampling location and corresponding importance weight
                            neuron_idx_I.add(idx);
                            neuron_wgt_I.add((ovlp_map[idx] & 0xff)/255f); // to make it [0,1]

                        }
                    }

                    neuron_idx.add(neuron_idx_I);
                    neuron_wgt.add(neuron_wgt_I);

                    IJ.log("         "+IJ.d2s(count_background,0)+" background kernels\t"+r.size()+" neuron annots");

                } // if it was a match mXY.tif
            } // it is a file
        } // go through the list of files

//        int samples_per_annot_neuron = Math.round(count_background / (float)count_neuron_annots);
        IJ.log("TOTAL:\t"+ count_background + " background kernels\t"+count_neuron_annots+" neuron annots"); // => " + samples_per_annot_neuron + " samples per annotated neuron");


        // go through the mosaic images again knowing the number of background patches and samples per annotated neuron
        int count_neuron = 0;
        for (int i = 0; i < mosaic_paths.size(); i++) {

            int D_in = D.get(i);

            ImagePlus mXY = new ImagePlus(mosaic_paths.get(i));
            byte[] mXYarray = (byte[])mXY.getProcessor().getPixels();

            // patch readout
            byte[] patch_readout = new byte[D_in * D_in];
            ByteProcessor patch_in  = new ByteProcessor(D_in,  D_in,  new byte[D_in * D_in]);

            Overlay ov_neurons = new Overlay();
            RoiManager rm_neurons = new RoiManager();

            // random sample positives
            int nr_samples = Math.round(count_background / (float)count_neuron_annots) * neuron_cnt.get(i);

            IJ.log(mosaic_paths.get(i) + " " + nr_samples + " samples");
            
            ArrayList<Float> neuron_wgt_cws = new ArrayList<Float>(); // cumulative weight sum used for sampling

            for (int j = 0; j < neuron_wgt.get(i).size(); j++) {
                neuron_wgt_cws.add(neuron_wgt.get(i).get(j) + ( (j==0)?0: neuron_wgt_cws.get(j-1)) );
            }

            float wmass = neuron_wgt_cws.get(neuron_wgt_cws.size()-1);
            float u1 = (wmass/nr_samples)  * new Random().nextFloat();

            int idx = 0;

            for (int j = 0; j < nr_samples; j++) { // pick nr_samples

                float uj = u1 + j * (wmass/nr_samples);

                while (uj>neuron_wgt_cws.get(idx) && idx<neuron_wgt_cws.size()-1) idx++;

                int x = neuron_idx.get(i).get(idx) % mXY.getWidth();
                int y = neuron_idx.get(i).get(idx) / mXY.getWidth();

                // extract the neuron patch and save to the destination folder
                Roi neuron_patch = new Roi(new Rectangle(x, y, D.get(i), D.get(i)));
                neuron_patch.setStrokeColor(Color.GREEN);
                ov_neurons.add(neuron_patch);
                rm_neurons.addRoi(neuron_patch);

                count_neuron++;

                // save the patch
                for (int k = 0; k < D_in*D_in; k++) {
                    patch_x = k % D_in;
                    patch_y = k / D_in;
                    patch_readout[patch_y * D_in + patch_x] = mXYarray[(y+patch_y)*mXY.getWidth()+(x+patch_x)];
                }

                patch_in.setPixels(patch_readout);
                patch_out_ip.setProcessor(patch_in.resize(D_out, D_out, true));
                IJ.saveAs(patch_out_ip, "Tiff",
                        mosaic_dir.getParent()+File.separator+outdir+File.separator+"NEURONS"+File.separator+
                                (i+1)+"_"+String.format("%010d",count_neuron)+".tif");

            }

            // export neuron sampling
            File ff= new File(mosaic_paths.get(i));
            rm_neurons.runCommand("Save", mosaic_dir.getParent()+File.separator+outdir+File.separator+"sampling_neuron"+File.separator+ff.getName()+".zip");
            rm_neurons.reset();
            rm_neurons.close();

        } // go through mosaics

        IJ.log("FINISHED");
    }

    public static void createAndCleanDir(String dirpath) {
        File f1 = new File(dirpath);
        if (!f1.exists()) {f1.mkdirs();}
        else { // clean it, reset with each exec
            File[] files = f1.listFiles();
            for (File file : files)
                if (!file.delete()) System.out.println("Failed to delete " + file);
        }
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

    public static String getFileName(String file_path)
    {
        String name = "";

        int i = file_path.lastIndexOf('.');
        int j = file_path.lastIndexOf(File.separator);

        if (i > 0) {
            if (j>=0) {
                name = file_path.substring(j+1, i);
            }
            else {
                name = file_path.substring(0, i);
            }
        }

        return name;
    }

    private static float overlap(int x1, int y1, int w1, int h1, int x2, int y2, int w2, int h2) {

        float x11 = x1; // r1 left up
        float y11 = y1;

        float x12 = x1 + w1; // r1 right down
        float y12 = y1 + h1;

        float x21 = x2;
        float y21 = y2;

        float x22 = x2 + w2;
        float y22 = y2 + h2;

        float xDiff;
        if (x11 < x21) {
            xDiff = x12 - x21;
        } else {
            xDiff = x22 - x11;
        }

        float yDiff;
        if (y11 < y21) {
            yDiff = y12 - y21;
        } else {
            yDiff = y22 - y11;
        }
        xDiff = (xDiff < 0) ? 0 : xDiff;
        yDiff = (yDiff < 0) ? 0 : yDiff;

        return xDiff * yDiff;

    }

    private static double overlap(Rectangle rect1, Rectangle rect2) {

        double x11 = rect1.getX();//left -up
        double y11 = rect1.getY();//left -up

        double x12 = rect1.getX() + rect1.getWidth();//right -up
        double y12 = rect1.getY() + rect1.getHeight();//left-down

        double x21 = rect2.getX();//left -up
        double y21 = rect2.getY();//left -up

        double x22 = rect2.getX() + rect2.getWidth();//right -up
        double y22 = rect2.getY() + rect2.getHeight();//left-down

        double xDiff;
        if (x11 < x21) {
            xDiff = x12 - x21;
        } else {
            xDiff = x22 - x11;
        }

        double yDiff;
        if (y11 < y21) {
            yDiff = y12 - y21;
        } else {
            yDiff = y22 - y11;
        }
        xDiff = (xDiff < 0) ? 0 : xDiff;
        yDiff = (yDiff < 0) ? 0 : yDiff;

        return xDiff * yDiff;
    }


    public static void createDir(String dirpath) {
        // create directory without cleaning it up
        File f1 = new File(dirpath);
        if (!f1.exists()) {f1.mkdirs();}
//        IJ.log("createDir " + dirpath);
    }



    public static int median(int[] a) {
        int n = a.length;
        int i, j, l, m, k;
        int x;
        if (n % 2 == 0) k = (n/2)-1;
        else k = (n/2);
        l=0 ; m=n-1 ;
        while (l < m)
        {
            x=a[k] ;
            i = l ;
            j = m ;
            do
            {
                while (a[i] < x) i++ ;
                while (x < a[j]) j-- ;
                if (i <= j) {
                    int temp = a[i];
                    a[i] = a[j];
                    a[j] = temp;
                    i++ ; j-- ;
                }
            } while (i <= j) ;
            if (j < k) l = i ;
            if (k < i) m = j ;
        }
        return a[k] ;
    }


}
