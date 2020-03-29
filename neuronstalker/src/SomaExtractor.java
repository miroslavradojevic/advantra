import ij.*;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.measure.Measurements;
import ij.measure.ResultsTable;
import ij.plugin.Filters3D;
import ij.plugin.Stack_Statistics;
import ij.plugin.filter.Analyzer;
import ij.process.ImageProcessor;
import ij.process.StackStatistics;

import java.awt.*;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * class that will do soma extraction in 2d/3d
 * input is image or image stack and the output is list of locations with radiuses: soma_list
 *
 * 3d soma
 * soma_list(0): x,y,z,r
 * soma_list(1): x,y,z,r
 * soma_list(2): x,y,z,r
 * ...
 *
 * 2d soma
 * soma_list(0): x,y,r
 * soma_list(1): x,y,r
 * soma_list(2): x,y,r
 * ...
 *
 * Created by miroslav on 1-3-15.
 */
public class SomaExtractor extends Filters3D { // Filters3D because of the erosion and dilatation

    public ArrayList<double[]> soma_list = new ArrayList<double[]>();

    private void exportSoma(String swc_path) {

        PrintWriter logWriter = null;
        try {
            logWriter = new PrintWriter(swc_path);
            logWriter.print("");
            logWriter.close();
        } catch (FileNotFoundException ex) {}

        try {
            logWriter = new PrintWriter(new BufferedWriter(new FileWriter(swc_path, true)));
            logWriter.println("# " + soma_list.size() + " somas found!");
        } catch (IOException e) {}

        int cnt = 1;
        for (int i = 0; i < soma_list.size(); i++) {

            logWriter.println(String.format(
                    "%-4d %-4d %-6.2f %-6.2f %-6.2f %-3.2f %-4d",
                    cnt, 2, soma_list.get(i)[0], soma_list.get(i)[1], (soma_list.get(i).length==3)?0f:soma_list.get(i)[2], (soma_list.get(i).length==3)?soma_list.get(i)[2]:soma_list.get(i)[3], -1));

        }

        logWriter.close();


    }


    public ArrayList<double[]> work(        // [x, y, z, r], z=0 if 2d image
            ImagePlus   input_img,
            float       strel_size,
            int         min_size,
            int         max_size,
            String      _midresults_dir
    )
    {

        soma_list.clear();

        double[] out = StackStats.getStats(input_img); // [pix_count, mean, stdDev, min, max, mode]

        ImageStack img3d;

        //System.out.println("eroding...");
        img3d = filter(input_img.getImageStack(),   Filters3D.MIN, strel_size, strel_size, strel_size); // erode

        //System.out.println("dilating...");
        img3d = filter(img3d,                       Filters3D.MAX, strel_size, strel_size, strel_size); // dilate

        int threshold = (int) (.5f*(out[3]+out[4]));  // coeff*(max - min)

        ArrayList<double[]> soma_lst = new ArrayList<double[]>();

        if (input_img.getImageStack().getSize()>1) { // 3D stack
            
            soma_lst = new ConnComponents3D().getObjectList(new ImagePlus("", img3d), threshold, min_size, max_size);

            img3d = threshGrey8(img3d, threshold);

        }
        else { // 2D image

            img3d = threshGrey8(img3d, threshold);

            ConnComponents2D cc2d = new ConnComponents2D(
                    new ImagePlus("", img3d),
                    ParticleAnalyzer.SHOW_NONE,
                    Measurements.CENTROID|Measurements.CENTER_OF_MASS|Measurements.AREA| Measurements.PERIMETER,
                    min_size,
                    max_size
                    );

            soma_lst = cc2d.getObjectList();

        }

        soma_list = soma_lst;


        if (!_midresults_dir.equalsIgnoreCase("")) {
            exportSoma(_midresults_dir + "soma.swc");
            IJ.saveAs(new ImagePlus("after_opening", img3d),   "Tiff", _midresults_dir + "after_opening.tif");
        }

        return soma_lst;

    }

    public Overlay getSomaOverlay() {

        Overlay ov = new Overlay();

        for (int i = 0; i < soma_list.size(); i++) {

            double x = soma_list.get(i)[0];
            double y = soma_list.get(i)[1];

            OvalRoi oroi = null;

            if (soma_list.get(i).length==4) { // 3d
                double z = soma_list.get(i)[2];
                double r = soma_list.get(i)[3];
                oroi = new OvalRoi(x+.5 - r, y+.5 - r, 2*r, 2*r);
                oroi.setPosition((int) (Math.round(z)) +1);
            }

            if (soma_list.get(i).length==3) { // 2d
                double r = soma_list.get(i)[2];
                oroi = new OvalRoi(x+.5 - r, y+.5 - r, 2*r, 2*r);
            }

            if (oroi!=null) {
                oroi.setFillColor(new Color(1,1,0,.5f));
                ov.add(oroi);
            }

        }

        return ov;

    }

    public ImageStack threshGrey8(ImageStack in_stack, int th_value) {


        ImageStack isout = in_stack.duplicate();

        if (isout.getBitDepth()!=8) {
            System.out.println("error: image needs to be 8bit");
            return null;
        }

        // go by slices and do the thresholding
        for (int i = 0; i < isout.getSize(); i++) {
            byte[] vals = (byte[]) isout.getPixels(i+1);
            // loop throuhg to threshold the values
            for (int j = 0; j < vals.length; j++) {
                if ( (vals[j]&0xff) < th_value )
                    vals[j] = (byte) 0;
                else
                    vals[j] = (byte) 255;
            }
        }
        return isout;
    }

}

// for measuring the basic stack statistics... works also for 2d images
class StackStats extends Stack_Statistics {

    public static double[] getStats(ImagePlus imp) { // [pix_count, mean, stdDev, min, max, mode]
//        imp = IJ.getImage();
        double histMax = imp.getBitDepth() != 8 && imp.getBitDepth() != 24?0.0D:256.0D;
        int measurements = Analyzer.getMeasurements();
        Analyzer.setMeasurements(measurements | 256);
        StackStatistics stats = new StackStatistics(imp, 256, 0.0D, histMax);
        Analyzer.setMeasurements(measurements);
//        ResultsTable rt = Analyzer.getResultsTable();
//        rt.incrementCounter();
        Roi roi = imp.getRoi();
        if(roi != null && !roi.isArea()) {
            imp.deleteRoi();
            roi = null;
        }

//        double stackVoxels = 0.0D;
        double images = (double)imp.getStackSize();
        if(roi == null) {
//            stackVoxels = (double)(imp.getWidth() * imp.getHeight()) * images;
        } else if(roi.getType() == 0) {
            Rectangle cal = roi.getBounds();
//            stackVoxels = (double)(cal.width * cal.height) * images;
        } else {
            Analyzer.setMeasurements(measurements & -257);
            StackStatistics cal1 = new StackStatistics(imp, 256, 0.0D, histMax);
            Analyzer.setMeasurements(measurements);
//            stackVoxels = (double)cal1.longPixelCount;
        }

//        Calibration cal2 = imp.getCalibration();
//        String units = cal2.getUnits();
//        double scale = cal2.pixelWidth * cal2.pixelHeight * cal2.pixelDepth;
//        rt.addValue("Voxels", (double)stats.longPixelCount);
//        if(scale != 1.0D) {
//            rt.addValue("Volume(" + units + "^3)", (double)stats.longPixelCount * scale);
//        }
//        rt.addValue("%Volume", (double)stats.longPixelCount * 100.0D / stackVoxels);
//        rt.addValue("Mean", stats.mean);
//        rt.addValue("StdDev", stats.stdDev);
//        rt.addValue("Min", stats.min);
//        rt.addValue("Max", stats.max);
//        rt.addValue("Mode", stats.dmode);
//        rt.show("Results");

        return new double[]{(double)stats.longPixelCount, stats.mean, stats.stdDev, stats.min, stats.max, stats.dmode};

    }

}

// for extracting connected components in 3D
class ConnComponents3D extends Object_Counter3D {

    public void run(String arg) {
        System.out.println("alternative run...");

        //if (! setupGUI(arg)) return;
        //analyze();
    }

    void analyze_plain() {
        IJ.showStatus("3D Objects Counter");
        long start = System.currentTimeMillis();
        debug = IJ.debugMode;
        int x, y, z;
        int xn, yn, zn;
        int i, j, k, arrayIndex, offset;
        int voisX = -1, voisY = -1, voisZ = -1;
        int maxX = Width - 1, maxY = Height - 1;

        int index;
        int val;
        double col;

        int minTag;
        int minTagOld;

        cal = img.getCalibration();
        if (cal == null) {
            cal = new Calibration(img);
        }
        double pixelDepth = cal.pixelDepth;
        double pixelWidth = cal.pixelWidth;
        double pixelHeight = cal.pixelHeight;
        double zOrigin = cal.zOrigin;
        double yOrigin = cal.yOrigin;
        double xOrigin = cal.xOrigin;
        double voxelSize = pixelDepth * pixelWidth * pixelHeight;

        pict = new int[Height * Width * NbSlices];
        thr = new boolean[Height * Width * NbSlices];
        tag = new int[Height * Width * NbSlices];
        surf = new boolean[Height * Width * NbSlices];
        Arrays.fill(thr, false);
        Arrays.fill(surf, false);

        //Load the image in a one dimension array
        ImageStack stack = img.getStack();
        arrayIndex = 0;
        for (z = 1; z <= NbSlices; z++) {
            ip = stack.getProcessor(z);
            for (y = 0; y < Height; y++) {
                for (x = 0; x < Width; x++) {
                    PixVal = ip.getPixel(x, y);
                    pict[arrayIndex] = PixVal;
                    if (PixVal > ThrVal) {
                        thr[arrayIndex] = true;
                    }
                    arrayIndex++;
                }
            }
        }

        //First ID attribution
        int tagvois;
        ID = 1;
        arrayIndex = 0;
        for (z = 1; z <= NbSlices; z++) {
            for (y = 0; y < Height; y++) {
                for (x = 0; x < Width; x++) {
                    if (thr[arrayIndex]) {
                        tag[arrayIndex] = ID;
                        minTag = ID;
                        i = 0;
                        //Find the minimum tag in the neighbours pixels
                        for (voisZ = z - 1; voisZ <= z + 1; voisZ++) {
                            for (voisY = y - 1; voisY <= y + 1; voisY++) {
                                for (voisX = x - 1; voisX <= x + 1; voisX++) {
                                    if (withinBounds(voisX, voisY, voisZ)) {
                                        offset = offset(voisX, voisY, voisZ);
                                        if (thr[offset]) {
                                            i++;
                                            tagvois = tag[offset];
                                            if (tagvois != 0 && tagvois < minTag) minTag = tagvois;
                                        }
                                    }
                                }
                            }
                        }
                        if (i != 27) surf[arrayIndex] = true;
                        tag[arrayIndex] = minTag;
                        if (minTag == ID) {
                            ID++;
                        }
                    }
                    arrayIndex++;
                }
            }
            IJ.showStatus("Finding structures");
            IJ.showProgress(z, NbSlices);
        }
        ID++;

        //Minimization of IDs=connection of structures
        arrayIndex = 0;
        for (z = 1; z <= NbSlices; z++) {
            for (y = 0; y < Height; y++) {
                for (x = 0; x < Width; x++) {
                    if (thr[arrayIndex]) {
                        minTag = tag[arrayIndex];
                        //Find the minimum tag in the neighbours pixels
                        for (voisZ = z - 1; voisZ <= z + 1; voisZ++) {
                            for (voisY = y - 1; voisY <= y + 1; voisY++) {
                                for (voisX = x - 1; voisX <= x + 1; voisX++) {
                                    if (withinBounds(voisX, voisY, voisZ)) {
                                        offset = offset(voisX, voisY, voisZ);
                                        if (thr[offset]) {
                                            tagvois = tag[offset];
                                            if (tagvois != 0 && tagvois < minTag) minTag = tagvois;
                                        }
                                    }
                                }
                            }
                        }
                        //Replacing tag by the minimum tag found
                        for (voisZ = z - 1; voisZ <= z + 1; voisZ++) {
                            for (voisY = y - 1; voisY <= y + 1; voisY++) {
                                for (voisX = x - 1; voisX <= x + 1; voisX++) {
                                    if (withinBounds(voisX, voisY, voisZ)) {
                                        offset = offset(voisX, voisY, voisZ);
                                        if (thr[offset]) {
                                            tagvois = tag[offset];
                                            if (tagvois != 0 && tagvois != minTag) replacetag(tagvois, minTag);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    arrayIndex++;
                }
            }
            IJ.showStatus("Connecting structures");
            IJ.showProgress(z, NbSlices);
        }

        //Parameters determination 0:volume; 1:surface; 2:intensity; 3:barycenter x; 4:barycenter y; 5:barycenter z; 6:barycenter x int; 7:barycenter y int; 8:barycenter z int
        arrayIndex = 0;
        Paramarray = new double[ID][9];
        for (z = 1; z <= NbSlices; z++) {
            for (y = 0; y < Height; y++) {
                for (x = 0; x < Width; x++) {
                    index = tag[arrayIndex];
                    val = pict[arrayIndex];
                    Paramarray[index][0]++;
                    if (surf[arrayIndex]) Paramarray[index][1]++;
                    Paramarray[index][2] += val;
                    Paramarray[index][3] += x;
                    Paramarray[index][4] += y;
                    Paramarray[index][5] += z;
                    Paramarray[index][6] += x * val;
                    Paramarray[index][7] += y * val;
                    Paramarray[index][8] += z * val;
                    arrayIndex++;
                }
            }
            IJ.showStatus("Retrieving structures' parameters");
            IJ.showProgress(z, NbSlices);
        }

        double voxCount, intensity;
        for (i = 0; i < ID; i++) {
            voxCount = Paramarray[i][0];
            intensity = Paramarray[i][2]; // sum over all intensity values
            if (voxCount >= minSize && voxCount <= maxSize) {
                if (voxCount != 0) {
                    Paramarray[i][2] /= voxCount;
                    Paramarray[i][3] /= voxCount;
                    Paramarray[i][4] /= voxCount;
                    Paramarray[i][5] /= voxCount;
                }
                if (intensity != 0) {
                    Paramarray[i][6] /= intensity;
                    Paramarray[i][7] /= intensity;
                    Paramarray[i][8] /= intensity;
                }
            } else {
                for (j = 0; j < 9; j++) Paramarray[i][j] = 0;
            }
            IJ.showStatus("Calculating barycenters' coordinates");
            IJ.showProgress(i, ID);
        }

//        //Log data
//        if (new_results) {
//            rt = new ResultsTable();
//        } else {
//            rt = ResultsTable.getResultsTable();
//        }
//        IDarray = new int[ID];
//
//        String[] head = {"Volume", "Surface", "Intensity", "Centre X", "Centre Y", "Centre Z", "Centre int X", "Centre int Y", "Centre int Z"};
//        for (i = 0; i < head.length; i++) rt.setHeading(i, head[i]);
//
//        k = 1;
//        for (i = 1; i < ID; i++) {
//            if (Paramarray[i][0] != 0) {
//                rt.incrementCounter();
//                IDarray[i] = k;
//                voxCount = Paramarray[i][0];
//                rt.addValue(0, voxCount * voxelSize);
//                rt.addValue(1, Paramarray[i][1] / voxCount);
//                rt.addValue(2, Paramarray[i][2]);
//                rt.addValue(3, cal.getX(Paramarray[i][3]));
//                rt.addValue(4, cal.getY(Paramarray[i][4]));
//                rt.addValue(5, cal.getZ(Paramarray[i][5] - 1));
//                rt.addValue(6, cal.getX(Paramarray[i][6]));
//                rt.addValue(7, cal.getY(Paramarray[i][7]));
//                rt.addValue(8, cal.getZ(Paramarray[i][8] - 1));
//                k++;
//            }
//        }
//        if (new_results) {
//            rt.show("Results from " + imgtitle);
//        } else {
//            //if (! IJ.isResultsWindow()) IJ.showResults();
//            rt.show("Results");
//        }
//        int nParticles = rt.getCounter();
    }

    public boolean setup_plain(ImagePlus input_img, int threshold, int min_sz, int max_sz) {
        img = input_img;//WindowManager.getCurrentImage();
        if (img==null){
            IJ.noImage();
            return false;
        } else if (img.getStackSize() == 1) {
            IJ.error("Stack required");
            return false;
        } else if (img.getType() != ImagePlus.GRAY8) { // && img.getType() != ImagePlus.GRAY16
            // In order to support 32bit images, pict[] must be changed to float[], and  getPixel(x, y); requires a Float.intBitsToFloat() conversion
            IJ.error("8 bit greyscale image required");
            return false;
        }
        Width=img.getWidth();
        Height=img.getHeight();
        NbSlices=img.getStackSize();
        arrayLength=Width*Height*NbSlices;
        imgtitle = img.getTitle();

        ip=img.getProcessor();
        ThrVal=ip.getAutoThreshold();
        ip.setThreshold(ThrVal,Math.pow(2,16), ImageProcessor.RED_LUT);
        img.setSlice((int)NbSlices/2);
        img.updateAndDraw();

//        GenericDialog gd=new GenericDialog("3D objects counter");
//        gd.addSlider("Threshold: ",ip.getMin(), ip.getMax(),ThrVal);
//        gd.addSlider("Slice: ",1, NbSlices,(int) NbSlices/2);
//        sliders=gd.getSliders();
//        ((Scrollbar)sliders.elementAt(0)).addAdjustmentListener(this);
//        ((Scrollbar)sliders.elementAt(1)).addAdjustmentListener(this);
//        value = gd.getNumericFields();
//        ((TextField)value.elementAt(0)).addTextListener(this);
//        ((TextField)value.elementAt(1)).addTextListener(this);
//        gd.addNumericField("Min number of voxels: ",minSize_default,0);
//        gd.addNumericField("Max number of voxels: ",Math.min(maxSize_default, Height*Width*NbSlices),0);
//        gd.addCheckbox("New_Results Table", new_results_default);
//        gd.addMessage("Show:");
//        gd.addCheckbox("Particles",showParticles_default);
//        gd.addCheckbox("Edges",showNumbers_default);
//        gd.addCheckbox("Geometrical centre", showCentres_default);
//        gd.addCheckbox("Intensity based centre", showCentresInt_default);
//        gd.addNumericField("Dot size",DotSize_default,0);
//        gd.addCheckbox("Numbers",showNumbers_default);
//        gd.addNumericField("Font size",FontSize_default,0);
//        gd.addMessage("");
//        gd.addCheckbox("Summary", summary_default);
//        gd.showDialog();

//        if (gd.wasCanceled()){
//            ip.resetThreshold();
//            img.updateAndDraw();
//            return false;
//        }

        ThrVal=threshold;//(int) gd.getNextNumber();

//        gd.getNextNumber();

        minSize=min_sz;//(int) gd.getNextNumber();
        minSize_default = minSize;

        maxSize=max_sz;//(int) gd.getNextNumber();
        maxSize_default = maxSize;

        new_results=false;//gd.getNextBoolean();
        new_results_default = new_results;

        showParticles=false;//gd.getNextBoolean();
        showParticles_default = showParticles;

        showEdges=false;//gd.getNextBoolean();
        showEdges_default = showEdges;

        showCentres=false;//gd.getNextBoolean();
        showCentres_default = showCentres;

        showCentresInt=false;//gd.getNextBoolean();
        showCentresInt_default = showCentresInt;

        DotSize=12;//(int)gd.getNextNumber();
        DotSize_default = DotSize;

        showNumbers=false;//gd.getNextBoolean();
        showNumbers_default = showNumbers;

        FontSize=12;//(int)gd.getNextNumber();
        FontSize_default = FontSize;

        summary=false;//gd.getNextBoolean();
        summary_default = summary;

//        IJ.register(Object_Counter3D.class); // static fields preserved when plugin is restarted
        //Reset the threshold
        ip.resetThreshold();
        img.updateAndDraw();
        return true;
    }

    public ArrayList<double[]> getObjectList(ImagePlus input_image, int threshold, int min_obj_size, int max_obj_size) {

        setup_plain(input_image, threshold, min_obj_size, max_obj_size);// initialization

        analyze_plain(); // Paramarray is filled up

        // extract the object list from Paramarray
        ArrayList<double[]> soma_list = new ArrayList<double[]>(ID);
        double voxCount, intensity;
        for (int i = 1; i < ID; i++) {
            voxCount = Paramarray[i][0];
            intensity = Paramarray[i][2]; // sum over all intensity values
            if (voxCount >= minSize && voxCount <= maxSize) {
                if (voxCount != 0 && intensity != 0) {
                    // String[] head = {"Volume", "Surface", "Intensity", "Centre X", "Centre Y", "Centre Z", "Centre int X", "Centre int Y", "Centre int Z"};
                    double[] circle_details = new double[4];
                    circle_details[0] = Paramarray[i][3]; // x
                    circle_details[1] = Paramarray[i][4]; // y
                    circle_details[2] = Paramarray[i][5]; // z
                    circle_details[3] = Math.pow((3*Paramarray[i][0])/(4*Math.PI), 1f/3); // r (radius based on the voxel count, volume of the sphere)
                    soma_list.add(circle_details);
                }
            }
            IJ.showStatus("Calculating barycenters' coordinates");
            IJ.showProgress(i, ID);
        }

        return soma_list;

    }

}

// for extracting connected components in 2D
class ConnComponents2D extends ParticleAnalyzer {

    public ConnComponents2D(
            ImagePlus binary_image,
            int options,
            int measurements,
            double minSize,
            double maxSize
    ){
        super(options, measurements, new ResultsTable(), minSize, maxSize);
        setup("", binary_image);
        run(binary_image.getProcessor());
    }

    public ArrayList<double[]> getObjectList() { // additional method to loop through the ResultsTable and extract the scores

        ArrayList<double[]> out = new ArrayList<double[]>();

        // adding soma traces to the list and to the tag_map
        for (int i = 0; i < rt.getCounter(); i++) {

            double[] x_y_r = new double[3];

            x_y_r[0] = rt.getColumn(ResultsTable.X_CENTER_OF_MASS)[i];
            x_y_r[1] = rt.getColumn(ResultsTable.Y_CENTER_OF_MASS)[i];
            x_y_r[2] = (float) Math.sqrt(rt.getColumn(ResultsTable.AREA)[i]/Math.PI);

            out.add(x_y_r); // here trace has only one element

        }

        return out;

    }

//    opts = ParticleAnalyzer.SHOW_NONE;
//    int meas = Measurements.CENTROID|Measurements.CENTER_OF_MASS|Measurements.AREA|Measurements.PERIMETER;
//    int maxSz = Math.round(W*H*0.25f);
//    ResultsTable rt = new ResultsTable();
//    ParticleAnalyzer pa = new ParticleAnalyzer(opts, meas, rt, minSz, maxSz);
//    pa.setup("", blur_temp); // binary image as an argument
//    pa.run(blur_temp.getProcessor());
//    System.out.print("->" + rt.getCounter() + " soma(s)... ");

}
