package detection3d;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;

import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.File;
import java.util.ArrayList;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 11/16/13
 * Time: 4:15 PM
 */
public class Sphere3DTest implements PlugInFilter, MouseListener {


    // fixed sphere parameters
    static float neuronDiameter = 4f;
    static float scale = 1.5f;
    static Sphere3D    sphere = new Sphere3D(neuronDiameter, scale);

    float[][][] img3d_zxy;
    float zDist;
    ImageCanvas canv;

    /*
        terminal command
     */

	public static void main(String[] args) {

		/*
			demo for the Sphere3D class
		 */

		String dir = System.getProperty("user.home")+File.separator+"Sphere3D_Demo";
		File theDir = new File(dir);

		// if the directory does not exist, create it
		if (!theDir.exists()) {
			System.out.println("creating directory: " + dir);
			boolean result = theDir.mkdir();

			if(result) {
				System.out.println("dir created");
			}
		}

		// this is the demo on Sphere3D class
		System.out.println("Sphere3D with neuron diameter "+neuronDiameter+" and scale "+scale+ " formed.");

		/*
			export sphere sampling
		 */
		String fileName = dir+File.separator+"S3D_direction_sampling.swc";
        sphere.printFilterDirections(fileName);
		System.out.println(fileName);

		/*
			export sampling used to filtering the profiles (swcs with sampling points)
		 */
		System.out.println("\nexport filter offsets...\n");
        sphere.printFilterOffsets(dir);

		/*
			export masks (swcs with sampling points)
		 */
		System.out.println("\nexport filter neighbourhood(s)...\n");
        sphere.printMasks(dir);

        /*
            scheme for finding neighbours on the sphere in planar - every stack corresponds to a sample direction in 3d
         */
        ImageStack isMasks = sphere.visualizeMasks();
        ImagePlus  imMasks = new ImagePlus("neighbour_masks", isMasks);
        imMasks.show();
        imMasks.getCanvas().zoomIn(0, 0);
        imMasks.getCanvas().zoomIn(0, 0);
        imMasks.getCanvas().zoomIn(0, 0);
        imMasks.getCanvas().zoomIn(0, 0);
        imMasks.getCanvas().zoomIn(0, 0);
        imMasks.getCanvas().zoomIn(0, 0);
        imMasks.getCanvas().zoomIn(0, 0);
        imMasks.getCanvas().zoomIn(0, 0);
        imMasks.getCanvas().zoomIn(0, 0);
        imMasks.getCanvas().zoomIn(0, 0);

        /*
            3d stack with filter weights used to filter out one direction in 3d
         */

//        sphere.

	}

    public void mouseClicked(MouseEvent e) {

        int clickX = canv.offScreenX(e.getX());
        int clickY = canv.offScreenY(e.getY());
        int clickZ = canv.getImage().getCurrentSlice()-1;

        showProfile(clickX, clickY, clickZ); // detect peaks on selected profile

    }

    @Override
    public void mousePressed(MouseEvent e) {}

    @Override
    public void mouseReleased(MouseEvent e) {}

    @Override
    public void mouseEntered(MouseEvent e) {}

    @Override
    public void mouseExited(MouseEvent e) {}

    public int setup(String s, ImagePlus imagePlus) {

        if (imagePlus==null || imagePlus.getStackSize()<=1) return DONE;

        img3d_zxy = stackToZxyArray(imagePlus.getStack());

        canv = imagePlus.getCanvas();
        canv.addMouseListener(this);

        GenericDialog gd = new GenericDialog("");
        gd.addNumericField("z dist: ", 3.03, 2);
        gd.showDialog();
        if (gd.wasCanceled()) return DONE;
        zDist = (float) gd.getNextNumber();

        /*
            take into account the fact that the structures are elongated in z 3x
         */
        zDist = zDist/3;


        return DOES_32+DOES_8G+NO_CHANGES;

    }

    public void run(ImageProcessor imageProcessor) {
        IJ.log("Sphere3D demo...");
    }

    private void showProfile(int clickX, int clickY, int clickZ) {

        // find indexes
        float[] profile = sphere.extractProfile(clickX, clickY, clickZ, img3d_zxy, zDist);

        // find planisphere profile
        ByteProcessor profileProc = sphere.drawProfile(profile);
        ImagePlus profileImage = new ImagePlus("profile", profileProc);

        // mark the peaks
        ArrayList<int[]> pks = sphere.profilePeaksXY(profile);
        Overlay ov = new Overlay();
        for (int tt=0; tt<pks.size(); tt++) {
            ov.add(new PointRoi(pks.get(tt)[0], pks.get(tt)[1]));
        }
        profileImage.setOverlay(ov);
        profileImage.show();

        // plot neighbouring [x,y,z] coordinates
        ArrayList<int[]> pks_xyz = sphere.profilePeaksXYZ(profile, zDist);
        for (int tt=0; tt<pks_xyz.size(); tt++) {
            IJ.log("peak "+tt+" : x="+pks_xyz.get(tt)[0]+"[px], y="+pks_xyz.get(tt)[1]+"[px], z="+pks_xyz.get(tt)[2]+"[lay]");
        }

    }

    private static float[][][] stackToZxyArray(ImageStack inis) {

        int W = inis.getWidth();
        int H = inis.getHeight();
        int L = inis.getSize();

        float[][][] img3d_zxy = new float[L][][];

        for (int l=0; l<L; l++) {

            img3d_zxy[l] = new float[W][H];
            float[] readSlice = (float[]) inis.getProcessor(l+1).convertToFloat().getPixels();

            for (int ww=0; ww<W; ww++) {
                for (int hh=0; hh<H; hh++) {
                    img3d_zxy[l][ww][hh] = readSlice[hh*W+ww];
                }
            }

        }

        return img3d_zxy;
    }

}
