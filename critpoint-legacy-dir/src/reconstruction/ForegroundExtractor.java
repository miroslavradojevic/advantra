package reconstruction;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ByteProcessor;

import java.util.ArrayList;

/**
 * will use Masker2D and Masker3D to separate out foreground points
 * the essential usage of the Masker2D/Masker3D class
 * output is the list of foreground locations and corresponding image
 *
 * Created by miroslav on 3-3-15.
 */

public class ForegroundExtractor {

    public ArrayList<boolean[][]>   mask_xy;   // mask at each layer - location x,y
    public ArrayList<int[][]>		i2xy;       // mapping index -> loc, index -> xy,xyz
    public ArrayList<int[][]>		xy2i;       // mapping loc -> index, xy,xyz -> index

    public ForegroundExtractor(){
        mask_xy    = new ArrayList<boolean[][]>();
        i2xy        = new ArrayList<int[][]>();
        xy2i        = new ArrayList<int[][]>();
    }

    public void work2D(
            float[][]   inimg_xy,
            int         margin,
            float       radius_check,
            float       percentile,
            String      _midresults_dir
    )
    {

        Masker2D.loadTemplate(
                    inimg_xy,
                    margin, // (int)Math.ceil(radius_check),
                    radius_check,
                    percentile); //image, margin, check, percentile

        int totalLocs = inimg_xy.length*inimg_xy[0].length;

        int CPU_NR = Runtime.getRuntime().availableProcessors() + 1;

        Masker2D ms_jobs[] = new Masker2D[CPU_NR];

        for (int i = 0; i < ms_jobs.length; i++) {
                ms_jobs[i] = new Masker2D(i*totalLocs/CPU_NR,  (i+1)*totalLocs/CPU_NR);
                ms_jobs[i].start();
        }

        for (int i = 0; i < ms_jobs.length; i++) {
                try {
                    ms_jobs[i].join();
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
        }

        Masker2D.defineThreshold();
        Masker2D.formRemainingOutputs();

        // append outputs to the output variables
//        mask_xy.clear();
        mask_xy.add(Masker2D.mask_xy.clone());

//        i2xy.clear();
        i2xy.add(Masker2D.i2xy.clone());

//        xy2i.clear();
        xy2i.add(Masker2D.xy2i.clone());

        if (!_midresults_dir.equalsIgnoreCase("")) {
            IJ.saveAs(getMask(), "Tiff", _midresults_dir + "mask.tif");
        }

        Masker2D.clean();

    }

    public void work3D(
            ArrayList<float[][]>    inimg_xyz,
            int                     margin,
            float                   radius,
            float                   percentile,
            String                  _midresults_dir
    )
    {

        for (int layIdx = 0; layIdx < inimg_xyz.size(); layIdx++)
            work2D(inimg_xyz.get(layIdx), margin, radius, percentile, ""); // will append each execution, so inimg_xyz.size()

        if (!_midresults_dir.equalsIgnoreCase("")) {
            IJ.saveAs(getMask(),   "Tiff", _midresults_dir + "mask.tif");
        }

//        Masker3D.loadTemplate(
//                inimg_xyz,
//                margin_xy,
//                margin_z,
//                radius_check,
//                percentile,
//                zDist);
//        int totalLocs = inimg_xyz.length*inimg_xyz[0].length*inimg_xyz[0][0].length;
//        int CPU_NR = Runtime.getRuntime().availableProcessors() + 1;
//        Masker3D ms_jobs[] = new Masker3D[CPU_NR];
//        for (int i = 0; i < ms_jobs.length; i++) {
//            ms_jobs[i] = new Masker3D(i*totalLocs/CPU_NR,  (i+1)*totalLocs/CPU_NR);
//            ms_jobs[i].start();
//        }
//        for (int i = 0; i < ms_jobs.length; i++) {
//            try {
//                ms_jobs[i].join();
//            } catch (InterruptedException e) {
//                e.printStackTrace();
//            }
//        }
//        Masker3D.defineThreshold();
//        Masker3D.formRemainingOutputs();
//        mask_loc.clear();
//        for (int i = 0; i < Masker3D.mask_xyz.size(); i++) mask_loc.add(Masker3D.mask_xyz.get(i).clone());
//        i2loc = Masker3D.i2xyz.clone();
//        loc2i.clear();
//        for (int i = 0; i < Masker3D.mask_xyz.size(); i++) loc2i.add(Masker3D.xyz2i.get(i).clone());
//        Masker3D.clean();

    }

    public ImagePlus getMask()
    {

        int W = mask_xy.get(0).length;
        int H = mask_xy.get(0)[0].length;
        int L = mask_xy.size();

        ImageStack is_out = new ImageStack(W,H);
//        System.out.println("exporting " + mask_xy.size());
        for (int i = 0; i < mask_xy.size(); i++) {
            byte[] out = new byte[W*H];
            for (int j=0; j<out.length; j++) out[j] = mask_xy.get(i)[j % W][j / W] ? (byte) 255 : (byte) 0;
            is_out.addSlice(new ByteProcessor(W, H, out));
        }

//        if (dim==2) { // there is 1 slice to add in 2d  case
//            byte[] out = new byte[W*H];
//            for (int i=0; i<out.length; i++) out[i] = mask_loc.get(0)[i % W][i / W] ? (byte) 255 : (byte) 0;
//            is_out.addSlice(new ByteProcessor(W, H, out));
//        }
//        else { // 3d there is sequence of slices
//            for (int i = 0; i < loc2i.size(); i++) {
//                byte[] out = new byte[W*H];
//                for (int j=0; j<out.length; j++) out[j] = mask_loc.get(i)[j % W][j / W] ? (byte) 255 : (byte) 0;
//                is_out.addSlice(new ByteProcessor(W, H, out));
//            }
//        }

        return new ImagePlus("", is_out);

    }

//    public void save_midresults(String tif_path){
//
//        ImagePlus fgrd = new ImagePlus("", getMask());
//        IJ.saveAs(fgrd, "Tiff", tif_path);
//
//
//    }

}
