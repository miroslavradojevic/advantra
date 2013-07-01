import ij.ImagePlus;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.io.File;
import java.io.FilenameFilter;

/**
 * Created with IntelliJ IDEA.
 * User: miroslav
 * Date: 6/22/13
 * Time: 12:43 PM
 */
public class Tools {

	public static ImagePlus convertToFloatImage(
		ImagePlus inim
	)
	{

		int W = inim.getWidth();
		int H = inim.getHeight();

		ImageProcessor ip = new FloatProcessor(W, H);
		for (int i = 0; i < H*W; i++ ) {
			ip.setf(i, inim.getProcessor().getPixelValue(i%W, i/W));
		}

		ImagePlus outim = new ImagePlus(inim.getTitle()+"ToFloat", ip);
		return outim;

	}

	public static String getExtension(String filePath)
	{
		String extension = "";
		int dotIdx = filePath.lastIndexOf('.');
		if (dotIdx > 0)
			extension = filePath.substring(dotIdx+1);
		return  extension;
	}

	public static String removeExtension(String filePath)
	{
		String extension = "";
		int dotIdx = filePath.lastIndexOf('.');
		if (dotIdx > 0)
			extension = filePath.substring(0, dotIdx);
		return  extension;
	}

    public static File[] listFilesEndingWith(
            File dir,
            String suffix
    )
    {
        final String sfx = suffix;
        File[] tif_train = dir.listFiles(
                new FilenameFilter() {
                    public boolean accept(File dir, String name) {
                        return name.toLowerCase().endsWith(sfx);
                    }
                }
        );
        return tif_train;
    }

    public static void extractMoments2D(
            float[][]   patchin,
            double[]    center,
            double[]    theta
    )
    {

        // moments[0] = centroid(x)
        // moments[1] = centroid(y)
        // moments[2] = central_moment(x,x)
        // moments[3] = central_moment(y,y)
        // moments[4] = central_moment(x,y)

        double M00  = 0;
        double M10  = 0;
        double M01  = 0;
        double M11  = 0;
        double M20  = 0;
        double M02 	= 0;

        // 1st order... - M10, M01,...  2nd order mu11, mu20, mu02 etc...
        // image_coordinates has to be 2xN where rows 1..2 are coords
        // image_values      has to be 1xN and contains image intensities

        for (int i = 0; i < patchin.length; i++) {
            for (int j = 0; j < patchin[0].length; j++) {
                M00 += patchin[i][j];
                M10 += patchin[i][j] * i;
                M01 += patchin[i][j] * j;
                M11 += patchin[i][j] * i * j;
                M20 += patchin[i][j] * i * i;
                M02 += patchin[i][j] * j * j;
            }
        }
        // first centroids
        center[0] = M10 / M00; // centroid(x), mi10
        center[1] = M01 / M00; // centroid(y), mi01
        // second central moments
        double mu11 = M11/M00 - center[0]*center[1]; // central_moment(x,y)
        double mu20 = M20/M00 - center[0]*center[0]; // central_moment(x,x)
        double mu02 = M02/M00 - center[1]*center[1]; // central_moment(y,y)

        theta[0] = 0.5 * Math.atan((2*mu11)/(mu20-mu02));

    }

}
