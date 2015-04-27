package cuda;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import jcuda.Pointer;
import jcuda.Sizeof;
import jcuda.utils.KernelLauncher;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;

import static jcuda.runtime.JCuda.cudaFree;
import static jcuda.runtime.JCuda.cudaMalloc;
import static jcuda.runtime.JCuda.cudaMemcpy;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyDeviceToHost;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyHostToDevice;

/**
 * Created by miroslav on 19-3-15.
 */
public class CudaKernelFiltering implements PlugIn {


    public void run(String s) {

        OpenDialog dlg = new OpenDialog("Select Text File");
        String path = dlg.getPath();
        if (path == null) return;

        ImagePlus imp = IJ.openImage(path);
        imp.show();

        float[][] image = imp.getProcessor().getFloatArray();
        int dimsx = image.length;
        int dimsy = image[0].length;

        // create 1D array for input and output images for GPU
        float[] inimg_cuda =  new float[dimsx * dimsy];
        float[] outimg_cuda = new float[dimsx * dimsy];
        for (int j = 0; j < dimsy; j++) {
            for (int i = 0; i < dimsx; i++) {
                inimg_cuda[i + j * dimsx] = (float)image[i][j];
            }
        }

        // upload images to GPU
        Pointer pointer_inimg = new Pointer();
        cudaMalloc(pointer_inimg, inimg_cuda.length * Sizeof.FLOAT);
        cudaMemcpy(pointer_inimg, Pointer.to(inimg_cuda), inimg_cuda.length * Sizeof.FLOAT, cudaMemcpyHostToDevice);

        Pointer pointer_outimg = new Pointer();
        cudaMalloc(pointer_outimg, outimg_cuda.length * Sizeof.FLOAT);
        cudaMemcpy(pointer_outimg, Pointer.to(outimg_cuda), outimg_cuda.length * Sizeof.FLOAT, cudaMemcpyHostToDevice);

        System.out.println("done");

        int nrOfOrientations = 60;
        int maskSize = 31;
        int nrOfScales = 3;
        float sigma_min = 1.f;

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // create kernels
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        System.out.println("creating kernels");
        // Allocate memory on the device, and copy the host data to the device
        float[] kernels1D = new float[maskSize * maskSize * nrOfOrientations * nrOfScales]; // in one row, [first all x, then y, then orientations, then scales]
        Pointer pointer_kernels = new Pointer();
        cudaMalloc(pointer_kernels, kernels1D.length * Sizeof.FLOAT);
        cudaMemcpy(pointer_kernels, Pointer.to(kernels1D), kernels1D.length * Sizeof.FLOAT, cudaMemcpyHostToDevice);




        int blockSize = nrOfOrientations; // nr. of cuda threads per block
        int gridSize = nrOfScales; // nr. of blocks



        KernelLauncher kernelLauncher = KernelLauncher.compile(getSourceCodeFromFile("myCUDA.cu"), "createKernels");
        kernelLauncher.setGridSize(gridSize, 1);
        kernelLauncher.setBlockSize(blockSize, 1, 1);
        kernelLauncher.call(
                pointer_kernels,
                maskSize,
                nrOfOrientations,
                nrOfScales,
                sigma_min,
                nrOfOrientations * nrOfScales);

        // read the output from the GPU
        cudaMemcpy(Pointer.to(kernels1D), pointer_kernels, kernels1D.length * Sizeof.FLOAT, cudaMemcpyDeviceToHost);
        System.out.println("done creating kernels");
        visualizeKernels(kernels1D, maskSize, nrOfOrientations * nrOfScales, "Orientation filters");

//        if (true) return;

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // filter the input image
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // create an array with coordinates
        ArrayList<int[]> pointsArray = new ArrayList<int[]>();
        for (int j = maskSize; j < dimsy - maskSize; j++) {
            for (int i = maskSize; i < dimsx - maskSize; i++) {
                pointsArray.add(new int[]{i, j});
            }
        }

        // create a 1D array of (x,y) locations for which filtering will be done (in the example below it is a whole image
        int[] points1D = new int[pointsArray.size() * 2];
        int counter = 0;
        for (int[] p : pointsArray) {
            points1D[counter++] = p[0];
            points1D[counter++] = p[1];
        }

        // push to GPU
        Pointer pointer_points = new Pointer();
        cudaMalloc(pointer_points, points1D.length * Sizeof.INT);
        cudaMemcpy(pointer_points, Pointer.to(points1D), points1D.length * Sizeof.INT, cudaMemcpyHostToDevice);

        int n_threads_max = 8192;//pointsArray.size();
        blockSize = 256; // nr. of cuda threads
        gridSize = n_threads_max / blockSize;

        // Create the kernelLauncher that will execute the kernel
        kernelLauncher = KernelLauncher.compile(getSourceCodeFromFile("myCUDA.cu"), "applyKernels");
        kernelLauncher.setGridSize(gridSize, 1);
        kernelLauncher.setBlockSize(blockSize, 1, 1);

        // because of the timeout for kernel execution of about 5-7 secs, each scale is computed separately
        for (int si = 0; si < nrOfScales; si++)
        {
            kernelLauncher.call(
                    pointer_kernels,
                    pointer_inimg,
                    pointer_outimg,
                    pointer_points,
                    maskSize,
                    nrOfOrientations,
                    si,
                    dimsx,
                    dimsy,
                    n_threads_max,
                    pointsArray.size());
            System.out.println("Computing scale = " + si);
        }

        cudaMemcpy(Pointer.to(outimg_cuda), pointer_outimg, outimg_cuda.length * Sizeof.FLOAT, cudaMemcpyDeviceToHost);
        ImageProcessor ip = new FloatProcessor(dimsx, dimsy, outimg_cuda);
        ImagePlus imp_out = new ImagePlus("Output", ip);
        imp_out.show();

        cudaFree(pointer_kernels);
        cudaFree(pointer_points);
        cudaFree(pointer_inimg);
        cudaFree(pointer_outimg);

    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////// EXTRA METHODS ///////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    public String getSourceCodeFromFile(String cuFileName)
    {
        // Obtain the CUDA source code from the .cu file
        InputStream cuInputStream = getClass().getResourceAsStream(cuFileName);
        String sourceCode = "";
        try
        {
            BufferedReader br = new BufferedReader(new InputStreamReader(cuInputStream));
            StringBuilder sb = new StringBuilder();
            while (true)
            {
                String line = br.readLine();
                if (line == null)
                {
                    break;
                }
                sb.append(line).append("\n");
//                System.out.println(sb.toString());
            }
            sourceCode = sb.toString();

        }
        catch (IOException e)  {
            e.printStackTrace();
        }
        finally {
            try {cuInputStream.close();}
            catch (IOException e)
            { e.printStackTrace();}
        }
        return sourceCode;
    }

    public static void visualizeKernels(float[] kernel, int size, int N, String name)
    {
        ImageStack ims = new ImageStack(size, size);

        for (int i = 0; i < N; i++) {
            float[][] image = new float[size][size];
            for (int ii = 0; ii < size; ii++) {
                for (int jj = 0; jj < size; jj++) {
                    image[ii][jj] = (float)kernel[i * size * size + ii + size * jj];
                }
            }
            FloatProcessor ip_temp = new FloatProcessor(image);
            ims.addSlice(ip_temp);
        }
        ImagePlus imp_temp = new ImagePlus(name,ims);
        imp_temp.show();
        for (int i = 0; i < 256 / size; i++) {
            imp_temp.getWindow().getCanvas().zoomIn(0, 0);
        }
    }

}
