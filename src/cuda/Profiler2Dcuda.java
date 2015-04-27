package cuda;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
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
import java.util.Arrays;

import static jcuda.runtime.JCuda.cudaFree;
import static jcuda.runtime.JCuda.cudaMalloc;
import static jcuda.runtime.JCuda.cudaMemcpy;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyDeviceToHost;
import static jcuda.runtime.cudaMemcpyKind.cudaMemcpyHostToDevice;

/**
 * Created by miroslav on 26-3-15.
 */
public class Profiler2Dcuda {

    public float radius;

    public float[][]            kernels;          // Nscales x Ndirections x L x L
    public float[]              kernels_avg;      // average
    public float[][]            kernels_hat;      // kernel-average

    public float[]              kernels1D_hat;
    public float[]              kernels_hat_sum_2;// sum(kernel-average)^2

    public float[]              inimg_cuda;       // input image
    public int[]                i2xy_cuda;        // locations mapping

    // output
    public float[]   i2zncc_cuda;                 // correlations per indexed foreground location
    public float[]   i2sigma_cuda;                // cross-section gaussian sigma per foreground location (the one with largest zncc score)
    public float[]   i2vx_cuda;                   // local orientation per foreground location (the one with largest zncc score)
    public float[]   i2vy_cuda;                   // local orientation per foreground location (the one with largest zncc score)

    public float[] sigmas;
    private float sigma_min = 1f;
    private float sigma_step = .5f;
    private float PI = (float) (Math.PI);
    private float arcRes = 1.0f;

    public int L, N;

    private int blockSize;
    private int n_threads_max;
    private int gridSize;

    public int imgW, imgH;

    public Profiler2Dcuda(
            float _radius,
            int[][] _i2xy,
            float[][] _inimg_xy,
            int _block_size,
            int _max_threads_nr
    )
    {

        n_threads_max = _max_threads_nr;// 1024*20; // pointsArray.size();
        blockSize = _block_size;//1024; // nr. cuda threads (per block)
        gridSize = n_threads_max / blockSize; // nr. blocks

        imgW = _inimg_xy.length;
        imgH = _inimg_xy[0].length;

        i2xy_cuda = new int[2*_i2xy.length]; // 1d storage of foreground locations

        int counter = 0;

        for (int i = 0; i < _i2xy.length; i++) {
            i2xy_cuda[counter++] = _i2xy[i][0];
            i2xy_cuda[counter++] = _i2xy[i][1];
        }

        inimg_cuda = new float[_inimg_xy.length*_inimg_xy[0].length];
        for (int j = 0; j < _inimg_xy[0].length; j++) {
            for (int i = 0; i < _inimg_xy.length; i++) {
                inimg_cuda[i + j * _inimg_xy.length] = _inimg_xy[i][j];
            }
        }

        // allocate outputs
        i2zncc_cuda     = new float[_i2xy.length];
        i2sigma_cuda    = new float[_i2xy.length];
        i2vx_cuda       = new float[_i2xy.length];
        i2vy_cuda       = new float[_i2xy.length];

        // generate kernels based on the radius

        radius = _radius;

        int L2 = (int) Math.ceil(radius);    // transversal sampling limits
        L = 2 * L2 + 1; // how many to take radially with given sampling step

        // sigmas define
        int cnt = 0;
        for (float sg = sigma_min; sg <= .4f*L2; sg+=sigma_step) cnt++;
        sigmas = new float[cnt];
        cnt = 0;
        for (float sg = sigma_min; sg <= .4f*L2; sg+=sigma_step) sigmas[cnt++] = sg;

        // form N theta (N directions)
        N 	= (int) Math.ceil(((PI * radius)/arcRes));

        kernels             = new float[N*sigmas.length][L*L];
        kernels_avg         = new float[N*sigmas.length];
        kernels_hat         = new float[N*sigmas.length][L*L];
        kernels_hat_sum_2   = new float[N*sigmas.length];

        // kernel formation
        for (int i = 0; i < kernels.length; i++) {

            int direc_idx = i % N;
            int scale_idx = i / N;

            float sigx = sigmas[scale_idx];
            float sigy = L2;                      // broader than sigmax
            float ang = direc_idx * (PI / N);

            float vx =  ((float) Math.cos(ang));
            float vy =  ((float) Math.sin(ang));

            kernels_avg[i] = 0; // average

            for (int j = 0; j < kernels[i].length; j++) {

                int xx = j % L;
                int yy = j / L;

                float currx = (xx - L2) *   vx  + (yy - L2) * vy;
                float curry = (xx - L2) * (-vy) + (yy - L2) * vx;

                kernels[i][j] = (float) Math.exp( -(  ((Math.pow(currx,2)/(2*Math.pow(sigx,2))) + (Math.pow(curry,2)/(2*Math.pow(sigy,2)))   ) )   );
                kernels_avg[i] += kernels[i][j];

            }

            kernels_avg[i] /= kernels[i].length;

            kernels_hat_sum_2[i] = 0;

            for (int j = 0; j < kernels[i].length; j++) {

                kernels_hat[i][j] = kernels[i][j] - kernels_avg[i];
                kernels_hat_sum_2[i] += Math.pow(kernels_hat[i][j], 2);

            }
        }

        // 1d format for CUDA specially
        kernels1D_hat   = new float[L * L * N * sigmas.length];

        cnt = 0;

        for (int scale_idx = 0; scale_idx < sigmas.length; scale_idx++) {
            for (int dir_idx = 0; dir_idx < N; dir_idx++) {
                for (int y_idx = 0; y_idx < L; y_idx++) {
                    for (int x_idx = 0; x_idx < L; x_idx++) {
                        kernels1D_hat[cnt]  = kernels_hat[ scale_idx * N + dir_idx][y_idx * L + x_idx];
                        cnt++;
                    }
                }
            }
        }

    }

    public ImagePlus getKernels() {

        ImageStack iskernels = new ImageStack(L, L);

        for (int i = 0; i <kernels.length; i++) {
            int direc_idx = i % N;
            int scale_idx = i / N;
            float ang = direc_idx * (PI / N);
            iskernels.addSlice(sigmas[scale_idx]+","+ IJ.d2s(ang, 1), new FloatProcessor(L, L, kernels[i]));
        }

        return new ImagePlus("", iskernels);

    }

    public void run() { // use cuda to calculate output values in parallel

//        System.out.print("run CUDA... ");

        // inimg_cuda upload on GPU
        Pointer p_inimg_cuda = new Pointer();
        cudaMalloc(p_inimg_cuda, inimg_cuda.length * Sizeof.FLOAT);
        cudaMemcpy(p_inimg_cuda, Pointer.to(inimg_cuda), inimg_cuda.length * Sizeof.FLOAT, cudaMemcpyHostToDevice);

        // i2xy_cuda upload on GPU
        Pointer p_i2xy_cuda = new Pointer();
        cudaMalloc(p_i2xy_cuda, i2xy_cuda.length * Sizeof.INT);
        cudaMemcpy(p_i2xy_cuda, Pointer.to(i2xy_cuda), i2xy_cuda.length * Sizeof.INT, cudaMemcpyHostToDevice);

        // kernels upload on GPU - make it 1d vector before uploading
        Pointer p_kernels1D_hat = new Pointer();
        cudaMalloc(p_kernels1D_hat, kernels1D_hat.length * Sizeof.FLOAT);
        cudaMemcpy(p_kernels1D_hat, Pointer.to(kernels1D_hat), kernels1D_hat.length * Sizeof.FLOAT, cudaMemcpyHostToDevice);

        // kernels_hat_sum_2 upload
        Pointer p_kernels_hat_sum_2 = new Pointer();
        cudaMalloc(p_kernels_hat_sum_2, kernels_hat_sum_2.length * Sizeof.FLOAT);
        cudaMemcpy(p_kernels_hat_sum_2, Pointer.to(kernels_hat_sum_2), kernels_hat_sum_2.length * Sizeof.FLOAT, cudaMemcpyHostToDevice);

        // i2zncc_cuda upload on GPU
        Pointer p_i2zncc_cuda = new Pointer();
        cudaMalloc(p_i2zncc_cuda, i2zncc_cuda.length * Sizeof.FLOAT);
        cudaMemcpy(p_i2zncc_cuda, Pointer.to(i2zncc_cuda), i2zncc_cuda.length * Sizeof.FLOAT, cudaMemcpyHostToDevice);

        // i2sigma_cuda upload on GPU
        Pointer p_i2sigma_cuda = new Pointer();
        cudaMalloc(p_i2sigma_cuda, i2sigma_cuda.length * Sizeof.FLOAT);
        cudaMemcpy(p_i2sigma_cuda, Pointer.to(i2sigma_cuda), i2sigma_cuda.length * Sizeof.FLOAT, cudaMemcpyHostToDevice);

        // i2vx_cuda
        Pointer p_i2vx_cuda = new Pointer();
        cudaMalloc(p_i2vx_cuda, i2vx_cuda.length * Sizeof.FLOAT);
        cudaMemcpy(p_i2vx_cuda, Pointer.to(i2vx_cuda), i2vx_cuda.length * Sizeof.FLOAT, cudaMemcpyHostToDevice);

        // i2vy_cuda
        Pointer p_i2vy_cuda = new Pointer();
        cudaMalloc(p_i2vy_cuda, i2vy_cuda.length * Sizeof.FLOAT);
        cudaMemcpy(p_i2vy_cuda, Pointer.to(i2vy_cuda), i2vy_cuda.length * Sizeof.FLOAT, cudaMemcpyHostToDevice);

        // Create the kernelLauncher that will execute the kernel
        KernelLauncher kernelLauncher = KernelLauncher.compile(getSourceCodeFromFile("zncc.cu"), "applyKernels");
        kernelLauncher.setGridSize(gridSize, 1);
        kernelLauncher.setBlockSize(blockSize, 1, 1);

        // because of the timeout for kernel execution of about 5-7 secs, each scale is computed separately
        for (int si = 0; si < sigmas.length; si++) {
            kernelLauncher.call(
                    p_kernels1D_hat,
                    p_kernels_hat_sum_2,
                    p_inimg_cuda,
                    p_i2zncc_cuda,
                    p_i2sigma_cuda,
                    p_i2vx_cuda,
                    p_i2vy_cuda,
                    p_i2xy_cuda,
                    L,
                    N,
                    si,
                    imgW,
                    imgH,
                    n_threads_max,
                    i2zncc_cuda.length
            );
        }

        cudaMemcpy(Pointer.to(i2zncc_cuda), p_i2zncc_cuda,  i2zncc_cuda.length * Sizeof.FLOAT, cudaMemcpyDeviceToHost);
        cudaMemcpy(Pointer.to(i2sigma_cuda),p_i2sigma_cuda, i2sigma_cuda.length * Sizeof.FLOAT, cudaMemcpyDeviceToHost);
        cudaMemcpy(Pointer.to(i2vx_cuda),   p_i2vx_cuda,    i2vx_cuda.length * Sizeof.FLOAT, cudaMemcpyDeviceToHost);
        cudaMemcpy(Pointer.to(i2vy_cuda), p_i2vy_cuda, i2vy_cuda.length * Sizeof.FLOAT, cudaMemcpyDeviceToHost);

        cudaFree(p_inimg_cuda);
        cudaFree(p_i2xy_cuda);
        cudaFree(p_kernels1D_hat);
        cudaFree(p_kernels_hat_sum_2);
        cudaFree(p_i2zncc_cuda);
        cudaFree(p_i2sigma_cuda);
        cudaFree(p_i2vx_cuda);
        cudaFree(p_i2vy_cuda);

//        System.out.println("done.");
    }

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

    public void clean() {

        for (int i = 0; i < kernels.length; i++) kernels[i] = null;
        kernels = null;

        kernels_avg = null;

        for (int i = 0; i < kernels_hat.length; i++) kernels_hat[i] = null;
        kernels_hat = null;

        kernels1D_hat = null;

        kernels_hat_sum_2 = null;

        inimg_cuda = null;

        i2xy_cuda = null;

        i2zncc_cuda = null;

        i2sigma_cuda = null;

        i2vx_cuda = null;

        i2vy_cuda = null;

    }

}
