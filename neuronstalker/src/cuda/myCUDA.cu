extern "C"
__global__ void createKernels(
float* kernels,
int size,
int nrOfOrientations,
int nrOfScales,
float sigma_min,
int N)
{
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    int orientation = threadIdx.x;
    int scale = blockIdx.x;

    if (index < N)
    {
        int s2 = size / 2;
        int nn = 11;
        float gamma = 20;
        // see: http://en.wikipedia.org/wiki/Gabor_filter

        float alpha = (3.141592654 * orientation) / nrOfOrientations;
        float sigma_adjusted = 0.6 * sigma_min * scale + sigma_min;
        float s2d_min = 2 * sigma_adjusted * sigma_adjusted;
        float s2d_max = 200;

        float lambda = 5 * sigma_adjusted;

        float totalSum = 0;
        for (int j = -s2; j <= s2; j++) {
            for (int i = -s2; i <= s2; i++) {
                float sum = 0;
                for (int ii = 0; ii < nn; ii++) {
                    float xx = i - 0.5 + (1 + 2 * ii) / (2.0 * nn);
                    for (int jj = 0; jj < nn; jj++)
                    {
                        float yy = j - 0.5 + (1 + 2 * jj) / (2.0 * nn);

                        float xx_ = yy * sinf(alpha) + xx * cosf(alpha);
                        float yy_ = yy * cosf(alpha) - xx * sinf(alpha);

                        //sum += expf(-(xx_ * xx_) / s2d_max - (yy_ * yy_) / s2d_min) * cosf(2 * 3.141592654 * yy_ / (size));
                        // Gabor filter
                        sum += expf(-(xx_ * xx_) / s2d_min / gamma  - (yy_ * yy_) / s2d_min) * cosf(2 * 3.141592654 * yy_ / lambda);
                    }
                }

                kernels[(i + s2) + (j + s2) * size + size * size * orientation + size * size * nrOfOrientations * scale] = sum;
                totalSum += sum;
            }
        }
        for (int j = -s2; j <= s2; j++) {
            for (int i = -s2; i <= s2; i++) {
                //kernels[(i + s2) + (j + s2) * size + size * size * orientation + size * size * nrOfOrientations * scale] /= totalSum;
            }
        }
    }
}

extern "C"
__global__ void applyKernels(
float* kernels,
float* inimg,
float* outimg,
int* positions,
int size,
int nrOfOrientations,
int scale,
int dimsx,
int dimsy,
int N_threads,
int NN)
{
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    int tid = index;

    while (tid < NN)
    {
        int x = positions[2 * tid];
        int y = positions[2 * tid + 1];
        int s2 = size / 2;

        int bestOrientation = 0;
        float bestResponse = 0;
        for (int r = 0; r < nrOfOrientations; r++) {
            float sum = 0;
            for (int j = -s2; j <= s2; j++) {
                for (int i = -s2; i <= s2; i++) {
                    sum += inimg[x + i + dimsx * (y + j)] * kernels[i + s2 + size * (j + s2) + size * size * r + size * size * nrOfOrientations * scale];
                }
            }
            if (sum > bestResponse) {
                bestResponse = sum;
                bestOrientation = r;
            }
        }
        if (outimg[x + dimsx * y] < bestResponse)
        {
            outimg[x + dimsx * y] = bestResponse;
        }

        tid += N_threads;
    }

 }