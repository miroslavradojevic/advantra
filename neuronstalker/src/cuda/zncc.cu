extern "C"
__global__ void applyKernels(   float* kernels_hat,
                                float* kernels_hat_sum_2,
                                float* inimg,
                                float* pos2zncc,
                                float* pos2sigma,
                                float* pos2vx,
                                float* pos2vy,
                                int*   pos,
                                int L,
                                int N,
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

        int x = pos[2 * tid];
        int y = pos[2 * tid + 1];
        int L2 = L / 2;

        float imgVal;
        float imgValsAvg = 0;

        for (int pidx = 0; pidx < L*L; pidx++) {

            int xx = x + (pidx % L) - L2;
            int yy = y + (pidx / L) - L2;

            imgVal = (xx>=0 && xx<dimsx && yy>=0 && yy<dimsy)?inimg[xx + yy * dimsx]:0;

            imgValsAvg += imgVal;

        }

        imgValsAvg /= (float) (L*L);

        for (int didx = 0; didx < N; didx++) {

            float num = 0;
            float den = 0;
            float zncc;

            for (int pidx = 0; pidx < L*L; pidx++) {

                // sample from the image
                int xx = x + (pidx % L) - L2;
                int yy = y + (pidx / L) - L2;

                imgVal = (xx>=0 && xx<dimsx && yy>=0 && yy<dimsy)?inimg[xx + yy * dimsx]:0;
 
                num += (imgVal - imgValsAvg) * kernels_hat[scale*N*L*L + didx*L*L + pidx];
                den += (imgVal - imgValsAvg) * (imgVal - imgValsAvg);

            }

            zncc = num / (float) sqrtf(den * kernels_hat_sum_2[didx + scale*N]);

            if (zncc > pos2zncc[tid]) {
                
                pos2zncc[tid] = zncc;
                
                pos2sigma[tid] = scale;

                float ang = didx * (3.141592654 / N);

                float vx = -sinf(ang);
                float vy =  cosf(ang);

                pos2vx[tid] = vx;
                pos2vy[tid] = vy; 

            }

        }

        tid += N_threads;

    }

 }