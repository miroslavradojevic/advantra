package MLdetection;

import ij.ImagePlus;
import imagescience.image.Axes;
import imagescience.image.Coordinates;
import imagescience.image.Dimensions;
import imagescience.image.Image;
import java.util.ArrayList;

public class FeaturePool {

int N; //size of the window

ArrayList<HaarLikeFeature> featurepool;

FeaturePool(int N) {
    this.N = N;
    featurepool = new ArrayList<HaarLikeFeature>();
    createPool();
}

public HaarLikeFeature getFeature(int i) {
    return featurepool.get(i);
}
// creates 511225 possible combinations 

private void createPoolAllPossible() {
    for (int y0 = 0; y0 < N; y0++) {
        for (int x0 = 0; x0 < N; x0++) {

            for (int ny = 1; ny <= N - y0; ny++) {
                for (int nx = 1; nx <= N - x0; nx++) {

                    for (int y1 = y0; y1 < ny + y0; y1++) {
                        for (int x1 = x0; x1 < nx + x0; x1++) {

                            for (int j = 1; j <= ny + y0 - y1; j++) {
                                for (int i = 1; i <= nx + x0 - x1; i++) {

                                    featurepool.add(new HaarLikeFeature(x0, y0, nx, ny, x1, y1, i, j));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

private void createPool() {

    
    // half - stripes
    for (int s = 2; s <= N; s += 2) {
        for (int y0 = 0; y0 < N - s + 1; y0++) {
            for (int x0 = 0; x0 < N - s + 1; x0++) {
               // featurepool.add(new HaarLikeFeature(x0, y0, s, s, x0 + s / 2, y0, s / 2, s, -1, 2));
//                featurepool.add(new HaarLikeFeature(x0, y0, s, s, x0, y0, s / 2, s, -1, 2));
                featurepool.add(new HaarLikeFeature(x0, y0, s, s, x0, y0 + s / 2, s, s / 2, -1, 2));
//                featurepool.add(new HaarLikeFeature(x0, y0, s, s, x0, y0, s, s / 2, -1, 2));
            }
        }
    }

    // 1/3 stripes
    for (int s = 3; s <= N; s += 3) {
        for (int y0 = 0; y0 < N - s + 1; y0++) {
            for (int x0 = 0; x0 < N - s + 1; x0++) {
               // featurepool.add(new HaarLikeFeature(x0, y0, s, s, x0 + s / 3, y0, s / 3, s, -1, 3));
                featurepool.add(new HaarLikeFeature(x0, y0, s, s, x0, y0 + s / 3, s, s / 3, -1, 3));
            }
        }
    }

//    // 2/4 stripes 
    for (int s = 4; s <= N; s += 4) {
        for (int y0 = 0; y0 < N - s + 1; y0++) {
            for (int x0 = 0; x0 < N - s + 1; x0++) {
                double w = s / (s / 2.0);
               // featurepool.add(new HaarLikeFeature(x0, y0, s, s, x0 + s / 4, y0, s / 2, s, -1, w));
                featurepool.add(new HaarLikeFeature(x0, y0, s, s, x0, y0 + s / 4, s, s / 2, -1, w));
            }
        }
    }

    // squares 
//    for (int k = 1; k <= (N - 1) / 2; k += 1) {
//        for (int s = 2 * k + 1 ; s <= N; s += 1) {
//            for (int y0 = 0; y0 < N - s + 1; y0++) {
//                for (int x0 = 0; x0 < N - s + 1; x0++) {
//                    double w = s * s / (1.0 * (s - 2 * k) * (s - 2 * k));
//                    featurepool.add(new HaarLikeFeature(x0, y0, s, s, x0 + k, y0 + k, s - 2 * k, s - 2 * k, -1, w));
//                }
//            }
//        }
//    }

}

private void createPoolx2() {

    // half - stripes
    for (int s = 4; s <= N; s += 2) {
        for (int y0 = 0; y0 < N - s + 1; y0++) {
            for (int x0 = 0; x0 < N - s + 1; x0++) {
                featurepool.add(new HaarLikeFeature(x0, y0, s, s, x0 + s / 2, y0, s / 2, s, -1, 2));
                featurepool.add(new HaarLikeFeature(x0, y0, s, s, x0, y0 + s / 2, s, s / 2, -1, 2));
            }
        }
    }

    // 1/3 stripes
    for (int s = 6; s <= N; s += 3) {
        for (int y0 = 0; y0 < N - s + 1; y0++) {
            for (int x0 = 0; x0 < N - s + 1; x0++) {
                featurepool.add(new HaarLikeFeature(x0, y0, s, s, x0 + s / 3, y0, s / 3, s, -1, 3));
                featurepool.add(new HaarLikeFeature(x0, y0, s, s, x0, y0 + s / 3, s, s / 3, -1, 3));
            }
        }
    }
////
    // 2/4 stripes 
    for (int s = 8; s <= N; s += 4) {
        for (int y0 = 0; y0 < N - s + 1; y0++) {
            for (int x0 = 0; x0 < N - s + 1; x0++) {
                featurepool.add(new HaarLikeFeature(x0, y0, s, s, x0 + s / 4, y0, s / 2, s, -1, 2));
                featurepool.add(new HaarLikeFeature(x0, y0, s, s, x0, y0 + s / 4, s, s / 2, -1, 2));
            }
        }
    }

    // squares 
    for (int k = 2; k <= (N - 1) / 2; k += 1) {
        for (int s = 2 * (k + 1); s <= N; s += 1) {
            for (int y0 = 0; y0 < N - s + 1; y0++) {
                for (int x0 = 0; x0 < N - s + 1; x0++) {
                    double w = s * s / (1.0 * (s - 2 * k) * (s - 2 * k));
                    featurepool.add(new HaarLikeFeature(x0, y0, s, s, x0 + k, y0 + k, s - 2 * k, s - 2 * k, -1, w));
                }
            }
        }
    }
}

public void visualizeFeaturePool() {
    int size = featurepool.size();

    Dimensions dims = new Dimensions(N, N, 1, size, 1);
    Coordinates cout = new Coordinates();
    Image outimg = Image.create(dims, "imagescience.image.ByteImage");
    outimg.axes(Axes.X + Axes.Y);


    for (int k = 0; k < size; k++) {
        int[] r0 = featurepool.get(k).getWhiteSquare();
        int[] r1 = featurepool.get(k).getBlackSquare();
        double[][] imageZ = new double[N][N];
        for (int j = 1; j <= r0[3]; j++) {
            for (int i = 1; i <= r0[2]; i++) {
                imageZ[r0[1] + j - 1][r0[0] + i - 1] = 128;
            }

        }
        for (int j = 1; j <= r1[3]; j++) {
            for (int i = 1; i <= r1[2]; i++) {
                imageZ[r1[1] + j - 1][r1[0] + i - 1] = 255;
            }

        }
        cout.t = k;
        outimg.set(cout, imageZ);
    }

    ImagePlus imp = outimg.imageplus();
    imp.show();
    imp.getCanvas().zoomIn(0, 0);
    imp.getCanvas().zoomIn(0, 0);
    imp.getCanvas().zoomIn(0, 0);
    imp.getCanvas().zoomIn(0, 0);
    imp.getCanvas().zoomIn(0, 0);
    imp.getCanvas().zoomIn(0, 0);
    imp.getCanvas().zoomIn(0, 0);

}

public int getSize() {
    return featurepool.size();
}
}
