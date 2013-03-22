package MLdetection;

/**
 *
 * Haar like feature, contains r0=(x0, y0, u0, v0) for the "white" (-1) rectangle 
 * and r1=(x1, y1, u1, v1) for the "black" square (+1)
 */
public class HaarLikeFeature {

private int x0,  y0,  u0,  v0,  x1,  y1,  u1,  v1;
private double w0,  w1;

HaarLikeFeature() {
    this.x0 = 0;
    this.y0 = 0;
    this.u0 = 0;
    this.v0 = 0;

    this.x1 = 0;
    this.y1 = 0;
    this.u1 = 0;
    this.v1 = 0;
}

HaarLikeFeature(int x0, int y0, int u0, int v0, int x1, int y1, int u1, int v1) {
    this.x0 = x0;
    this.y0 = y0;
    this.u0 = u0;
    this.v0 = v0;

    this.x1 = x1;
    this.y1 = y1;
    this.u1 = u1;
    this.v1 = v1;
}

HaarLikeFeature(int x0, int y0, int u0, int v0, int x1, int y1, int u1, int v1, 
        double w0, double w1) {
    this.x0 = x0;
    this.y0 = y0;
    this.u0 = u0;
    this.v0 = v0;

    this.x1 = x1;
    this.y1 = y1;
    this.u1 = u1;
    this.v1 = v1;

    this.w0 = w0;
    this.w1 = w1;

}

public void setWeigths(double w0, double w1) {
    this.w0 = w0;
    this.w1 = w1;
}

public void setX0Y0(int x0, int y0) {
    this.x0 = x0;
    this.y0 = y0;
}

public void setX1Y1(int x1, int y1) {
    this.x1 = x1;
    this.y1 = y1;
}

public void setU0V0(int u0, int v0) {
    this.u0 = u0;
    this.v0 = v0;
}

public void setU1V1(int u1, int v1) {
    this.u1 = u1;
    this.v1 = v1;
}

public int[] getWhiteSquare() {
    return new int[]{x0, y0, u0, v0};
}

public int[] getBlackSquare() {
    return new int[]{x1, y1, u1, v1};
}
public double[] getWeights() {
    return new double[]{w0, w1};
}

}
