package detection2d;

/**
 * Created by miroslav on 11-11-14.
 */
public class SemiCircle {

    public float[][]     p; // 2d locations, x,y
    public float[][]     v; // 2d orientation,vx,vy
    public float[]       w; // prior: weight assigned wrt the direction (geomrtry used as prior in tracing)
    public int           NN;
    public float         radius;
    private static float arcRes 	        = 1.0f; // default sampling resolution along arc

    public SemiCircle()
    {
        radius = Float.NaN;
        NN = 0;
        p = null;
        v = null;
        w = null;
    }

    public SemiCircle(float radius)
    {

        this.radius = radius;
        // radius will define the number of points tha will cover pi angle
        NN = (int) Math.ceil( ( (1 * Math.PI * radius) / arcRes) );
        p = new float[NN][2];
        v = new float[NN][2];
        w = new float[NN];
    }

    public void reset(float radius)
    {
        this.radius = radius;
        NN = (int) Math.ceil( ( (1 * Math.PI * radius) / arcRes) );
        p = new float[NN][2];
        v = new float[NN][2];
        w = new float[NN];
    }

    public void reset(float radius, int Nset)
    {
        this.radius = radius;
        NN = Nset;
        p = new float[NN][2];
        v = new float[NN][2];
        w = new float[NN];
    }

    public void set(float px, float py, float vx, float vy, float sigma_deg)
    {

        float sigma_rad = (float) ((sigma_deg/180f)*Math.PI);

        // p, v, w
        for (int i = 0; i < NN; i++) {

            float alfa = (float) (-Math.PI/2 + Math.PI * ((float)i/(NN-1))); // -pi/2 -- +pi/2

            v[i][0] = (float) (vx * Math.cos(alfa) - vy * Math.sin(alfa)); // vx
            v[i][1] = (float) (vx * Math.sin(alfa) + vy * Math.cos(alfa)); // vy

            p[i][0] = px + radius * v[i][0]; // x
            p[i][1] = py + radius * v[i][1]; // y


            w[i] = (float) Math.exp(-Math.pow(alfa,2)/(2*sigma_rad*sigma_rad));

        }


    }

}
