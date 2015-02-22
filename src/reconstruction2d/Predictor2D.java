package reconstruction2d;

/**
 * will make predictions in 2D and assign prior probabilities to each predicted state
 * state = [x, y, vx, vy, sigma]
 * state -> state' = [x', y', vx', vy', sigma']
 * states are predicted using semi circle geometry, semi circle has radius defined with 'step'
 * prior probabilities are assigned with respect to the divergence from the root direction [vx, vy]
 * and with respect to how sigma' changed with respect to the sigma of the root state
 * Created by miroslav on 22-2-15.
 */

public class Predictor2D {

    public float[][]            p; // 2d locations x,y
    public float[][]            v; // 2d direction vx, vy
    public float[]              s; // sigma gaussian
    //----------------------------------------------------
    public float[]              w; // priorweight assigned wrt. the direction divergence and radius change
    public int                  NN; // number of states capturing the circle
    public float                radius;
    private static float        arcRes = 1f; // sampling resolution along arc

    public Predictor2D(float _radius) {
        radius = _radius;
    }

}
