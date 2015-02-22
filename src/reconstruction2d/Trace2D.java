package reconstruction2d;

import ij.ImagePlus;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.util.ArrayList;

/**
 * format to store traces in an image
 * Created by miroslav on 12-2-15.
 */
public class Trace2D {

    public ArrayList<Float> locs_x;
    public ArrayList<Float> locs_y;
    public ArrayList<Float> rads;

    public int tag_init; // will be the first available tag
    public int tag_prev; // will be the tag of the region where it ended

    // typical usage
    // set_init(tag_init)   first available tag, using the counting machine
    // add(x,y,r)           as many times as it needs - depending on the tracer2D, it stops at some existing tag
    // set_prev(tag_prev)   set that tag at which the trace stopped

    public Trace2D(int _tag_init){

        tag_init   = _tag_init; // first available tag given by the tracer
        tag_prev   = Integer.MIN_VALUE; // yet to be set by a separate method

        locs_x      = new ArrayList<Float>();
        locs_y      = new ArrayList<Float>();
        rads        = new ArrayList<Float>();

    }

    public void add(float x, float y, float r) {
        locs_x.add(x);
        locs_y.add(y);
        rads.add(r);
    }

    public void add(float[] x, float[] y, float[] r) {
        for (int i = 0; i < x.length; i++) {
            add(x[i], y[i], r[i]);
        }
    }

    public void set_prev(int _tag_prev) {
        tag_prev = _tag_prev;
    }

    public boolean check_trace(){
        return (tag_prev!=Integer.MIN_VALUE)&&(locs_x.size()>0);
    }

}
