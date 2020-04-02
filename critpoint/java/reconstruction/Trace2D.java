package reconstruction;

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

    public Trace2D(int _tag_init, int _tag_prev){

        tag_init   = _tag_init; // first available tag given by the tracer
        tag_prev   = _tag_prev; // tag that was obtained while tracing

        // last tag will be determined with the number of elements
        // so:
        // tag_init    ...(0)... tag_prev
        // tag_init+1  ...(1)... tag_init
        // tag_init+2  ...(2)... tag_init+1
        // tag_init+3  ...(3)... tag_init+2
        // ...
        // ...
        // tag_init+(L-1) (L-1)  tag_init+(L-1)-1

        // follow up tag will therefore always be tag_init+locs_x.length() - that's set outside
        // for each individual trace it becomes important to know the initial tag and predcessor
        // initial tag and predcessor tag is known at the moment of creating the trace
        // last tag is not necessary to book keep here as it can be inferred from the list length

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

//    public boolean check_trace(){
//        return (tag_prev!=Integer.MIN_VALUE)&&(locs_x.size()>0);
//    }

}
