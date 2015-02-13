package reconstruction2d;

import java.util.ArrayList;

/**
 * format to store traces in an image
 * Created by miroslav on 12-2-15.
 */
public class Trace2D {

    public ArrayList<Float> locs_x;
    public ArrayList<Float> locs_y;
    public ArrayList<Float> rads;

    public int tag_start;
    public int tag_end;

    public Trace2D(int _tag_start){

        tag_start   = _tag_start;
        tag_end     = -1;

        locs_x      = new ArrayList<Float>();
        locs_y      = new ArrayList<Float>();
        rads        = new ArrayList<Float>();

    }

    public void add(float x, float y, float r) {

        locs_x.add(x);
        locs_y.add(y);
        rads.add(r);

    }

    public void terminate(int _tag_end) {
        tag_end = _tag_end;
    }

    public void add_to_map(int[] _map) {
        // will update tag map over the trace
        // update with tag indexes incremented starting from the start tag
        // TODO

    }

}
