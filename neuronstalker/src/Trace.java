import java.util.ArrayList;

/**
 * Created by miroslav on 8-3-15.
 */
public class Trace {

    public ArrayList<float[]>   locs; // [x y] or [x y z]
    public ArrayList<Float>     rads;
    public ArrayList<Integer>   prev;
    public ArrayList<Integer>   curr;

    public int type;

//    public int tag_init;
//    public int tag_prev;

    static int UNDEFINED    = 0;
    static int SOMA         = 1;
    static int AXON         = 2;
    static int LOOP         = 6;
    static int LOOSE        = 7; // those that don't reach either of soma(s)

    public Trace() {

        type = UNDEFINED;

//        tag_init = Integer.MIN_VALUE;
//        tag_prev = Integer.MIN_VALUE;

        locs = new ArrayList<float[]>();
        rads = new ArrayList<Float>();
        prev = new ArrayList<Integer>();
        curr = new ArrayList<Integer>();

    }

    public Trace(int _type) { //int _tag_init, int _tag_prev, int _type

        type = _type;

//        tag_init = _tag_init;
//        tag_prev = _tag_prev;

        locs = new ArrayList<float[]>();
        rads = new ArrayList<Float>();
        prev = new ArrayList<Integer>();
        curr = new ArrayList<Integer>();
    }

    public void add(int prev_label, int curr_label,     float x, float y,               float r) {
        locs.add(new float[]{x, y});
        rads.add(r);
        prev.add(prev_label);
        curr.add(curr_label);
    }

    public void add(int prev_label, int curr_label,     float x, float y, float z,      float r) {
        locs.add(new float[]{x, y, z});
        rads.add(r);
        prev.add(prev_label);
        curr.add(curr_label);
    }

}
