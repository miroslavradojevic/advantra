import java.util.ArrayList;

/**
 * Created by miroslav on 13-3-15.
 */
public class Node {

    public float[]  loc;
    public float    r;
    public ArrayList<Integer> neighbors;

    // dummy initialization
    public Node() {
        loc         = null;
        r           = -1;
        neighbors   = null;
    }

    // 2d nodes
    public Node(float xn, float yn, float rn) {
        loc = new float[2];
        loc[0] = xn;
        loc[1] = yn;
        r      = rn;
        neighbors = new ArrayList<Integer>();
    }

    public Node(float xn, float yn, float rn, int nbr1){
        this(xn, yn, rn);
        neighbors.add(nbr1);
    }

    public Node(float xn, float yn, float rn, int nbr1, int nbr2){
        this(xn, yn, rn, nbr1);
        neighbors.add(nbr2);
    }

    // 3d nodes
    public Node(float xn, float yn, float zn, float rn) { // no neighbours at the moment
        loc = new float[3];
        loc[0]  = xn;
        loc[1]  = yn;
        loc[2]  = zn;
        r       = rn;
        neighbors = new ArrayList<Integer>();
    }

    public Node(float xn, float yn, float zn, float rn, int nbr1) { // know only prev
        this(xn, yn, zn, rn);
        neighbors.add(nbr1);
    }

    public Node(float xn, float yn, float zn, float rn, int nbr1, int nbr2) {// knew prev, know the next
        this(xn, yn, zn, rn, nbr1);
        neighbors.add(nbr2);
    }

}
