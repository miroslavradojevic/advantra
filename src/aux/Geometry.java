package aux;

/**
 * Created by miroslav on 3/23/14.
 * some common geometry operations in 2d and 3d
 */
public class Geometry {

    /*
    2d
     */

    // euclidean distance point (p) to line segment (l1, l2)
    public static final float point_to_line(float px, float py, float l1x, float l1y, float l2x, float l2y)
    {
        float ll = (l2x-l1x)*(l2x-l1x)+(l2y-l1y)*(l2y-l1y); // length squared
        if (ll<=Float.MIN_VALUE) return point_to_point(px, py, l1x, l1y); // l1 == l2
        // projection of point onto the line
        // it falls where t = [(p-l1) . (l2-l1)] / |l2-l1|^2
        float t = ((px-l1x)*(l2x-l1x)+(py-l1y)*(l2y-l1y)) / ll;
        if (t<0)        return point_to_point(px, py, l1x, l1y); // beyond l1 end of the segment
        else if (t>1.0) return point_to_point(px, py, l2x, l2y); // beyond l2 end of the segment
        // falls on the segment
        float proj_x = l1x + t * (l2x-l1x);
        float proj_y = l1y + t * (l2y-l1y);
        return point_to_point(px, py, proj_x, proj_y);
    }

    // euclidean distance point (p) to point (l)
    public static final float point_to_point(float px, float py, float lx, float ly)
    {
        return (float) Math.sqrt((lx-px)*(lx-px) + (ly-py)*(ly-py));
    }

    public static final void line_seg_intersec(float px, float py, float rx, float ry, float qx, float qy, float sx, float sy, float[] intersec)
    {
        // reference: http://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
        // points are: p, p+r, q, q+s
        // p + tr = q + us
        // t = (q-p) x s / (r x s)
        // u = (q-p) x r / (r x s)
        float r_cross_s = rx*sy-ry*sx;
        float tt = (qx-px)*sy-(qy-py)*sx;
        float uu = (qx-px)*ry-(qy-py)*rx;
        float t = tt / r_cross_s;
        float u = uu / r_cross_s;

        // cases
        if (Math.abs(r_cross_s)<=Float.MIN_VALUE && Math.abs(uu)<=Float.MIN_VALUE) {

            // two lines are  collinear
            float check_r = dot(qx-px, qy-py, rx, ry); // (q-p).r
            float check_s = dot(qx-px, qy-py, sx, sy); // (q-p).s

            if ( (check_r>=0 && check_r<=dot(rx, ry, rx, ry)) || (check_s>=0 && check_s<=dot(sx, sy, sx, sy)) ) {
                // overlapping and collinear
                intersec = null;
            }
            else {
                // disjoint and collinear
                intersec = null;
            }



        }
        else if (Math.abs(r_cross_s)<=Float.MIN_VALUE && !(Math.abs(uu)<=Float.MIN_VALUE)) {
            // parallel and non-intersecting
            intersec = null;
        }
        else if (!(Math.abs(r_cross_s)<=Float.MIN_VALUE)) {

        }



    }

    public static final void line_intersec()
    {
        // reference: http://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
        // points are: p, p+r, q, q+s


    }

    // cross product a x b
    public static final float cross(float ax, float ay, float bx, float by)
    {
        return ax*by-ay*bx;
    }

    // dot product a . b
    public static final float dot(float ax, float ay, float bx, float by)
    {
        return ax*bx+ay*by;
    }



}
