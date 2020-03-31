package aux;

/**
 * Created by miroslav on 3/23/14.
 * some common euclidean geometry operations in 2d and 3d
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
        float r_cross_s = cross(rx, ry, sx, sy); //r x s
        float u = cross(qx-px, qy-py, rx, ry); // u = (q-p) x r

        // cases
        if (Math.abs(r_cross_s)<=Float.MIN_VALUE && Math.abs(u)<=Float.MIN_VALUE) {

            // two lines are  collinear
            float check_r = dot(qx-px, qy-py, rx, ry); // (q-p).r
            float check_s = dot(qx-px, qy-py, sx, sy); // (q-p).s

            if ( (check_r>=0 && check_r<=dot(rx, ry, rx, ry)) || (check_s>=0 && check_s<=dot(sx, sy, sx, sy)) ) {
                // overlapping and collinear
                intersec[0] = Float.NaN;
                intersec[1] = Float.NaN;
            }
            else {
                // disjoint and collinear
                intersec[0] = Float.NaN;
                intersec[1] = Float.NaN;
            }

        }
        else if (Math.abs(r_cross_s)<=Float.MIN_VALUE && !(Math.abs(u)<=Float.MIN_VALUE)) {
            // parallel and non-intersecting
            intersec[0] = Float.NaN;
            intersec[1] = Float.NaN;
        }
        else { //Math.abs(r_cross_s)!=0

        	// they are not parallel, and they intersect
			float t = cross(qx-px, qy-py, sx, sy) / r_cross_s;
			u = u / r_cross_s;

			if ((t>=0 && t<=1) && (u>=0 && u<=1)) {
				// segments intersect, store the point p + tr = q + us
				intersec[0] = px + t * rx;
				intersec[1] = py + u * ry;
			}
			else {
				// they are not parallel, but do not intersect
				intersec[0] = Float.NaN;
				intersec[1] = Float.NaN;
			}

        }

    }

    public static final void line_intersec(float px, float py, float rx, float ry, float qx, float qy, float sx, float sy, float[] intersec)
    {
        // reference: http://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
        // points are: p, p+r, q, q+s
		// p + tr = q + us
		// t = (q-p) x s / (r x s)
		// u = (q-p) x r / (r x s)
		float r_cross_s = cross(rx, ry, sx, sy); //r x s
		float u = cross(qx-px, qy-py, rx, ry); // u = (q-p) x r

		// cases
		if (Math.abs(r_cross_s)<=Float.MIN_VALUE && Math.abs(u)<=Float.MIN_VALUE) {
			// two lines are  collinear
			intersec[0] = Float.NaN;
			intersec[1] = Float.NaN;
		}
		else if (Math.abs(r_cross_s)<=Float.MIN_VALUE && !(Math.abs(u)<=Float.MIN_VALUE)) {
			// parallel and non-intersecting
			intersec[0] = Float.NaN;
			intersec[1] = Float.NaN;
		}
		else { //Math.abs(r_cross_s)!=0

			// segments intersect
			float t = cross(qx-px, qy-py, sx, sy) / r_cross_s;
			u = u / r_cross_s;

			// store the point p + tr = q + us
			intersec[0] = px + t * rx;
			intersec[1] = py + u * ry;

//			if ((t>=0 && t<=1) && (u>=0 && u<=1)) {
//			}
//			else {
//				intersec[0] = Float.NaN;
//				intersec[1] = Float.NaN;
//			}

		}

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

    public static final void orthogonal(float x1, float x2, float y1, float y2, float[] ort)
    {
        //float l = (float) Math.sqrt(Math.pow(x2-x1, 2)+Math.pow(y2-y1, 2));
        float vx = (x2-x1);///l;
        float vy = (y2-y1);///l;
//        float wx = vy;
        ort[0] = vy;
//        float wy = -vx;
        ort[1] = -vx;
    }



}
