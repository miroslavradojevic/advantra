package detection2d;

/**
 * Created by miroslav on 6/8/14.
 */
public class CritpointRegion {

	public float radius;
	public float[] centroid;
	public float score;
	public RegionType type;
	public float[][] outward_directions;


	public CritpointRegion(RegionType _type, float _centroidX, float _centroidY, float _radius, float _score, float[][] _outward_directions) {

		this.type = _type;
		this.centroid = new float[2];
		this.centroid[0] = _centroidX;
		this.centroid[1] = _centroidY;
		this.radius = _radius;
		this.score = _score;
		this.outward_directions = new float[_outward_directions.length][_outward_directions[0].length];
		for (int i = 0; i < _outward_directions.length; i++) {
			for (int j = 0; j < _outward_directions[0].length; j++) {
				this.outward_directions[i][j] = _outward_directions[i][j];
			}
		}

	}

	public boolean isOverlapping(float[] _other_centroid, float _other_radius) {

		float dist = (float) Math.sqrt(Math.pow(_other_centroid[0]-this.centroid[0],2) + Math.pow(_other_centroid[1]-this.centroid[1],2));
		if (dist <= this.radius + _other_radius) return true;
		else return false;

	}

	public boolean isOverlaping(CritpointRegion _other_critpoint_region) {

		float dist = (float) Math.sqrt(Math.pow(_other_critpoint_region.centroid[0]-this.centroid[0],2) + Math.pow(_other_critpoint_region.centroid[1],2));
		if (dist <= this.radius + _other_critpoint_region.radius) return true;
		else return false;

	}

	public enum RegionType {
		END, BIF, CROSS, BIF_CROSS; // BIF_CROSS represents both, need it for some function calls clearness, appears as a pseudo-category here
	}


}

