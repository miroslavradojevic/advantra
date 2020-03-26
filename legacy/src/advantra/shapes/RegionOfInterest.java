package advantra.shapes;

public class RegionOfInterest {

	// class that generalizes all the shapes taken from the image stack
	// center of the region
	protected 	double x;
	protected 	double y;
	protected 	double z;
	
	protected 	RoiType roi_type;	// = RoiType.GENERAL_ROI;
	
	public static enum RoiType{GENERAL_ROI, SPHERE, CYLINDER, CONE, CONE_CUT, POINT, ORIENTED_PROJECTIVE_PLANE};
	
	public RegionOfInterest(double x, double y, double z){
		
		roi_type = RoiType.GENERAL_ROI;
		
		this.x = x;
		this.y = y;
		this.z = z;
		
	}

	public void setCenter(double x, double y, double z){
		this.x = x;
		this.y = y;
		this.z = z;
	}

	public double getCenterX(){
		return this.x;
	}	

	public double getCenterY(){
		return this.y;
	}	

	public double getCenterZ(){
		return this.z;
	}	
	
	public RoiType getRoiType(){
		return roi_type;
	}
	
	
	
}
