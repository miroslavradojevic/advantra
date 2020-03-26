package advantra.trace;

public class SeedPoint3D {
	
	// all that is necessary to know in order to start the new trace
	double[] 	position3D;
	double[] 	orientation3D;
	double		radius;

	public SeedPoint3D(){
		position3D 		= new double[3];
		orientation3D 	= new double[3];
		radius			= 1.0;
	}
	
	public SeedPoint3D(double[] pos3D, double[] ort3D, double r){
		
		if(pos3D.length!=3 || ort3D.length!=3){
			System.err.println("SeedPoint3D:SeedPoint3D(): postion & orientation vectors have to have length 3");
			System.exit(1);
		}
		
		position3D = new double[3];
		orientation3D = new double[3];
		
		position3D[0] = pos3D[0];
		position3D[1] = pos3D[1];
		position3D[2] = pos3D[2];
		
		orientation3D[0] = ort3D[0];
		orientation3D[1] = ort3D[1];
		orientation3D[2] = ort3D[2];
		
		radius 			= r;
		
	}
	
	public SeedPoint3D(double posX, double posY, double posZ, double ortX, double ortY, double ortZ, double r){
		
		position3D = new double[3];
		orientation3D = new double[3];
		
		position3D[0] = posX;
		position3D[1] = posY;
		position3D[2] = posZ;
		
		orientation3D[0] = ortX;
		orientation3D[1] = ortY;
		orientation3D[2] = ortZ;
		
		radius  	= r;
		
	}
	
	public void setSeedPoint(double posX, double posY, double posZ, double ortX, double ortY, double ortZ, double r){
		
		position3D[0] = posX;
		position3D[1] = posY;
		position3D[2] = posZ;
		
		orientation3D[0] = ortX;
		orientation3D[1] = ortY;
		orientation3D[2] = ortZ;
		
		radius		= r;
		
	}
	
	public void setSeedPoint(double[] pos3D, double[] ort3D, double r){
		
		if(pos3D.length!=3 || ort3D.length!=3){
			System.err.println("SeedPoint3D:setSeedPoint(): postion & orientation vectors have to have length 3");
			System.exit(1);
		}
		
		position3D[0] = pos3D[0];
		position3D[1] = pos3D[1];
		position3D[2] = pos3D[2];
		
		orientation3D[0] = ort3D[0];
		orientation3D[1] = ort3D[1];
		orientation3D[2] = ort3D[2];
		
		radius 	= r;
		
	}
	
	public void getSeedPointPosition(double[] position3D){
		
		if(position3D.length!=3){
			System.err.println("SeedPoint3D:getSeedPointPosition():storage arrays for position and orientation have to be 3 dimensional.");
			System.exit(1);
		}
		
		position3D[0] = this.position3D[0];
		position3D[1] = this.position3D[1];
		position3D[2] = this.position3D[2];
		
	}
	
	public void getSeedPointOrientation(double[] orientation3D){
		
		if(orientation3D.length!=3){
			System.err.println("SeedPoint3D:getSeedPointOrientation():storage arrays for position and orientation have to be 3 dimensional.");
			System.exit(1);
		}
		
		orientation3D[0] = this.orientation3D[0];
		orientation3D[1] = this.orientation3D[1];
		orientation3D[2] = this.orientation3D[2];
		
	}
	
	public double[] getSeedPointPosition(){
		double[] out = new double[3];
		out[0] = position3D[0];
		out[1] = position3D[1];
		out[2] = position3D[2];
		return out;
	}

	public double getSeedPointPositionX(){
		return this.position3D[0];
	}
	
	public double getSeedPointPositionY(){
		return this.position3D[1];
	}
	
	public double getSeedPointPositionZ(){
		return this.position3D[2];
	}
	
	public double getSeedPointOrientationX(){
		return this.orientation3D[0];
	}
	
	public double getSeedPointOrientationY(){
		return this.orientation3D[1];
	}
	
	public double getSeedPointOrientationZ(){
		return this.orientation3D[2];
	}
	
 	public double[] getSeedPointOrientation(){
		double[] out = new double[3];
		out[0] = orientation3D[0];
		out[1] = orientation3D[1];
		out[2] = orientation3D[2];
		return out;
	}
	
	public double getRadius(){
		return radius;
	}
}
