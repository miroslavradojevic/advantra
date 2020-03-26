package advantra.trace.component;

public class NeuronNode {
	// cartesian coordinates of the node
	private double 	x;
	private double 	y;
	private double 	z;
	// radius at that node
	private double 	r;
	// node type
	private Type 	typ;
	
	
	
	// direction 1
	private double	vx1;
	private double	vy1;
	private double	vz1;
	// direction 2
	private double	vx2;
	private double	vy2;
	private double	vz2;	
	// direction 3
	private double	vx3;
	private double	vy3;
	private double	vz3;
	
	// types of the node - depends on the number of neighbouring nodes
	public static enum Type{
		UNDEFINED_POINT, END_POINT, BODY_POINT, BIFURCATION_POINT 
	}
	
	public NeuronNode(double x, double y, double z, double r){
		this.x = x;
		this.y = y;
		this.z = z;
		this.r = r;
		
		this.typ = Type.UNDEFINED_POINT;
		
		this.vx1 = 0;
		this.vy1 = 0;
		this.vz1 = 0;
		
		this.vx2 = 0;
		this.vy2 = 0;
		this.vz2 = 0;	
		
		this.vx3 = 0;
		this.vy3 = 0;
		this.vz3 = 0;
	}
	
	public NeuronNode(){
		this.x 		= 0;
		this.y 		= 0;
		this.z 		= 0;
		this.r 		= 0;
		
		this.typ 	= Type.UNDEFINED_POINT;
		
		this.vx1 	= 0;
		this.vy1 	= 0;
		this.vz1 	= 0;
		
		this.vx2 	= 0;
		this.vy2 	= 0;
		this.vz2 	= 0;	
		
		this.vx3 	= 0;
		this.vy3 	= 0;
		this.vz3 	= 0;
	}
	
	public void set(double x, double y, double z, double r){
		this.x = x;
		this.y = y;
		this.z = z;
		this.r = r;	
		this.typ = Type.UNDEFINED_POINT;
		
		this.vx1 = 0;
		this.vy1 = 0;
		this.vz1 = 0;
		
		this.vx2 = 0;
		this.vy2 = 0;
		this.vz2 = 0;	
		
		this.vx3 = 0;
		this.vy3 = 0;
		this.vz3 = 0;
	}
	
	public void addNeighbour(double point_x, double point_y, double point_z){
		
		double dx 	= point_x-this.x;
		double dy 	= point_y-this.y;
		double dz 	= point_z-this.z;
		double V 	= Math.sqrt(dx*dx+dy*dy+dz*dz);
		
		if(V>0.0005){
			switch(this.typ){
				case UNDEFINED_POINT: // make it an end-point now that it has at least one neighbour
					this.typ 	= Type.END_POINT;
					this.vx1 	= dx/V; // direction 1
					this.vy1 	= dy/V;
					this.vz1 	= dz/V;
					break;
				case END_POINT: // make it a body-point
					this.typ 	= Type.BODY_POINT;
					this.vx2 	= dx/V; // direction 2
					this.vy2 	= dy/V;
					this.vz2 	= dz/V;
					break;
				case BODY_POINT: // make it a bifurcation
					this.typ 	= Type.BIFURCATION_POINT;
					this.vx3 	= dx/V; // direction 3
					this.vy3 	= dy/V;
					this.vz3 	= dz/V;
					break;
				case BIFURCATION_POINT: // make it undefined again - it has more than 3 neighbours ?!
					this.typ 	= Type.UNDEFINED_POINT;
					this.vx3 	= 0;	this.vy3 	= 0;	this.vz3 	= 0;
					this.vx2 	= 0;	this.vy2 	= 0;	this.vz2 	= 0;
					this.vx1 	= 0;	this.vy1 	= 0;	this.vz1 	= 0;
					break;
				default:
					System.err.println("NeuronNode:addNeighbour():"+this.typ+" type of node does not exist.");
					System.exit(1);
					break;	
			}
		}
		else{
			System.out.format("adding neighbour (%5.2f, %5.2f, %5.2f) to the structure was skipped... distance was too small...\n" +
					"starting the new structure", point_x, point_y, point_z);
		}
	}
	
	/*
	 * get values
	 */
	
	public double getX(){
		return this.x;
	}

	public double getY(){
		return this.y;
	}
	
	public double getZ(){
		return this.z;
	}

	public double getR(){
		return this.r;
	}	
	
	public Type getType(){
		return this.typ;
	}
	
	public double[] getV1(){
		double[] out = new double[3];
		out[0] = this.vx1;
		out[1] = this.vy1;
		out[2] = this.vz1;
		return out;
	}
	
	public double[] getV2(){
		double[] out = new double[3];
		out[0] = this.vx2;
		out[1] = this.vy2;
		out[2] = this.vz2;
		return out;
	}

	public double[] getV3(){
		double[] out = new double[3];
		out[0] = this.vx3;
		out[1] = this.vy3;
		out[2] = this.vz3;
		return out;
	}

}
