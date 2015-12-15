/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package lambertmrev;

import org.apache.commons.math3.geometry.euclidean.threed.*;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.*;
import org.orekit.frames.*;
import org.orekit.time.*;
import org.orekit.utils.*;
import org.orekit.orbits.*;              
import org.orekit.propagation.analytical.*;


public class Lambert{
    Vector3D default_r1 = new Vector3D(1.0,0.0,0.0);
    Vector3D default_r2 = new Vector3D(1.0,0.0,0.0);
    private RealMatrix m_v1;

    /// Constructor
    /** Constructs and solves a Lambert problem.
     *
     * \param[in] R1 first cartesian position
     * \param[in] R2 second cartesian position
     * \param[in] tof time of flight
     * \param[in] mu gravity parameter
     * \param[in] cw when 1 a retrograde orbit is assumed
     * \param[in] multi_revs maximum number of multirevolutions to compute
     */

    public void lambert_problem(Vector3D r1, Vector3D r2, double tof, double mu, Boolean cw, int multi_revs)
{
	// sanity checks
	if(tof <= 0){
	    System.out.println("ToF is negative! \n");
	}
	if(mu<=0){
	    System.out.println("mu is below zero");
	}

	// 1 - getting lambda and T
	double m_c = FastMath.sqrt((r2.getX() - r1.getX())*(r2.getX() - r1.getX())+ (r2.getY()-r1.getY())*(r2.getY()-r1.getY()) + (r2.getZ()-r1.getZ())*(r2.getZ()-r1.getZ()));
      	double R1 = r1.getNorm();
	double R2 = r2.getNorm();
	double m_s = (m_c + R1 + R2) / 2.0;

	Vector3D ir1 = r1.normalize();
	Vector3D ir2 = r2.normalize();
	Vector3D ih = Vector3D.crossProduct(ir1,ir2);
	ih = ih.normalize();

	if(ih.getZ() == 0){
	    System.out.println("angular momentum vector has no z component \n");
	}
	double lambda2 = 1.0 - m_c / m_s;
	double m_lambda = FastMath.sqrt(lambda2);

	Vector3D it1 = new Vector3D(0.0,0.0,0.0);
	Vector3D it2 = new Vector3D(0.0,0.0,0.0);

	if(ih.getZ() < 0.0){ // Transfer angle is larger than 180 degrees as seen from abive the z axis
	    m_lambda = -m_lambda;
	    it1 = Vector3D.crossProduct(ir1,ih);
	    it2 = Vector3D.crossProduct(ir2,ih);
	} else {
	    it1 = Vector3D.crossProduct(ih,ir1);
	    it2 = Vector3D.crossProduct(ih,ir2);
	}
	it1.normalize();
	it2.normalize();

	if(cw){ // Retrograde motion
	    m_lambda = -m_lambda;
	    it1.negate();
	    it2.negate();
	}
	double lambda3 = m_lambda * lambda2;
	double T = FastMath.sqrt(2.0*mu/m_s/m_s/m_s)*tof;

	// 2 - We now hava lambda, T and we will find all x
	// 2.1 - let us first detect the maximum number of revolutions for which there exists a solution
	int m_Nmax = FastMath.toIntExact(FastMath.round(T/FastMath.PI)); 
	double T00 = FastMath.acos(m_lambda) + m_lambda*FastMath.sqrt(1.0-lambda2);
	double T0 = (T00 + m_Nmax * FastMath.PI);
	double T1 = 2.0/3.0 * (1.0 - lambda3);
	double DT = 0.0; 
	double DDT = 0.0;
	double DDDT = 0.0;

	if(m_Nmax > 0){
	    if(T < T0){ // We use Halley iterations to find xM and TM
		int it = 0;
		double err = 1.0;
		double T_min = T0;
		double x_old = 0.0, x_new = 0.0;
		while(true){
                    ArrayRealVector deriv = dTdx(x_old, T_min, m_lambda);
		    DT = deriv.getEntry(0);
		    DDT = deriv.getEntry(1);
		    DDDT = deriv.getEntry(2);
		    if(DT != 0.0){
			x_new = x_old - DT * DDT / (DDT * DDT - DT * DDDT / 2.0);
		    }
		    err = FastMath.abs(x_old-x_new);
		    if((err < 1e-13) || (it > 12) ){
			break;
		    }
		    tof = x2tof(x_new,m_Nmax, m_lambda);
		    x_old = x_new;
		    it++;
		}
		if(T_min > T){
		    m_Nmax -= 1;
		} 
	    }
	}
	// We exit this if clause with Mmax being the maximum number of revolutions
	// for which there exists a solution. We crop it to multi_revs
	m_Nmax = FastMath.min(multi_revs,m_Nmax);

	// 2.2 We now allocate the memory for the output variables
	m_v1 = MatrixUtils.createRealMatrix(m_Nmax*2+1, 3);
	RealMatrix m_v2 = MatrixUtils.createRealMatrix(m_Nmax*2+1, 3);
	RealMatrix m_iters = MatrixUtils.createRealMatrix(m_Nmax*2+1, 3);
	//RealMatrix m_x = MatrixUtils.createRealMatrix(m_Nmax*2+1, 3);
        ArrayRealVector m_x = new ArrayRealVector(m_Nmax*2+1);
               
	// 3 - We may now find all solution in x,y
	// 3.1 0 rev solution
	// 3.1.1 initial guess
	if(T >= T00){
	    m_x.setEntry(0, -(T - T00)/(T - T00 + 4));
	}else if(T <= T1){
	    m_x.setEntry(0, T1*(T1-T) / ( 2.0/5.0*(1-lambda2*lambda3) * T ) + 1);
	}else {
	    m_x.setEntry(0, FastMath.pow((T/T00),0.69314718055994529 / FastMath.log(T1/T00)) - 1.0);
	}
	// 3.1.2 Householder iterations
	//m_iters.setEntry(0, 0, housOutTmp.getEntry(0));
	m_x.setEntry(0, householder(T, m_x.getEntry(0), 0, 1e-5, 15, m_lambda));
      
	// 3.2 multi rev solutions
	double tmp;
        double x0;

	for(int i = 1; i< m_Nmax+1; i++){
	    // 3.2.1 left householder iterations
	    tmp = FastMath.pow((i*FastMath.PI+FastMath.PI) / (8.0*T), 2.0/3.0);
	    m_x.setEntry(2*i-1, (tmp-1)/(tmp+1));
	    x0 = householder(T, m_x.getEntry(2*i-1), i, 1e-8, 15, m_lambda);
	    m_x.setEntry(2*i-1, x0);
	    //m_iters.setEntry(2*i-1, 0, housOutTmp.getEntry(0));

	    //3.2.1 right Householder iterations
	    tmp = FastMath.pow((8.0*T)/(i*FastMath.PI), 2.0/3.0);
	    m_x.setEntry(2*i, (tmp-1)/(tmp+1));
	    x0 = householder(T, m_x.getEntry(2*i), i, 1e-8, 15, m_lambda);
	    m_x.setEntry(2*i, x0);
	    //m_iters.setEntry(2*i, 0, housOutTmp.getEntry(0));
	}
	
	// 4 - For each found x value we recontruct the terminal velocities
	double gamma = FastMath.sqrt(mu*m_s/2.0);
	double rho = (R1-R2) / m_c;
        
	double sigma = FastMath.sqrt(1-rho*rho);
	double vr1,vt1,vr2,vt2,y;

	ArrayRealVector ir1_vec = new ArrayRealVector(3);
	ArrayRealVector ir2_vec = new ArrayRealVector(3);
	ArrayRealVector it1_vec = new ArrayRealVector(3);
	ArrayRealVector it2_vec = new ArrayRealVector(3);

	// set Vector3D values to a mutable type
	ir1_vec.setEntry(0, ir1.getX());
        ir1_vec.setEntry(1, ir1.getY());
	ir1_vec.setEntry(2, ir1.getZ());
	ir2_vec.setEntry(0, ir2.getX());
	ir2_vec.setEntry(1, ir2.getY());
	ir2_vec.setEntry(2, ir2.getZ());
	it1_vec.setEntry(0, it1.getX());
	it1_vec.setEntry(1, it1.getY());
	it1_vec.setEntry(2, it1.getZ());
	it2_vec.setEntry(0, it2.getX());
	it2_vec.setEntry(1, it2.getY());
	it2_vec.setEntry(2, it2.getZ());

	for(int i = 0; i<m_x.getDimension(); i++){
	    y = FastMath.sqrt(1.0-lambda2+lambda2*m_x.getEntry(i)*m_x.getEntry(i));
            vr1 = gamma *((m_lambda*y-m_x.getEntry(i))-rho*(m_lambda*y+m_x.getEntry(i)))/R1;           
            vr2 = -gamma*((m_lambda*y-m_x.getEntry(i))+rho*(m_lambda*y+m_x.getEntry(i)))/R2;
         
            double vt = gamma*sigma*(y+m_lambda*m_x.getEntry(i));
            
            vt1 = vt/R1;
            vt2 = vt/R2;
            
            for (int j=0; j<3;++j) m_v1.setEntry(i,j, vr1 * ir1_vec.getEntry(j) + vt1 * it1_vec.getEntry(j));             
            for (int j=0; j<3;++j) m_v2.setEntry(i,j, vr2 * ir2_vec.getEntry(j) + vt2 * it2_vec.getEntry(j));
	}


    }

    // define the methods
    public double householder(double T, double x0, int N, double eps, int iter_max, double m_lambda){
	int it = 0;
	double err = 1.0;
	double xnew=0.0;
	double tof=0.0, delta=0.0,DT=0.0,DDT=0.0,DDDT=0.0;
	ArrayRealVector out = new ArrayRealVector();
	while ( (err>eps) && (it < iter_max) )
	{
	    tof = x2tof(x0,N, m_lambda);
            ArrayRealVector deriv = dTdx(x0, tof, m_lambda);
	    DT = deriv.getEntry(0);
	    DDT = deriv.getEntry(1);
	    DDDT = deriv.getEntry(2);
	    delta = tof-T;
	    double DT2 = DT*DT;
	    xnew = x0 - delta * (DT2-delta*DDT/2.0) / (DT*(DT2-delta*DDT) + DDDT*delta*delta/6.0);
	    err=FastMath.abs(x0-xnew);
	    x0=xnew;
	    it++;
	}

	return x0;
    }


    public ArrayRealVector dTdx(double x, double T, double m_lambda){
	double l2 = m_lambda*m_lambda;
	double l3 = l2 * m_lambda;
	double umx2 = 1.0 - x*x;
	double y = FastMath.sqrt(1.0-l2*umx2);
	double y2 = y*y;
	double y3 = y2*y;
        
        ArrayRealVector out = new ArrayRealVector(3);
	out.setEntry(0, 1.0/umx2 * (3.0*T*x-2.0+2.0*l3*x/y));
        out.setEntry(1,1.0 / umx2 * (3.0*T+5.0*x*out.getEntry(0)+2.0*(1.0-l2)*l3/y3));
        out.setEntry(2, 1.0 / umx2 * (7.0*x*out.getEntry(1) +8.0*out.getEntry(0)-6.0*(1.0-l2)*l2*l3*x/y3/y2));
	
	return out;
    }
    

    public double x2tof2(double x, int N, double m_lambda){
	double a = 1.0 / (1.0-x*x);
	double tof;
	if (a>0)	//ellipse
	    {
		double alfa = 2.0*FastMath.acos(x);
		double beta = 2.0 * FastMath.asin (FastMath.sqrt(m_lambda*m_lambda/a));
		if (m_lambda<0.0) beta = -beta;
		tof =  ((a * FastMath.sqrt (a)* ( (alfa - FastMath.sin(alfa)) - (beta - FastMath.sin(beta)) + 2.0*FastMath.PI*N)) / 2.0);
	    }
	else
	    {
		double alfa = 2.0*FastMath.acosh(x);
		double beta = 2.0 * FastMath.asinh(FastMath.sqrt(-m_lambda*m_lambda/a));
		if (m_lambda<0.0) beta = -beta;
		tof =  ( -a * FastMath.sqrt (-a)* ( (beta - FastMath.sinh(beta)) - (alfa - FastMath.sinh(alfa)) ) / 2.0);
	    }
	return tof;
    }

    public double hypergeometricF(double z, double tol){
	double Sj=1.0;
	double Cj=1.0;
	double err=1.0;
	double Cj1=0.0;
	double Sj1=0.0;
	int j=0;
	while (err > tol){
		Cj1 = Cj*(3.0+j)*(1.0+j)/(2.5+j)*z/(j+1);
		Sj1 = Sj + Cj1;
		err=FastMath.abs(Cj1);
		Sj = Sj1;
		Cj=Cj1;
		j=j+1;
	}
	return Sj;
    }

    public double x2tof(double x, int N, double m_lambda){
	double battin = 0.01;
	double lagrange = 0.2;
	double dist = FastMath.abs(x-1);
	double tof;
	if (dist < lagrange && dist > battin) { // We use Lagrange tof expression
	    tof = x2tof2(x,N,m_lambda);
		return tof;
	}
	double K = m_lambda*m_lambda;
	double E = x*x-1.0;
	double rho = FastMath.abs(E);
	double z = FastMath.sqrt(1+K*E);
	if (dist < battin) { // We use Battin series tof expression
		double eta = z-m_lambda*x;
		double S1 = 0.5*(1.0-m_lambda-x*eta);
		double Q = hypergeometricF(S1,1e-11);
		Q = 4.0/3.0*Q;
		tof = (eta*eta*eta*Q+4.0*m_lambda*eta)/2.0 + N*FastMath.PI / FastMath.pow(rho,1.5);
		return tof;
	} else { // We use Lancaster tof expresion
		double y=FastMath.sqrt(rho);
		double g = x*z - m_lambda*E;
		double d = 0.0;
		if (E<0) {
			double l = FastMath.acos(g);
			d=N*FastMath.PI+l;
		} else {
			double f = y*(z-m_lambda*x);
			d=FastMath.log(f+g);
		}
		tof = (x-m_lambda*z-d/y)/E;
		return tof;
	}
	//return tof;
    }

    // define methods to get values
    public RealMatrix get_v1(){
	return m_v1;
    }
      
}

