/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package lambertmrev;

import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.FastMath;
import org.orekit.errors.OrekitException;
import org.orekit.frames.Frame;
import org.orekit.frames.TopocentricFrame;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.DateComponents;
import org.orekit.time.TimeComponents;
import org.orekit.time.TimeScale;
import org.orekit.utils.PVCoordinates;
import org.apache.commons.math3.stat.regression.*;
import org.orekit.orbits.CartesianOrbit;
import org.orekit.utils.Constants;

/**
 *
 * @author zittersteijn
 */
public class conversions {

    
    public static double[] geo2radec(PVCoordinates obj, TopocentricFrame staF, Frame inertialFrame, AbsoluteDate epoch) {

        Vector3D rho = new Vector3D(0, 0, 0);

        try {
            rho = obj.getPosition().subtract(staF.getPVCoordinates(epoch, inertialFrame).getPosition());
        } catch (OrekitException ex) {
            Logger.getLogger(conversions.class.getName()).log(Level.SEVERE, null, ex);
        }

        double rho_mag = rho.getNorm();
        double DEC = FastMath.asin(rho.getZ() / rho_mag);
        double cosRA = 0.0;
        double sinRA = 0.0;
        double RA = 0.0;

        Vector3D v_site = new Vector3D(0, 0, 0);
        try {
            v_site = staF.getPVCoordinates(epoch, inertialFrame).getVelocity();
        } catch (OrekitException ex) {
            Logger.getLogger(conversions.class.getName()).log(Level.SEVERE, null, ex);
        }

        Vector3D rhoDot = obj.getVelocity().subtract(v_site);

        if (FastMath.sqrt(FastMath.pow(rho.getX(), 2) + FastMath.pow(rho.getY(), 2)) != 0) {
            cosRA = rho.getX() / FastMath.sqrt(FastMath.pow(rho.getX(), 2) + FastMath.pow(rho.getY(), 2));
            sinRA = rho.getY() / FastMath.sqrt(FastMath.pow(rho.getX(), 2) + FastMath.pow(rho.getY(), 2));
            RA = FastMath.atan2(sinRA, cosRA);
            if (RA <= 0) {
                RA = RA + 2 * FastMath.PI;
            }
        } else {
            sinRA = rhoDot.getY() / FastMath.sqrt(FastMath.pow(rhoDot.getX(), 2) + FastMath.pow(rhoDot.getY(), 2));
            cosRA = rhoDot.getX() / FastMath.sqrt(FastMath.pow(rhoDot.getX(), 2) + FastMath.pow(rhoDot.getY(), 2));
            RA = FastMath.atan2(sinRA, cosRA);
            if (RA <= 0) {
                RA = RA + 2 * FastMath.PI;
            }
        }

        double rhoDot_mag = rho.dotProduct(rhoDot) / rho_mag;
        double RAdot = (rhoDot.getX() * rho.getY() - rhoDot.getY() * rho.getX()) / (-1 * FastMath.pow(rho.getY(), 2) - FastMath.pow(rho.getX(), 2));
        double DECdot = (rhoDot.getZ() - rhoDot_mag * FastMath.sin(DEC)) / FastMath.sqrt(FastMath.pow(rho.getX(), 2) + FastMath.pow(rho.getY(), 2));

        double[] out = {RA, RAdot, DEC, DECdot, rho_mag, rhoDot_mag};

        return out;
    }

    public static Vector3D radec2geo(Tracklet tracklet, double rho, TopocentricFrame staF, Frame inertialFrame, TimeScale utc) {

        double rho_k = rho * FastMath.sin(tracklet.getDEC());
        double rho_j = rho * FastMath.cos(tracklet.getDEC()) * FastMath.sin(tracklet.getRA());
        double rho_i = rho * FastMath.cos(tracklet.getDEC()) * FastMath.cos(tracklet.getRA());

        Vector3D rho_vec = new Vector3D(rho_i, rho_j, rho_k);
//        System.out.println("rho used: " +  rho_vec);

        AbsoluteDate year = new AbsoluteDate(Main.YEAR, utc);
        AbsoluteDate epoch = new AbsoluteDate(year, tracklet.getDOY() * 24 * 3600);

        Vector3D out = new Vector3D(0, 0, 0);

        try {
//            System.out.println("at epoch " + epoch + "\t the station is at: " + staF.getPVCoordinates(epoch, inertialFrame).getPosition());
            out = rho_vec.add(staF.getPVCoordinates(epoch, inertialFrame).getPosition());
        } catch (OrekitException e) {
            System.out.println("station coordinates not computed");
        }

        return out;
    }

    public static CartesianOrbit hyp2ele(double[] hyp, Vector<Tracklet> subSet, constraints set, TopocentricFrame staF, Frame inertialFrame, TimeScale utc) {
        int N = subSet.size();
        AbsoluteDate year = new AbsoluteDate(Main.YEAR, utc);
        double tof = FastMath.abs(subSet.elementAt(0).getDOY() - subSet.elementAt(N - 1).getDOY()) * 24 * 3600;

//        System.out.println("RA\t:" + subSet.elementAt(0).getRA() + "\tDEC:\t" + subSet.elementAt(0).getDEC() + "\thyp:\t" + hyp[0]);
        
        Vector3D geoPos1 = conversions.radec2geo(subSet.elementAt(0), hyp[0], staF, inertialFrame, utc);
        Vector3D geoPos2 = conversions.radec2geo(subSet.elementAt(N - 1), hyp[1], staF, inertialFrame, utc);

//        System.out.println("x1:\t" + geoPos1.getX() + "\ty1:\t" + geoPos1.getY() + "\tz1:\t" + geoPos1.getZ());

        // define epoch of first tracklet
        AbsoluteDate epoch0 = new AbsoluteDate(year, subSet.elementAt(0).getDOY() * 24 * 3600);

        // perform lambert solver
        Lambert solve = new Lambert();
        solve.lambert_problem(geoPos1, geoPos2, tof, Constants.EIGEN5C_EARTH_MU, Boolean.FALSE, 1);
        RealMatrix v1 = solve.get_v1();

        // get the orbital elements at epoch 0
        Vector3D v1vec = new Vector3D(v1.getEntry(0, 0), v1.getEntry(0, 1), v1.getEntry(0, 2));
        PVCoordinates PVgeo1 = new PVCoordinates(geoPos1, v1vec);
        CartesianOrbit orbit = new CartesianOrbit(PVgeo1, inertialFrame, epoch0, Constants.EIGEN5C_EARTH_MU);

//        System.out.println("a:\t" + orbit.getA() + "\ti:\t" + orbit.getI() + "\te:\t" + orbit.getE());
        return orbit;
    }

    public static RealVector kMat2Vec(int[][] in) {
        int dim = in.length;
        RealVector vec = new ArrayRealVector(dim);
        for (int ii = 0; ii < dim; ii++) {
            for (int jj = 0; jj < dim; jj++) {
                if (in[ii][jj] == 1) {
                    double tmp = (double) jj;
//                    System.out.println("there is a one on row: " + ii + " and col: " + tmp);
                    vec.setEntry(ii, tmp);
//                    System.out.println("vector at: " + ii + " is " + vec.getEntry(ii));
                }
            }
        }
        return vec;
    }

    public static int[][] vec2kMat(RealVector in) {
        int dim = in.getDimension();
        int[][] mat = new int[dim][dim];
        for (int ii = 0; ii < dim; ii++) {
            int tmp = (int) in.getEntry(ii);
            mat[ii][tmp] = 1;
        }
        return mat;
    }
}
