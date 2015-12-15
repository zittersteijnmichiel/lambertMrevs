/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package lambertmrev;

import org.orekit.frames.*;
import org.orekit.time.*;
import org.orekit.utils.*;
import org.orekit.orbits.*;
import java.io.PrintWriter;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import org.orekit.orbits.KeplerianOrbit;
import org.orekit.propagation.analytical.*;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.analytical.KeplerianPropagator;
import org.apache.commons.math3.geometry.euclidean.threed.*;
import org.orekit.data.DataProvidersManager;
import org.orekit.time.TimeScalesFactory;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.random.*;
import org.orekit.orbits.PositionAngle;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author zittersteijn
 */
public class LambertMRev {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // Want to test the Lambert class so you can specify the number of revs for which to compute
        //System.out.print("this is the frames tutorial \n");
        try {
            Frame inertialFrame = FramesFactory.getEME2000();
            TimeScale utc = TimeScalesFactory.getTAI();
            AbsoluteDate initialDate = new AbsoluteDate(2004, 01, 01, 23, 30, 00.000, utc);
            double mu = 3.986004415e+14;

            double a = 24396159;                 // semi major axis in meters
            double e = 0.72831215;               // eccentricity
            double i = Math.toRadians(7);        // inclination
            double omega = Math.toRadians(180);  // perigee argument
            double raan = Math.toRadians(261);   // right ascension of ascending node
            double lM = 0;                       // mean anomaly

            Orbit initialOrbit = new KeplerianOrbit(a, e, i, omega, raan, lM, PositionAngle.MEAN, inertialFrame, initialDate, mu);
        //KeplerianPropagator kepler = new KeplerianPropagator(initialOrbit);

            // set geocentric positions
            Vector3D r1 = new Vector3D(-6.88999e3, 3.92763e4, 2.67053e3);
            Vector3D r2 = new Vector3D(-3.41458e4, 2.05328e4, 3.44315e3);
            Vector3D r1_site = new Vector3D(4.72599e3, 1.26633e3, 4.07799e3);
            Vector3D r2_site = new Vector3D(4.70819e3, 1.33099e3, 4.07799e3);

            // get the topocentric positions
            Vector3D top1 = Transform.geo2radec(r1.scalarMultiply(1000), r1_site.scalarMultiply(1000));
            Vector3D top2 = Transform.geo2radec(r2.scalarMultiply(1000), r2_site.scalarMultiply(1000));

            // time of flight in seconds
            double tof = 3 * 3600;

            // propagate to 0 and tof
            Lambert test = new Lambert();

            boolean cw = false;
            int multi_revs = 1;
            RealMatrix v1_mat;
            Random randomGenerator = new Random();
            
            PrintWriter out_a = new PrintWriter("out_java_a.txt");
            PrintWriter out_e = new PrintWriter("out_java_e.txt");
            PrintWriter out_rho1 = new PrintWriter("out_java_rho1.txt");
            PrintWriter out_rho2 = new PrintWriter("out_java_rho2.txt");

            // start the loop
            double A, Ecc, rho1, rho2, tof_hyp;

            long time1 = System.nanoTime();
            for (int ll = 0; ll < 1e6; ll++) {

                rho1 = top1.getZ() / 1000 + 1e-3 * randomGenerator.nextGaussian() * top1.getZ() / 1000;
                rho2 = top2.getZ() / 1000 + 1e-3 * randomGenerator.nextGaussian() * top2.getZ() / 1000;
            //tof_hyp = FastMath.abs(tof + 0.1*3600 * randomGenerator.nextGaussian());
                // from topo to geo
                Vector3D r1_hyp = Transform.radec2geo(top1.getX(), top1.getY(), rho1, r1_site);
                Vector3D r2_hyp = Transform.radec2geo(top2.getX(), top2.getY(), rho2, r2_site);
//            System.out.println(r1_hyp.scalarMultiply(1000).getNorm());
//            System.out.println(r2_hyp.scalarMultiply(1000).getNorm());
//            System.out.println(tof/3600);
                test.lambert_problem(r1_hyp.scalarMultiply(1000), r2_hyp.scalarMultiply(1000), tof, mu, cw, multi_revs);

                v1_mat = test.get_v1();

                Vector3D v1 = new Vector3D(v1_mat.getEntry(0, 0), v1_mat.getEntry(0, 1), v1_mat.getEntry(0, 2));
//            System.out.println(v1);
                PVCoordinates rv1 = new PVCoordinates(r1_hyp.scalarMultiply(1000), v1);
                Orbit orbit_out = new KeplerianOrbit(rv1, inertialFrame, initialDate, mu);
                A = orbit_out.getA();
                Ecc = orbit_out.getE();

//            System.out.println(ll + " - " +A);
                out_a.println(A);
                out_e.println(Ecc);
                out_rho1.println(rho1);
                out_rho2.println(rho2);
            }
            long time2 = System.nanoTime();
            long timeTaken = time2 - time1;

            out_a.close();
            out_e.close();
            out_rho1.close();
            out_rho2.close();

            System.out.println("Time taken " + timeTaken / 1000 / 1000 + " milli secs");

            // get the truth
            test.lambert_problem(r1.scalarMultiply(1000), r2.scalarMultiply(1000), tof, mu, cw, multi_revs);
            v1_mat = test.get_v1();
            Vector3D v1 = new Vector3D(v1_mat.getEntry(0, 0), v1_mat.getEntry(0, 1), v1_mat.getEntry(0, 2));
            PVCoordinates rv1 = new PVCoordinates(r1.scalarMultiply(1000), v1);
            Orbit orbit_out = new KeplerianOrbit(rv1, inertialFrame, initialDate, mu);
            //System.out.println(orbit_out.getA());
        } catch (FileNotFoundException ex) {
            Logger.getLogger(LambertMRev.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

}


