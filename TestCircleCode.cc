/* 
 * File:   main.cpp
 * Author: Sandip
 *
 * 13/10/13
 */
 
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <math.h>
 
#ifndef max
  #define max(a,b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
  #define min(a,b) (((a) < (b)) ? (a) : (b))
#endif
using namespace std;  

// Let EPS (epsilon) be a small value
const double EPS = 0.0000001f;
const double PI  =3.141592653589793238463f;
const float  PI_F=3.14159265358979f;
 

#define earthRadiusKm 6371.0
#define earthRadiusMeter  6371000
// This function converts decimal degrees to radians
double deg2rad(double deg) {
  return (deg * M_PI / 180);
}

//  This function converts radians to decimal degrees
double rad2deg(double rad) {
  return (rad * 180 / M_PI);
}

/**
 * Returns the distance between two points on the Earth.
 * Direct translation from http://en.wikipedia.org/wiki/Haversine_formula
 * @param lat1d Latitude of the first point in degrees
 * @param lon1d Longitude of the first point in degrees
 * @param lat2d Latitude of the second point in degrees
 * @param lon2d Longitude of the second point in degrees
 * @return The distance between the two points in kilometers
 */
double distanceEarth(double lat1d, double lon1d, double lat2d, double lon2d) {
  double lat1r, lon1r, lat2r, lon2r, u, v;
  lat1r = deg2rad(lat1d);
  lon1r = deg2rad(lon1d);
  lat2r = deg2rad(lat2d);
  lon2r = deg2rad(lon2d);
  u = sin((lat2r - lat1r)/2);
  v = sin((lon2r - lon1r)/2);
  return 2.0 * earthRadiusKm * asin(sqrt(u * u + cos(lat1r) * cos(lat2r) * v * v));
}


double distanceEuclideanCalculate(double x1, double y1, double x2, double y2)
{
    double x = x1 - x2;
    double y = y1 - y2;
    double dist;

    dist = pow(x,2)+pow(y,2);           //calculating distance by euclidean formula
    dist = sqrt(dist);                  //sqrt is function in math.h

    return dist;
}

 //  public domain function by Darel Rex Finley, 2006



//  Determines the intersection point of the line segment defined by points A and B
//  with the line segment defined by points C and D.
//
//  Returns YES if the intersection point was found, and stores that point in X,Y.
//  Returns NO if there is no determinable intersection point, in which case X,Y will
//  be unmodified.

bool lineSegmentIntersection(
    double Ax, double Ay,
    double Bx, double By,
    double Cx, double Cy,
    double Dx, double Dy,
    double *X, double *Y) {
    
      double  distAB, theCos, theSin, newX, ABpos ;
    
      //  Fail if either line segment is zero-length.
      if (Ax==Bx && Ay==By || Cx==Dx && Cy==Dy) return false;
    
      //  Fail if the segments share an end-point.
      if (Ax==Cx && Ay==Cy || Bx==Cx && By==Cy
      ||  Ax==Dx && Ay==Dy || Bx==Dx && By==Dy) {
        return false; }
    
      //  (1) Translate the system so that point A is on the origin.
      Bx-=Ax; By-=Ay;
      Cx-=Ax; Cy-=Ay;
      Dx-=Ax; Dy-=Ay;
    
      //  Discover the length of segment A-B.
      distAB=sqrt(Bx*Bx+By*By);
    
      //  (2) Rotate the system so that point B is on the positive X axis.
      theCos=Bx/distAB;
      theSin=By/distAB;
      newX=Cx*theCos+Cy*theSin;
      Cy  =Cy*theCos-Cx*theSin; Cx=newX;
      newX=Dx*theCos+Dy*theSin;
      Dy  =Dy*theCos-Dx*theSin; Dx=newX;
    
      //  Fail if segment C-D doesn't cross line A-B.
      if (Cy<0. && Dy<0. || Cy>=0. && Dy>=0.) return false;
    
      //  (3) Discover the position of the intersection point along line A-B.
      ABpos=Dx+(Cx-Dx)*Dy/(Dy-Cy);
    
      //  Fail if segment C-D crosses line A-B outside of segment A-B.
      if (ABpos<0. || ABpos>distAB) return false;
    
      //  (4) Apply the discovered position to line A-B in the original coordinate system.
      *X=Ax+ABpos*theCos;
      *Y=Ay+ABpos*theSin;
    
      //  Success.
      return true; }
    


/*
 * Find the intersection point(s) of two circles,
 * when their centers and radiuses are given (2D).
 */
 
class Point2d{
public:
    Point2d() {}
    Point2d(double x, double y)
        : X(x), Y(y) {}
     
    double x() const { return X; }
    double y() const { return Y; }
     
    /**
     * Returns the norm of this vector.
     * @return the norm
    */
    double norm() const {
        return sqrt( X * X + Y * Y );
    }
     
    void setCoords(double x, double y) {
        X = x; Y = y;
    }
     
    // Print point
    friend std::ostream& operator << ( std::ostream& s, const Point2d& p )  {
      s << p.x() << " " << p.y();
      return s;
    }
public:
    double X;
    double Y;
};
 
class Circle{
public:
    /**
     * @param R - radius
     * @param C - center
     */
    Circle(double R, Point2d& C) 
        : r(R), c(C) {}
         
    /**
     * @param R - radius
     * @param X - center's x coordinate
     * @param Y - center's y coordinate
     */
    Circle(double R, double X, double Y) 
        : r(R), c(X, Y) {}    
     
    Point2d getC() const { return c; }
    double getR() const { return r; }
     
    size_t intersect(const Circle& C2, Point2d& i1, Point2d& i2) {
        // distance between the centers
        double d = Point2d(c.x() - C2.c.x(), 
                c.y() - C2.c.y()).norm();
         
        // find number of solutions
        if(d > r + C2.r) // circles are too far apart, no solution(s)
        {
            std::cout << "Circles are too far apart\n";
            return 0;
        }
        else if(d == 0 && r == C2.r) // circles coincide
        {
            std::cout << "Circles coincide\n";
            return 0;
        }
        // one circle contains the other
        else if(d + min(r, C2.r) < max(r, C2.r))
        {
            std::cout << "One circle contains the other\n";
            return 0;
        }
        else
        {
            double a = (r*r - C2.r*C2.r + d*d)/ (2.0*d);
            double h = sqrt(r*r - a*a);
             
            // find p2
            Point2d p2( c.x() + (a * (C2.c.x() - c.x())) / d,
                    c.y() + (a * (C2.c.y() - c.y())) / d);
             
            // find intersection points p3
            i1.setCoords( p2.x() + (h * (C2.c.y() - c.y())/ d),
                    p2.y() - (h * (C2.c.x() - c.x())/ d)
            );
            i2.setCoords( p2.x() - (h * (C2.c.y() - c.y())/ d),
                    p2.y() + (h * (C2.c.x() - c.x())/ d)
            );
             
            if(d == r + C2.r)
                return 1;
            return 2;
        }
    }
     
    // Print circle
    friend std::ostream& operator << ( std::ostream& s, const Circle& C )  {
      s << "Center: " << C.getC() << ", r = " << C.getR();
      return s;
    }
private:
    // radius
    double r;
    // center
    Point2d c;
     
};
int main(void)
{
    // radius and center of circles
    // Circle c1(2, 0, 0);
    // Circle c2(1, 0, 2);
    // Point2d i1, i2;
    
    Point2d i1(-9999, -9999);
    Point2d i2(-9999, -9999);

    //point intersect	 
    Point2d l1(37.396346133189255, -122.56210327148439);
    Point2d l2(37.39307301476918, -121.77520751953125);
    Point2d pt1(37.48684571271661, -122.34649658203126);
    //point failed intersect
    // Point2d l1(37.396346133189255, -122.56210327148439);
    // Point2d l2(37.39307301476918, -121.77520751953125);
    // Point2d pt1(37.497741887143576,-122.67059326171876);

    // double c1Radius =  distanceEarth (l1.X,l1.Y, pt1.X,pt1.Y);
    // double c2Radius =  distanceEarth (l2.X, l2.Y, pt1.X,pt1.Y);
    double c1Radius =  distanceEuclideanCalculate (l1.X,l1.Y, pt1.X,pt1.Y);
    double c2Radius =  distanceEuclideanCalculate (l2.X, l2.Y, pt1.X,pt1.Y);
    cout<<"c1Radius : "<<c1Radius <<" c2Radius : " <<c2Radius <<endl;

    Circle c1(l1.X,l1.Y, c1Radius);
    Circle c2(l2.X,l2.Y, c2Radius);


    std::cout << c1 << "\n" << c2 << "\n";
    // intersections point(s)
    size_t i_points = c1.intersect(c2, i1, i2);

     
    std::cout << "Circle Intersection point(s)\n";
    if(i_points == 2)
    {
        std::cout << i1 << "\n";
        std::cout << i2 << "\n";
    }
    else if(i_points)
        std::cout << i1 << "\n";
     
    double LnX = -9999.00, LnY = -9999.00;
    bool ret =  lineSegmentIntersection(l1.X,l1.Y,l2.X, l2.Y, \
        i1.X, i1.Y, i2.X, i2.Y,&LnX,&LnY);

    std::cout << "Line Intersection point(s) " <<ret <<"  " << LnX <<"  :  "<<LnY<<endl;
        
    return 0;
}
