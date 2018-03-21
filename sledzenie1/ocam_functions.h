/*------------------------------------------------------------------------------
   Example code that shows the use of the 'cam2world" and 'world2cam" functions
   Shows also how to undistort images into perspective or panoramic images
   Copyright (C) 2008 DAVIDE SCARAMUZZA, ETH Zurich
   Author: Davide Scaramuzza - email: davide.scaramuzza@ieee.org
------------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include "opencv2/video/tracking.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include <iostream>
#include <ctype.h>

using namespace cv;
using namespace std;


#define CMV_MAX_BUF 1024
#define MAX_POL_LENGTH 64

struct ocam_model
{
  double pol[MAX_POL_LENGTH];    // the polynomial coefficients: pol[0] + x"pol[1] + x^2*pol[2] + ... + x^(N-1)*pol[N-1]
  int length_pol;                // length of polynomial
  double invpol[MAX_POL_LENGTH]; // the coefficients of the inverse polynomial
  int length_invpol;             // length of inverse polynomial
  double xc;         // row coordinate of the center
  double yc;         // column coordinate of the center
  double c;          // affine parameter
  double d;          // affine parameter
  double e;          // affine parameter
  int width;         // image width
  int height;        // image height
};


/*------------------------------------------------------------------------------
 This function reads the parameters of the omnidirectional camera model from
 a given TXT file
------------------------------------------------------------------------------*/
int get_ocam_model(struct ocam_model *myocam_model, char *filename);

/*------------------------------------------------------------------------------
WORLD2CAM projects a 3D point on to the image
    WORLD2CAM(POINT2D, POINT3D, OCAM_MODEL)
    projects a 3D point (point3D) on to the image and returns the pixel coordinates (point2D).

    POINT3D = [X;Y;Z] are the coordinates of the 3D point.
    OCAM_MODEL is the model of the calibrated camera.
    POINT2D = [rows;cols] are the pixel coordinates of the reprojected point

    Copyright (C) 2009 DAVIDE SCARAMUZZA
    Author: Davide Scaramuzza - email: davide.scaramuzza@ieee.org

    NOTE: the coordinates of "point2D" and "center" are already according to the C
    convention, that is, start from 0 instead than from 1.
------------------------------------------------------------------------------*/
void world2cam(double point2D[2], double point3D[3], struct ocam_model *myocam_model);

/*------------------------------------------------------------------------------
 CAM2WORLD projects a 2D point onto the unit sphere
    CAM2WORLD(POINT3D, POINT2D, OCAM_MODEL)
    back-projects a 2D point (point2D), in pixels coordinates,
    onto the unit sphere returns the normalized coordinates point3D = [x;y;z]
    where (x^2 + y^2 + z^2) = 1.

    POINT3D = [X;Y;Z] are the coordinates of the 3D points, such that (x^2 + y^2 + z^2) = 1.
    OCAM_MODEL is the model of the calibrated camera.
    POINT2D = [rows;cols] are the pixel coordinates of the point in pixels

    Copyright (C) 2009 DAVIDE SCARAMUZZA
    Author: Davide Scaramuzza - email: davide.scaramuzza@ieee.org

    NOTE: the coordinates of "point2D" and "center" are already according to the C
    convention, that is, start from 0 instead than from 1.
------------------------------------------------------------------------------*/
void cam2world(double point3D[3], double point2D[2], struct ocam_model *myocam_model);


/*------------------------------------------------------------------------------
 Create Look Up Table for undistorting the image into a panoramic image
 It computes a trasformation from cartesian to polar coordinates
 Therefore it does not need the calibration parameters
 The region to undistorted in contained between Rmin and Rmax
 xc, yc are the row and column coordinates of the image center
------------------------------------------------------------------------------*/
void create_panoramic_undistortion_LUT ( Mat mapx, Mat mapy, float Rmin, float Rmax, float xc, float yc );


// Funkcja Ocam dodane
void Ovu2point2d(double point2d[2], double pointVU[2], struct ocam_model *myocam_model); // przeksztalcenie z VU do 2D
void Opoint2d2vu(double pointVU[2], double point2d[2], struct ocam_model *myocam_model); // przeksztalcenie z 2D do VU
void O3d2vu(double pointVU[2], double point3d[3], struct ocam_model *myocam_model); // przeksztalcenie z 3D do VU
void Ovu2d3(double point3d[3], double pointVU[2], struct ocam_model *myocam_model); // przeksztalcenie z VU do 3D
void Ocreate_panoramic_central(Mat mapx, Mat mapy, float Rprocent, float Rmax_new, ocam_model *myocam_model); // tworzy panorame z centraln� lini� poziom� na �rodku panoramy
void Ocreate_panoramic(Mat mapx, Mat mapy, float Rmin, float Rmax, struct ocam_model *myocam_model); // tworzy panorame z parametrami Rmin,Rmax
void Ocreate_projection(Mat mapx, Mat mapy, double f, double fi, double theta, struct ocam_model *myocam_model); // tworzy projection plane
																												 // 3 uklady osi odniesienia (3D, VU, 2D)