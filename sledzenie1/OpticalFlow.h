#pragma once
#ifndef OPTICALFLOW
#define OPTICALFLOW

#include "libraries.h"
#include "Frame.h"


class OpticalFlow
{
	//____________________________________________________________________________________________________________________________________________________
	//---Variables - "Tracking unknown moving targets on omnidirectional vision"----------------------
	/*
	Mat Ir;		///gradient image in r-axis (uchar)
	Mat Itheta;	///gradient image in angle axis (uchar)
	Mat It;		///gradient image with respect to time (uchar)
	
	Mat Ir_int;	///gradient image in r-axis (int)
	Mat Itheta_int;	///gradient image in angle axis (int)
	Mat It_int;	///gradient image with respect to time (int)

	
	int alpha;	///hyperparameter chosen in article as alpha=500
	int magnitude_threshold;	///threshold above which change in optical flow is taken as movemnent
	*/
	//_____________________________________________________________________________________________________________________________________________________
	
	Mat mask;			//maska binarna umo¿liwiaj¹ca utworzenie okr¹g³ego obrazu 
	int R;				//promieñ obrazu omni w pikselach  
	int perimeter;		//obwód obrazu omni w pikselach (k¹t pe³ny)
	int frame_number;	//numer klatki
	
	Mat bsMask(int x, int y, int width, int height);
	
	
public:
	Frame* previous_frame_obj;
	Frame* current_frame_obj;


	OpticalFlow::OpticalFlow();
	~OpticalFlow();
	Mat cart_2_polar_matrix;	///3-channel matrix in which every element M_{ij} contains coordinates of pixel with r = i and theta = j on omni image (3rd channel not used)

	Mat transformOmniImage(Mat src_image);			///omni image transformed to coordinates r and theta
	Mat cart2PolarTransformation(Mat src_image);	///method for computation of cart_2_polar_matrix
	int getFrameNumber() { return frame_number; }
	void getPreviousAndCurrentFrame();

	int r_coordinate(int x, int y);
	int theta_coordinate(int x, int y);

	//____________________________________________________________________________________________________________________________________________________
	//---Variables & Functions - "Tracking unknown moving targets on omnidirectional vision"----------------------
	/*
	Mat saved_u;	///wektory przechowuj¹ce macierze u i v
	Mat saved_v;
	Mat u, u_previous;	///matrix containing values of optical flow in x axis
	Mat v, v_previous;	///matrix containing values of optical flow in y axis

	void preprocessing(Mat src_image);		///initial operations on image
	Mat calculate();	///method for calculation of optical flow
	*/
	//_____________________________________________________________________________________________________________________________________________________
	//---Variables & Functions - Farneback Optical Flow----------------------------------------------------------


	struct movingPixel
	{
		int x;
		int y;
		int r;
		int theta;
		double r_motion;
		double theta_motion;

	};

	vector<movingPixel> all_moving_pixels;


	Mat showFarnebackOpticalFlow();
	Mat r_flow, theta_flow;
	void calculateFarnebackOpticalFlow();
	void searchForMovingPixels();
	vector<movingPixel> getAllMovingPixels() { return all_moving_pixels; }
	
	void nextframe() { frame_number++; }

	//_____________________________________________________________________________________________________________________________________________________
	//---Unused Functions----------------------------------------------------------------------------------------
	//Mat calculateMagnitude();		///method for calculation of magnitude of optical flow

};



/*
Mat calculateIt(Mat frame_previous, Mat frame_current);
void performOpticalFlow(Mat frame_previous, Mat frame_current, Mat frame_next, OpticalFlow object_frame_previous, OpticalFlow object_frame_current);
Mat polar2CartTraformation(Mat cart_source, Mat polar_transformed);
Mat sobel1(Mat image);
*/
#endif