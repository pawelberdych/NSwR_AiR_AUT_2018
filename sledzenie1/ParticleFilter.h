#pragma once
#ifndef PARTICLEFILTER
#define PARTICLEFILTER

#include "libraries.h"
#include "OpticalFlow.h"
#include "Frame.h"


class ParticleFilter: public OpticalFlow
{
	int frame_width;
	int frame_height;
	float w1, w2, w3;
	int number_of_particles;
	Rect picked_found_object;
	Rect found_object_1;
	Rect found_object_2;
	Rect found_object_3;

	Mat pol;
	struct trackedObject 
	{
		Rect area;
		vector<movingPixel> moving_pixels;
		float avg_flow_r;
		float avg_flow_theta;

	};

	trackedObject target;

	struct particle
	{
		int x;
		int y;
		float weight;
		float avg;
	};

	vector<particle> random_particles;

	//vector<trackedObject> f;
	void fillTarget( Mat polar, vector<movingPixel> vec_all_moving_pixels);
	

public:

	float minimum_weight_1 = 200;
	float minimum_weight_2 = 200;
	float minimum_weight_3 = 200;

	ParticleFilter(int x, int y, Rect area,Mat polar, vector<movingPixel> vec_all_moving_pixels);
	~ParticleFilter();
	float calculateWeight(int x, int y, vector<movingPixel> vec_all_moving_pixels);
	void createParticles(Mat matrix,int f_height, vector<movingPixel> vec_all_moving_pixels);
	void createParticles2(Mat matrix, int f_height, vector<movingPixel> vec_all_moving_pixels);
	void createParticles3(Mat matrix, int f_height, vector<movingPixel> vec_all_moving_pixels);
	Rect getFoundObject();
	Mat showTracking(Mat current_frame);
	void clearTargetVector();
	//Mat ComputeDescriptors(const Mat &image, vector<KeyPoint> &keyPoints);
	//Rect mouseCallback(int event, int x, int y, int flags, void* userdata);
	//int getX() { return x; }
	//int getY() { return y; }
};

#endif
