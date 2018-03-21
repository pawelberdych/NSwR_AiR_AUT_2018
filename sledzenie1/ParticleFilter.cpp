#include "ParticleFilter.h"

ParticleFilter::ParticleFilter(int x, int y, Rect area, Mat polar, vector<movingPixel> vec_all_moving_pixels) :OpticalFlow()
{
	frame_width = x;
	frame_height = y;

	target.area = area;
	fillTarget(polar, vec_all_moving_pixels);

	pol = polar;

	w1 = 0.5; //flow weight
	w2 = 0.5; //location weight
	w3 = 0.002; //number of pixels weight

	number_of_particles = 30;

	found_object_1.x = 0;
	found_object_1.y = 0;

	found_object_2.x = 0;
	found_object_2.y = 0;

	found_object_3.x = 0;
	found_object_3.y = 0;


	cout << "TARGET  x:" << target.area.x << "   y:" << target.area.y << endl;
}

ParticleFilter::~ParticleFilter() 
{

}

void ParticleFilter::fillTarget(Mat polar, vector<movingPixel> vec_all_moving_pixels)
{
	target.moving_pixels.push_back(movingPixel());
	int vector_index = 0;
	Mat a,b;
	polar.copyTo(a);
	previous_frame_obj->getMatrix().copyTo(b);

	int x_center, y_center;
	x_center = current_frame_obj->getMatrix().cols /2;
	y_center = current_frame_obj->getMatrix().rows / 2;


	for (int i = 0; i < vec_all_moving_pixels.size(); i++)
	{
		if (vec_all_moving_pixels[i].x >= target.area.x && vec_all_moving_pixels[i].x <= (target.area.x + target.area.width) && vec_all_moving_pixels[i].y >= target.area.y && vec_all_moving_pixels[i].y <= (target.area.y + target.area.height))
		{

			//Assign variables to vector
			target.moving_pixels[vector_index].r = vec_all_moving_pixels[i].r;
			target.moving_pixels[vector_index].theta = vec_all_moving_pixels[i].theta;
			target.moving_pixels[vector_index].r_motion = vec_all_moving_pixels[i].r_motion;
			target.moving_pixels[vector_index].theta_motion = vec_all_moving_pixels[i].theta_motion;

			//Add all values of flow
			target.avg_flow_r = target.avg_flow_r + target.moving_pixels[vector_index].r_motion;
			target.avg_flow_theta = target.avg_flow_theta + target.moving_pixels[vector_index].theta_motion;

			vector_index++;
			target.moving_pixels.push_back(movingPixel());
		}

	}

	//Calculate average flow
	//target.avg_flow_r = target.avg_flow_r / vector_index ;
	//target.avg_flow_theta = target.avg_flow_theta / vector_index ;

	target.avg_flow_r = target.avg_flow_r / vector_index / 2;
	target.avg_flow_theta = target.avg_flow_theta / vector_index / 2;

	/*
	for (int i = 0; i < target.moving_pixels.size(); i++)
	{
		circle(a, Point(target.moving_pixels[i].theta, target.moving_pixels[i].r), 1, Scalar(0, 255, 255), -1);
		int x = target.moving_pixels[i].r*cos(CV_PI*target.moving_pixels[i].theta / 180) + (previous_frame_obj->getFrameCenterX());
		int y = target.moving_pixels[i].r*sin(CV_PI*target.moving_pixels[i].theta / 180) + (previous_frame_obj->getFrameCenterY());
		circle(b, Point(x,y), 1, Scalar(0, 255, 255), -1);
	}
	
	imshow("polar", a);
	waitKey(0);
	imshow("cart", b);
	waitKey(0);
	*/
}

void ParticleFilter::createParticles(Mat matrix,int f_height, vector<movingPixel> vec_all_moving_pixels)
{
	Mat b;
	matrix.copyTo(b);

	float min = 0;
	srand(time(NULL));
	random_particles.push_back(particle());

	random_device rd;
	mt19937 e2(rd());

	normal_distribution<> x(target.area.x, 100);
	normal_distribution<> y(target.area.y, 100);

	for (int i = 0; i < number_of_particles; i++)
	{
		
		//cout << "\n Particle nr " << i << endl;
		//random_particles[i].x = rand() % (f_height - target.area.width) + 280;
		//random_particles[i].y = rand() % (f_height - target.area.height);
		random_particles[i].x = round(x(e2));
		random_particles[i].y = round(y(e2));

		random_particles[i].weight = calculateWeight(random_particles[i].x, random_particles[i].y, vec_all_moving_pixels);
		//---calculateWage---------------------------
		//cout << "x: " << random_particles[i].x << "  y: " << random_particles[i].y << endl;
		//rectangle(b, Point(random_particles[i].x, random_particles[i].y), Point(random_particles[i].x + target.area.width, random_particles[i].y + target.area.height), Scalar(255, 0, 255), 1, 8, 0);
		
		if (i == 0)
		{
			min = random_particles[i].weight;
		}
		
		if (random_particles[i].weight < min)
		{
			min = random_particles[i].weight;
			//cout << "ZMIANA MINIMUM" << endl;
			//cout << "min: " << min << endl;
		}
		
		random_particles.push_back(particle());
	}

	for (int i = 0; i < number_of_particles; i++)
	{
		if (random_particles[i].weight == min)
		{
			//rectangle(b, Point(random_particles[i].x , random_particles[i].y ), Point(random_particles[i].x + target.area.width , random_particles[i].y + target.area.height ), Scalar(0, 255, 0), 1);
			cout << "min wage = " << random_particles[i].weight << "  dla x = " << random_particles[i].x << "  y = " << random_particles[i].y << endl;
			found_object_1 = Rect(Point(random_particles[i].x, random_particles[i].y), Point(random_particles[i].x + target.area.width, random_particles[i].y + target.area.height));
			minimum_weight_1 = min;
		}
	}

	random_particles.clear();
	//imshow("1", b);
	//waitKey(0);
	//destroyWindow("1");
}

void ParticleFilter::createParticles2(Mat matrix,int f_height, vector<movingPixel> vec_all_moving_pixels)
{
	Mat b;
	matrix.copyTo(b);

	float min = 0;
	srand(time(NULL));
	random_particles.push_back(particle());

	random_device rd;
	mt19937 e2(rd());

	normal_distribution<> x(found_object_1.x, 60);
	normal_distribution<> y(found_object_1.y, 60);

	for (int i = 0; i < number_of_particles; i++)
	{
		//cout << "\n Particle nr " << i << endl;
		
		//random_particles[i].x = rand() % (f_height - target.area.width) + (found_object_1.x - (f_height / 2));
		//random_particles[i].y = rand() % (f_height - target.area.height) + (found_object_1.y - (f_height / 2));

		random_particles[i].x = round(x(e2));
		random_particles[i].y = round(y(e2));


		random_particles[i].weight = calculateWeight(random_particles[i].x, random_particles[i].y, vec_all_moving_pixels);
		//---calculateWage---------------------------
		//cout << "x: " << random_particles[i].x << "  y: " << random_particles[i].y << endl;
		//rectangle(b, Point(random_particles[i].x, random_particles[i].y), Point(random_particles[i].x + target.area.width, random_particles[i].y + target.area.height), Scalar(255, 0, 255), 1, 8, 0);

		if (i == 0)
		{
			min = random_particles[i].weight;
		}

		if (random_particles[i].weight < min)
		{
			min = random_particles[i].weight;
			//cout << "ZMIANA MINIMUM" << endl;
			//cout << "min: " << min << endl;
		}

		random_particles.push_back(particle());
	}

	for (int i = 0; i < number_of_particles; i++)
	{
		if (random_particles[i].weight == min)
		{
			//rectangle(b, Point(random_particles[i].x , random_particles[i].y ), Point(random_particles[i].x + target.area.width , random_particles[i].y + target.area.height ), Scalar(0, 255, 0), 1);
			cout << "min wage = " << random_particles[i].weight << "  dla x = " << random_particles[i].x << "  y = " << random_particles[i].y << endl;
			found_object_2 = Rect(Point(random_particles[i].x, random_particles[i].y), Point(random_particles[i].x + target.area.width, random_particles[i].y + target.area.height));
			minimum_weight_2 = min;
		}
	}

	random_particles.clear();
	//imshow("2", b);
	//waitKey(0);
	//destroyWindow("2");
}

void ParticleFilter::createParticles3(Mat matrix, int f_height, vector<movingPixel> vec_all_moving_pixels)
{
	Mat b;
	matrix.copyTo(b);

	float min = 0;
	srand(time(NULL));
	random_particles.push_back(particle());

	random_device rd;
	mt19937 e2(rd());

	normal_distribution<> x(found_object_2.x, 15);
	normal_distribution<> y(found_object_2.y, 15);

	for (int i = 0; i < number_of_particles; i++)
	{
		//cout << "\n Particle nr " << i << endl;
		
		//random_particles[i].x = rand() % (f_height - target.area.width) + (found_object_2.x - (f_height / 2));
		//random_particles[i].y = rand() % (f_height - target.area.height) + (found_object_2.y - (f_height / 2));

		random_particles[i].x = round(x(e2));
		random_particles[i].y = round(y(e2));

		random_particles[i].weight = calculateWeight(random_particles[i].x, random_particles[i].y, vec_all_moving_pixels);
		//---calculateWage---------------------------
		//cout << "x: " << random_particles[i].x << "  y: " << random_particles[i].y << endl;
		//rectangle(b, Point(random_particles[i].x, random_particles[i].y), Point(random_particles[i].x + target.area.width, random_particles[i].y + target.area.height), Scalar(255, 0, 255), 1, 8, 0);

		if (i == 0)
		{
			min = random_particles[i].weight;
		}

		if (random_particles[i].weight < min)
		{
			min = random_particles[i].weight;
			//cout << "ZMIANA MINIMUM" << endl;
			//cout << "min: " << min << endl;
		}

		random_particles.push_back(particle());
	}

	for (int i = 0; i < number_of_particles; i++)
	{
		if (random_particles[i].weight == min)
		{
		//	rectangle(b, Point(random_particles[i].x , random_particles[i].y ), Point(random_particles[i].x + target.area.width , random_particles[i].y + target.area.height ), Scalar(0, 255, 0), 1);
			cout << "min wage = " << random_particles[i].weight << "  dla x = " << random_particles[i].x << "  y = " << random_particles[i].y << endl;
			found_object_3 = Rect(Point(random_particles[i].x, random_particles[i].y), Point(random_particles[i].x + target.area.width, random_particles[i].y + target.area.height));
			minimum_weight_3 = min;
			
		}
	}
	
	random_particles.clear();
	//imshow("3", b);
	//waitKey(0);
	//destroyWindow("3");
}

float ParticleFilter::calculateWeight(int particle_x, int particle_y, vector<movingPixel> vec_all_moving_pixels)
{
	
	int x_center, y_center;
	x_center = frame_width / 2;
	y_center = frame_height / 2;
	float avg_flow_r = 0, avg_flow_theta = 0;
	float delta_x = 0, delta_y = 0, delta_flow_r = 0, delta_flow_theta = 0, delta_n = 0;
	int counter = 0;
	float total_weight = 1000;

	for (int i = 0; i < vec_all_moving_pixels.size(); i++)
	{
		if (vec_all_moving_pixels[i].x >= particle_x && vec_all_moving_pixels[i].x <= (particle_x + target.area.width) && vec_all_moving_pixels[i].y >= particle_y && vec_all_moving_pixels[i].y <= (particle_y + target.area.height))
		{
			avg_flow_r = avg_flow_r + vec_all_moving_pixels[i].r_motion;
			avg_flow_theta = avg_flow_theta + vec_all_moving_pixels[i].theta_motion;
			counter++;
		}
			
	}
	
	if (counter != 0)
	{
		
		

		avg_flow_r = avg_flow_r / counter;
		avg_flow_theta = avg_flow_theta / counter;

		//Delta flow r
		delta_flow_r = (abs(avg_flow_r - target.avg_flow_r));

		//Delta flow theta
		delta_flow_theta = (abs(avg_flow_theta - target.avg_flow_theta));

		//Delta x
		int r, theta;
		r = r_coordinate(target.area.x - x_center, target.area.y - y_center);
		theta = theta_coordinate(target.area.x - x_center, target.area.y - y_center);
	
		r = cvRound(r + (target.avg_flow_r));
		theta = cvRound(theta + (target.avg_flow_theta));

		int est_x, est_y;
		est_x = r*cos(CV_PI*theta/180) + (previous_frame_obj->getFrameCenterX());
		est_y = r*sin(CV_PI*theta / 180) + (previous_frame_obj->getFrameCenterY());

		
		//Mat b;
		//current_frame_obj->getMatrix().copyTo(b);
		//rectangle(b, target.area, Scalar(0, 255, 0)); //target
		//rectangle(b, Point(target.area.x, target.area.y), Point(target.area.x + target.area.width, target.area.y + target.area.height), Scalar(255, 0, 0));
		//rectangle(b, Point(particle_x, particle_y), Point(particle_x + target.area.width, particle_y + target.area.height), Scalar(255, 0, 255));
		//rectangle(b, Point(est_x, est_y), Point(est_x + target.area.width, est_y + target.area.height), Scalar(255, 0, 0));
		//imshow("f", b);
		//waitKey(0);
		
		

		delta_x = abs(particle_x - est_x);
		delta_y = abs(particle_y - est_y);


		total_weight = w1*(delta_flow_r + delta_flow_theta) + w2*(delta_x + delta_y);  // +w3*delta_n;

		
	}
	
	return total_weight;
}

Mat ParticleFilter::showTracking(Mat current_frame)
{
	
	Mat frame_with_object;
	current_frame.copyTo(frame_with_object);
	rectangle(frame_with_object, picked_found_object, Scalar(0, 0, 255), 1);
	rectangle(frame_with_object,Point( picked_found_object.x-1, picked_found_object.y-1),Point(picked_found_object.x+picked_found_object.width, picked_found_object.y+picked_found_object.height), Scalar(0, 0, 255), 1);
	//rectangle(frame_with_object, target.area, Scalar(255, 0, 255));
	//imshow("l", frame_with_object);
	//waitKey(0);

	return frame_with_object;
}

void ParticleFilter::clearTargetVector()
{
	target.moving_pixels.clear();
}

Rect ParticleFilter::getFoundObject()
{
	//cout << "min1: " << minimum_weight_1 << "  min2: " << minimum_weight_2 << "  min3: " << minimum_weight_3 << endl;
	if (minimum_weight_1 < minimum_weight_2 && minimum_weight_1 < minimum_weight_3)
		picked_found_object = found_object_1;

	if (minimum_weight_2 < minimum_weight_1 && minimum_weight_2 < minimum_weight_3)
		picked_found_object = found_object_2;

	if (minimum_weight_3 < minimum_weight_2 && minimum_weight_3 < minimum_weight_1)
		picked_found_object = found_object_3;

	return picked_found_object;
}