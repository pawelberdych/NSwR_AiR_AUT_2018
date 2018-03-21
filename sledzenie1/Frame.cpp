#include "OpticalFlow.h"

Frame::Frame() :Film() 
{
	this->frame_path = image_path[0];
	this->setMatrix();
}


Frame::~Frame()
{

}

void Frame::changePath( int i)
{
	this->frame_path = image_path[i];
}

void Frame::setMatrix()
{
	this->frame_matrix = imread(frame_path, CV_LOAD_IMAGE_COLOR);
}

void Frame::setPolarMatrix(Mat matrix)
{
	this->polar_frame_matrix = matrix;

}

void Frame::setDimensions()
{
	this->frame_size_x = frame_matrix.cols;
	this->frame_size_y = frame_matrix.rows;
	this->frame_center_x = frame_matrix.cols / 2;
	this->frame_center_y = frame_matrix.rows / 2;
}