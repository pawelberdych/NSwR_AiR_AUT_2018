#pragma once
#ifndef FRAME
#define FRAME

#include "libraries.h"
#include "Film.h"

class Frame : public Film
{
private:
	string frame_path;		//�cie�ka do danej klatki
	Mat frame_matrix;		//macierz kartezja�ska [x,y] obrazu klatki
	Mat polar_frame_matrix;	//macierz biegunowa [360,R] klatki
	int frame_size_x;		//szeroko�� klatki
	int frame_size_y;		//wysoko�� klatki
	int frame_center_x;		//centrum obrazu - x
	int frame_center_y;		//centrum obrazu - y
	

public:
	
	Frame();
	~Frame();

	void changePath(int i);

	void setMatrix();
	Mat getMatrix() { return frame_matrix; }			//zwraca macierz kartezja�sk�

	void setPolarMatrix(Mat matrix);
	Mat getPolarMatrix() { return polar_frame_matrix; }	//zwraca macierz biegunow� 
	
	void setDimensions();
	int getFrameSizeX() { return frame_size_x; }		//zwraca szeroko�� klatki
	int getFrameSizeY() { return frame_size_y; }		//zwraca wysoko�c klatki
	int getFrameCenterX() { return frame_center_x; }	//zwraca �rodek klatki w osi X
	int getFrameCenterY() { return frame_center_y; }	//zwraca �rodek klatki w osi Y
};

#endif
