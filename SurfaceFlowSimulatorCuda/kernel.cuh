#pragma once

#include <thrust/device_vector.h>

#include "stdafx.h"
#include "ComBase.h"

struct colorTable;

template <typename T>
void copyToDevice(T* data, int size, T* gpuData);

void runKernel();

template <typename T>
void copyToHost(T* data, int size, T* hostData);

void copyFlowTrack(std::vector<polyline3D>& flowTracks, point3D*& flowTrack, int*& lineSize, int& size);

void copyRainIdx(int endTime, int rainIdxSize, int*& rainIdx);

void copyFlowVal(int endTime, int plineNum, float*& flowVal);

template<typename T>
void copyDataHtD(T* cpuData, int size, T*& gpuData);

extern "C" void runMarkRainIdx(point3D* flowTracks, int* flowTracksSize, int flowTracksSizeSize, double xmin, double ymax, double dx, double dy,
	float* pRVal, int imgDemWidth, int* rainIdx, float* flowVal, int plineNum, int endTime, int* rainIdxLen);

extern "C" void runDrawPoints(point3D* flowTracks, int* flowTracksSize, int flowTracksSizeSize, int* rainIdx, float* flowVal, int plineNum, int endTime,
	int* rainIdxLen, int imgWidth, int imgHeight, double dx, double dy, double xmin, double ymax, BYTE* outColor, colorTable* colorTables, int colorTableSize, float* pOutNum, int outNumSize);