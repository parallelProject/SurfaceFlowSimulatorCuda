
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <device_atomic_functions.h>

#include "kernel.cuh"

#include <vector>
#include "SurfaceFlowSimulator.h"


template <typename T>
void copyToDevice(T* data, int size, T* gpuData)
{
	if (gpuData == nullptr) {
		cudaMalloc(&gpuData, sizeof(T) * size);
	}
	cudaMemcpy(gpuData, data, sizeof(T) * size, cudaMemcpyHostToDevice);
}

void copyFlowTrack(std::vector<polyline3D>& flowTracks, point3D*& flowTrack, int*& lineSize, int& size)
{
	int sum = 0;
	int* cpuLineSize = new int[flowTracks.size()];
	std::vector<point3D> cpuLines;
	for (size_t i = 0; i < flowTracks.size(); ++i) {
		sum += flowTracks[i].size();
		cpuLineSize[i] = flowTracks[i].size();
		if (i > 0) {
			cpuLineSize[i] += cpuLineSize[i - 1];
		}
		for (size_t j = 0; j < flowTracks[i].size(); ++j) {
			cpuLines.emplace_back(flowTracks[i][j]);
		}
	}

	cudaMalloc(&flowTrack, sizeof(point3D) * sum);
	cudaMemcpy(flowTrack, cpuLines.data(), sizeof(point3D) * sum, cudaMemcpyHostToDevice);
	cudaMalloc(&lineSize, sizeof(int) * flowTracks.size());
	cudaMemcpy(lineSize, cpuLineSize, sizeof(int) * flowTracks.size(), cudaMemcpyHostToDevice);
	size = flowTracks.size();
}

void copyRainIdx(int endTime, int plineNum, int*& rainIdx)
{
	cudaMalloc(&rainIdx, sizeof(int) * endTime * plineNum);
	cudaMemset(rainIdx, 0, sizeof(int) * endTime * plineNum);
}

void copyFlowVal(int endTime, int plineNum, float*& flowVal)
{
	cudaMalloc(&flowVal, sizeof(float) * endTime * plineNum);
	cudaMemset(flowVal, 0, sizeof(float) * endTime * plineNum);
}

//void runKernel();

template <typename T>
void copyToHost(T* data, int size, T* hostData)
{
	if (hostData == nullptr) {
		hostData = new T[size];
	}
	cudaMemcpy(hostData, data, sizeof(T) * size, cudaMemcpyDeviceToHost);
}

template<typename T>
void copyDataHtD(T* cpuData, int size, T*& gpuData)
{
	if (cpuData == nullptr) {
		return;
	}
	cudaMalloc(&gpuData, sizeof(T) * size);
	cudaMemcpy(gpuData, cpuData, sizeof(T) * size);
}

__global__ void markRainIdx(point3D* flowTracks, int* flowTracksSize, int flowTracksSizeSize, double xmin, double ymax, double dx, double dy,
	float* pRVal, int imgDemWidth, int* rainIdx, float* flowVal, int plineNum, int endTime, int* rainIdxLen)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int idy = threadIdx.y;
	int id = idy * (gridDim.x * blockDim.x) + idx;
	if (id >= plineNum)
		return;

	int plineStart = 0;
	if (id > 0)
		plineStart = flowTracksSize[id - 1];
	point3D runPnt = flowTracks[plineStart];
	
	int col = (runPnt.x - xmin) / dx;
	int row = (ymax - runPnt.y) / dy;

	if (col < 0) col = 0;
	if (row < 0) row = 0;

	float rValue = pRVal[row * imgDemWidth + col];
	if (rValue > 0) {
		int rainIdxEnd = id * endTime + rainIdxLen[id];
		rainIdx[rainIdxEnd] = 0;
		flowVal[rainIdxEnd] = rValue;
		rainIdxLen[id]++;
	}
}

void runMarkRainIdx(point3D* flowTracks, int* flowTracksSize, int flowTracksSizeSize, double xmin, double ymax, double dx, double dy,
	float* pRVal, int imgDemWidth, int* rainIdx, float* flowVal, int plineNum, int endTime, int* rainIdxLen)
{
	dim3 tpb(32, 32);
	int blockNum = (plineNum + 32 * 32 - 1) / (32 * 32);
	markRainIdx << <blockNum, tpb >> >(flowTracks, flowTracksSize, flowTracksSizeSize, xmin, ymax, dx, dy,
		pRVal, imgDemWidth, rainIdx, flowVal, plineNum, endTime, rainIdxLen);
}

__device__ void DeviceGetOneColor(colorTable* colors, int colorTableSize, float colDate, BYTE& colorR, BYTE& colorG, BYTE& colorB)
{
	for (int i = 0; i < colorTableSize; i++) {
		colorTable   aCol = colors[i];
		float starNum = aCol.starData;
		float endNum = aCol.endData;
		int colR = aCol.colorR;
		int colG = aCol.colorG;
		int colB = aCol.colorB;

		if (colDate >= starNum && colDate <= endNum) {
			colorR = (BYTE)colR;
			colorG = (BYTE)colG;
			colorB = (BYTE)colB;
			return;
		}
	}
}

__device__  inline void atomicFloatAdd(float *address, float val)
{
	int i_val = __float_as_int(val);
	int tmp0 = 0;
	int tmp1;

	while ((tmp1 = atomicCAS((int *)address, tmp0, i_val)) != tmp0) {
		tmp0 = tmp1;
		i_val = __float_as_int(val + __int_as_float(tmp1));
	}
}

//绘制一个雨滴单元点runPnt在BMP位图中的表现pOutColor，并修改对应栅格的计数pOutNum。(以雨滴实际的流量进行累加)
__device__ void DeviceDrawOnePoint(point3D& runPnt, int imgWidth, int imgHeight, 
	double dx, double dy, double Xmin, double Ymax, float RVal,
	BYTE* pOutColor, colorTable* colors, int colorTableSize, float* pOutNum)
{
	if (runPnt.z <= 0)
		return;

	int col = (runPnt.x - Xmin) / dx;
	int row = (Ymax - runPnt.y) / dy;

	if (row <= 2 || row >= imgHeight - 2) return;
	if (col <= 2 || col >= imgWidth - 2) return;

	//pOutNum[row*imgWidth + col] = pOutNum[row*imgWidth + col] + RVal;
	atomicFloatAdd(&pOutNum[row*imgWidth + col], RVal);
	row = imgHeight - row;

	BYTE colorR = 0, colorG = 0, colorB = 255;
	DeviceGetOneColor(colors, colorTableSize, pOutNum[(imgHeight - row)*imgWidth + col], colorR, colorG, colorB);

	pOutColor[3 * ((row)*imgWidth + col) + 0] = colorG;   //R颜色  也许后面是GRB
	pOutColor[3 * ((row)*imgWidth + col) + 1] = colorR;   //G颜色
	pOutColor[3 * ((row)*imgWidth + col) + 2] = colorB;   //B颜色
}

__global__ void drawPoints(point3D* flowTracks, int* flowTracksSize, int flowTracksSizeSize, 
	int* rainIdx, float* flowVal, int plineNum, int endTime,
	int* rainIdxLen, int imgWidth, int imgHeight, double dx, double dy, double xmin, double ymax, 
	BYTE* outColor, colorTable* colorTables, int colorTableSize, float* pOutNum, int outNumSize)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int idy = threadIdx.y;
	int id = idy * (gridDim.x * blockDim.x) + idx;
	if (id >= plineNum)
		return;

	int plineStart = 0;
	int plineSize = flowTracksSize[id];
	if (id > 0) {
		plineStart = flowTracksSize[id - 1];
		plineSize = plineSize - flowTracksSize[id - 1];
	}

	if (rainIdxLen[id] < 1 || plineSize < 1)
		return;

	for (int n = 0; n < rainIdxLen[id]; n++) {
		if (rainIdx[id * endTime + n] >= plineSize || rainIdx[id * endTime + n] < 0)
			continue;
		point3D runPnt = flowTracks[plineStart + rainIdx[id * endTime + n]];//flowTracks[rainIdx[id * endTime + n]];
		DeviceDrawOnePoint(runPnt, imgWidth, imgHeight, dx, dy, xmin, ymax, 
			flowVal[id * endTime + n], outColor, colorTables, colorTableSize, pOutNum);
		rainIdx[id * endTime + n] = rainIdx[id * endTime + n] + 1;
	}
}

void runDrawPoints(point3D* flowTracks, int* flowTracksSize, int flowTracksSizeSize, int* rainIdx, float* flowVal, int plineNum, int endTime,
	int* rainIdxLen, int imgWidth, int imgHeight, double dx, double dy, double xmin, double ymax, BYTE* outColor, colorTable* colorTables, int colorTableSize, float* pOutNum, int outNumSize)
{
	dim3 tpb(32, 32);
	int blockNum = (plineNum + 32 * 32 - 1) / (32 * 32);
	drawPoints << <blockNum, tpb >> >(flowTracks, flowTracksSize, flowTracksSizeSize, rainIdx, flowVal, plineNum, endTime,
		rainIdxLen, imgWidth, imgHeight, dx, dy, xmin, ymax, outColor, colorTables, colorTableSize, pOutNum, outNumSize);
}