// SurfaceFlowSimulator.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "SurfaceFlowSimulator.h"

#include <iostream>
#include <time.h>
#include "include\include\gdal.h"
#include "ComBase.h"
#include "TIN2Flow.h"
#include <fstream>
#include <sstream>
#include <map>
#include <memory>
#include <thread>

#include <time.h>

#include "kernel.cuh"

#include "ImgCacher.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// 唯一的应用程序对象

CWinApp theApp;

//从颜色表中查找colDate的颜色值。
void GetOneColor(std::vector<colorTable>& colors, float colDate, BYTE& colorR, BYTE& colorG, BYTE& colorB)
{
	for (int i = 0; i < colors.size(); i++)
	{
		colorTable   aCol = colors[i];
		float starNum = aCol.starData;
		float endNum = aCol.endData;
		int colR = aCol.colorR;
		int colG = aCol.colorG;
		int colB = aCol.colorB;

		if (colDate >= starNum && colDate <= endNum)
		{
			colorR = (BYTE)colR;
			colorG = (BYTE)colG;
			colorB = (BYTE)colB;
			return;
		}
	}
}

//绘制一个雨滴单元点runPnt在BMP位图中的表现pOutColor，并修改对应栅格的计数pOutNum。(以雨滴实际的流量进行累加)
void   DrawOnePoint(point3D& runPnt, int imgWidth, int imgHeight, double dx, double dy, double Xmin, double Ymax, float RVal,
	BYTE* pOutColor, std::vector<colorTable>& colors, float* pOutNum)
{
	if (runPnt.z <= 0)
		return;

	int col = (runPnt.x - Xmin) / dx;
	int row = (Ymax - runPnt.y) / dy;

	if (row <= 2 || row >= imgHeight - 2) return;
	if (col <= 2 || col >= imgWidth - 2) return;

	pOutNum[row*imgWidth + col] = pOutNum[row*imgWidth + col] + RVal;
	row = imgHeight - row;

	BYTE colorR = 0, colorG = 0, colorB = 255;
	GetOneColor(colors, pOutNum[(imgHeight - row)*imgWidth + col], colorR, colorG, colorB);

	pOutColor[3 * ((row)*imgWidth + col) + 0] = colorG;   //R颜色  也许后面是GRB
	pOutColor[3 * ((row)*imgWidth + col) + 1] = colorR;   //G颜色
	pOutColor[3 * ((row)*imgWidth + col) + 2] = colorB;   //B颜色

}

//计算雨滴的颜色变化表
void SetColorTable(BSTR colorTableName, std::vector<colorTable>& colors)
{
	CString strcolTableName(colorTableName);
	if (strcolTableName.GetLength() < 1)
		return;

	CStdioFile infile;
	if (!infile.Open(strcolTableName, CFile::modeRead))
		return;

	//获得颜色数目
	CString strLine, out[5];
	infile.ReadString(strLine);

	int num = atol(strLine);
	for (int i = 0; i < num; i++)
	{
		infile.ReadString(strLine);
		if (strLine.GetLength() < 1)
		{
			infile.Close();
			return;
		}

		for (int j = 0; j < 5; j++)
		{
			AfxExtractSubString(out[j], strLine, j, ',');
			if (out[j].GetLength() < 1)
			{
				infile.Close();
				return;
			}
		}

		colorTable col;
		col.starData = atof(out[0]);
		col.endData = atof(out[1]);
		col.colorR = atoi(out[2]);
		col.colorG = atoi(out[3]);
		col.colorB = atoi(out[4]);

		colors.push_back(col);
	}

	infile.Close();
}

void writeImg(std::string fileName, std::shared_ptr<float> img, int width, int height, double* geoTransform)
{
	//std::string format = fileName.substr(fileName.find_last_of('.'), fileName.size());
	GDALDriverH hDriver = GDALGetDriverByName("HFA");
	//if (format == "img") {
	//	hDriver = GDALGetDriverByName("HFA");
	//} else if (format == "tif") {
	//	hDriver = GDALGetDriverByName("GTIFF");
	//} else if (format == "bmp") {
	//	hDriver = GDALGetDriverByName("BMP");
	//} else {
	//	return;
	//}

	if (hDriver == NULL || GDALGetMetadataItem(hDriver, GDAL_DCAP_CREATE, NULL) == NULL)
		return;

	char **papszOptions = NULL;
	GDALDataset *pDataset = (GDALDataset *)GDALCreate(hDriver, fileName.c_str(), width, height, 1, GDT_Float32, papszOptions);
	pDataset->SetGeoTransform(geoTransform);
	GDALRasterBand *pBand = pDataset->GetRasterBand(1);

	//pBand->SetNoDataValue((double)-9999.0);  //设置波段的无效值，这里设置影像=invalidVal为无效值。
	pBand->RasterIO(GF_Write, 0, 0, width, height, img.get(), width, height, GDT_Float32, 0, 0);
	double min = 0, max = 0, mean = 0, dev = 0; //设置的影像的最大值和最小值。
	pBand->ComputeStatistics(FALSE, &min, &max, &mean, &dev, NULL, NULL);

	GDALClose(pDataset);
}

//保存BMP影像数据。返回True说明保存成功，否则保存失败。
//ImgName —— 输入影像的路径；
//imgWidth,imgHeight —— 影像的宽、高；
//BmpHeader —— 影像的文件头信息；
//BmpInfo —— 影像的文件信息；
//pBuffer —— 影像的各个像素值。BYTE[3*imgWidth*imgHeight]
bool SaveBMP(CString strImgName, int imgWidth, int imgHeight, BITMAPFILEHEADER& BmpHeader, BITMAPINFO& BmpInfo, BYTE *pBuffer)
{
	CFile outfile;
	outfile.Open(strImgName, CFile::modeCreate | CFile::modeWrite);
	outfile.Write((LPSTR)&BmpHeader, sizeof(BmpHeader)); //首先读文件头部分到BmpHeader中。
	outfile.Write(&BmpInfo, sizeof(BITMAPINFO));
	outfile.Write(pBuffer, 3 * imgWidth*imgHeight);

	outfile.Flush();
	outfile.Close();

	return true;
}


//绘制一个时间段的已知产流数据的汇流过程模拟，在TIN的基础上绘制(已知区域的一个水深文件和流域出口点的坐标)
//inTINName—输入TIN 文件
//inDEMName—输入用于获取降雨源点的DEM文件
//inRunoffDirName— 输入每一时间段的产流数据*.img影像目录
//inBackBMPName—为指定的背景BMP图片，生成的FlowTrack是在该背景图片上绘制
//inColorTableName—指定颜色表的文本文件，根据雨滴的多少绘制不同的颜色
//inFirstDayDepthName一曼宁公式中用到的水力半径water flow depth
//inManingIndData一曼宁公式中用到曼宁系数，如N=0.04
//intelTimeData—产流数据的间隔时间，以秒为单位，如86400秒（一天时间）
//endTimeData—指定多长时间结束整个流域绘制成图像过程，以天为单位
//(outPntXData, outPntYData,outPntZData)—指定的流域出口点的坐标，该流域出口点将加入到构建的每一条流水线的最后一点
//outFlowBMPName—生成每一时刻具有颜色分类的背景图
//outFlowDischargeName—生成文件保存每个栅格在指定时刻的汇流量
int FlowDischarge(BSTR inTINName, BSTR inDEMName, BSTR inRunoffDirName, BSTR inBackBMPName, BSTR inColorTableName,
	BSTR inFirstDayDepthName, double inManingIndData, BSTR inTWIName, BSTR inLName, BSTR inCName, BSTR inSName, double intelTimeData, double endTimeData, double outPntXData,
	double outPntYData, double outPntZData, BSTR outFlowBMPName, BSTR outFlowDischargeName)
  {
	clock_t t1 = clock();

	CString strTINName(inTINName);
	CString strOutBMPName(outFlowBMPName);
	CString strOutDischarge(outFlowDischargeName);
	CString strFirstDayDepth(inFirstDayDepthName);
	if (strTINName.GetLength() < 1 || strOutBMPName.GetLength() < 1 || strOutDischarge.GetLength() < 1) {
		return 0;
	}

	//打开采样降雨源点的DEM影像
	CComBase comBase;
	int imgDEMWidth = 0, imgDEMHeight = 0;
	double dx = 0, dy = 0, Xmin = 0, Ymax = 0;
	if (!comBase.OpenImg(inDEMName, imgDEMWidth, imgDEMHeight, dx, dy, Xmin, Ymax)) {
		return 0;
	}
	double geoTransform[6]{Xmin, dx, 0, Ymax, 0, -dy};

	//打开TWI影像
	float *pTWI = new float[imgDEMWidth*imgDEMHeight];
	if (!comBase.OpenImg(inTWIName, imgDEMWidth, imgDEMHeight, dx, dy, Xmin, Ymax, pTWI))
	{
		return FALSE;
	}

	//打开L影像
	float *pL = new float[imgDEMWidth*imgDEMHeight];
	if (!comBase.OpenImg(inLName, imgDEMWidth, imgDEMHeight, dx, dy, Xmin, Ymax, pL))
		return FALSE;

	//打开C影像
	float *pC = new float[imgDEMWidth*imgDEMHeight];
	if (!comBase.OpenImg(inCName, imgDEMWidth, imgDEMHeight, dx, dy, Xmin, Ymax, pC))
		return FALSE;

	//打开S影像
	float *pS = new float[imgDEMWidth*imgDEMHeight];
	if (!comBase.OpenImg(inSName, imgDEMWidth, imgDEMHeight, dx, dy, Xmin, Ymax, pS))
		return FALSE;


	//读取TIN和降雨源点的DEM生成静态的FlowTrack(polylines)
	TIN2Flow tin;
	if (!tin.CreateFlow(strTINName, inDEMName)) {
		return 0;
	}
 	std::vector<polyline3D> polylines; //这里每一个起始点形成的流域线看作一个polyline3D,其中polyline3D的坐标串是在inteltime时间间隔时生成的坐标点。
 	tin.CalFlowPnts(intelTimeData, polylines, strFirstDayDepth, inManingIndData, pTWI, pL, pC, pS, outPntXData, outPntYData, outPntZData);

	//打开BMP背景影像
	BITMAPFILEHEADER BmpHeader; //记录"Bmp"格式图象的头信息
	BITMAPINFO BmpInfo;			//记录"Bmp"格式图象的信息
	int imgWidth = 0, imgHeight = 0;
	if (!comBase.OpenBMP(inBackBMPName, imgWidth, imgHeight, BmpHeader, BmpInfo)) {
		return 0;
	}
	BYTE* pInColor = new BYTE[3 * imgWidth*imgHeight];
	if (!comBase.OpenBMP(inBackBMPName, imgWidth, imgHeight, pInColor)) {
		delete[]pInColor;  pInColor = NULL;
		return 0;
	}

	//设置雨滴颜色表
	std::vector<colorTable>  colTables;
	SetColorTable(inColorTableName, colTables);

	//计算生成的BMP影像值
	const char* cEvenTime = "0000";
	CString eveTime(cEvenTime);

	//记录每一时刻，各个点有几滴雨rainIdx[i].size()，各个雨滴的位移量(rainIdx[i])[j]，即polylines[i]中的idx值
	int rainEndtime = comBase.OpenDirfile(inRunoffDirName);
	if (rainEndtime < 1) return 0;
	int plineNum = polylines.size();
	intArray*    rainIdx = new intArray[plineNum];    //rainIdx是一个二维变量，即rainIdx[plineNum][n]。其中(rainIdx[plineNum])表示每一条流水线上记录的多个雨滴串。(rainIdx[plineNum])[n]记录一条流水线上的第n个雨滴流到哪一个节点的下标，该下标从0开始(流水线中第一个节点开始)，每走一个时间间隔inteltime，下标数+1。
	floatArray*  pflowVal = new floatArray[plineNum]; //pflowVal是一个二维变量，即pflowVal[plineNum][n]，与rainIdx一一对应，记录各个雨滴的产流量信息。

	std::vector<std::thread> bmpThreads;
	std::vector<std::thread> imgThreads;

	point3D* gpuFlowTracks = nullptr;
	int* gpuFlowTracksSizes = nullptr;
	int gpuFlowTrackSize;
	copyFlowTrack(polylines, gpuFlowTracks, gpuFlowTracksSizes, gpuFlowTrackSize);

	int endTime = static_cast<int>(endTimeData);
	int* cpuRainIdx = new int[endTime * plineNum];
	int* gpuRainIdx = nullptr;
	cudaMalloc((void**)&gpuRainIdx, sizeof(int) * endTime * plineNum);
	
	float* cpuFlowVal = new float[endTime * plineNum];
	float* gpuFlowVal = nullptr;
	cudaMalloc((void**)&gpuFlowVal, sizeof(float) * endTime * plineNum);

	int* cpuRainIdxLen = new int[plineNum];
	int* gpuRainIdxLen = nullptr;
	cudaMalloc((void**)&gpuRainIdxLen, sizeof(int) * plineNum);
	cudaMemset(gpuRainIdxLen, 0, sizeof(int) * plineNum);

	BYTE* gpuOutColor = nullptr;
	cudaMalloc((void**)&gpuOutColor, sizeof(BYTE) * 3 * imgWidth * imgHeight);

	colorTable* gpuColorTables = nullptr;
	cudaMalloc((void**)&gpuColorTables, sizeof(colorTable) * colTables.size());
	cudaMemcpy(gpuColorTables, colTables.data(), sizeof(colorTable) * colTables.size(), cudaMemcpyHostToDevice);
	int gpuColorTablsSize = colTables.size();

	float* gpuOutNum = nullptr;
	cudaMalloc((void**)&gpuOutNum, sizeof(float) * imgWidth * imgHeight);
	int gpuOutNumSize = imgWidth * imgHeight;

	float* gpuRVal = nullptr;
	cudaMalloc((void**)&gpuRVal, sizeof(float) * imgDEMWidth * imgDEMHeight);

	ImgCacher cacher(imgDEMWidth, imgDEMHeight, (int)endTimeData, inRunoffDirName);

	for (int t = 0; t < endTimeData; t += 1) {
		//绘制背景
		cudaMemcpy(gpuOutColor, pInColor, 3 * imgWidth * imgHeight * sizeof(BYTE), cudaMemcpyHostToDevice);

		//将新产生的雨滴加入到rainidx[i]中
		if (t < rainEndtime) {
			//得到产流信息
			std::shared_ptr<float> prval = cacher.getImages(t);
			cudaMemcpy(gpuRVal, prval.get(), sizeof(float) * imgDEMWidth * imgDEMHeight, cudaMemcpyHostToDevice);

			runMarkRainIdx(gpuFlowTracks, gpuFlowTracksSizes, gpuFlowTrackSize, Xmin, Ymax, dx, dy,
				gpuRVal, imgDEMWidth, gpuRainIdx, gpuFlowVal, plineNum, endTime, gpuRainIdxLen);
		}

		//绘制雨滴
		cudaMemset(gpuOutNum, 0, sizeof(float) * imgWidth * imgHeight);

		runDrawPoints(gpuFlowTracks, gpuFlowTracksSizes, gpuFlowTrackSize, gpuRainIdx, gpuFlowVal, plineNum,
			endTime, gpuRainIdxLen, imgWidth, imgHeight, dx, dy, Xmin, Ymax,
			gpuOutColor, gpuColorTables, gpuColorTablsSize, gpuOutNum, gpuOutNumSize);

		BYTE* pOutColor = new BYTE[3 * imgWidth*imgHeight];
		cudaMemcpy(pOutColor, gpuOutColor, sizeof(BYTE) * 3 * imgWidth * imgHeight, cudaMemcpyDeviceToHost);

		float* pOutNum = new float[imgWidth*imgHeight];
		cudaMemcpy(pOutNum, gpuOutNum, imgWidth * imgHeight * sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy(cpuRainIdx, gpuRainIdx, sizeof(int) * plineNum * endTime, cudaMemcpyDeviceToHost);

		eveTime.Format("%.4d", t + 1);
		bmpThreads.emplace_back(std::thread(SaveBMP, strOutBMPName + eveTime + ".bmp", imgWidth, imgHeight, BmpHeader, BmpInfo, pOutColor));

		std::string rainFileName(strOutDischarge + eveTime + ".img");
		std::shared_ptr<float> poutnum(pOutNum);
		imgThreads.emplace_back(std::thread(writeImg, rainFileName, poutnum, imgWidth, imgHeight, geoTransform));
	}
	cacher.stop();
	for (auto& item : bmpThreads) {
		item.join();
	}
	for (auto& item : imgThreads) {
		item.join();
	}

	delete[]pInColor;   pInColor = NULL;
	delete[]rainIdx;    rainIdx = NULL;
	delete[]pflowVal;   pflowVal = NULL;

	return 1;
}



int _tmain(int argc, TCHAR* argv[], TCHAR* envp[])
{
	GDALAllRegister();

	clock_t begintime, endtime;
	double time;

	double inManingInd = 0.04;

	//模拟日流量
	double intelTime = 3600; 
	double endTime = 2;
	double outPntX = 596234.836;
	double outPntY = 5216225.614;
	double outPntZ = 151.591003;

	BSTR inTIN = L".\\TIN\\TIN.shp";
	BSTR inDEM = L".\\RainfallSourcePoints\\DEM30.img";
	BSTR inRunoffDir = L".\\R2001\\";
	BSTR inBackBMP = L".\\back.bmp";
	BSTR inColorTable = L".\\aColorTable.txt";
	BSTR inFirstDayDepth = L".\\FlowDepthFirstDay\\0120_IDW.img"; 
	BSTR inTWI = L".\\TWI1.img";
	BSTR inL = L".\\L1.img";
	BSTR inC = L".\\C1.img";
	BSTR inS = L".\\Slope.img";
	BSTR outFlowBMP = L".\\Flow30M\\";
	BSTR outFlowDischarge = L".\\Flow30N\\";
	std::ofstream timeTxt;
	timeTxt.open(".\\time.txt");

	begintime = clock();
	FlowDischarge(inTIN, inDEM, inRunoffDir, inBackBMP, inColorTable, inFirstDayDepth,
		inManingInd, inTWI, inL, inC, inS, intelTime, endTime, outPntX, outPntY, outPntZ, outFlowBMP, outFlowDischarge);
	endtime = clock();

	time = double(endtime - begintime) / CLOCKS_PER_SEC;

	timeTxt << "time:" << time << std::endl;
	timeTxt << std::flush;
	timeTxt.close();

	std::cout << "Simulating is finished!" << std::endl;
	std::cout << "The sum time is " << time << std::endl;
	system("pause");
	return 0;
}
