// SurfaceFlowSimulator.cpp : �������̨Ӧ�ó������ڵ㡣
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


// Ψһ��Ӧ�ó������

CWinApp theApp;

//����ɫ���в���colDate����ɫֵ��
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

//����һ����ε�Ԫ��runPnt��BMPλͼ�еı���pOutColor�����޸Ķ�Ӧդ��ļ���pOutNum��(�����ʵ�ʵ����������ۼ�)
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

	pOutColor[3 * ((row)*imgWidth + col) + 0] = colorG;   //R��ɫ  Ҳ�������GRB
	pOutColor[3 * ((row)*imgWidth + col) + 1] = colorR;   //G��ɫ
	pOutColor[3 * ((row)*imgWidth + col) + 2] = colorB;   //B��ɫ

}

//������ε���ɫ�仯��
void SetColorTable(BSTR colorTableName, std::vector<colorTable>& colors)
{
	CString strcolTableName(colorTableName);
	if (strcolTableName.GetLength() < 1)
		return;

	CStdioFile infile;
	if (!infile.Open(strcolTableName, CFile::modeRead))
		return;

	//�����ɫ��Ŀ
	CString strLine, out[5];
	infile.ReadString(strLine);

	const char* cStrLine = CStringA(strLine);
	int num = atol(cStrLine);
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
		
		const char* cOut0 = CStringA(out[0]);
		col.starData = atof(cOut0);
		const char* cOut1 = CStringA(out[1]);
		col.endData = atof(cOut1);
		const char* cOut2 = CStringA(out[2]);
		col.colorR = atoi(cOut2);
		const char* cOut3 = CStringA(out[3]);
		col.colorG = atoi(cOut3);
		const char* cOut4 = CStringA(out[4]);
		col.colorB = atoi(cOut4);

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

	//pBand->SetNoDataValue((double)-9999.0);  //���ò��ε���Чֵ����������Ӱ��=invalidValΪ��Чֵ��
	pBand->RasterIO(GF_Write, 0, 0, width, height, img.get(), width, height, GDT_Float32, 0, 0);
	double min = 0, max = 0, mean = 0, dev = 0; //���õ�Ӱ������ֵ����Сֵ��
	pBand->ComputeStatistics(FALSE, &min, &max, &mean, &dev, NULL, NULL);

	GDALClose(pDataset);
}

//����BMPӰ�����ݡ�����True˵������ɹ������򱣴�ʧ�ܡ�
//ImgName ���� ����Ӱ���·����
//imgWidth,imgHeight ���� Ӱ��Ŀ��ߣ�
//BmpHeader ���� Ӱ����ļ�ͷ��Ϣ��
//BmpInfo ���� Ӱ����ļ���Ϣ��
//pBuffer ���� Ӱ��ĸ�������ֵ��BYTE[3*imgWidth*imgHeight]
bool SaveBMP(CString strImgName, int imgWidth, int imgHeight, BITMAPFILEHEADER& BmpHeader, BITMAPINFO& BmpInfo, BYTE *pBuffer)
{
	CFile outfile;
	outfile.Open(strImgName, CFile::modeCreate | CFile::modeWrite);
	outfile.Write((LPSTR)&BmpHeader, sizeof(BmpHeader)); //���ȶ��ļ�ͷ���ֵ�BmpHeader�С�
	outfile.Write(&BmpInfo, sizeof(BITMAPINFO));
	outfile.Write(pBuffer, 3 * imgWidth*imgHeight);

	outfile.Flush();
	outfile.Close();

	return true;
}

//����һ��ʱ��ε���֪�������ݵĻ�������ģ�⣬��TIN�Ļ����ϻ���(��֪�����һ��ˮ���ļ���������ڵ������)
//inTINName������TIN �ļ�
//inDEMName���������ڻ�ȡ����Դ���DEM�ļ�
//inRunoffDirName�� ����ÿһʱ��εĲ�������*.imgӰ��Ŀ¼
//inBackBMPName��Ϊָ���ı���BMPͼƬ�����ɵ�FlowTrack���ڸñ���ͼƬ�ϻ���
//inColorTableName��ָ����ɫ����ı��ļ���������εĶ��ٻ��Ʋ�ͬ����ɫ
//inFirstDayDepthNameһ������ʽ���õ���ˮ���뾶water flow depth
//inManingIndDataһ������ʽ���õ�����ϵ������N=0.04
//inLName-�³��Ĵ洢·��
//inCName-�������ʵĴ洢·��
//inSName-�¶ȵĴ洢·��
//intelTimeData���������ݵļ��ʱ�䣬����Ϊ��λ����86400�루һ��ʱ�䣩
//endTimeData��ָ���೤ʱ���������������Ƴ�ͼ����̣�����Ϊ��λ
//(outPntXData, outPntYData,outPntZData)��ָ����������ڵ�����꣬��������ڵ㽫���뵽������ÿһ����ˮ�ߵ����һ��
//outFlowBMPName������ÿһʱ�̾�����ɫ����ı���ͼ
//outFlowDischargeName�������ļ�����ÿ��դ����ָ��ʱ�̵Ļ�����
int FlowDischarge(BSTR inTINName, BSTR inDEMName, BSTR inRunoffDirName, BSTR inBackBMPName, BSTR inColorTableName,
	BSTR inFirstDayDepthName, double inManingIndData, BSTR inLName, BSTR inCName, BSTR inSName, double intelTimeData, double endTimeData, double outPntXData,
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

	//�򿪲�������Դ���DEMӰ�� 
	CComBase comBase;
	int imgDEMWidth = 0, imgDEMHeight = 0;
	double dx = 0, dy = 0, Xmin = 0, Ymax = 0;
	if (!comBase.OpenImg(inDEMName, imgDEMWidth, imgDEMHeight, dx, dy, Xmin, Ymax)) {
		return 0;
	}
	double geoTransform[6]{Xmin, dx, 0, Ymax, 0, -dy};

	//��LӰ��
	float *pL = new float[imgDEMWidth*imgDEMHeight];
	if (!comBase.OpenImg(inLName, imgDEMWidth, imgDEMHeight, dx, dy, Xmin, Ymax, pL))
		return FALSE;

	//��CӰ��
	float *pC = new float[imgDEMWidth*imgDEMHeight];
	if (!comBase.OpenImg(inCName, imgDEMWidth, imgDEMHeight, dx, dy, Xmin, Ymax, pC))
		return FALSE;

	//��SӰ��
	float *pS = new float[imgDEMWidth*imgDEMHeight];
	if (!comBase.OpenImg(inSName, imgDEMWidth, imgDEMHeight, dx, dy, Xmin, Ymax, pS))
		return FALSE;

	//��ȡTIN�ͽ���Դ���DEM���ɾ�̬��FlowTrack(polylines)
	TIN2Flow tin;
	if (!tin.CreateFlow(strTINName, inDEMName)) {
		return 0;
	}
 	std::vector<polyline3D> polylines; //����ÿһ����ʼ���γɵ������߿���һ��polyline3D,����polyline3D�����괮����inteltimeʱ����ʱ���ɵ�����㡣
	tin.CalFlowPnts(intelTimeData, polylines, strFirstDayDepth, inManingIndData, pL, pC, pS, outPntXData, outPntYData, outPntZData);

	//��BMP����Ӱ��
	BITMAPFILEHEADER BmpHeader; //��¼"Bmp"��ʽͼ���ͷ��Ϣ
	BITMAPINFO BmpInfo;			//��¼"Bmp"��ʽͼ�����Ϣ
	int imgWidth = 0, imgHeight = 0;
	if (!comBase.OpenBMP(inBackBMPName, imgWidth, imgHeight, BmpHeader, BmpInfo)) {
		return 0;
	}
	BYTE* pInColor = new BYTE[3 * imgWidth*imgHeight];
	if (!comBase.OpenBMP(inBackBMPName, imgWidth, imgHeight, pInColor)) {
		delete[]pInColor;  pInColor = NULL;
		return 0;
	}

	//���������ɫ��
	std::vector<colorTable>  colTables;
	SetColorTable(inColorTableName, colTables);

	//�������ɵ�BMPӰ��ֵ
	const char* cEvenTime = "0000";
	CString eveTime(cEvenTime);

	//��¼ÿһʱ�̣��������м�����rainIdx[i].size()��������ε�λ����(rainIdx[i])[j]����polylines[i]�е�idxֵ
	int rainEndtime = comBase.OpenDirfile(inRunoffDirName);
	if (rainEndtime < 1) return 0;
	int plineNum = polylines.size();
	intArray*    rainIdx = new intArray[plineNum];    //rainIdx��һ����ά��������rainIdx[plineNum][n]������(rainIdx[plineNum])��ʾÿһ����ˮ���ϼ�¼�Ķ����δ���(rainIdx[plineNum])[n]��¼һ����ˮ���ϵĵ�n�����������һ���ڵ���±꣬���±��0��ʼ(��ˮ���е�һ���ڵ㿪ʼ)��ÿ��һ��ʱ����inteltime���±���+1��
	floatArray*  pflowVal = new floatArray[plineNum]; //pflowVal��һ����ά��������pflowVal[plineNum][n]����rainIdxһһ��Ӧ����¼������εĲ�������Ϣ��

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
		//���Ʊ���
		cudaMemcpy(gpuOutColor, pInColor, 3 * imgWidth * imgHeight * sizeof(BYTE), cudaMemcpyHostToDevice);

		//���²�������μ��뵽rainidx[i]��
		if (t < rainEndtime) {
			//�õ�������Ϣ
			std::shared_ptr<float> prval = cacher.getImages(t);
			cudaMemcpy(gpuRVal, prval.get(), sizeof(float) * imgDEMWidth * imgDEMHeight, cudaMemcpyHostToDevice);

			runMarkRainIdx(gpuFlowTracks, gpuFlowTracksSizes, gpuFlowTrackSize, Xmin, Ymax, dx, dy,
				gpuRVal, imgDEMWidth, gpuRainIdx, gpuFlowVal, plineNum, endTime, gpuRainIdxLen);
		}

		//�������
		cudaMemset(gpuOutNum, 0, sizeof(float) * imgWidth * imgHeight);

		runDrawPoints(gpuFlowTracks, gpuFlowTracksSizes, gpuFlowTrackSize, gpuRainIdx, gpuFlowVal, plineNum,
			endTime, gpuRainIdxLen, imgWidth, imgHeight, dx, dy, Xmin, Ymax,
			gpuOutColor, gpuColorTables, gpuColorTablsSize, gpuOutNum, gpuOutNumSize);

		BYTE* pOutColor = new BYTE[3 * imgWidth*imgHeight];
		cudaMemcpy(pOutColor, gpuOutColor, sizeof(BYTE) * 3 * imgWidth * imgHeight, cudaMemcpyDeviceToHost);

		float* pOutNum = new float[imgWidth*imgHeight];
		cudaMemcpy(pOutNum, gpuOutNum, imgWidth * imgHeight * sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy(cpuRainIdx, gpuRainIdx, sizeof(int) * plineNum * endTime, cudaMemcpyDeviceToHost);

 		eveTime.Format(_T("%.4d"), t + 1);
		bmpThreads.emplace_back(std::thread(SaveBMP, strOutBMPName + eveTime + L".bmp", imgWidth, imgHeight, BmpHeader, BmpInfo, pOutColor));

		CString cstrFileName = strOutDischarge + eveTime + L".img";
		std::string strFileName = CT2A(cstrFileName);
		std::string rainFileName(strFileName);
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

	//ģ��������
	double intelTime = 3600; 
	double endTime = 365;
	double outPntX = 596257.336;
	double outPntY = 5216215.614;
	double outPntZ = 151.590988;

	BSTR inTIN = L".\\TIN\\TIN.shp";
	BSTR inDEM = L"RainfallSourcePoints\\DEM30.img";
	BSTR inRunoffDir = L".\\R2001\\";
	BSTR inBackBMP = L".\\back.bmp";
	BSTR inColorTable = L".\\aColorTable.txt";
	BSTR inFirstDayDepth = L".\\FlowDepthFirstDay\\0120_IDW.img"; 
	BSTR inL = L".\\USL3_Navi.img";
	BSTR inC = L".\\C1.img";
	BSTR inS = L".\\Slope.img";
	BSTR outFlowBMP = L".\\Flow30M\\";
	BSTR outFlowDischarge = L".\\Flow30N\\";
	std::ofstream timeTxt;
	timeTxt.open(".\\time.txt");

	begintime = clock();
	FlowDischarge(inTIN, inDEM, inRunoffDir, inBackBMP, inColorTable, inFirstDayDepth,
		inManingInd, inL, inC, inS, intelTime, endTime, outPntX, outPntY, outPntZ, outFlowBMP, outFlowDischarge);
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