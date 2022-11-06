// ComBase.cpp: implementation of the CComBase class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "resource.h"
#include "ComBase.h"

#include <string>
#include <cstringt.h>
using namespace std;

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CComBase::CComBase()
{

}

CComBase::~CComBase()
{

}

//获得影像的后缀名。
CString CComBase::GetImgFormat(LPCTSTR lpstrFileName)
{
	//string saveFileName=lpstrFileName;
	string saveFileName=string((char*)lpstrFileName);

	string suffix = "";		 

	if (saveFileName.length() != 0 )
	{
		const char* charName ;
		charName = saveFileName.data();
		charName = strrchr(charName,'.');
		if (charName) 
			suffix = charName+1;		
	}
	//CString result = suffix.c_str();
	CString result(suffix.c_str());
	
	return result;
}

//通过影像的后缀名得到对应GDAL能处理的类型名
CString CComBase::GetGDALFormat(CString imgFormat)
{
	//CString strFormat="";
	CString strFormat;

	if (imgFormat.GetLength() == 0 ) 
	{
		return strFormat;
	}

	if (imgFormat=="bmp")
		strFormat = "BMP";
	else if (imgFormat=="bt")
		strFormat = "BT";
	else if (imgFormat=="gif")
		strFormat = "GIF";
	else if (imgFormat=="img")
		strFormat = "HFA";
	else if (imgFormat=="jpg")
		strFormat = "JPEG";
	else if (imgFormat=="png")
		strFormat = "PNG";
	else if (imgFormat=="tif")
		strFormat = "GTiff";
	else if (imgFormat=="vrt")
		strFormat = "VRT";

	return strFormat;
}

//释放GDALDataset数据集。
void  CComBase::RELEASE(GDALDataset* pData)
{
	if (pData!=NULL)
	{
		GDALClose((GDALDatasetH)pData);
		pData = NULL;
	}
}

//判断输入的strImgName是否有效，此处的imgName支持*.tif & *.img格式。
BOOL  CComBase::bImgNameVerf(CString strImgName,CString& imgFormat,CString& strFormat)
{
    if (strImgName.GetLength()<1)	return FALSE;
    
    imgFormat = GetImgFormat(strImgName);
    strFormat = GetGDALFormat(imgFormat);
	if(imgFormat!="tif" && imgFormat!="img")
	{
		printf("仅支持tif img格式。");
		return FALSE;
	}
    return TRUE;
}

//创建一个新的影像。返回True说明创建成功，否则创建失败。
//ImgName —— 输入影像的路径；
//imgWidth、imgHeight —— 设置影像的宽和高；
//Xmin、Ymax —— 设置影像中X、Y的范围；
//dx、dy —— 设置影像的格网分辨率；
//invalidVal —— 设置影像中的无效值；
//demZ —— 设置影像中每个像素的值。
BOOL  CComBase::CreateNewImg(BSTR ImgName, int imgWidth,int imgHeight,double Xmin, double Ymax,  
                             double dx, double dy, double invalidVal, float *demZ)
{
    //检查输出影像格式是否正确。
    CString strImgName(ImgName);	
	BOOL flag = CreateNewImg(strImgName,imgWidth,imgHeight,
                             Xmin,Ymax,dx,dy,invalidVal,demZ);

    return flag;
}

//创建一个新的影像。返回True说明创建成功，否则创建失败。
//ImgName —— 输入影像的路径；
//imgWidth、imgHeight —— 设置影像的宽和高；
//Xmin、Ymax —— 设置影像中X、Y的范围；
//dx、dy —— 设置影像的格网分辨率；
//invalidVal —— 设置影像中的无效值；
//demZ —— 设置影像中每个像素的值。
BOOL  CComBase::CreateNewImg(CString strImgName, int imgWidth,int imgHeight,double Xmin, double Ymax,  
                             double dx, double dy, double invalidVal, float *demZ)
{
    //检查输出影像格式是否正确
    CString imgFormat,strFormat;
	if (bImgNameVerf(strImgName,imgFormat,strFormat)==FALSE)
        return FALSE;

    GDALDriverH hDriver = NULL;
	const char* cStrFormat = CStringA(strFormat);
	hDriver = GDALGetDriverByName(cStrFormat);
	if( hDriver == NULL || GDALGetMetadataItem( hDriver, GDAL_DCAP_CREATE, NULL ) == NULL )
		return FALSE;
    
	char **papszOptions = NULL;
	const char* cStrImgName = CStringA(strImgName);
	GDALDataset *pDataset = (GDALDataset *)GDALCreate(hDriver, cStrImgName, imgWidth, imgHeight, 1, GDT_Float32, papszOptions);
	double adfGeoTransform[6] = { Xmin, dx, 0, Ymax , 0, -dy };
	pDataset->SetGeoTransform( adfGeoTransform );        
	GDALRasterBand *pBand = pDataset->GetRasterBand(1);
   
    pBand->SetNoDataValue((double)invalidVal);  //设置波段的无效值，这里设置影像=invalidVal为无效值。
	pBand->RasterIO( GF_Write,0, 0, imgWidth, imgHeight, demZ ,imgWidth ,imgHeight ,GDT_Float32, 0, 0);  
    double min=0,max=0,mean=0,dev=0; //设置的影像的最大值和最小值。
	pBand->ComputeStatistics(FALSE,&min,&max,&mean,&dev,NULL,NULL);

	GDALDeleteDataset(hDriver, cStrImgName);
	RELEASE(pDataset); 

    return TRUE;
}


//打开一影像数据，获得该影像的信息。返回True说明创建成功，否则创建失败。
//ImgName —— 输入影像的路径；
//imgWidth,imgHeight —— 返回影像的宽、高；
BOOL  CComBase::OpenImg(BSTR ImgName, int& imgWidth, int& imgHeight)
{
    CString strImgName(ImgName);
    if (strImgName.GetLength()<1)
	{
		printf("length小于1，返回FALSE \n");
		return FALSE;
	}

    //打开DEM影像    
    GDALDataset *pImgDataset=NULL;
	//此处有改动
	USES_CONVERSION;
	LPCSTR charImgName=T2A(strImgName);

	int len = strlen(charImgName);
	//printf("长度为：%d\n",len);
	
	//printf(strImgName);
	pImgDataset = (GDALDataset *)GDALOpen(charImgName, GA_ReadOnly);
	if( pImgDataset == NULL ) 
	{
		printf("Image cannot open! here");
		return FALSE;
	}	
    
    imgHeight = pImgDataset->GetRasterYSize();
	imgWidth = pImgDataset->GetRasterXSize();

	//printf("已经打开，高为%d,宽为%d\n",imgHeight,imgWidth);

    RELEASE(pImgDataset); 

    if (imgHeight<=0||imgWidth<=0)
        return FALSE;

    return TRUE;
}

//打开一影像数据，获得该影像的信息。返回True说明创建成功，否则创建失败。
//ImgName —— 输入影像的路径；
//imgWidth,imgHeight —— 返回影像的宽、高；
BOOL  CComBase::OpenImg(CString strImgName, int& imgWidth, int& imgHeight)
{
    if (strImgName.GetLength()<1)	
		return FALSE;

    //打开DEM影像    
    GDALDataset *pImgDataset=NULL;
	//此处有改动
	//const char* charImgName=(char*)(LPCTSTR)strImgName;
	USES_CONVERSION;
	LPCSTR charImgName=T2A(strImgName);

	pImgDataset = (GDALDataset *)GDALOpen(charImgName, GA_ReadOnly);
	if( pImgDataset == NULL ) 
	{
		printf("Img cannot open!");
		return FALSE;
	}	
    
    imgHeight = pImgDataset->GetRasterYSize();
	imgWidth = pImgDataset->GetRasterXSize();

    RELEASE(pImgDataset); 

    if (imgHeight<=0||imgWidth<=0)
        return FALSE;

    return TRUE;
}


//打开一影像数据，获得该影像的信息。返回True说明创建成功，否则创建失败。
//ImgName —— 输入影像的路径；
//imgWidth,imgHeight —— 影像的宽、高；
//dx,dy —— 返回影像的分辨率；
//Xmin,Ymax —— 返回影像的坐标范围
BOOL  CComBase::OpenImg(BSTR ImgName,int& imgWidth, int& imgHeight,double& dx, double& dy, double& Xmin, double& Ymax)
{
	CString strImgName(ImgName);
    if (strImgName.GetLength()<1)	
		return FALSE;

    //打开DEM影像    
    GDALDataset *pImgDataset=NULL;
	//此处有改动
	//const char* charImgName=(char*)(LPCTSTR)strImgName;
	USES_CONVERSION;
	LPCSTR charImgName=T2A(strImgName);

	pImgDataset = (GDALDataset *)GDALOpen(charImgName, GA_ReadOnly);
	if( pImgDataset == NULL ) 
	{
		printf("Img cannot open!");
		return FALSE;
	}	
    imgHeight = pImgDataset->GetRasterYSize();
	imgWidth = pImgDataset->GetRasterXSize();

    if (imgHeight<=0||imgWidth<=0)
        return FALSE;

	double geoTransform[6];
	pImgDataset->GetGeoTransform(geoTransform);  //获得DEM的坐标范围，以及像素分辨率。
	dx = geoTransform[1], dy = fabs(geoTransform[5]);
    Xmin = geoTransform[0], Ymax = geoTransform[3];  
	//printf("dx= %lf,dy= %lf,Xmin= %lf,Ymax= %lf",dx,dy,Xmin,Ymax);
	RELEASE(pImgDataset); 

    return TRUE;
}
    
//打开一影像数据，获得该影像的信息。返回True说明创建成功，否则创建失败。
//ImgName —— 输入影像的路径；
//imgWidth,imgHeight —— 影像的宽、高；
//pBuffer —— 返回影像的各个像素值。
BOOL  CComBase::OpenImg(CString strImgName,int imgWidth, int imgHeight, float *pBuffer)
{
    if (strImgName.GetLength()<1||pBuffer==NULL)	
		return FALSE;

    //打开DEM影像    
    GDALDataset *pImgDataset=NULL;
	//此处有改动
	//const char* charImgName=(char*)(LPCTSTR)strImgName;
	USES_CONVERSION;
	LPCSTR charImgName=T2A(strImgName);

	pImgDataset = (GDALDataset *)GDALOpen(charImgName, GA_ReadOnly);
	if( pImgDataset == NULL ) 
	{
		printf("Img cannot open!");
		return FALSE;
	}	

    int BandNum = pImgDataset->GetRasterCount();
    GDALRasterBand *pBand = pImgDataset->GetRasterBand(BandNum);	
	if (pBand==NULL)  
    {
        RELEASE(pImgDataset);
        return FALSE;
    }
	pBand ->RasterIO(GF_Read, 0, 0, imgWidth, imgHeight, pBuffer, imgWidth, imgHeight, GDT_Float32, 0, 0);

    RELEASE(pImgDataset); 

    if (pBuffer==NULL)
        return FALSE;

    return TRUE;
}

//打开一影像数据，获得该影像的信息。返回True说明创建成功，否则创建失败。
//ImgName —— 输入影像的路径；
//imgWidth,imgHeight —— 影像的宽、高；
//dx,dy —— 返回影像的分辨率；
//Xmin,Ymax —— 返回影像的坐标范围；
//pBuffer —— 返回影像的各个像素值。
BOOL  CComBase::OpenImg(BSTR ImgName,int imgWidth, int imgHeight,double& dx, double& dy, double& Xmin, double& Ymax, float *pBuffer)
{
    CString strImgName(ImgName);
    if (strImgName.GetLength()<1||pBuffer==NULL)	
		return FALSE;

    //打开DEM影像    
    GDALDataset *pImgDataset=NULL;
	//此处有改动
	//const char* charImgName=(char*)(LPCTSTR)strImgName;
	USES_CONVERSION;
	LPCSTR charImgName=T2A(strImgName);

	pImgDataset = (GDALDataset *)GDALOpen(charImgName, GA_ReadOnly);
	if( pImgDataset == NULL ) 
	{
		printf("Img cannot open!");
		return FALSE;
	}	

    double geoTransform[6];
	pImgDataset->GetGeoTransform(geoTransform);  //获得DEM的坐标范围，以及像素分辨率。
	dx = geoTransform[1], dy = fabs(geoTransform[5]);
    Xmin = geoTransform[0], Ymax = geoTransform[3];  

	//此处有改动
    //CString prj = pImgDataset->GetProjectionRef();
	CString prj(pImgDataset->GetProjectionRef());

    int BandNum = pImgDataset->GetRasterCount();
    GDALRasterBand *pBand = pImgDataset->GetRasterBand(BandNum);	
	if (pBand==NULL)  
    {
        RELEASE(pImgDataset);
        return FALSE;
    }
	pBand ->RasterIO(GF_Read, 0, 0, imgWidth, imgHeight, pBuffer, imgWidth, imgHeight, GDT_Float32, 0, 0);

    RELEASE(pImgDataset); 

    if (pBuffer==NULL||imgWidth<=0||imgHeight<=0||dx<=0||dy<=0)
        return FALSE;

    return TRUE;
}

//打开一影像数据，获得该影像的信息。返回True说明创建成功，否则创建失败。
//ImgName —— 输入影像的路径；
//imgWidth,imgHeight —— 影像的宽、高；
//dx,dy —— 返回影像的分辨率；
//Xmin,Ymax —— 返回影像的坐标范围；
//pBuffer —— 返回影像的各个像素值。
BOOL  CComBase::OpenImg(CString strImgName,int imgWidth, int imgHeight,double& dx, double& dy, double& Xmin, double& Ymax, float *pBuffer)
{
    if (strImgName.GetLength()<1||pBuffer==NULL)	
		return FALSE;

    //打开DEM影像    
    GDALDataset *pImgDataset=NULL;
	//此处有改动
	//const char* charImgName=(char*)(LPCTSTR)strImgName;
	USES_CONVERSION;
	LPCSTR charImgName=T2A(strImgName);

	pImgDataset = (GDALDataset *)GDALOpen(charImgName, GA_ReadOnly);
	if( pImgDataset == NULL ) 
	{
		printf("Img cannot open!");
		return FALSE;
	}	

    double geoTransform[6];
	pImgDataset->GetGeoTransform(geoTransform);  //获得DEM的坐标范围，以及像素分辨率。
	dx = geoTransform[1], dy = fabs(geoTransform[5]);
    Xmin = geoTransform[0], Ymax = geoTransform[3];  

	//此处有改动
    CString prj(pImgDataset->GetProjectionRef());

    int BandNum = pImgDataset->GetRasterCount();
    GDALRasterBand *pBand = pImgDataset->GetRasterBand(BandNum);	
	if (pBand==NULL)  
    {
        RELEASE(pImgDataset);
        return FALSE;
    }
	pBand ->RasterIO(GF_Read, 0, 0, imgWidth, imgHeight, pBuffer, imgWidth, imgHeight, GDT_Float32, 0, 0);

    RELEASE(pImgDataset); 

    if (pBuffer==NULL||imgWidth<=0||imgHeight<=0||dx<=0||dy<=0)
        return FALSE;

    return TRUE;
}

//打开一影像数据，获得该影像的信息(包括投影信息)。返回True说明创建成功，否则创建失败。
//ImgName —— 输入影像的路径；
//imgWidth,imgHeight —— 影像的宽、高；
//dx,dy —— 返回影像的分辨率；
//Xmin,Ymax —— 返回影像的坐标范围；
//pBuffer —— 返回影像的各个像素值。
BOOL  CComBase::OpenImg(BSTR ImgName,int imgWidth, int imgHeight,double& dx, double& dy, double& Xmin, double& Ymax, CString& projRef, float *pBuffer)
{
    CString strImgName(ImgName);
    if (strImgName.GetLength()<1||pBuffer==NULL)
	{
		printf("strImgName.GetLength()<1||pBuffer==NULL,返回FALSE");
		return FALSE;
	}

    //打开DEM影像    
    GDALDataset *pImgDataset=NULL;
	//此处有改动
	//const char* charImgName=(char*)(LPCTSTR)strImgName;
	USES_CONVERSION;
	LPCSTR charImgName=T2A(strImgName);

	pImgDataset = (GDALDataset *)GDALOpen(charImgName, GA_ReadOnly);
	if( pImgDataset == NULL ) 
	{
		printf("img cannot open!");
		return FALSE;
	}	

    double geoTransform[6];
	pImgDataset->GetGeoTransform(geoTransform);  //获得DEM的坐标范围，以及像素分辨率。
	dx = geoTransform[1], dy = fabs(geoTransform[5]);
    Xmin = geoTransform[0], Ymax = geoTransform[3];  

    projRef = pImgDataset->GetProjectionRef();

    int BandNum = pImgDataset->GetRasterCount();
    GDALRasterBand *pBand = pImgDataset->GetRasterBand(BandNum);	
	if (pBand==NULL)  
    {
        RELEASE(pImgDataset);
        return FALSE;
    }
	pBand ->RasterIO(GF_Read, 0, 0, imgWidth, imgHeight, pBuffer, imgWidth, imgHeight, GDT_Float32, 0, 0);
	//
	//
	
    RELEASE(pImgDataset); 

    if (pBuffer==NULL||imgWidth<=0||imgHeight<=0||dx<=0||dy<=0)
        return FALSE;
	printf("return TRUE.\n");
	for(int i=0;i<30;i++)
		printf("pBuffer[%d]=%f\n",i,pBuffer[i]);
    return TRUE;
}


//已知三顶点坐标求三角形的面积。
double CComBase::AreaTrig(double x1, double y1, double x2, double y2, double x3, double y3)
{
	double a = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
	double b = sqrt((x1-x3)*(x1-x3)+(y1-y3)*(y1-y3));
	double c = sqrt((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2));
	
	double p = (a+b+c)/2;

	double area= sqrt(p*(p-a)*(p-b)*(p-c));
	return area;
}

//判断点(cenX,cenY,cenZ)是否在三角形(xt,yt,zt)的内(包括边界)。判断一个点到三角形各个顶点的面积和是否等于三角形面积。
BOOL CComBase::IsInTriangle(double cenX,double cenY,double* xt,double* yt)
{
	double c01=AreaTrig(cenX,cenY,xt[0],yt[0],xt[1],yt[1]);
	double c02=AreaTrig(cenX,cenY,xt[0],yt[0],xt[2],yt[2]);
	double c21=AreaTrig(cenX,cenY,xt[2],yt[2],xt[1],yt[1]);
	
	double trig=AreaTrig(xt[0],yt[0],xt[1],yt[1],xt[2],yt[2]);

	if (fabs(trig-(c01+c02+c21))<0.00001)
		return true;

	return false;
}

//判断点(cenX,cenY,cenZ)是否在三角形(xt,yt,zt)的边界上。判断一个点到线两个端点距离和是否等于线的距离。
BOOL CComBase::IsInEdge(double cenX,double cenY,double* xt,double* yt)
{
	double dis01= sqrt((xt[0]-xt[1])*(xt[0]-xt[1])+(yt[0]-yt[1])*(yt[0]-yt[1]));
	double dis02= sqrt((xt[0]-xt[2])*(xt[0]-xt[2])+(yt[0]-yt[2])*(yt[0]-yt[2]));
	double dis12= sqrt((xt[2]-xt[1])*(xt[2]-xt[1])+(yt[2]-yt[1])*(yt[2]-yt[1]));

	double disc0 = sqrt((cenX-xt[0])*(cenX-xt[0])+(cenY-yt[0])*(cenY-yt[0]));
	double disc1 = sqrt((cenX-xt[1])*(cenX-xt[1])+(cenY-yt[1])*(cenY-yt[1]));
	double disc2 = sqrt((cenX-xt[2])*(cenX-xt[2])+(cenY-yt[2])*(cenY-yt[2]));

	if (fabs((disc0+disc1)-dis01)<0.00001 || fabs((disc0+disc2)-dis02)<0.00001 || fabs((disc2+disc1)-dis12)<0.00001)
		return true;

	return false;
}

//打开BMP影像数据，获得该影像的信息。返回True说明创建成功，否则创建失败。
//ImgName —— 输入影像的路径；
//imgWidth,imgHeight —— 影像的宽、高；
BOOL  CComBase::OpenBMP(BSTR ImgName, int& imgWidth, int& imgHeight, BITMAPFILEHEADER& BmpHeader, BITMAPINFO& BmpInfo)
{
	CString strImgName(ImgName);
	if (strImgName.GetLength()<1)
		return FALSE;

	CFile orgfile;
	orgfile.Open(strImgName, CFile::modeRead);
	orgfile.Read((LPSTR)&BmpHeader, sizeof(BmpHeader)); //首先读文件头部分到BmpHeader中。
	orgfile.Read(&BmpInfo, sizeof(BITMAPINFO));  //接着读Bmp图象的信息到m_BmpInfo中。
	int InfoWidth = BmpInfo.bmiHeader.biWidth;  //获得图像的宽
	imgHeight = BmpInfo.bmiHeader.biHeight;  //获得图像的高
	imgWidth = ((InfoWidth + 3) >> 2) << 2; //bmp在存放数据的时候，每一行宽度值必须是4的整倍数，比如影像宽度是15，那么实际存放的时候会多加个1凑够16。

	orgfile.Close();

	return TRUE;
}

//打开BMP影像数据，获得该影像的信息。返回True说明创建成功，否则创建失败。
//ImgName —— 输入影像的路径；
//imgWidth,imgHeight —— 影像的宽、高；
//BmpHeader —— 影像的文件头信息；
//BmpInfo —— 影像的文件信息；
//pBuffer —— 返回影像的各个像素值。BYTE[3*imgWidth*imgHeight]
BOOL  CComBase::OpenBMP(BSTR ImgName, int imgWidth, int imgHeight, BYTE *pBuffer)
{
	CString strImgName(ImgName);
	if (strImgName.GetLength()<1)
		return FALSE;

	CFile orgfile;
	orgfile.Open(strImgName, CFile::modeRead);

	BITMAPFILEHEADER BmpHeader; //记录"Bmp"格式图象的头信息
	BITMAPINFO       BmpInfo;   //记录"Bmp"格式图象的信息

	orgfile.Read((LPSTR)&BmpHeader, sizeof(BmpHeader)); //首先读文件头部分到BmpHeader中。
	orgfile.Read(&BmpInfo, sizeof(BITMAPINFO));  //接着读Bmp图象的信息到m_BmpInfo中。
	int InfoWidth = BmpInfo.bmiHeader.biWidth;  //获得图像的宽
	imgHeight = BmpInfo.bmiHeader.biHeight;  //获得图像的高
	imgWidth = ((InfoWidth + 3) >> 2) << 2; //bmp在存放数据的时候，每一行宽度值必须是4的整倍数，比如影像宽度是15，那么实际存放的时候会多加个1凑够16。

	orgfile.Read(pBuffer, 3 * imgWidth*imgHeight);

	orgfile.Close();

	return TRUE;
}

//打开一个目录InDir,返回该目录下的文件总数目。
int CComBase::OpenDirfile(BSTR InDir)
{
	int nCount = 0;

	CString strDir(InDir);
	CFileFind ff;

	const char* charLink = "*.*";
	CString strLink(charLink);


	BOOL bFind = ff.FindFile(strDir + strLink);

	while (bFind)
	{
		bFind = ff.FindNextFile();

		if (ff.IsDots() || ff.IsSystem() || ff.IsHidden())
			continue;

		nCount++;
	}

	return nCount;
}

//打开一个目录InDir，获得其中第fileNum幅影像的数据信息，放入pVal中。imgWidth, imgHeight为该影像的宽度和高度。
//nCount—为该目录下影像的总数目。
BOOL CComBase::OpenDirfile(BSTR InDir, int fileNum, int imgWidth, int imgHeight, float* pVal)
{
	int nCount = 0;

	CString strDir(InDir);
	CFileFind ff;

	const char* charLink = "*.*";
	CString strLink(charLink);

	BOOL bFind = ff.FindFile(strDir + strLink);

	while (bFind)
	{
		bFind = ff.FindNextFile();

		if (ff.IsDots() || ff.IsSystem() || ff.IsHidden())
			continue;

		nCount++;

		if (nCount == fileNum + 1)
		{
			CString strFileName = ff.GetFilePath();
			OpenImg(strFileName, imgWidth, imgHeight, pVal);

			break;
		}
	}

	return TRUE;
}

//保存BMP影像数据。返回True说明保存成功，否则保存失败。
//ImgName —— 输入影像的路径；
//imgWidth,imgHeight —— 影像的宽、高；
//BmpHeader —— 影像的文件头信息；
//BmpInfo —— 影像的文件信息；
//pBuffer —— 影像的各个像素值。BYTE[3*imgWidth*imgHeight]
BOOL  CComBase::SaveBMP(CString strImgName, int imgWidth, int imgHeight, BITMAPFILEHEADER& BmpHeader, BITMAPINFO& BmpInfo, BYTE *pBuffer)
{
	CFile outfile;
	outfile.Open(strImgName, CFile::modeCreate | CFile::modeWrite);
	outfile.Write((LPSTR)&BmpHeader, sizeof(BmpHeader)); //首先读文件头部分到BmpHeader中。
	outfile.Write(&BmpInfo, sizeof(BITMAPINFO));
	outfile.Write(pBuffer, 3 * imgWidth*imgHeight);

	outfile.Flush();
	outfile.Close();

	return TRUE;
}