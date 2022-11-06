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

//���Ӱ��ĺ�׺����
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

//ͨ��Ӱ��ĺ�׺���õ���ӦGDAL�ܴ����������
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

//�ͷ�GDALDataset���ݼ���
void  CComBase::RELEASE(GDALDataset* pData)
{
	if (pData!=NULL)
	{
		GDALClose((GDALDatasetH)pData);
		pData = NULL;
	}
}

//�ж������strImgName�Ƿ���Ч���˴���imgName֧��*.tif & *.img��ʽ��
BOOL  CComBase::bImgNameVerf(CString strImgName,CString& imgFormat,CString& strFormat)
{
    if (strImgName.GetLength()<1)	return FALSE;
    
    imgFormat = GetImgFormat(strImgName);
    strFormat = GetGDALFormat(imgFormat);
	if(imgFormat!="tif" && imgFormat!="img")
	{
		printf("��֧��tif img��ʽ��");
		return FALSE;
	}
    return TRUE;
}

//����һ���µ�Ӱ�񡣷���True˵�������ɹ������򴴽�ʧ�ܡ�
//ImgName ���� ����Ӱ���·����
//imgWidth��imgHeight ���� ����Ӱ��Ŀ�͸ߣ�
//Xmin��Ymax ���� ����Ӱ����X��Y�ķ�Χ��
//dx��dy ���� ����Ӱ��ĸ����ֱ��ʣ�
//invalidVal ���� ����Ӱ���е���Чֵ��
//demZ ���� ����Ӱ����ÿ�����ص�ֵ��
BOOL  CComBase::CreateNewImg(BSTR ImgName, int imgWidth,int imgHeight,double Xmin, double Ymax,  
                             double dx, double dy, double invalidVal, float *demZ)
{
    //������Ӱ���ʽ�Ƿ���ȷ��
    CString strImgName(ImgName);	
	BOOL flag = CreateNewImg(strImgName,imgWidth,imgHeight,
                             Xmin,Ymax,dx,dy,invalidVal,demZ);

    return flag;
}

//����һ���µ�Ӱ�񡣷���True˵�������ɹ������򴴽�ʧ�ܡ�
//ImgName ���� ����Ӱ���·����
//imgWidth��imgHeight ���� ����Ӱ��Ŀ�͸ߣ�
//Xmin��Ymax ���� ����Ӱ����X��Y�ķ�Χ��
//dx��dy ���� ����Ӱ��ĸ����ֱ��ʣ�
//invalidVal ���� ����Ӱ���е���Чֵ��
//demZ ���� ����Ӱ����ÿ�����ص�ֵ��
BOOL  CComBase::CreateNewImg(CString strImgName, int imgWidth,int imgHeight,double Xmin, double Ymax,  
                             double dx, double dy, double invalidVal, float *demZ)
{
    //������Ӱ���ʽ�Ƿ���ȷ
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
   
    pBand->SetNoDataValue((double)invalidVal);  //���ò��ε���Чֵ����������Ӱ��=invalidValΪ��Чֵ��
	pBand->RasterIO( GF_Write,0, 0, imgWidth, imgHeight, demZ ,imgWidth ,imgHeight ,GDT_Float32, 0, 0);  
    double min=0,max=0,mean=0,dev=0; //���õ�Ӱ������ֵ����Сֵ��
	pBand->ComputeStatistics(FALSE,&min,&max,&mean,&dev,NULL,NULL);

	GDALDeleteDataset(hDriver, cStrImgName);
	RELEASE(pDataset); 

    return TRUE;
}


//��һӰ�����ݣ���ø�Ӱ�����Ϣ������True˵�������ɹ������򴴽�ʧ�ܡ�
//ImgName ���� ����Ӱ���·����
//imgWidth,imgHeight ���� ����Ӱ��Ŀ��ߣ�
BOOL  CComBase::OpenImg(BSTR ImgName, int& imgWidth, int& imgHeight)
{
    CString strImgName(ImgName);
    if (strImgName.GetLength()<1)
	{
		printf("lengthС��1������FALSE \n");
		return FALSE;
	}

    //��DEMӰ��    
    GDALDataset *pImgDataset=NULL;
	//�˴��иĶ�
	USES_CONVERSION;
	LPCSTR charImgName=T2A(strImgName);

	int len = strlen(charImgName);
	//printf("����Ϊ��%d\n",len);
	
	//printf(strImgName);
	pImgDataset = (GDALDataset *)GDALOpen(charImgName, GA_ReadOnly);
	if( pImgDataset == NULL ) 
	{
		printf("Image cannot open! here");
		return FALSE;
	}	
    
    imgHeight = pImgDataset->GetRasterYSize();
	imgWidth = pImgDataset->GetRasterXSize();

	//printf("�Ѿ��򿪣���Ϊ%d,��Ϊ%d\n",imgHeight,imgWidth);

    RELEASE(pImgDataset); 

    if (imgHeight<=0||imgWidth<=0)
        return FALSE;

    return TRUE;
}

//��һӰ�����ݣ���ø�Ӱ�����Ϣ������True˵�������ɹ������򴴽�ʧ�ܡ�
//ImgName ���� ����Ӱ���·����
//imgWidth,imgHeight ���� ����Ӱ��Ŀ��ߣ�
BOOL  CComBase::OpenImg(CString strImgName, int& imgWidth, int& imgHeight)
{
    if (strImgName.GetLength()<1)	
		return FALSE;

    //��DEMӰ��    
    GDALDataset *pImgDataset=NULL;
	//�˴��иĶ�
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


//��һӰ�����ݣ���ø�Ӱ�����Ϣ������True˵�������ɹ������򴴽�ʧ�ܡ�
//ImgName ���� ����Ӱ���·����
//imgWidth,imgHeight ���� Ӱ��Ŀ��ߣ�
//dx,dy ���� ����Ӱ��ķֱ��ʣ�
//Xmin,Ymax ���� ����Ӱ������귶Χ
BOOL  CComBase::OpenImg(BSTR ImgName,int& imgWidth, int& imgHeight,double& dx, double& dy, double& Xmin, double& Ymax)
{
	CString strImgName(ImgName);
    if (strImgName.GetLength()<1)	
		return FALSE;

    //��DEMӰ��    
    GDALDataset *pImgDataset=NULL;
	//�˴��иĶ�
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
	pImgDataset->GetGeoTransform(geoTransform);  //���DEM�����귶Χ���Լ����طֱ��ʡ�
	dx = geoTransform[1], dy = fabs(geoTransform[5]);
    Xmin = geoTransform[0], Ymax = geoTransform[3];  
	//printf("dx= %lf,dy= %lf,Xmin= %lf,Ymax= %lf",dx,dy,Xmin,Ymax);
	RELEASE(pImgDataset); 

    return TRUE;
}
    
//��һӰ�����ݣ���ø�Ӱ�����Ϣ������True˵�������ɹ������򴴽�ʧ�ܡ�
//ImgName ���� ����Ӱ���·����
//imgWidth,imgHeight ���� Ӱ��Ŀ��ߣ�
//pBuffer ���� ����Ӱ��ĸ�������ֵ��
BOOL  CComBase::OpenImg(CString strImgName,int imgWidth, int imgHeight, float *pBuffer)
{
    if (strImgName.GetLength()<1||pBuffer==NULL)	
		return FALSE;

    //��DEMӰ��    
    GDALDataset *pImgDataset=NULL;
	//�˴��иĶ�
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

//��һӰ�����ݣ���ø�Ӱ�����Ϣ������True˵�������ɹ������򴴽�ʧ�ܡ�
//ImgName ���� ����Ӱ���·����
//imgWidth,imgHeight ���� Ӱ��Ŀ��ߣ�
//dx,dy ���� ����Ӱ��ķֱ��ʣ�
//Xmin,Ymax ���� ����Ӱ������귶Χ��
//pBuffer ���� ����Ӱ��ĸ�������ֵ��
BOOL  CComBase::OpenImg(BSTR ImgName,int imgWidth, int imgHeight,double& dx, double& dy, double& Xmin, double& Ymax, float *pBuffer)
{
    CString strImgName(ImgName);
    if (strImgName.GetLength()<1||pBuffer==NULL)	
		return FALSE;

    //��DEMӰ��    
    GDALDataset *pImgDataset=NULL;
	//�˴��иĶ�
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
	pImgDataset->GetGeoTransform(geoTransform);  //���DEM�����귶Χ���Լ����طֱ��ʡ�
	dx = geoTransform[1], dy = fabs(geoTransform[5]);
    Xmin = geoTransform[0], Ymax = geoTransform[3];  

	//�˴��иĶ�
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

//��һӰ�����ݣ���ø�Ӱ�����Ϣ������True˵�������ɹ������򴴽�ʧ�ܡ�
//ImgName ���� ����Ӱ���·����
//imgWidth,imgHeight ���� Ӱ��Ŀ��ߣ�
//dx,dy ���� ����Ӱ��ķֱ��ʣ�
//Xmin,Ymax ���� ����Ӱ������귶Χ��
//pBuffer ���� ����Ӱ��ĸ�������ֵ��
BOOL  CComBase::OpenImg(CString strImgName,int imgWidth, int imgHeight,double& dx, double& dy, double& Xmin, double& Ymax, float *pBuffer)
{
    if (strImgName.GetLength()<1||pBuffer==NULL)	
		return FALSE;

    //��DEMӰ��    
    GDALDataset *pImgDataset=NULL;
	//�˴��иĶ�
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
	pImgDataset->GetGeoTransform(geoTransform);  //���DEM�����귶Χ���Լ����طֱ��ʡ�
	dx = geoTransform[1], dy = fabs(geoTransform[5]);
    Xmin = geoTransform[0], Ymax = geoTransform[3];  

	//�˴��иĶ�
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

//��һӰ�����ݣ���ø�Ӱ�����Ϣ(����ͶӰ��Ϣ)������True˵�������ɹ������򴴽�ʧ�ܡ�
//ImgName ���� ����Ӱ���·����
//imgWidth,imgHeight ���� Ӱ��Ŀ��ߣ�
//dx,dy ���� ����Ӱ��ķֱ��ʣ�
//Xmin,Ymax ���� ����Ӱ������귶Χ��
//pBuffer ���� ����Ӱ��ĸ�������ֵ��
BOOL  CComBase::OpenImg(BSTR ImgName,int imgWidth, int imgHeight,double& dx, double& dy, double& Xmin, double& Ymax, CString& projRef, float *pBuffer)
{
    CString strImgName(ImgName);
    if (strImgName.GetLength()<1||pBuffer==NULL)
	{
		printf("strImgName.GetLength()<1||pBuffer==NULL,����FALSE");
		return FALSE;
	}

    //��DEMӰ��    
    GDALDataset *pImgDataset=NULL;
	//�˴��иĶ�
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
	pImgDataset->GetGeoTransform(geoTransform);  //���DEM�����귶Χ���Լ����طֱ��ʡ�
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


//��֪�����������������ε������
double CComBase::AreaTrig(double x1, double y1, double x2, double y2, double x3, double y3)
{
	double a = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
	double b = sqrt((x1-x3)*(x1-x3)+(y1-y3)*(y1-y3));
	double c = sqrt((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2));
	
	double p = (a+b+c)/2;

	double area= sqrt(p*(p-a)*(p-b)*(p-c));
	return area;
}

//�жϵ�(cenX,cenY,cenZ)�Ƿ���������(xt,yt,zt)����(�����߽�)���ж�һ���㵽�����θ��������������Ƿ���������������
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

//�жϵ�(cenX,cenY,cenZ)�Ƿ���������(xt,yt,zt)�ı߽��ϡ��ж�һ���㵽�������˵������Ƿ�����ߵľ��롣
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

//��BMPӰ�����ݣ���ø�Ӱ�����Ϣ������True˵�������ɹ������򴴽�ʧ�ܡ�
//ImgName ���� ����Ӱ���·����
//imgWidth,imgHeight ���� Ӱ��Ŀ��ߣ�
BOOL  CComBase::OpenBMP(BSTR ImgName, int& imgWidth, int& imgHeight, BITMAPFILEHEADER& BmpHeader, BITMAPINFO& BmpInfo)
{
	CString strImgName(ImgName);
	if (strImgName.GetLength()<1)
		return FALSE;

	CFile orgfile;
	orgfile.Open(strImgName, CFile::modeRead);
	orgfile.Read((LPSTR)&BmpHeader, sizeof(BmpHeader)); //���ȶ��ļ�ͷ���ֵ�BmpHeader�С�
	orgfile.Read(&BmpInfo, sizeof(BITMAPINFO));  //���Ŷ�Bmpͼ�����Ϣ��m_BmpInfo�С�
	int InfoWidth = BmpInfo.bmiHeader.biWidth;  //���ͼ��Ŀ�
	imgHeight = BmpInfo.bmiHeader.biHeight;  //���ͼ��ĸ�
	imgWidth = ((InfoWidth + 3) >> 2) << 2; //bmp�ڴ�����ݵ�ʱ��ÿһ�п��ֵ������4��������������Ӱ������15����ôʵ�ʴ�ŵ�ʱ����Ӹ�1�չ�16��

	orgfile.Close();

	return TRUE;
}

//��BMPӰ�����ݣ���ø�Ӱ�����Ϣ������True˵�������ɹ������򴴽�ʧ�ܡ�
//ImgName ���� ����Ӱ���·����
//imgWidth,imgHeight ���� Ӱ��Ŀ��ߣ�
//BmpHeader ���� Ӱ����ļ�ͷ��Ϣ��
//BmpInfo ���� Ӱ����ļ���Ϣ��
//pBuffer ���� ����Ӱ��ĸ�������ֵ��BYTE[3*imgWidth*imgHeight]
BOOL  CComBase::OpenBMP(BSTR ImgName, int imgWidth, int imgHeight, BYTE *pBuffer)
{
	CString strImgName(ImgName);
	if (strImgName.GetLength()<1)
		return FALSE;

	CFile orgfile;
	orgfile.Open(strImgName, CFile::modeRead);

	BITMAPFILEHEADER BmpHeader; //��¼"Bmp"��ʽͼ���ͷ��Ϣ
	BITMAPINFO       BmpInfo;   //��¼"Bmp"��ʽͼ�����Ϣ

	orgfile.Read((LPSTR)&BmpHeader, sizeof(BmpHeader)); //���ȶ��ļ�ͷ���ֵ�BmpHeader�С�
	orgfile.Read(&BmpInfo, sizeof(BITMAPINFO));  //���Ŷ�Bmpͼ�����Ϣ��m_BmpInfo�С�
	int InfoWidth = BmpInfo.bmiHeader.biWidth;  //���ͼ��Ŀ�
	imgHeight = BmpInfo.bmiHeader.biHeight;  //���ͼ��ĸ�
	imgWidth = ((InfoWidth + 3) >> 2) << 2; //bmp�ڴ�����ݵ�ʱ��ÿһ�п��ֵ������4��������������Ӱ������15����ôʵ�ʴ�ŵ�ʱ����Ӹ�1�չ�16��

	orgfile.Read(pBuffer, 3 * imgWidth*imgHeight);

	orgfile.Close();

	return TRUE;
}

//��һ��Ŀ¼InDir,���ظ�Ŀ¼�µ��ļ�����Ŀ��
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

//��һ��Ŀ¼InDir��������е�fileNum��Ӱ���������Ϣ������pVal�С�imgWidth, imgHeightΪ��Ӱ��Ŀ�Ⱥ͸߶ȡ�
//nCount��Ϊ��Ŀ¼��Ӱ�������Ŀ��
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

//����BMPӰ�����ݡ�����True˵������ɹ������򱣴�ʧ�ܡ�
//ImgName ���� ����Ӱ���·����
//imgWidth,imgHeight ���� Ӱ��Ŀ��ߣ�
//BmpHeader ���� Ӱ����ļ�ͷ��Ϣ��
//BmpInfo ���� Ӱ����ļ���Ϣ��
//pBuffer ���� Ӱ��ĸ�������ֵ��BYTE[3*imgWidth*imgHeight]
BOOL  CComBase::SaveBMP(CString strImgName, int imgWidth, int imgHeight, BITMAPFILEHEADER& BmpHeader, BITMAPINFO& BmpInfo, BYTE *pBuffer)
{
	CFile outfile;
	outfile.Open(strImgName, CFile::modeCreate | CFile::modeWrite);
	outfile.Write((LPSTR)&BmpHeader, sizeof(BmpHeader)); //���ȶ��ļ�ͷ���ֵ�BmpHeader�С�
	outfile.Write(&BmpInfo, sizeof(BITMAPINFO));
	outfile.Write(pBuffer, 3 * imgWidth*imgHeight);

	outfile.Flush();
	outfile.Close();

	return TRUE;
}