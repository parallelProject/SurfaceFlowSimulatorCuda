// TIN2Flow.h: interface for the TIN2Flow class
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_TIN2FLOW_H__57620B20_6CB9_48DC_9571_E64855CAE6CE__INCLUDED_)
#define AFX_TIN2FLOW_H__57620B20_6CB9_48DC_9571_E64855CAE6CE__INCLUDED_

#include <vector>

#if _MSC_VER > 1000
#pragma once
#endif

#include "shapefil.h"
#include "ComBase.h"

struct node;
struct edge;
struct tin;
struct topoGrid;
class GridCounter;

//��TIN����Flow����Ϣ
class TIN2Flow  
{
    std::vector<tin>    m_tinA;   //��¼����TIN
    std::vector<edge>   m_edgeA;  //��¼���еı�
	std::vector<node*>  m_nodeA;  //��¼���еĽڵ㣬m_flow��ʾ����m_TIN!=NULL��ʾ���ĵ㡣flow����ʼ������m_TIN!=NULL

    //��shp�ļ�������Ϣ������m_tinA��m_nodeA
    BOOL ReadTin(CString fin,std::vector<tin>& tinA,std::vector<node*>& nodeA);

    //��������ͬ�Ľڵ��Ϊͬһ��node����ָ��
    void MapNodes(std::vector<tin>& tinA,std::vector<node*>& nodeA);

    //�����ߣ�����m_edgeA
    void CreateEdges(std::vector<tin>& tinA,std::vector<edge>& edgeA);

    //����ƽ̹���m_degrees
    void SetDegrees(std::vector<node*>& nodeA);

	//��DEM���������ĵ���뵽m_nodeA������ÿһ�����ĵ��ҵ���Ӧ��TIN
	BOOL AddStartNode(std::vector<tin>& tinA, BSTR InDEM, std::vector<node*>& nodeA);

	//��polylines������һ���߶�(������)֮������е�
	//void AddOneToPline(polyline3D& pline, double X0, double Y0, double Z0, double X1, double Y1, double Z1, double& velocity, double& lefttime, double intelTime, double depth, double inManingInd, double Xmin, double Ymax, int imgWidth, int dx, int dy, float *pTWI, float *pL, float *pC, float *pS);
	void AddOneToPline(polyline3D& pline, double X0, double Y0, double Z0, double X1, double Y1, double Z1, double& velocity, double& lefttime, double intelTime, double depth, double inManingInd, double Xmin, double Ymax, int imgWidth, int dx, int dy, float *pL, float *pC, float *pS);

	//����һ���߶�(������)֮���ͨ��ʱ��
	//double CalOneLineTime(double X0, double Y0, double Z0, double X1, double Y1, double Z1, double& velocity, double depth, double inManingInd, double Xmin, double Ymax, int imgWidth, int dx, int dy, float *pTWI, float *pL, float *pC, float *pS);
	double CalOneLineTime(double X0, double Y0, double Z0, double X1, double Y1, double Z1, double& velocity, double depth, double inManingInd, double Xmin, double Ymax, int imgWidth, int dx, int dy, float *pL, float *pC, float *pS);

public:
    TIN2Flow();

	virtual ~TIN2Flow();

	//���������TIN����֮��Ӧ��DEM����Flow������m_nodeA��
	BOOL CreateFlow(CString fin, BSTR InDEM);

	//�����starPnts��ʼ�����е�FlowTrack�ģ���inteltimeʱ������Track������plines������plines�����괮��inteltime���ʱ������������ֵ
	//BOOL CalFlowPnts(double intelTimeData, std::vector<polyline3D>& polylines, CString strFirstDayDepth, double inManingIndData, float *pTWI, float *pL, float *pC, float *pS, double outPntXData, double outPntYData, double outPntZData);
	BOOL CalFlowPnts(double intelTimeData, std::vector<polyline3D>& polylines, CString strFirstDayDepth, double inManingIndData, float *pL, float *pC, float *pS, double outPntXData, double outPntYData, double outPntZData);
};

#endif
