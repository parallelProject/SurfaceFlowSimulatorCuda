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

//由TIN计算Flow的信息
class TIN2Flow  
{
    std::vector<tin>    m_tinA;   //记录所有TIN
    std::vector<edge>   m_edgeA;  //记录所有的边
	std::vector<node*>  m_nodeA;  //记录所有的节点，m_flow表示流向。m_TIN!=NULL表示重心点。flow的起始点来自m_TIN!=NULL

    //从shp文件读出信息，放入m_tinA和m_nodeA
    BOOL ReadTin(CString fin,std::vector<tin>& tinA,std::vector<node*>& nodeA);

    //将坐标相同的节点归为同一个node对象指针
    void MapNodes(std::vector<tin>& tinA,std::vector<node*>& nodeA);

    //创建边，放入m_edgeA
    void CreateEdges(std::vector<tin>& tinA,std::vector<edge>& edgeA);

    //设置平坦点的m_degrees
    void SetDegrees(std::vector<node*>& nodeA);

	//将DEM的所有中心点加入到m_nodeA，并给每一个中心点找到对应的TIN
	BOOL AddStartNode(std::vector<tin>& tinA, BSTR InDEM, std::vector<node*>& nodeA);

	//向polylines中增加一条线段(两个点)之间的所有点
	//void AddOneToPline(polyline3D& pline, double X0, double Y0, double Z0, double X1, double Y1, double Z1, double& velocity, double& lefttime, double intelTime, double depth, double inManingInd, double Xmin, double Ymax, int imgWidth, int dx, int dy, float *pTWI, float *pL, float *pC, float *pS);
	void AddOneToPline(polyline3D& pline, double X0, double Y0, double Z0, double X1, double Y1, double Z1, double& velocity, double& lefttime, double intelTime, double depth, double inManingInd, double Xmin, double Ymax, int imgWidth, int dx, int dy, float *pL, float *pC, float *pS);

	//计算一条线段(两个点)之间的通过时间
	//double CalOneLineTime(double X0, double Y0, double Z0, double X1, double Y1, double Z1, double& velocity, double depth, double inManingInd, double Xmin, double Ymax, int imgWidth, int dx, int dy, float *pTWI, float *pL, float *pC, float *pS);
	double CalOneLineTime(double X0, double Y0, double Z0, double X1, double Y1, double Z1, double& velocity, double depth, double inManingInd, double Xmin, double Ymax, int imgWidth, int dx, int dy, float *pL, float *pC, float *pS);

public:
    TIN2Flow();

	virtual ~TIN2Flow();

	//根据输入的TIN和与之对应的DEM创建Flow，放入m_nodeA中
	BOOL CreateFlow(CString fin, BSTR InDEM);

	//计算从starPnts开始的所有的FlowTrack的，在inteltime时间间隔的Track流域线plines。其中plines的坐标串是inteltime间隔时间计算出的坐标值
	//BOOL CalFlowPnts(double intelTimeData, std::vector<polyline3D>& polylines, CString strFirstDayDepth, double inManingIndData, float *pTWI, float *pL, float *pC, float *pS, double outPntXData, double outPntYData, double outPntZData);
	BOOL CalFlowPnts(double intelTimeData, std::vector<polyline3D>& polylines, CString strFirstDayDepth, double inManingIndData, float *pL, float *pC, float *pS, double outPntXData, double outPntYData, double outPntZData);
};

#endif
