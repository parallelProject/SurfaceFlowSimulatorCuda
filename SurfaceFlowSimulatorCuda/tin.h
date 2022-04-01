// tin.h: interface for the tin class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_TIN_H__B60ED367_CD92_4543_8CC7_F96B8F656F31__INCLUDED_)
#define AFX_TIN_H__B60ED367_CD92_4543_8CC7_F96B8F656F31__INCLUDED_

#include <vector>
#include <math.h>

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

struct node;
struct edge;
struct tin;

struct node
{
    double              m_x, m_y, m_z;

    BOOL                m_isbound;  //是否是最外围的边界点
    unsigned int        m_degree;   //该节点到非平坦点的深度

	std::vector<edge*>  m_edges;    //和该节点相连的所有边。(该点是三角形的角点)
	edge*               m_edge;     //如果this是某条边上的点，
	//则m_edge是此边，否则为NULL。(该点在三角形边上)
    tin*                m_tin;      //该点所在的三角形。
	//如果this是某个三角形内的点，则m_tin是此三角形，否则为NULL。(该点在三角形内)

    node*               m_flow;     //this流向的节点m_to。
	//如果this有流向的节点则m_flow!=null, 否则m_flow=null。

	node(double x,double y,double z):m_degree(0),m_flow(0),
        m_x(x),m_y(y),m_z(z),m_tin(0),m_edge(0),m_isbound(0){}

    void  neighbors(std::vector<node*>& A);  //返回和该节点相连的所有节点A

    edge* find_edge(node* n); //在m_edges中找到一个edge对象，
	//from和to必须是n和this对象，找不到返回NULL。

	BOOL  find_degree();      //判断neighbors相连点中是否有非平.坦点或已设置degree点。
    
    node* find_flow();        //找出m_to。
	//如果m_to是已经存在的节点则返回NULL，否则返回m_to，便于删除。
	//找到的结果m_to放入m_flow中。

};


struct edge
{
    node*               m_from;
    node*               m_to;

    tin*                m_left;  //遇到的第一个TIN（与方向无关）
    tin*                m_right; //遇到的第二个TIN（与方向无关）

    edge(){ memset(this,0,sizeof(*this)); }  

    edge*   create(tin* t,node* a,node* b);   //根据三角形t，和三角形的节点a、b初始化edge

    double  intersect(double x,double y,double a);  //根据x,y和aspect找到和此边交点的参数。
	//返回值=0表示from点，=1表示to点，在(0,1)表示内部点，否则无交点。

    node*   createnode(double t)      //intersect()函数后，需要new一个点。0<t<1
    {
        double x,y,z; node* n;
        x = m_from->m_x+(m_to->m_x-m_from->m_x)*t;
        y = m_from->m_y+(m_to->m_y-m_from->m_y)*t;
        z = m_from->m_z+(m_to->m_z-m_from->m_z)*t;
        n = new node(x,y,z);
        n->m_edge = this;
        return n;
    }

	tin*    find_flow();       //返回为NULL，表示在该边流动，
	//否则流向返回的三角形。

};


struct tin
{
    node*               m_nodes[3];
    double              m_slope;    //坡度
    double              m_aspect;   //坡向

	tin(){ memset(this,0,sizeof(*this));}

    node* center()     //得到tin的中心点
    {
        double x,y,z;
        x = m_nodes[0]->m_x+m_nodes[1]->m_x+m_nodes[2]->m_x;
        y = m_nodes[0]->m_y+m_nodes[1]->m_y+m_nodes[2]->m_y;
        z = m_nodes[0]->m_z+m_nodes[1]->m_z+m_nodes[2]->m_z;
        node* p = new node(x/3,y/3,z/3);
        p->m_tin = this;
        return p;
    }

	node* InPnt(double x, double y, double z)
	{
		node* p = new node(x,y,z);
		p->m_tin = this;
		return p;
	}

	double projectArea()  //tin的投影面积(二维)
	{
		double x0,y0,z0,x1,y1,z1,x2,y2,z2;
		double areaVal=0.0;

		x0 = m_nodes[0]->m_x;
		y0 = m_nodes[0]->m_y;
		z0 = m_nodes[0]->m_z;

		x1 = m_nodes[1]->m_x;
		y1 = m_nodes[1]->m_y;
		z1 = m_nodes[1]->m_z;

		x2 = m_nodes[2]->m_x;
		y2 = m_nodes[2]->m_y;
		z2 = m_nodes[2]->m_z;

		double a = sqrt((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1));
		double b = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
		double c = sqrt((x2-x0)*(x2-x0)+(y2-y0)*(y2-y0));
		double p = (a+b+c)/2;

		areaVal = sqrt(p*(p-a)*(p-b)*(p-c));

		return areaVal;
	}

    node* find_flow(node* n);  //在三角形内部寻找点n的流向点。如果成功则n的m_flow将被赋予其流向，否则n->m_flow仍然为NULL。

};


struct topoGrid   //记录flowtrack的流水线节点所对应的DEM格网关联的topoEdge边。
{
	std::vector<long>    m_topoEdgeid;   //记录该节点所关联的所有topoEdge边。
	std::vector<double>  m_x;            //记录该节点的x坐标。
	std::vector<double>  m_y;            //记录该节点的y坐标。
	std::vector<double>  m_z;            //记录该节点的z坐标。
};


#endif // !defined(AFX_TIN_H__B60ED367_CD92_4543_8CC7_F96B8F656F31__INCLUDED_)
