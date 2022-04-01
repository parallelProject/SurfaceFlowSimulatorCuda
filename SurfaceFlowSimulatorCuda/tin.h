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

    BOOL                m_isbound;  //�Ƿ�������Χ�ı߽��
    unsigned int        m_degree;   //�ýڵ㵽��ƽ̹������

	std::vector<edge*>  m_edges;    //�͸ýڵ����������бߡ�(�õ��������εĽǵ�)
	edge*               m_edge;     //���this��ĳ�����ϵĵ㣬
	//��m_edge�Ǵ˱ߣ�����ΪNULL��(�õ��������α���)
    tin*                m_tin;      //�õ����ڵ������Ρ�
	//���this��ĳ���������ڵĵ㣬��m_tin�Ǵ������Σ�����ΪNULL��(�õ�����������)

    node*               m_flow;     //this����Ľڵ�m_to��
	//���this������Ľڵ���m_flow!=null, ����m_flow=null��

	node(double x,double y,double z):m_degree(0),m_flow(0),
        m_x(x),m_y(y),m_z(z),m_tin(0),m_edge(0),m_isbound(0){}

    void  neighbors(std::vector<node*>& A);  //���غ͸ýڵ����������нڵ�A

    edge* find_edge(node* n); //��m_edges���ҵ�һ��edge����
	//from��to������n��this�����Ҳ�������NULL��

	BOOL  find_degree();      //�ж�neighbors���������Ƿ��з�ƽ.̹���������degree�㡣
    
    node* find_flow();        //�ҳ�m_to��
	//���m_to���Ѿ����ڵĽڵ��򷵻�NULL�����򷵻�m_to������ɾ����
	//�ҵ��Ľ��m_to����m_flow�С�

};


struct edge
{
    node*               m_from;
    node*               m_to;

    tin*                m_left;  //�����ĵ�һ��TIN���뷽���޹أ�
    tin*                m_right; //�����ĵڶ���TIN���뷽���޹أ�

    edge(){ memset(this,0,sizeof(*this)); }  

    edge*   create(tin* t,node* a,node* b);   //����������t���������εĽڵ�a��b��ʼ��edge

    double  intersect(double x,double y,double a);  //����x,y��aspect�ҵ��ʹ˱߽���Ĳ�����
	//����ֵ=0��ʾfrom�㣬=1��ʾto�㣬��(0,1)��ʾ�ڲ��㣬�����޽��㡣

    node*   createnode(double t)      //intersect()��������Ҫnewһ���㡣0<t<1
    {
        double x,y,z; node* n;
        x = m_from->m_x+(m_to->m_x-m_from->m_x)*t;
        y = m_from->m_y+(m_to->m_y-m_from->m_y)*t;
        z = m_from->m_z+(m_to->m_z-m_from->m_z)*t;
        n = new node(x,y,z);
        n->m_edge = this;
        return n;
    }

	tin*    find_flow();       //����ΪNULL����ʾ�ڸñ�������
	//�������򷵻ص������Ρ�

};


struct tin
{
    node*               m_nodes[3];
    double              m_slope;    //�¶�
    double              m_aspect;   //����

	tin(){ memset(this,0,sizeof(*this));}

    node* center()     //�õ�tin�����ĵ�
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

	double projectArea()  //tin��ͶӰ���(��ά)
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

    node* find_flow(node* n);  //���������ڲ�Ѱ�ҵ�n������㡣����ɹ���n��m_flow�������������򣬷���n->m_flow��ȻΪNULL��

};


struct topoGrid   //��¼flowtrack����ˮ�߽ڵ�����Ӧ��DEM����������topoEdge�ߡ�
{
	std::vector<long>    m_topoEdgeid;   //��¼�ýڵ�������������topoEdge�ߡ�
	std::vector<double>  m_x;            //��¼�ýڵ��x���ꡣ
	std::vector<double>  m_y;            //��¼�ýڵ��y���ꡣ
	std::vector<double>  m_z;            //��¼�ýڵ��z���ꡣ
};


#endif // !defined(AFX_TIN_H__B60ED367_CD92_4543_8CC7_F96B8F656F31__INCLUDED_)
