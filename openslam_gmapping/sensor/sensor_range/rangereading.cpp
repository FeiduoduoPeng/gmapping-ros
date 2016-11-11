#include <limits>
#include <iostream>
#include <assert.h>
#include <sys/types.h>
#include <gmapping/utils/gvalues.h>
#include <gmapping/sensor/sensor_range/rangereading.h>

namespace GMapping{

using namespace std;


/*���캯��*/
RangeReading::RangeReading(const RangeSensor* rs, double time):
	SensorReading(rs,time){}

/*���캯��
n_beams  ָ���˼�����������
d        ��ʾ��������ľ���
*/
RangeReading::RangeReading(unsigned int n_beams, const double* d, const RangeSensor* rs, double time):
	SensorReading(rs,time)
{
	assert(n_beams==rs->beams().size());
	resize(n_beams);
	for (unsigned int i=0; i<size(); i++)
		(*this)[i]=d[i];
}

RangeReading::~RangeReading(){
//	cerr << __PRETTY_FUNCTION__ << ": CAZZZZZZZZZZZZZZZZZZZZOOOOOOOOOOO" << endl;
}

/*
���ؾ���density����֮��ļ�������ֵ
���б����˵��ļ����ֵ��������ֵΪMAX
*/
unsigned int RangeReading::rawView(double* v, double density) const{
	//���density����0 �����м�����ȫ���Ϸ�
	if (density==0)
	{
		for (unsigned int i=0; i<size(); i++)
			v[i]=(*this)[i];
	}
	else 
	{
		Point lastPoint(0,0);
		uint suppressed=0;
		for (unsigned int i=0; i<size(); i++)
		{
			const RangeSensor* rs=dynamic_cast<const RangeSensor*>(getSensor());
			assert(rs);
			Point lp(
				cos(rs->beams()[i].pose.theta)*(*this)[i],
				sin(rs->beams()[i].pose.theta)*(*this)[i]);
			Point dp=lastPoint-lp;
			double distance=sqrt(dp*dp);
			/*ֻ����������ľ������ĳһ��ֵ��֮�󣬲��ǺϷ���*/
			if (distance<density)
			{
				v[i]=std::numeric_limits<double>::max();
				suppressed++;
			}
			else
			{
				lastPoint=lp;
				v[i]=(*this)[i];
			}
		}
	}
	//	return size();
	return static_cast<unsigned int>(size());

};

/*
������ָ�������ܶȵ�ʱ���ж��ٵ���Ч������
density��ʾ��������Ķ˵����������density
������Ч�ļ��������(�������������ڽӽ�������Ϊ��Ч)
*/
unsigned int RangeReading::activeBeams(double density) const{
	if (density==0.)
		return size();
		int ab=0;
	Point lastPoint(0,0);
	uint suppressed=0;
	for (unsigned int i=0; i<size(); i++)
	{
		const RangeSensor* rs=dynamic_cast<const RangeSensor*>(getSensor());
		assert(rs);
		Point lp(
			cos(rs->beams()[i].pose.theta)*(*this)[i],
			sin(rs->beams()[i].pose.theta)*(*this)[i]);
		Point dp=lastPoint-lp;
		double distance=sqrt(dp*dp);
		/*����������ڽӽ� ����Ϊ��Ч*/
		if (distance<density)
		{
			suppressed++;
		}
		/*������Ϊ��Ч��������+1*/
		else
		{
			lastPoint=lp;
			ab++;
		}
		//std::cerr<< __PRETTY_FUNCTION__ << std::endl;
		//std::cerr<< "suppressed " << suppressed <<"/"<<size() << std::endl;
	}
	return ab;
}

/*�Ѽ�������ת��Ϊ�ڵѿ�������ϵ�µ�����*/
std::vector<Point> RangeReading::cartesianForm(double maxRange) const{
	const RangeSensor* rangeSensor=dynamic_cast<const RangeSensor*>(getSensor());
	assert(rangeSensor && rangeSensor->beams().size());
	//	uint m_beams=rangeSensor->beams().size();
	uint m_beams=static_cast<unsigned int>(rangeSensor->beams().size());
	std::vector<Point> cartesianPoints(m_beams);
	double px,py,ps,pc;
	px=rangeSensor->getPose().x;
	py=rangeSensor->getPose().y;
	ps=sin(rangeSensor->getPose().theta);
	pc=cos(rangeSensor->getPose().theta);
	for (unsigned int i=0; i<m_beams; i++)
	{
		const double& rho=(*this)[i];
		const double& s=rangeSensor->beams()[i].s;
		const double& c=rangeSensor->beams()[i].c;
		if (rho>=maxRange)
		{
			cartesianPoints[i]=Point(0,0);
		} 
		else 
		{
			Point p=Point(rangeSensor->beams()[i].pose.x+c*rho, rangeSensor->beams()[i].pose.y+s*rho);
			cartesianPoints[i].x=px+pc*p.x-ps*p.y;
			cartesianPoints[i].y=py+ps*p.x+pc*p.y;
		}
	}
	return cartesianPoints;
}

};

