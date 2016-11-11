#ifndef RANGESENSOR_H
#define RANGESENSOR_H

#include <vector>
#include <gmapping/sensor/sensor_base/sensor.h>
#include <gmapping/utils/point.h>

namespace GMapping{

/*���⴫������*/
class RangeSensor: public Sensor{
	friend class Configuration;
	friend class CarmenConfiguration;
	friend class CarmenWrapper;
	public:
		/*������*/
		struct Beam{
			OrientedPoint pose;	//pose relative to the center of the sensor
			double span;	//spam=0 indicates a line-like beam
			double maxRange;	//maximum range of the sensor
			double s,c;		//sinus and cosinus of the beam (optimization);
		};	
		/*���캯�� ����������*/
		RangeSensor(std::string name);
		
		/**/
		RangeSensor(std::string name, unsigned int beams, double res, const OrientedPoint& position=OrientedPoint(0,0,0), double span=0, double maxrange=89.0);
		
		/*�������еļ�������(����)*/
		inline const std::vector<Beam>& beams() const {return m_beams;}
		
		/*�������еļ�������*/
		inline std::vector<Beam>& beams() {return m_beams;}
		
		/*���ؼ��⴫������λ��*/
		inline OrientedPoint getPose() const {return m_pose;}
		
		/*���¼����ѯ��*/
		void updateBeamsLookup();
		bool newFormat;
	protected:
		OrientedPoint m_pose;		//���⴫������λ��
		std::vector<Beam> m_beams;  //��������
};

};

#endif
