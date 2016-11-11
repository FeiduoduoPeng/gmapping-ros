#ifndef SENSORREADING_H
#define SENSORREADING_H

#include "sensor.h"
namespace GMapping{

/*��������ȡ�Ļ���*/
class SensorReading{
	public:
		SensorReading(const Sensor* s=0, double time=0);
		inline double getTime() const {return m_time;}
		inline const Sensor* getSensor() const {return m_sensor;}
	protected:
		double m_time;				//ʱ��
		const Sensor* m_sensor;		//������
};

}; //end namespace
#endif
