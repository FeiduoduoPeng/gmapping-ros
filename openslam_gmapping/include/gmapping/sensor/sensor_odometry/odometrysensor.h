#ifndef ODOMETRYSENSOR_H
#define ODOMETRYSENSOR_H

#include <string>
#include <gmapping/sensor/sensor_base/sensor.h>

namespace GMapping{

//��̼ƴ�����  
class OdometrySensor: public Sensor{
	public:
		OdometrySensor(const std::string& name, bool ideal=false);
		inline bool isIdeal() const { return m_ideal; }
	protected:
		bool m_ideal;			//�����Ƿ�������Ĵ�����
};

};

#endif

