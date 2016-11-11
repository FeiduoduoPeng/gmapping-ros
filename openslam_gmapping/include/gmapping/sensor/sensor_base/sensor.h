#ifndef SENSOR_H
#define SENSOR_H

#include <string>
#include <map>

namespace GMapping{

//�������Ļ���  ���봫����������
class Sensor{
	public:
		Sensor(const std::string& name="");
		virtual ~Sensor();
		inline std::string getName() const {return m_name;}
		inline void setName(const std::string& name) {m_name=name;}
	protected:
		std::string m_name;
};

//�������ĵ�ͼ�����������ж��ٴ����������Դ����ͼ����
typedef std::map<std::string, Sensor*> SensorMap;

}; //end namespace

#endif

