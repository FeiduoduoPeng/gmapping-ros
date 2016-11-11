#ifndef RANGESENSOR_H
#define RANGESENSOR_H

#include <vector>
#include <gmapping/sensor/sensor_base/sensor.h>
#include <gmapping/utils/point.h>

namespace GMapping{

/*激光传感器类*/
class RangeSensor: public Sensor{
	friend class Configuration;
	friend class CarmenConfiguration;
	friend class CarmenWrapper;
	public:
		/*激光束*/
		struct Beam{
			OrientedPoint pose;	//pose relative to the center of the sensor
			double span;	//spam=0 indicates a line-like beam
			double maxRange;	//maximum range of the sensor
			double s,c;		//sinus and cosinus of the beam (optimization);
		};	
		/*构造函数 给激光命令*/
		RangeSensor(std::string name);
		
		/**/
		RangeSensor(std::string name, unsigned int beams, double res, const OrientedPoint& position=OrientedPoint(0,0,0), double span=0, double maxrange=89.0);
		
		/*返回所有的激光数据(常量)*/
		inline const std::vector<Beam>& beams() const {return m_beams;}
		
		/*返回所有的激光数据*/
		inline std::vector<Beam>& beams() {return m_beams;}
		
		/*返回激光传感器的位姿*/
		inline OrientedPoint getPose() const {return m_pose;}
		
		/*更新激光查询表*/
		void updateBeamsLookup();
		bool newFormat;
	protected:
		OrientedPoint m_pose;		//激光传感器的位姿
		std::vector<Beam> m_beams;  //激光数据
};

};

#endif
