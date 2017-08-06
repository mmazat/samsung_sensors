#include "ros/ros.h"
#include "std_msgs/String.h"
#include <sensor_msgs/Imu.h>
#include <sensor_msgs/NavSatFix.h>

#include "stdUtils.h"
#include <sstream>
#include <vector>
#include <cstdlib>

#include <boost/asio.hpp>

#include "blocking_udp_client.cpp"

sensor_msgs::Imu imuMsg;
ros::Publisher imuPub; 

using boost::asio::ip::udp;

void HandleMSGLine(const std::string line)
{
	double dt;
		auto spl = split(line,',');
		double time = 0;
		for (auto it = spl.begin(); it != spl.end(); ++it)
		{		
			if (it == spl.begin())
			{
				time = std::stof(*it);
				++it;
			}
			if (it == spl.end()) break;

			int code = 0;
			if (has_only_digits(*it))
				code = std::stoi(*it++);
			else
				continue;

			switch (code)
			{
			case 1:
				
				/*gps->time = time;
				gps->rl[0] = std::stof(*it++);
				gps->rl[1] = std::stof(*it++);
				gps->rl[2] = std::stof(*it++);
				sensors->GetGPS(gps);*/
				break;

			case 3:
				imuMsg.header.stamp=ros::Time(time);										
				imuMsg.linear_acceleration.x  = std::stof(*it++);
				imuMsg.linear_acceleration.y  = std::stof(*it++);
				imuMsg.linear_acceleration.z  = std::stof(*it++);
				if (it == spl.end()) break;
				if (std::stoi(*it++) != 4)
					break;
				imuMsg.angular_velocity.x = std::stof(*it++);
				imuMsg.angular_velocity.y = std::stof(*it++);
				imuMsg.angular_velocity.z = std::stof(*it++);			
				imuPub.publish(imuMsg);
				break;
			default:
				break;
			}

			if (it == spl.end()) 
			break;

		}
	
}




int main(int argc, char *argv[])
{
	/**
	 * The ros::init() function needs to see argc and argv so that it can perform
	 * any ROS arguments and name remapping that were provided at the command line.
	 * For programmatic remappings you can use a different version of init() which takes
	 * remappings directly, but for most command-line programs, passing argc and argv is
	 * the easiest way to do it.  The third argument to init() is the name of the node.
	 *
	 * You must call one of the versions of ros::init() before using any other
	 * part of the ROS system.
	 */
	ros::init(argc, argv, "driver");

	/**
	 * NodeHandle is the main access point to communications with the ROS system.
	 * The first NodeHandle constructed will fully initialize this node, and the last
	 * NodeHandle destructed will close down the node.
	 */
	ros::NodeHandle n;

	/**
	 * The advertise() function is how you tell ROS that you want to
	 * publish on a given topic name. This invokes a call to the ROS
	 * master node, which keeps a registry of who is publishing and who
	 * is subscribing. After this advertise() call is made, the master
	 * node will notify anyone who is trying to subscribe to this topic name,
	 * and they will in turn negotiate a peer-to-peer connection with this
	 * node.  advertise() returns a Publisher object which allows you to
	 * publish messages on that topic through a call to publish().  Once
	 * all copies of the returned Publisher object are destroyed, the topic
	 * will be automatically unadvertised.
	 *
	 * The second parameter to advertise() is the size of the message queue
	 * used for publishing messages.  If messages are published more quickly
	 * than we can send them, the number here specifies how many messages to
	 * buffer up before throwing some away.
	 */
	imuPub = n.advertise<sensor_msgs::Imu>("imu", 1000);
	imuMsg.header.frame_id="imu_link";
	int port=5555; //default port for the GPS IMU app
	n.param<int>("port",port,5555);

	try
	{
	
    udp::endpoint local_endpoint = boost::asio::ip::udp::endpoint
	(udp::v4(),port);

	std::cout << "Local bind " << local_endpoint << std::endl;	
		
	client c(local_endpoint);
	int64_t count = 0;
	char data_[1024];
	boost::system::error_code ec;

	while (ros::ok())
	{			       
       std::size_t n = c.receive(boost::asio::buffer(data_), 
	   boost::posix_time::seconds(1), ec);
			
		if(n>0)
		{
			HandleMSGLine(std::string(data_,n));
			++count;
		}
		ros::spinOnce();		
	}
	std::cout<<"receieved messages "<<count<<std::endl;
	}
	catch (std::exception& e)
	{
		std::cerr << e.what() << std::endl;		
	}
	
	return 0;
}




