#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"

#include "spline.h" // for data interpolation

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
  return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

  double closestLen = 100000; //large number
  int closestWaypoint = 0;

  for(int i = 0; i < maps_x.size(); i++)
  {
    double map_x = maps_x[i];
    double map_y = maps_y[i];
    double dist = distance(x,y,map_x,map_y);
    if(dist < closestLen)
    {
      closestLen = dist;
      closestWaypoint = i;
    }

  }

  return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

  int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

  double map_x = maps_x[closestWaypoint];
  double map_y = maps_y[closestWaypoint];

  double heading = atan2((map_y-y),(map_x-x));

  double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
    if (closestWaypoint == maps_x.size())
    {
      closestWaypoint = 0;
    }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
  int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

  int prev_wp;
  prev_wp = next_wp-1;
  if(next_wp == 0)
  {
    prev_wp  = maps_x.size()-1;
  }

  double n_x = maps_x[next_wp]-maps_x[prev_wp];
  double n_y = maps_y[next_wp]-maps_y[prev_wp];
  double x_x = x - maps_x[prev_wp];
  double x_y = y - maps_y[prev_wp];

  // find the projection of x onto n
  double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
  double proj_x = proj_norm*n_x;
  double proj_y = proj_norm*n_y;

  double frenet_d = distance(x_x,x_y,proj_x,proj_y);

  //see if d value is positive or negative by comparing it to a center point

  double center_x = 1000-maps_x[prev_wp];
  double center_y = 2000-maps_y[prev_wp];
  double centerToPos = distance(center_x,center_y,x_x,x_y);
  double centerToRef = distance(center_x,center_y,proj_x,proj_y);

  if(centerToPos <= centerToRef)
  {
    frenet_d *= -1;
  }

  // calculate s value
  double frenet_s = 0;
  for(int i = 0; i < prev_wp; i++)
  {
    frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
  }

  frenet_s += distance(0,0,proj_x,proj_y);

  return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
  int prev_wp = -1;

  while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
  {
    prev_wp++;
  }

  int wp2 = (prev_wp+1)%maps_x.size();

  double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
  // the x,y,s along the segment
  double seg_s = (s-maps_s[prev_wp]);

  double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
  double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

  double perp_heading = heading-pi()/2;

  double x = seg_x + d*cos(perp_heading);
  double y = seg_y + d*sin(perp_heading);

  return {x,y};

}


//---------------------------------------------------------------------------
//[0,1,2] --> [left lane, middle lane, right lane]
int lane = 1; //start @ middle lane
int current_lane = lane; //used to track changes in lanes

//Defining Target Ref. Velocity, close to speed limit
double ref_vel =  0; // 49.5 mph



//---------------------------------------------------------------------------


int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
    istringstream iss(line);
    double x;
    double y;
    float s;
    float d_x;
    float d_y;
    iss >> x;
    iss >> y;
    iss >> s;
    iss >> d_x;
    iss >> d_y;
    map_waypoints_x.push_back(x);
    map_waypoints_y.push_back(y);
    map_waypoints_s.push_back(s);
    map_waypoints_dx.push_back(d_x);
    map_waypoints_dy.push_back(d_y);
  }



  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
    uWS::OpCode opCode) {
      // "42" at the start of the message means there's a websocket message event.
      // The 4 signifies a websocket message
      // The 2 signifies a websocket event
      //auto sdata = string(data).substr(0, length);
      //cout << sdata << endl;
      if (length && length > 2 && data[0] == '4' && data[1] == '2')
      {

        auto s = hasData(data);

        if (s != "") {
          auto j = json::parse(s);

          string event = j[0].get<string>();

          if (event == "telemetry") {
            // j[1] is the data JSON object

            // Main car's localization Data
            double car_x = j[1]["x"];
            double car_y = j[1]["y"];
            double car_s = j[1]["s"];
            double car_d = j[1]["d"];
            double car_yaw = j[1]["yaw"];
            double car_speed = j[1]["speed"];

            // Previous path data given to the Planner
            auto previous_path_x = j[1]["previous_path_x"];
            auto previous_path_y = j[1]["previous_path_y"];
            // Previous path's end s and d values
            double end_path_s = j[1]["end_path_s"];
            double end_path_d = j[1]["end_path_d"];

            // Sensor Fusion Data, a list of all other cars on the same side of the road.
            auto sensor_fusion = j[1]["sensor_fusion"];

            json msgJson;



            //----------------------------------------------------------------------------
            // define the acutal x,y points used in the Planner
            vector<double> next_x_vals;
            vector<double> next_y_vals;


            int prev_size = previous_path_x.size(); //the last path the car was following

            if(prev_size > 0){ car_s = end_path_s; }

            bool too_close = false;
            bool CONSIDER_CHANGE_LANE = false;


            //------------------------------------------------------------------------------



            // Collision Avoidance with the car infront of us by Sensor Fusion
            // check if the car is in the same lane & if it is close enough
            //------------------------------------------------------------------
            //possible lane boundries
            //-----------------------
            int center_Lane_B1 = 4; //2+4*lane -2;
            int center_Lane_B2 = 8; //2+4*lane +2;

            int right_Lane_B1 = 8;  //2+4*(lane+1) -2;
            int right_Lane_B2 = 12; //2+4*(lane+1)+2;

            int left_Lane_B1 = 0;   //2+4*(lane-1) -2;
            int left_Lane_B2 = 4;   //2+4*(lane-1)+2;

            // Parameters to make the code readable
            //---------------------------------------
            bool OTHER_CAR_IN_CENTER = false, OTHER_CAR_IN_RIGHT_LANE = false, OTHER_CAR_IN_LEFT_LANE=false ;
            bool CHANGE_RIGHT = false, CHANGE_LEFT = false;

            bool CLOSE_2_CAR_IN_FRONT = false, CLOSE_2_CAR_BEHIND = false  ;;

            int LOCATION_DIFFERENCE = 0; // detailed below



            for(int i=0; i<sensor_fusion.size(); i++)
            {
              // loop over all sensor data from nearby cars and estimate their location and speed

              //find ref_velcocity
              //get S and d for all other detected cars
              //---------------------------------------
              float d = sensor_fusion[i][6]; // check the lane of car[i]
              //every lane is bounded by [B1,B2] representing left and right boundries consequently
              //check in which lanes the detected cars are in. A car might be in my lane, right of left if exists.

              //get nearby car properties
              double vx_Next = sensor_fusion[i][3];
              double vy_Next = sensor_fusion[i][4];
              double check_speed_Next = sqrt(vx_Next*vx_Next + vy_Next*vy_Next); //veclocity magnigude of car[i]
              double check_car_s_Next = sensor_fusion[i][5]; //the s value of car[i]

              //check what this car will behave in future by using previous speed to project s value out
              check_car_s_Next += ((double)prev_size *0.02*check_speed_Next);

              //Next, determine where this car is relative to my position to take further action.

              //difference between my S and other detected cars S.
              int delta_S = check_car_s_Next - car_s;
              int distance_threshold = 30; //initial distance threshold to judge if the car is close by. increase if you warry a lot.


              int other_cars_lane;  //used to track the changes in my lanes

              bool CHECK_CENTER = d < center_Lane_B2 && d > center_Lane_B1;
              bool CHECK_RIGTH  = d < right_Lane_B2 && d > right_Lane_B1;
              bool CHECK_LEFT   = d < left_Lane_B2 && d > left_Lane_B1;

              if(CHECK_LEFT){other_cars_lane = 0; }
              if(CHECK_CENTER){other_cars_lane = 1; }
              if(CHECK_RIGTH){other_cars_lane = 2; }



              /*
              Differene between my Lane and other cars lane
                              |OTHER CARS LANE|
                              | 0  |  1 | 2 |
                          ======================
                          |0  | 0  | -1 | -2 |
              [MY LANE]   |1  | 1  |  0 | -1 |
                          |2  | 2  |  1 |  0 |

              if(diff == 0): other cars are in my lane
              if(diff == -1): {I am @ 0--> move(1) || I am @ 1--> move {0,2}  i.e. moving right  :)
              if(diff == 1):  {Iam @ 2--> move(1) || I am at 1--> move{0,2}}  i.e. moving left :)
              if(diff > |1|): don't care for the moment, but could be used as heuristic in future
              */

              CLOSE_2_CAR_IN_FRONT = delta_S  < distance_threshold ; //to account for situation where two cars are in front and next to each other
              CLOSE_2_CAR_BEHIND = delta_S  > -distance_threshold ;

              bool CLOSE_DIST = CLOSE_2_CAR_IN_FRONT && CLOSE_2_CAR_BEHIND;


              LOCATION_DIFFERENCE = lane - other_cars_lane;


              switch(LOCATION_DIFFERENCE) // considering the states given my position
              {
                case (-1):
                {
                  if(CLOSE_DIST){  OTHER_CAR_IN_RIGHT_LANE = true;  }
                }
                break;//end of case lane -1

              case (0):
              {

                if(delta_S > 0 && delta_S < 30)  //the car is close to the one infront
                {  CONSIDER_CHANGE_LANE = true;  }

              }
              break;// end case 0

              case (1):
              {
                if(CLOSE_DIST){ OTHER_CAR_IN_LEFT_LANE = true; }
              }
              break;//end case 1

              default:
              {
                continue; // for the time being not interested in case 2 or -2: coult be used for future heuristics
              }

            }//end of swtich cases over all Location difference cases

        } // end of looping over all sensor fusion data


        /*
        If lane to be changed 4 cases can be distinguished depending on the current lane

                OTHER_CAR_IN_LEFT_LANE   | OTHER_CAR_IN_RIGHT_LANE |         CHANGE           |
                ------------------------------------------------------------------------------|
                        false            |          false          |      Right or Left       |
                        false            |          true           |           Left           |
                        true             |          false          |           Right          |
                        true             |           true          |         Slow Down        |

        */

        switch(lane) // considering the states given my position
        {
          case (0):
          {
              //left most lane. consider changing to right if necessary ..
              //i.e. there is car infront and it is safe to change.
              // Considering all cases in the first row in the above matrix. (0,-1,-2)

              if(!OTHER_CAR_IN_RIGHT_LANE) // no other cars in center
              {
                  CHANGE_RIGHT = true;  CHANGE_LEFT = false;
              }

          }
          break; //end of case lane 0

          case (1):
          {

              if(!OTHER_CAR_IN_LEFT_LANE) //No cars in left lane
              {
                CHANGE_LEFT = true; CHANGE_RIGHT = false;
              }
              else if(!OTHER_CAR_IN_RIGHT_LANE ) // no cars in right lane
              {
                CHANGE_RIGHT = true; CHANGE_LEFT = false;
              }

          }
          break; //end of case lane 1

          case (2):
          {
            //right most lane: lane=2
            if(!OTHER_CAR_IN_LEFT_LANE) //No cars in center
            {
              CHANGE_LEFT = true; CHANGE_RIGHT = false;
            }

          }
          break; //end case of right most lane.

        }//end of swtich cases over all near by cars in different lanes


        /*for debugging*: */
        // cout << "CONSIDER_CHANGE_LANE...."<<CONSIDER_CHANGE_LANE
        // << "\n =======================================" << endl;
        // cout <<" OTHER_CAR_IN_LEFT_LANE: " << OTHER_CAR_IN_LEFT_LANE << " OTHER_CAR_IN_CENTER: " << OTHER_CAR_IN_CENTER <<
        // " OTHER_CAR_IN_RIGHT_LANE: " << OTHER_CAR_IN_RIGHT_LANE <<  "\nCHANGE_LEFT: " << CHANGE_LEFT << " CHANGE_RIGHT: " << CHANGE_RIGHT <<
        // "\nCLOSE_2_CAR_IN_FRONT " << CLOSE_2_CAR_IN_FRONT << " CLOSE_2_CAR_BEHIND " << CLOSE_2_CAR_BEHIND <<
        // "\nSAFE_DIST_TO_CHANGE_RIGTH: "<< SAFE_DIST_TO_CHANGE_RIGTH << " SAFE_DIST_TO_CHANGE_LEFT: " << SAFE_DIST_TO_CHANGE_LEFT <<
        // endl << "-------------------------------------------" << endl ;

        if(CONSIDER_CHANGE_LANE)
        {

          cout << "Hello... Now Decide: " << " CHANGE_RIGHT: " << CHANGE_RIGHT << ", CHANGE_LEFT: " << CHANGE_LEFT << endl;
          cout << "*******************************************************" << endl;

          if(lane >2 || lane <0)
            cout << lane << " Something is Wrong: " << endl;
          //decrease speed a little

          ref_vel -= .1; // might need to slow down before changing

          //OTHER_CAR_IN_RIGHT_LANE
          if(CHANGE_RIGHT) // just in case if our car decided to go further
          {
            cout << "CHANGE RIGHT" << endl;
            lane = lane +1;
            CHANGE_RIGHT = false;

          }
          else if(CHANGE_LEFT)
          {
              cout << "CHANGE LEFT" << endl;
            lane = lane-1;
            CHANGE_LEFT = false;

          }else
          {
            //lower ref. velocity so not to crash into the car infront,
            ref_vel -= .224; // around 5m/s2 under the 10 requirment for the annoying jerk :)

          }
        }
        else if(ref_vel< 49.5)
        {
          CHANGE_LEFT = false; CHANGE_RIGHT = false;
          ref_vel += .224;
        }


        // Generating a Smooth PATH
        //------------------------------------------------------------------
        // create a sparse list of points (x,y) view points. They are evenly spaced @30m
        // Points in between will be later interploated peise-wise fashion

        vector<double> ptsx;
        vector<double> ptsy;

        //Keep track of the reference state (x,y,yaw)
        //either ref. the starting point of the car or its previous path end point

        double ref_x = car_x;
        double ref_y = car_y;
        double ref_yaw = deg2rad(car_yaw);

        // 1) Define a starting reference
        //================================
        //If the previous size is empty, use the car as the starting ref.
        if(prev_size < 2){

          //use two points that make the path tangent to the angle of the car.
          //check where the car is and go back in time based on its angle

          double prev_car_x = car_x - cos(car_yaw);
          double prev_car_y = car_y - sin(car_yaw);

          ptsx.push_back(prev_car_x);
          ptsx.push_back(car_x);

          ptsy.push_back(prev_car_y);
          ptsy.push_back(car_y);

        }
        // utilize previous pathpoints as a starting ref.
        else{

          //use previous path end points and find the slope
          ref_x = previous_path_x[prev_size -1];
          ref_y = previous_path_y[prev_size -1];

          double ref_x_prev = previous_path_x[prev_size -2];
          double ref_y_prev = previous_path_y[prev_size -2];

          ref_yaw = atan2(ref_y- ref_y_prev, ref_x - ref_x_prev);

          // use two points that make the path tangent to the prevoius path's end points
          ptsx.push_back(ref_x_prev);
          ptsx.push_back(ref_x);

          ptsy.push_back(ref_y_prev);
          ptsy.push_back(ref_y);


        }

        //2)
        //=====================================
        // Add 30m evenly spaced points ahead of the starting ref. using Frent Representation, later
        // the spaceing between those points will be filled by spline interpolation.
        // d is set to 2+4*lane. if we are in the center, lane is 1 and d = 6. The reasoning is below

        vector<double> next_w0 = getXY(car_s+30, (2+4*lane), map_waypoints_s, map_waypoints_x,map_waypoints_y);
        vector<double> next_w1 = getXY(car_s+60, (2+4*lane), map_waypoints_s, map_waypoints_x,map_waypoints_y);
        vector<double> next_w2 = getXY(car_s+90, (2+4*lane), map_waypoints_s, map_waypoints_x,map_waypoints_y);


        ptsx.push_back(next_w0[0]); ptsy.push_back(next_w0[1]);
        ptsx.push_back(next_w1[0]); ptsy.push_back(next_w1[1]);
        ptsx.push_back(next_w2[0]); ptsy.push_back(next_w2[1]);

        for(int i=0; i< ptsx.size(); i++){

          //shift car ref. angle to zero degrees
          double shift_x = ptsx[i]-ref_x;
          double shift_y = ptsy[i]-ref_y;

          /* Rotation
          | cos  -sin  |   | shift_x |    | ptsx[i] |
          | sin   cos  | x | shift_y | =  | ptsy[i] |
          */
          ptsx[i] = (shift_x*cos(0-ref_yaw) - shift_y*sin(0-ref_yaw));
          ptsy[i] = (shift_x*sin(0-ref_yaw) + shift_y*cos(0-ref_yaw));

        }

        //create a spline
        tk::spline s;

        // set x,y points to the spline
        s.set_points(ptsx,ptsy);


        // start with all of the previous pathpoints from last times
        for(int i=0; i< previous_path_x.size(); i++){

          next_x_vals.push_back(previous_path_x[i]);
          next_y_vals.push_back(previous_path_y[i]);

        }

        //Calculate how spline points are broken up so that we match the desired velocity
        //choose a target x-point and get its corrsoponding y from the spline.
        //Calculate the distance between the target point and car location using pythogoras call it d.
        // split this d into N points and return those N points that would approximate their corrosponding points
        // on the car curve while maintaining original velocity. such that. d/N = Car_Velocity * time

        double target_x = 30.0;
        double target_y = s(target_x);
        double target_dist = sqrt((target_x)*(target_x) + (target_y)*(target_y));

        double x_add_on =0;

        //instead of creating the path everytime, just add points to it
        //fill up the rest of the path planner after filling it with previous points, output 50 pts
        for(int i =1; i<= 50 - previous_path_x.size();i++)
        {

          double N = target_dist/(0.02*ref_vel/2.24); // t*v = d/N see above: divide by 2.24 to convert from mile/h --> m/s
          double x_point = x_add_on + (target_x)/N;
          double y_point = s(x_point); // corrsoponding d value :)

          x_add_on = x_point; // move to the second point

          double x_ref = x_point;
          double y_ref = y_point;

          // Correct the rotation again to move to global coordinates
          x_point = (x_ref*cos(ref_yaw) - y_ref*sin(ref_yaw));
          y_point = (x_ref*sin(ref_yaw) + y_ref*cos(ref_yaw));

          x_point += ref_x;
          y_point += ref_y;

          next_x_vals.push_back(x_point);
          next_y_vals.push_back(y_point);


        }


        //------------------------------------------------------------------




        //------------------------------------------------------------------
        // TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
        //------------------------------------------------------------------

        //Since the car moves 50 times a second, a distance of 0.5m per move will create a velocity of 25 m/s.
        // 25 m/s is close to 50 MPH.
        // Acceleration = delta(avg(speed))/0.2s
        // jerk = (average acceleration)/ 1.s .
        //In order for the passenger to have an enjoyable ride both the jerk and the total acceleration should not exceed 10 m/s^2.
        //Total acceleration contains a normal component, (AccN) which measures the centripetal acceleration from turning.
        //The tighter and faster a turn is made, the higher the AccN value will be.

        // double dist_inc = 0.5; // set the points 0.5 m apart
        // for(int i = 0; i < 50; i++)
        // {
        //       //To stay on the lane, get Frenet s,d coordinates
        //       double next_s = car_s + (i+1)*dist_inc; // i+1 so the car will move and not set still
        //
        //       // We are in the middle Lane, points are measured from the yellow line in the middle
        //       //of the road. i.e. the care is 1.5 way from the points. Since lanes are 4m wide ==> d = 4*1.5 = 6
        //       double next_d = 6;
        //
        //       vector<double> xy = getXY(next_s, next_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);
        //
        //       next_x_vals.push_back(xy[0]);
        //       next_y_vals.push_back(xy[1]);
        //
        //       // // move with the car angle at a constant velocity in a straight line
        //       // next_x_vals.push_back(car_x+(dist_inc*i)*cos(deg2rad(car_yaw)));
        //       // next_y_vals.push_back(car_y+(dist_inc*i)*sin(deg2rad(car_yaw)));
        // }


        //------------------------------------------------------------------



        msgJson["next_x"] = next_x_vals;
        msgJson["next_y"] = next_y_vals;

        auto msg = "42[\"control\","+ msgJson.dump()+"]";

        //this_thread::sleep_for(chrono::milliseconds(1000));
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

      }
    } else {
      // Manual driving
      std::string msg = "42[\"manual\",{}]";
      ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
    }
  }
});

// We don't need this since we're not using HTTP but if it's removed the
// program
// doesn't compile :-(
h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
  size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
    char *message, size_t length) {
      ws.close();
      std::cout << "Disconnected" << std::endl;
    });

    int port = 4567;
    if (h.listen(port)) {
      std::cout << "Listening to port " << port << std::endl;
    } else {
      std::cerr << "Failed to listen to port" << std::endl;
      return -1;
    }
    h.run();
  }
