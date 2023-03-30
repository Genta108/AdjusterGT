/*************************************/
/*                                   */
/*          AdjusterGT.cpp           */
/*  Incomplement coordination file   */
/*                                   */
/*************************************/

/******************** Coding Style in this code *********************/
/*
  # Indent tab: 2 space
  # Naming style
  - PascalCase: class
  - snake_case: method, function, variable, file
  - camelCase: local constant
  - SNAKE_CASE: macro, global constant
*/
/********************************************************************/

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string.h>
#include <sys/stat.h>
#include <iomanip>

using std::cout;
using std::endl;
using std::string;
using std::array;
using std::stringstream;
using std::ifstream;
using std::ofstream;
using std::ostringstream;
using std::ios_base;

#define JUMP_SW 1
#define STAY_SW 1
#define PARTS_SW 1
#define ORIGIN_SW 1
#define DISPLAY 1

#define PARTS 8
#define COLUMN 25
#define FRAME_RATE 24 //15or24
#define JUMP_FIX_MAX 10  //dep. cycle
#define PINCH_FIX_MAX 100 //dep. cycle
#define PARTS_FIX_MAX 1
#define ROLL_FIX_MAX 10 //dep. cycle
#define FRAMEMAX 43155 //486001, 43155(hakataya), 43155(8prince)
#define JUMPDET_WD FRAME_RATE*3 //sec. the number of records in the window for jump detection
#define JUMPMOD_WD FRAME_RATE //half window for modification of jumping
#define STAYMOD_WD FRAME_RATE*3 //half window for modification of staying(*CYCLE = 45000)
#define STAYMOD_WD_MIN FRAME_RATE*100
#define PARTSMOD_WD FRAME_RATE
#define LIKELI_UP 0.1
#define LIKELI_DOWN 1
#define PRECISION 1

//define mode
#define JUMPING 1
#define STAYING_PINCH 2
#define STAYING_ROLL 3
#define ORIGIN_PATCH 4
#define OUTPARTS 5

#define VCLIM 30 //velocity coefficient
#define DRLIM 300

namespace{
  string TERM;
  string FILENAME;
  string GROUP; // = "panda";
  string DATE; // = "220112";
  string PID; // = "B2";
  int CYCLE; // = 1000;
  int SAVE_POINT; // = 0
  int SKIP = 10000;
  double COMP_THRE = 0.8;

  double velocity_overlimit[PARTS] = {300, 300, 300, 300, 300, 300, 300, 300}; //nose, head, earR, earL, neck, body1, body2, tailbase
  double velocity_underlimit[PARTS] = {3, 3, 3, 3, 3, 3, 3, 3}; //nose, head, earR, earL, neck, body1, body2, tailbase
  double outlier_area[4] = {10, 1910, 10, 1070};
  int parts_order[PARTS] = {6, 4, 3, 2, 1, 7, 0, 5}; //order of typical likelihood; body2, neck, earR, earL, head, tailbase, nose, body1
  double distance_uplimit = 300;
  double distance_lwlimit = 1;
}

//========================================== Output directry =========================================//
//string access_route = "/TEST/"; //my computer in domitory
string ACCESS_ROUTE = "/mnt/data_complex/syncbox/rat_analysis/";  //Ryzen PC
//string ACCESS_ROUTE = "/mnt/4TBSandisk/rat_analysis/panda18/";  //Ryzen PC/
//string access_route = "/Users/toyagenta/Dropbox/B03_PopEco/Rat/rat_analysis/";  //Ryzen PC
//====================================================================================================//



//================================ csv file input/output ================================//
class CSVio{
private:
  string filename;
  string comp_file_name;

public:
  void reading_csv(double **coords_data);
  void writing_csv(double **coords_data, int save_point);
};
//=====================================================================================================//



//============================== jumping detection and completion class ===============================//
class OutlierComp{
private:
  int fix_num, fix_max;
  double avelikeli;
  double vclim, drlim;

  double return_d;
  double distance[PARTS], distance_pre[PARTS], distance_post[PARTS];

  void jump_detection(int p, int f, double **coords_data);
  void jump_correction(int p, int f, double **coords_data, double **original_data);
  void stay_correction(int p, int f, double **coords_data);
  void roll_compensation(int p, int f, double **coords_data);
  void outparts_detection(int p, int f, double **coords_data);
  void outparts_correction(int p, int f, double **coords_data);

  void display_result(int p);

public:
  double complemented_x, complemented_y, complemented_l;
  double complemented_px, complemented_py, complemented_pl;
  double cycle;
  int complement, incomplement;
  void coords_mod(double **coords_data, double **original_data, CSVio csvio, int mode);
};
//=====================================================================================================//



//#################################################### main ###################################################//
int main(int argc, char *argv[]){
  //command line argument
  char cname[1024];
  int opc;
  strcpy(cname, argv[0]);
  for(opc = 1; opc < argc; ++opc){
    char *s = argv[opc];
    int lgopt = strlen(s);
    if(strcmp(s, "-pid") == 0){
      PID = argv[opc+1];
      opc++;
    }
    else if(strcmp(s, "-group") == 0){
      GROUP = argv[opc+1];
      opc++;
    }
    else if(strcmp(s, "-cycle") == 0){
      CYCLE = atoi(argv[opc+1]);
      opc++;
    }
    else if(strcmp(s, "-sp") == 0){
      SAVE_POINT = atoi(argv[opc+1]);
      opc++;
    }
    else if(strcmp(s, "-date") == 0){
      DATE = argv[opc+1];
      opc++;
    }
    else if(strcmp(s, "-term") == 0){
      TERM = argv[opc+1];
      opc++;
    }
    else if(strcmp(s, "-f") == 0){
      FILENAME = argv[opc+1];
      opc++;
    }
  }

  //using class
  CSVio csvio;
  OutlierComp comp;

  //coordination matrix
  double **coords_data;
  coords_data = new double*[COLUMN];
  for(int c = 0; c < COLUMN; c++) coords_data[c] = new double[FRAMEMAX];
  double **original_data;
  original_data = new double*[COLUMN];
  for(int c = 0; c < COLUMN; c++) original_data[c] = new double[FRAMEMAX];

  csvio.reading_csv(coords_data); //data file reading
  memcpy(original_data, coords_data, sizeof(coords_data));

  int i;
  for(i = SAVE_POINT; i < CYCLE; ++i){
    comp.cycle = i+1;
    if(JUMP_SW) comp.coords_mod(coords_data, original_data, csvio, JUMPING);  //modification of jumping outlier
    //if((comp.complement == 0) && (comp.incomplement == 0)) break;
    if(STAY_SW) comp.coords_mod(coords_data, original_data, csvio, STAYING_PINCH);  //modification of staying outlier
    if(comp.cycle > CYCLE/2 && PARTS_SW) comp.coords_mod(coords_data, original_data, csvio, OUTPARTS);  //modification of body parts outlier
    if(STAY_SW) comp.coords_mod(coords_data, original_data, csvio, STAYING_ROLL);  //modification of staying outlier
    if((i+1)%SKIP == 0) csvio.writing_csv(coords_data, i+1);
    //cout << "==========" << endl;
  }
  csvio.writing_csv(coords_data, i);  //data file writing

  delete coords_data;
  delete original_data;
}
//###############################################################################################################//


//============================== coordination modification class ===============================//
void OutlierComp::coords_mod(double **coords_data, double **original_data, CSVio csvio, int mode){
  for(int o = 0; o < PARTS; ++o){
    int p = parts_order[o];  //body1, body2, neck, earR, earL, head, tailbase, nose
    if(mode == JUMPING){
      fix_max = cycle;
      if(fix_max > JUMP_FIX_MAX) fix_max = JUMP_FIX_MAX;
    }else if(mode == STAYING_PINCH){
      fix_max = 10*cycle;
      if(fix_max > PINCH_FIX_MAX) fix_max = PINCH_FIX_MAX;
    }else if(mode == STAYING_ROLL){
      fix_max = cycle;
      if(fix_max > ROLL_FIX_MAX) fix_max = ROLL_FIX_MAX;
    }else if(mode == OUTPARTS){
      fix_max = PARTS_FIX_MAX;
      int ro = PARTS-1 - o;
      p = parts_order[ro];  //body1, body2, neck, earR, earL, head, tailbase, nose
    }

    fix_num = 0;
    complement = 0, incomplement = 1;
    avelikeli = 0;
    vclim = VCLIM; //if VCLIM = 30
    drlim = DRLIM-(200*cycle/CYCLE); //if DRLIM = 300
    while((avelikeli < PRECISION) && (fix_num < fix_max)){
      avelikeli = 0; complement = 0; incomplement = 0;
      for(int f = 0; f < FRAMEMAX; ++f){
        if(mode == JUMPING){
          jump_detection(p, f, coords_data);
        }else if(mode == OUTPARTS){
          outparts_detection(p, f, coords_data);
        }else{
          break;
        }
      }
      for(int f = 0; f < FRAMEMAX; ++f){
        complemented_x = coords_data[3*p+1][f], complemented_y = coords_data[3*p+2][f], complemented_l = coords_data[3*p+3][f]; //default value
        if(coords_data[3*p+3][f] < COMP_THRE){
          if(mode == JUMPING){
            complemented_px = coords_data[3*p+1][f-1], complemented_py = coords_data[3*p+2][f-1], complemented_pl = coords_data[3*p+3][f-1]; //default value
            jump_correction(p, f, coords_data, original_data);
            coords_data[3*p+1][f-1] = complemented_px;
            coords_data[3*p+2][f-1] = complemented_py;
            coords_data[3*p+3][f-1] = complemented_pl;
          }else if(mode == STAYING_PINCH){
            stay_correction(p, f, coords_data);
          }else if(mode == STAYING_ROLL){
            roll_compensation(p, f, coords_data);
          }else if(mode == OUTPARTS){
            outparts_correction(p, f, coords_data);
          }
          coords_data[3*p+1][f] = complemented_x;
          coords_data[3*p+2][f] = complemented_y;
          coords_data[3*p+3][f] = complemented_l;
        }
        avelikeli += coords_data[3*p+3][f];

        //outlier position
        if(coords_data[3*p+1][f] < outlier_area[0]) coords_data[3*p+3][f] = 0;
        if(coords_data[3*p+1][f] > outlier_area[1]) coords_data[3*p+3][f] = 0;
        if(coords_data[3*p+2][f] < outlier_area[2]) coords_data[3*p+3][f] = 0;
        if(coords_data[3*p+2][f] > outlier_area[3]) coords_data[3*p+3][f] = 0;
        if(coords_data[3*p+1][f] < 0) coords_data[3*p+1][f] *= -1;
        if(coords_data[3*p+2][f] < 0) coords_data[3*p+2][f] *= -1;
      }
      avelikeli /= FRAMEMAX;
      ++fix_num;

      if(((int(cycle)-1)%10 == 0)&&(fix_num%10 == 0)){
        if(mode == JUMPING){
          //cout << "#1. JUMPING MOD >>" << endl;
          if(DISPLAY) display_result(p);
        }else if(mode == STAYING_PINCH){
          //cout << "#2. PINCHING MOD >>" << endl;
        }else if(mode == STAYING_ROLL){
          //cout << "#3. ROLLING MOD >>" << endl;
        }else if(mode == OUTPARTS){
          //cout << "#4. OUTPARTS DET >>" << endl;
        }
      }
      if(complement == 0 && mode != OUTPARTS) break;
    }
  }
}

//--//

void OutlierComp::jump_detection(int p, int f, double **coords_data){
  double v = 0, integral_v = 0, step_v = 0; //velocity of rat

  //calc. average velocity
  double integral_n = 0, stable = 0, w = 0;

  int jump = 0;
  for(int t = JUMPDET_WD; t > 0; t--){
    w = (JUMPDET_WD-double(t)+1)/JUMPDET_WD;
    if(f-t-1 > 0){
      if((coords_data[3*p+3][f-t] > COMP_THRE) && (coords_data[3*p+3][f-t-1] > COMP_THRE)){
        v = sqrt(pow(coords_data[3*p+1][f-t] - coords_data[3*p+1][f-t-1], 2)+pow(coords_data[3*p+2][f-t]-coords_data[3*p+2][f-t-1], 2)); //euclidean distance
        if(velocity_underlimit[p] < v < velocity_overlimit[p]){
          integral_v += v*w;
          integral_n += w;
          ++stable;
        }else if(v < velocity_underlimit[p]){ //static state
          integral_v += 1*w;
          integral_n += w;
          ++stable;
        }else{
          jump = 1;
        }
      }else{
        break;
      }
    }else{
      continue;
    }
    if(jump) break;
  }
  if(integral_n){
    integral_v /= integral_n;
    if(integral_v > velocity_overlimit[p]) integral_v = velocity_overlimit[p];

    //calc. step velocity
    step_v = sqrt(pow(coords_data[3*p+1][f] - coords_data[3*p+1][f-1], 2)+pow(coords_data[3*p+2][f]-coords_data[3*p+2][f-1], 2)); //euclidean distance
    if((integral_n > FRAME_RATE) && (step_v > vclim*integral_v)){
      coords_data[3*p+3][f] *= (vclim*integral_v)/step_v;
      coords_data[3*p+3][f-1] -= LIKELI_DOWN*(stable/JUMPDET_WD);
      if(coords_data[3*p+3][f-1] < 0) coords_data[3*p+3][f-1] = 0;
    }else{
      return;
    }
    if(step_v > velocity_overlimit[p]){
      coords_data[3*p+3][f] *= velocity_overlimit[p]/step_v;
      coords_data[3*p+3][f-1] -= LIKELI_DOWN*(stable/JUMPDET_WD);
      if(coords_data[3*p+3][f-1] < 0) coords_data[3*p+3][f-1] = 0;
    }
  }
}

//--//

void OutlierComp::jump_correction(int p, int f, double **coords_data, double **original_data){
  //Origin patch
  if(int(cycle+1)%(CYCLE/5) == 0 && (original_data[3*p+3][f] > COMP_THRE)){
    double midpoint_x, midpoint_y, original_d;
    midpoint_x = (coords_data[3*p+1][f-1] + coords_data[3*p+1][f+1])/2;
    midpoint_y = (coords_data[3*p+2][f-1] + coords_data[3*p+2][f+1])/2;
    original_d = sqrt(pow(original_data[3*p+1][f] - midpoint_x, 2)+pow(original_data[3*p+2][f]-midpoint_y, 2)); //euclidean distance
    if(original_d < drlim/10){
      complemented_x = original_data[3*p+1][f];
      complemented_y = original_data[3*p+2][f];
      complemented_l = COMP_THRE;
      ++complement;
      return;
    }
  }

//## linear complement ##//
  for(int t = 1; t < JUMPMOD_WD+1; ++t){
    //pre information complement
    if(f-1 > 0 && f+t < FRAMEMAX){
      return_d = sqrt(pow(coords_data[3*p+1][f-1] - coords_data[3*p+1][f+t], 2)+pow(coords_data[3*p+2][f-1]-coords_data[3*p+2][f+t], 2)); //euclidean distance
      if(return_d < drlim){
        complemented_x = (coords_data[3*p+1][f-1]+coords_data[3*p+1][f+t])/2;
        complemented_y = (coords_data[3*p+2][f-1]+coords_data[3*p+2][f+t])/2;
        complemented_l = (coords_data[3*p+3][f-1]+coords_data[3*p+3][f+t])/2;
        complemented_px = (coords_data[3*p+1][f-1]+complemented_x)/2;
        complemented_py = (coords_data[3*p+2][f-1]+complemented_y)/2;
        complemented_pl = (coords_data[3*p+3][f-1]+complemented_l)/2;
        ++complement;
        return;
      }
    }
    if(t == JUMPMOD_WD){
      if(f-JUMPMOD_WD < 0){++incomplement; return;}

      //## velocity complement ##//
      double v, vx, vy, integral_vx = 0, integral_vy = 0; //velocity of rat
      //calc. average velocity
      for(int t = JUMPMOD_WD; t > 0; t--){
        if((coords_data[3*p+3][f-t] > COMP_THRE) && (coords_data[3*p+3][f-t-1] > COMP_THRE)){
          vx = coords_data[3*p+1][f-t]-coords_data[3*p+1][f-t-1];
          vy = coords_data[3*p+2][f-t]-coords_data[3*p+2][f-t-1];
          v = sqrt(pow(vx, 2)+pow(vy, 2)); //euclidean distance
          if(velocity_underlimit[p] < v < velocity_overlimit[p]){
            integral_vx += vx;
            integral_vy += vy;
          }else if(v < velocity_underlimit[p]){ //static state
            integral_vx += 0;
            integral_vy += 0;
          }else{
            ++incomplement;
            return;
          }
        }else{
          ++incomplement;
          return;
        }
      }
      integral_vx /= JUMPMOD_WD;
      integral_vy /= JUMPMOD_WD;

      //calc. step movement
      complemented_x = coords_data[3*p+1][f-1]+integral_vx;
      complemented_y = coords_data[3*p+2][f-1]+integral_vy;
      complemented_l = COMP_THRE;
    }
  }
}

//--//

void OutlierComp::stay_correction(int p, int f, double **coords_data){
  //long time range complement
  int success = 0;
  double t_max = STAYMOD_WD*(CYCLE-(cycle-1));
  if(t_max < STAYMOD_WD_MIN) t_max = STAYMOD_WD_MIN;
  double return_d;
  for(int t = 1; t < t_max+1; ++t){
    if(f+t < FRAMEMAX && f-1 >= 0){
      return_d = sqrt(pow(coords_data[3*p+1][f+t] - coords_data[3*p+1][f-1], 2)+pow(coords_data[3*p+2][f+t]-coords_data[3*p+2][f-1], 2)); //euclidean distance
      if((return_d < drlim) && (coords_data[3*p+3][f+t] >= COMP_THRE) && (coords_data[3*p+3][f-1] >= COMP_THRE)){
        complemented_x = (coords_data[3*p+1][f-1]+coords_data[3*p+1][f+t])/2;
        complemented_y = (coords_data[3*p+2][f-1]+coords_data[3*p+2][f+t])/2;
        complemented_l = (coords_data[3*p+3][f-1]+coords_data[3*p+3][f+t])/2;
        ++complement;
        return;
      }
    }
    if(t == STAYMOD_WD){
      ++incomplement;
      return;
    }
  }
}

//--//

void OutlierComp::outparts_detection(int p, int f, double **coords_data){
  //parts distance check
  for(int pair = 0; pair < 8; ++pair){ //max pd = 7
    if(p == pair) continue;
    distance[pair] = sqrt(pow(coords_data[3*pair+1][f] - coords_data[3*p+1][f], 2)+pow(coords_data[3*pair+2][f]-coords_data[3*p+2][f], 2)); //euclidean distance
    if(distance[pair] > distance_uplimit || distance[pair] < distance_lwlimit){
      if(coords_data[3*pair+3][f] > COMP_THRE){
        coords_data[3*p+3][f] -= LIKELI_DOWN*(cycle/CYCLE);
        if(coords_data[3*p+3][f] < 0) coords_data[3*p+3][f] = 0;
      }
    }
  }
}

//--//

void OutlierComp::outparts_correction(int p, int f, double **coords_data){
  //## linear complement ##//
  for(int t = 1; t < PARTSMOD_WD+1; ++t){
    int best_method = 0;
    double point_max = 0, point = 0;

    //pre information complement
    if(f+t < FRAMEMAX && f-1 > 0){
      return_d = sqrt(pow(coords_data[3*p+1][f+t] - coords_data[3*p+1][f-1], 2)+pow(coords_data[3*p+2][f+t]-coords_data[3*p+2][f-1], 2)); //euclidean distance

      for(int pair = 0; pair < PARTS; ++pair){
        distance_pre[pair] = sqrt(pow(coords_data[3*pair+1][f-1] - coords_data[3*p+1][f-1], 2)+pow(coords_data[3*pair+2][f-1]-coords_data[3*p+2][f-1], 2)); //euclidean distance
        distance_post[pair] = sqrt(pow(coords_data[3*pair+1][f+t] - coords_data[3*p+1][f+t], 2)+pow(coords_data[3*pair+2][f+t]-coords_data[3*p+2][f+t], 2)); //euclidean distance

        if((return_d < drlim) && (distance_pre[pair] < distance_uplimit) && (distance_pre[pair] > distance_lwlimit) && (distance_post[pair] < distance_uplimit) && (distance_post[pair] > distance_lwlimit) && (coords_data[3*p+3][f+t] > COMP_THRE) && (coords_data[3*p+3][f-1] > COMP_THRE)){
          ++point;
        }
      }
    }
    if(point > point_max){
      best_method = 1;
      point_max = point;
    }
    point = 0;

    //post information complement
    if(f-t > 0 && f+1 < FRAMEMAX){
      return_d = sqrt(pow(coords_data[3*p+1][f-t] - coords_data[3*p+1][f+1], 2)+pow(coords_data[3*p+2][f-t]-coords_data[3*p+2][f+1], 2));

      for(int pair = 0; pair < PARTS; ++pair){ //max pd = 7
        distance_pre[pair] = sqrt(pow(coords_data[3*pair+1][f-t] - coords_data[3*p+1][f-t], 2)+pow(coords_data[3*pair+2][f-t]-coords_data[3*p+2][f-t], 2));
        distance_post[pair] = sqrt(pow(coords_data[3*pair+1][f+1] - coords_data[3*p+1][f+1], 2)+pow(coords_data[3*pair+2][f+1]-coords_data[3*p+2][f+1], 2));
        if((return_d < drlim) && (distance_pre[pair] < distance_uplimit) && (distance_pre[pair] > distance_lwlimit) && (distance_post[pair] < distance_uplimit) && (distance_post[pair] > distance_lwlimit) && (coords_data[3*p+3][f+1] > COMP_THRE) && (coords_data[3*p+3][f-t] > COMP_THRE)){
          ++point;
        }
      }
    }
    if(point > point_max){
      best_method = 2;
      point_max = point;
    }
    point = 0;

    if(best_method == 0){
    }else if(best_method == 1){
      complemented_x = (coords_data[3*p+1][f-1]+coords_data[3*p+1][f+t])/2;
      complemented_y = (coords_data[3*p+2][f-1]+coords_data[3*p+2][f+t])/2;
      complemented_l = coords_data[3*p+3][f] + LIKELI_UP*point_max;
      if(complemented_l > 1) complemented_l = 1;
      complement +=1;
      return;
    }else if(best_method == 2){
      complemented_x = (coords_data[3*p+1][f-t]+coords_data[3*p+1][f+1])/2;
      complemented_y = (coords_data[3*p+2][f-t]+coords_data[3*p+2][f+1])/2;
      complemented_l = coords_data[3*p+3][f] + LIKELI_UP*point_max;
      if(complemented_l > 1) complemented_l = 1;
      complement +=1;
      return;
    }
    if(t == PARTSMOD_WD){
      incomplement += 1;
    }
  }
}

//--//

void OutlierComp::roll_compensation(int p, int f, double **coords_data){
  if(coords_data[3*p+3][f-1] >= COMP_THRE){
    complemented_x = coords_data[3*p+1][f-1];
    complemented_y = coords_data[3*p+2][f-1];
    complemented_l = COMP_THRE;// * (cycle / CYCLE);
    //if(complemented_l > COMP_THRE) complemented_l = COMP_THRE;
    complement += 1;
    return;
  }else if(coords_data[3*p+3][f+1] >= COMP_THRE){
    complemented_x = coords_data[3*p+1][f+1];
    complemented_y = coords_data[3*p+2][f+1];
    complemented_l = COMP_THRE;// * (cycle / CYCLE);
    //if(complemented_l > COMP_THRE) complemented_l = COMP_THRE;
    complement += 1;
    return;
  }else{
    incomplement += 1;
    return;
  }
}

//--//

void OutlierComp::display_result(int p){
  cout << "date: " << DATE << endl;
  cout << "pid: " << PID << endl;
  cout << "cycle: " << cycle << endl;
  cout << "part: " << p << endl;
  cout << "#fix_num: " << fix_num << endl;
  cout << "#complement: " << complement << endl;
  cout << "#incomplement: " << incomplement << endl;
  cout << "#ave of likelihood: " << avelikeli << endl;
  cout << endl;
}

//=====================================================================================================//


//================================ weited moving average for smoothing ================================//
void CSVio::reading_csv(double **coords_data){
  //file loading
  ostringstream oss_filename;
  if(GROUP == "8prince"){
    oss_filename << GROUP << "/" << TERM << "/" << FILENAME;
  }else{
    oss_filename << GROUP << "/" << GROUP << DATE << "/" << GROUP << "_" << DATE << "_2k9h15fDLC_resnet50_" << GROUP << "_" << PID << "Oct16shuffle1_500000";
  }


  filename = oss_filename.str();
  ostringstream file_location;
  if(SAVE_POINT == 0){
    file_location << ACCESS_ROUTE << filename << ".csv";  //name of data file
  }else{
    file_location << ACCESS_ROUTE << filename << "_fixed" << SAVE_POINT << ".csv";  //name of data file
  }
  ifstream csvfile(file_location.str());

  if(csvfile.fail()){
    cout << "file load failed..." << endl;
    cout << file_location.str() << endl;
    std::exit(0);
  }

  string gab;
  getline(csvfile.seekg(0, ios_base::cur), gab, '\n');
  getline(csvfile.seekg(0, ios_base::cur), gab, '\n');
  getline(csvfile.seekg(0, ios_base::cur), gab, '\n');

  string column_value[COLUMN] = {
    "frame", "nose_x", "nose_y", "nose_l", "head_x", "head_y", "head_l", "earR_x", "earR_y", "earR_l", "earL_x", "earL_y", "earL_l", "neck_x", "neck_y", "neck_l",
    "body1_x", "body1_y", "body1_l", "body2_x", "body2_y", "body2_l", "tailbase_x", "tailbase_y", "tailbase_l"
  };


  for(int f = 0; f < FRAMEMAX; ++f){
    stringstream ssCoords;
    //ssCoords.precision(10);
    for(int c = 0; c < COLUMN; ++c){
      //read Value
      if(c == COLUMN-1){
        getline(csvfile.seekg(0, ios_base::cur), column_value[COLUMN-1], '\n');
      }else{
        getline(csvfile.seekg(0, ios_base::cur), column_value[c], ',');
      }
      ssCoords << column_value[c];
      ssCoords >> coords_data[c][f];
      //cout << "coords_data: " << coords_data[c][f] << endl;
      ssCoords.str("");
      ssCoords.clear(stringstream::goodbit);
    }
  }

  csvfile.close();

}

//--//

void CSVio::writing_csv(double **coords_data, int save_point){
  ostringstream complemented_data_file;
  complemented_data_file << ACCESS_ROUTE << filename << "_fixed"<< save_point << ".csv";  //name of data file
  comp_file_name = complemented_data_file.str();
  ofstream fout(comp_file_name.c_str());
  fout << "frame," << "nose_x," << "nose_y," << "nose_l," << "head_x," << "head_y," << "head_l," << "earR_x," << "earR_y," << "earR_l," << "earL_x," << "earL_y," << "earL_l," << "neck_x," << "neck_y," << "neck_l," <<
    "body1_x," << "body1_y," << "body1_l," << "body2_x," << "body2_y," << "body2_l," << "tailbase_x," << "tailbase_y," << "tailbase_l" << endl;

  for(int f = 0; f < FRAMEMAX; ++f){
    fout << f << ",";
    for(int c = 1; c < COLUMN-1; ++c){
      fout << coords_data[c][f] << ",";
    }
    fout << coords_data[COLUMN-1][f];
    fout << endl;
  }
  fout.close();

}
//=====================================================================================================//
