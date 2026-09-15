/** \file
 * 
 * \brief   Header-only time measuring tools.
 * \details All declarations and implementations of tictoc are made here to allow for inlining.
 * 
 * \author  Andreas Sommer  <code@andreas-sommer.eu>
 * \date    2016-2019
 * 
 */



#ifndef TICTOC_HEADER
#define TICTOC_HEADER


#include <cstdio>
#include <chrono>



/**
 * \namespace tictoc
 * \brief   Namespace of tictoc
 * \details Contains everything that belongs to tictoc.
 */
namespace tictoc {



// macros for easy access
#define TICTOC_START(x)  int ticid_ ## x = tictoc::tic();
#define TICTOC_STOP(x) double time_ ## x = tictoc::toc( ticid_ ## x ); 
#define TICTOC_PRINT(x) printf(#x " took: %15.9f milliseconds.\n", time_ ## x);



// types for duration and time point
typedef std::chrono::high_resolution_clock clock;
typedef std::chrono::duration<double, std::milli> duration;
typedef std::chrono::time_point<clock> timepoint;


static const double zerotime = duration::zero().count();  ///< Zero time as error flag

// locals for circular timing with IDs
static const int capacity = 100;            ///< maximum number of concurrent tictocs
static int       index = 1;                 ///< index (=ID) of last started tictoc
static timepoint timestamps[capacity];      ///< array keeping the starting timestamps


// forward declarations
static inline void print_time_ms(const double time_in_ms, const char* message);
static inline int tic();
static inline int tic(int id);
static inline double toc();
static inline double toc(const int id);
static inline double toctoc();
static inline double toctoc(const int id);
static inline double toctoc(const char* message);
static inline double toctoc(const int id, const char* message);



// implementation (in header, due to inlining)

/**
*  \brief Prints a message plus a time to stdout.
*  \param message is the message.
*  \param time_in_ms is the time in milliseconds.
*/ 
static inline
void print_time_ms(const double time_in_ms, const char* message) {
   printf("%s: %g milliseconds.\n", message, time_in_ms);
}


/**
*  \brief   Stars a new timer.
*  \return  Returns the ID (=index) of the new started timer.
*/ 
static inline
int tic() {
   index++;
   if (index >= capacity) {
      index = 1;
   }
   timestamps[index] = clock::now();
   return index;
}


/**
*  \brief  Stars a new timer with specified ID.
*  \param  id is the ID of the timer to be newly started.
*          If the ID is invalid, timer ID 0 is used.
*  \return Returns the ID (=index) of the new started timer.
*          If the ID is invalid, timer ID 0 is used and returned.
*/ 
static inline
int tic(int id) {
   if (id < 0  ||  id >= capacity) {
      id = 0;
   }
   timestamps[id] = clock::now();
   return id;
}

/**
*  \brief  Returns the time since timer with specified ID was started.
*  \param  id is the ID of the timer. 
           If the ID is invalid, zero is returned.
*  \return Returns the time elapsed in milliseconds since the specified timer was startd.
*/ 
static inline
double toc(const int id) {
   timepoint now = clock::now();
   if (id < 0  ||  id >= capacity) {
      return zerotime;
   }
   duration dur = (now - timestamps[id]);
   return dur.count();
}


/**
*  \brief  Returns the time since the last timer was started.
*  \return Returns the time elapsed in milliseconds since the last timer was startd.
*/ 
static inline
double toc() {
   return toc(index);
}

/**
*  \brief  Emits the specified message followed by the time elapsed since the last timer was started.
*  \param  message is the message to be displayed.
*  \return Returns the time elapsed in milliseconds since the last timer was startd.
*/ 
static inline
double toctoc(const char* message) {
   return toctoc(index, message);
}


/**
*  \brief  Emits a default message followed by the time elapsed since the specified timer was started.
*  \param  id specifies the timer id. 
*  \return Returns the time elapsed in milliseconds since the last timer was startd.
*/ 
static inline
double toctoc(const int id) {
   const double t = toc(id);
   print_time_ms(t, "Elapsed time");
   return t;
}

/**
*  \brief  Emits a default message containing the time elapsed since last timer was started.
*  \param  id specifies the timer id. 
*  \return Returns the time elapsed in milliseconds since the last timer was startd.
*/ 
static inline
double toctoc() {
   return toctoc(index);
}


/**
*  \brief  Emits the specified message followed by the time elapsed since the specified timer was started.
*  \param  id specifies the timer id. 
*  \param  message to be displayed.
*  \return Returns the time elapsed in milliseconds since the last timer was startd.
*/
static inline 
double toctoc(const int id, const char* message) {
   const double t = toc(id);
   print_time_ms(t, message);
   return t;
}




} // end of namespace




#endif // include guard