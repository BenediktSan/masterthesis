#include<vector>
#include<string>
#include<chrono>
#include<iostream>
#include<iomanip>

typedef double RealType;
typedef std::chrono::time_point<std::chrono::system_clock> time_point;
typedef std::chrono::duration<RealType> duration;

// TIME MEASURE CLASS 
class TimeMeasure
{ 
public:
  // CONSTRUCTORS
  TimeMeasure();

  // PUBLIC METHODS
  void measure( const std::string what );
  void stop();

  // PUBLIC MEMBERS
  time_point start_time; // construction time 
  time_point measure_time{}; // time of the previous measurement
  duration diagonalization{};

  duration total_time{}; // total time of the job
  std::time_t date_and_time{};

private:
  // PRIVATE METHODS 
  void print_measure( const std::string what, const duration time );
};


// constructor : start time measurement
TimeMeasure::TimeMeasure()
{
    start_time = std::chrono::system_clock::now();
    measure_time = start_time;
}

// PUBLIC METHODS
// method : measure the time for 'what'
void TimeMeasure::measure( const std::string what )
{
    time_point current_time = std::chrono::system_clock::now();
    duration time = current_time - measure_time;

    if( what == "diagonalization" )
    {
        diagonalization = time;
    }
    else
    {
        //std::cout << "error : measurement to TimeMeasure unknown\n";
    }
    print_measure( what, time );
    measure_time = current_time;
}

// method : stop measuring, evaluate and print the measurement results at core 0
void TimeMeasure::stop()
{   
    // measure the total time:
    measure_time = std::chrono::system_clock::now();
    total_time = measure_time - start_time;
    std::cout << "\033[1;36m" << "everything done: \033[0m" << "in total it took me " << total_time.count() << " seconds or " ;
    std::cout << std::chrono::duration_cast<std::chrono::minutes>(total_time).count() << " minutes or ";
    std::cout << std::chrono::duration_cast<std::chrono::hours>(total_time).count() << " hours\n";
    
    // estimate current time and date 
    time_point end = std::chrono::system_clock::now();
    date_and_time = std::chrono::system_clock::to_time_t(end);
}

// printing method : print the measurement results
void TimeMeasure::print_measure( const std::string what, const duration time )
{
    std::cout << "\033[1;33m" << what << " done: \033[0m" << "it took me " << time.count() << " seconds or " ;
    std::cout << std::chrono::duration_cast<std::chrono::minutes>(time).count() << " minutes or ";
    std::cout << std::chrono::duration_cast<std::chrono::hours>(time).count() << " hours\n";
}

