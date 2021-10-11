/*     
 * File:   cel_timer.c -- see headers for details and function's comments
 * Project CRESTA (see details on https://cresta-project.eu)
 * CRESTA Exascale library
 * Send you email about the bugs to
 * The High Performance Computing Center Stuttgart (HLRS) of the University of Stuttgart
 * Dmitry Khabi,  (email: khabi@hlrs.de)
 * Created on March 25, 2013
*/
//#include "Timer.h"

#include <stdio.h>
#ifndef __USE_GNU
#define __USE_GNU
#endif
#include <sched.h>
#include <errno.h>
#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 600
#endif
#include <string.h>
#ifndef __USE_POSIX199309
#define __USE_POSIX199309
#include <time.h>
#undef __USE_POSIX199309
#else
#include <time.h>
#endif
//#include "cel_timer.h"
#include <stdlib.h>
#include <unistd.h>
#ifdef VTRACE
#include "vt_user.h"
#endif

#ifdef CRAYPAT
#include  <pat_api.h>
#endif





#ifdef CRAY_COMPILER
extern  long getrdtsc(void);
#else
long getrdtsc(void)
{
  unsigned hi, lo;
  __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
  return (long)(( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 ));
}
#endif

#ifdef CRAY_COMPILER
extern  long nops(void);
#else
inline long nops(void)
{
 #if defined (__i386__)
   __asm__ __volatile__
   (
          "nop \n\t"
          "nop \n\t"
          "nop \n\t"
          "nop \n\t"
          "nop \n\t"
          "nop \n\t"
          "nop \n\t"
          "nop \n\t"
          "nop \n\t"
          "nop \n\t"
   );
 #elif defined (__x86_64__)
   __asm__ __volatile__
   (
          "nop \n\t"
          "nop \n\t"
          "nop \n\t"
          "nop \n\t"
          "nop \n\t"
          "nop \n\t"
          "nop \n\t"
          "nop \n\t"
          "nop \n\t"
          "nop \n\t"
   );
 #else
 #warning "Architecture not yet supported in ASM"
 #endif
   return (long)10L;
}
#endif
  //gcc -O2 -S cel_timer_cray.c and than link your application: craycc main.c cel_timer_cray.s

double get_time(void)
{
    struct timespec time_now;
    double result_sec_now;
 
    clock_gettime(CLOCK_REALTIME,&time_now);
    result_sec_now=(time_now.tv_sec)+(time_now.tv_nsec)*1e-9;
    return result_sec_now;
}

void get_time_stamp(long long* seconds,long long* nanoseconds)
{
    struct timespec time_now;

    clock_gettime(CLOCK_REALTIME,&time_now);
    *seconds=time_now.tv_sec;
    *nanoseconds=time_now.tv_nsec;
}

double calcTscCycle(double accuracy,double maximumDuration, int output_on) {

  double time_start, time_end;
  long tick_start, tick_end;

  time_start = get_time();
  tick_start = getrdtsc();
  usleep(100000);
  tick_end = getrdtsc();
  time_end = get_time();

  return (time_end-time_start)/(double)(tick_end-tick_start);
}


double get_tsc_cycle(double *accuracy,double *maximumDuration,int* output_on)
{
  double tsc_cycle;
  tsc_cycle = calcTscCycle(*accuracy,*maximumDuration,*output_on);
  if(*output_on>0)
  {
      printf("get_tsc_cycle: accuracy %g, maximumDuration %g\n",*accuracy,*maximumDuration);
      printf("get_tsc_cycle: tsc_cycle %g\n",tsc_cycle);
  }
  return tsc_cycle;
}


long craypat_off(void)
{
#ifdef CRAYPAT
  PAT_record (PAT_STATE_OFF);
#endif
 return 0;
}

long craypat_on(void)
{
#ifdef CRAYPAT

  PAT_record(PAT_STATE_ON);
#endif
 return 0;
}


long vtrace_on(void)
{
#ifdef VTRACE
 //printf("VTRACE ON\n");
 //VT_ON();
#endif
 return 0;
}

long vtrace_off(void)
{
#ifdef VTRACE
 //printf("VTRACE OFF\n");
 //VT_OFF();
#endif
 return 0;
}

void vt_start_1(void)
{
#ifdef VTRACE
 //VT_USER_START("1");
#endif
}

void vt_end_1(void)
{
#ifdef VTRACE
 //VT_USER_END("1");
#endif
}

void micro_sleep(int micro_seconds)
{
  usleep(micro_seconds);
}

int cel_free_int_c(void* array)
{
  free(array);
  	return 0;
}

int cel_free_double_c(void* array)
{
  free(array);
  	return 0;
}


//set thread policy for  thread
//id of core is theThreadId
int setThreadPolicy_sched(int theThreadId)
{
    int err;
    cpu_set_t processor_mask;
    
    CPU_ZERO(&processor_mask);
    CPU_SET(theThreadId,&processor_mask);
    err = sched_setaffinity( 0, sizeof(cpu_set_t), &processor_mask );
    if(err<0)
    {
      printf("setThreadPolicy_sched error:%s\n",strerror(errno));
    }
    return err;
}

