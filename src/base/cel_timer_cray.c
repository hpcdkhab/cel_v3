 inline long getrdtsc(void)
 {
   unsigned long long x;
 #if defined (__i386__)
   __asm__ __volatile__ ("rdtsc" : "=A" (x));
 #elif defined (__x86_64__)
   unsigned int tickl, tickh;
   __asm__ __volatile__ ("rdtsc" : "=a" (tickl), "=d" (tickh));
   x = ((unsigned long long)tickh << 32) | tickl;
 #else
 #warning "Architecture not yet supported in ASM"
 #endif
   return (long)x;
 }

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
