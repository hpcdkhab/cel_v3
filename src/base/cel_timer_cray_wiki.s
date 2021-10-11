        .text
        .p2align 4,,15
  .globl getrdtsc
        .type   getrdtsc, @function
  getrdtsc:
  .LFB2:
        rdtsc
        movl    %eax, %ecx
        movq    %rdx, %rax
        salq    $32, %rax
        mov     %ecx, %edx
        orq     %rdx, %rax
        ret
  .LFE2:
        .size   getrdtsc, .-getrdtsc
        .section        .eh_frame,"a",@progbits
