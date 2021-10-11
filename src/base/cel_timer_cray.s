	.file	"cel_timer_cray.c"
	.text
	.p2align 4,,15
.globl getrdtsc
	.type	getrdtsc, @function
getrdtsc:
.LFB2:
#APP
# 8 "cel_timer_cray.c" 1
	rdtsc
# 0 "" 2
#NO_APP
	movl	%eax, %ecx
	movq	%rdx, %rax
	salq	$32, %rax
	mov	%ecx, %edx
	orq	%rdx, %rax
	ret
.LFE2:
	.size	getrdtsc, .-getrdtsc
	.section	.eh_frame,"a",@progbits
.Lframe1:
	.long	.LECIE1-.LSCIE1
.LSCIE1:
	.long	0x0
	.byte	0x1
	.string	"zR"
	.uleb128 0x1
	.sleb128 -8
	.byte	0x10
	.uleb128 0x1
	.byte	0x3
	.byte	0xc
	.uleb128 0x7
	.uleb128 0x8
	.byte	0x90
	.uleb128 0x1
	.align 8
.LECIE1:
.LSFDE1:
	.long	.LEFDE1-.LASFDE1
.LASFDE1:
	.long	.LASFDE1-.Lframe1
	.long	.LFB2
	.long	.LFE2-.LFB2
	.uleb128 0x0
	.align 8
.LEFDE1:
	.ident	"GCC: (SUSE Linux) 4.3.4 [gcc-4_3-branch revision 152973]"
	.section	.comment.SUSE.OPTs,"MS",@progbits,1
	.string	"Ospwg"
	.section	.note.GNU-stack,"",@progbits
