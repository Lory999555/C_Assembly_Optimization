; --------------------------------------------------
; Macro per il linking di file oggetto
; NASM con la C standard library e 
; alcune utility per l'uso di SSE
;
; Uso:
;	%include "sseutils.nasm"
;	(da specificare nel file sorgente NASM
;	prima della section .data)
; --------------------------------------------------
; Fabrizio Angiulli, 4/6/2011
;

extern	printf

section	.bss

dbuf:	resq	1

section	.data

dmask:	db		'%f ',0
cr:		db		10,0
br1:	db		'( ',0
br2:	db		')',10,0
;------------roba messa da me-----
imask:	db		'%d',0
align 16
xmmtemp: db 0.0, 0.0, 0.0, 0.0
;------------fine roba messa da me---

%macro	start	0
		push	ebp
		mov	ebp, esp
		pushad
%endmacro

%macro	stop	0
		popad
		mov	esp, ebp
		pop	ebp
		ret
%endmacro

%macro	prints	1
		pushad
		push	%1
		call	printf
		add	esp, 4
		popad
%endmacro

%macro	dprint	1
		pushad
		mov		eax,[%1+4]
		push	eax
		mov		eax,[%1]
		push	eax
		push	dmask
		call	printf
		add		esp, 12
		popad
%endmacro

%macro	sprint	1
		finit
		fld		dword [%1]
		fst		qword [dbuf]
		dprint	dbuf
%endmacro

%macro	printss	1
		sprint	%1
		prints	cr
%endmacro

%macro	printps	2
		prints	br1
		push	edx
		push	ecx
		mov	edx, %1
		mov	ecx, %2
%%loopps:
		sprint	edx
		sprint	edx+4
		sprint	edx+8
		sprint	edx+12
		add	edx, 16
		dec	ecx
		jnz	%%loopps
		pop	ecx
		pop	edx
		prints	br2
%endmacro

%macro	vprintps	2
		prints	br1
		push	edx
		push	ecx
		mov		edx, %1
		mov		ecx, %2
%%loopps:
		sprint	edx
		sprint	edx+4
		sprint	edx+8
		sprint	edx+12
		sprint	edx+16
		sprint	edx+20
		sprint	edx+24
		sprint	edx+28
		add		edx, 32
		dec		ecx
		jnz		%%loopps
		pop		ecx
		pop		edx
		prints	br2
%endmacro

%macro	printsd	1
		dprint	%1
		prints	cr
%endmacro

%macro	printpd	2
		prints	br1
		push	edx
		push	ecx
		mov		edx, %1
		mov		ecx, %2
%%looppd:
		dprint	edx
		dprint	edx+8
		add		edx, 16
		dec		ecx
		jnz		%%looppd
		pop		ecx
		pop		edx
		prints	br2
%endmacro

%macro	vprintpd	2
		prints	br1
		push	edx
		push	ecx
		mov		edx, %1
		mov		ecx, %2
%%looppd:
		dprint	edx
		dprint	edx+8
		dprint	edx+16
		dprint	edx+24
		add		edx, 32
		dec		ecx
		jnz		%%looppd
		pop		ecx
		pop		edx
		prints	br2
%endmacro
;------------roba messa da me-----

%macro printregps 1
		vpushax
        movaps [xmmtemp], %1
        printps xmmtemp, 1
		vpopax
%endmacro

%macro printregpd 1
		vpushax
        movapd [xmmtemp], %1
        printpd xmmtemp, 1
		vpopax
%endmacro

%macro	vpush	1
	sub	esp, 16
	movups	[esp], %1
%endmacro

%macro	vpop	1
	movups	%1, [esp]
	add	esp, 16
%endmacro

%macro vpushax 0
	vpush	xmm0
	vpush	xmm1
	vpush	xmm2
	vpush	xmm3
	vpush	xmm4
	vpush	xmm5
	vpush	xmm6
	vpush	xmm7
	
%endmacro

%macro vpopax 0
	
	vpop	xmm7
	vpop	xmm6
	vpop	xmm5
	vpop	xmm4
	vpop	xmm3
	vpop	xmm2
	vpop	xmm1
	vpop	xmm0
%endmacro

;------------fine roba messa da me-----