%include "sseutils64.nasm"

section .data			; Sezione contenente dati inizializzati


distance		equ		8


min_distance		equ		12


label		equ		16

j			equ		20
i			equ		24

dim		equ		4
p		equ		8
UNROLL		equ		4
BLOCKSIZE	equ		32


section .bss			; Sezione contenente dati non inizializzati

section .text			; Sezione contenente il codice macchina

global	distanceControl64

distanceControl64:

	push		rbp				; salva il Base Pointer
	mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
	pushaq						; salva i registri generali


	;mov 	r8,[rbp+i]		;i
	imul	r8,dim			;i*4

	mov		r13, 0			; k = 0

	

fork2:		
	;printregps xmm7
	;mov 		r11,[rbp+distance]		;distance
	

	vmovaps			ymm0, [r13+rdi]	; distance[k..k+7]
	;printregps  xmm0
	vextractf128	xmm4,ymm0,1



	;mov 		r11,[rbp+min_distance]		

	vmovaps			ymm1, [r13+rsi]			; min_distance[k..k+7]
	;printregps	xmm1
	vmovaps			ymm2,ymm0		;sarà la maschera
	vcmpnltps		ymm2,ymm1		;>= in modo da farmi dare 0 dove mi serve
	vextractf128	xmm3,ymm2,1


	extractps	rax,xmm2,0
	cmp			rax,0
	jne			if1

	mov			r11,r8			;i*4
	mov 		r10,r13			;k*4
	add 		r11,r10			;i*4+k*4
	;mov 		rax,[rbp+label]		
	;mov		rax,[r10]		;label
	;mov		r10,[rbp+j]			;j
	mov 		[rdx+r11],rcx		;label[i+k]=j
	extractps	rax,xmm0,0			;

	;mov			r10,[rbp+min_distance]		;min_distance
	mov			[rsi+r13],rax				;min_distance+4+k*4

if1:

	extractps	rax,xmm2,1
	cmp			rax,0
	jne			if2
	;printreg	rax	

	;printregps	xmm2
	mov			r11,r8			;i*4
	mov 		r10,r13			;k*4
	add 		r11,r10			;i*4+k*4
	;mov 		rax,[rbp+label]		
	;mov			rax,[r10]		;label
	;mov			r10,[rbp+j]			;j
	;printregps	xmm2
	add			r11,dim			;(i+k+1)*4
	mov 		[rdx+r11],rcx		
	extractps	rax,xmm0,1		
	;printregps	xmm2
		
	;mov			r10,[rbp+min_distance]		;min_distance
	;add			r10,dim						;min_distance+4
	mov			[rsi+4+r13],rax				;min_distance+4+k*4
	;printregps	xmm2


if2:

	extractps	rax,xmm2,2
	cmp			rax,0
	jne			if3
	
	;printregps	xmm1
	mov			r11,r8			;i*4
	mov 		r10,r13			;k*4
	add 		r11,r10			;i*4+k*4
	;mov 		rax,[rbp+label]		
	;mov			rax,[r10]		;label
	;mov			r10,[rbp+j]			;j
	add			r11,dim*2				;(i+k+2)*4
	mov 		[rdx+r11],rcx		
	extractps	rax,xmm0,2		

	;mov			r10,[rbp+min_distance]
	;add			r10,dim*2				;min_distance + (i+k+2)*4
	mov			[rsi+8+r13],rax

if3:
	extractps	rax,xmm2,3
	cmp			rax,0
	jne			if4
	;printregps	xmm1
	mov			r11,r8			;i*4
	mov 		r10,r13			;k*4
	add 		r11,r10			;i*4+k*4
	;mov 		rax,[rbp+label]		
	;mov			rax,[r10]		;label
	;mov			r10,[rbp+j]			;j

	add			r11,dim*3				;(i+k+2)*4
	mov 		[rdx+r11],rcx		
	extractps	rax,xmm0,3	

	;mov			r10,[rbp+min_distance]
	;add			r10,dim*3				;min_distance + (i+k+2)*4
	mov			[rsi+12+r13],rax

if4:

	extractps	rax,xmm3,0
	cmp			rax,0
	jne			if5

	mov			r11,r8			;i*4
	mov 		r10,r13			;k*4
	add 		r11,r10			;i*4+k*4
	;mov 		rax,[rbp+label]		
	;mov		rax,[r10]		;label
	;mov		r10,[rbp+j]			;j
	add			r11,dim*4				;(i+k+4)*4
	mov 		[rdx+r11],rcx		;label[i+k]=j
	extractps	rax,xmm4,0			;

	;mov			r10,[rbp+min_distance]		;min_distance
	mov			[rsi+16+r13],rax				;min_distance+4+k*4

if5:

	extractps	rax,xmm3,1
	cmp			rax,0
	jne			if6
	;printreg	rax	

	;printregps	xmm2
	mov			r11,r8			;i*4
	mov 		r10,r13			;k*4
	add 		r11,r10			;i*4+k*4
	;mov 		rax,[rbp+label]		
	;mov			rax,[r10]		;label
	;mov			r10,[rbp+j]			;j
	;printregps	xmm2
	add			r11,dim*5			;(i+k+1)*4
	mov 		[rdx+r11],rcx		
	extractps	rax,xmm4,1		
	;printregps	xmm2
		
	;mov			r10,[rbp+min_distance]		;min_distance
	;add			r10,dim						;min_distance+4
	mov			[rsi+20+r13],rax				;min_distance+4+k*4
	;printregps	xmm2


if6:

	extractps	rax,xmm3,2
	cmp			rax,0
	jne			if7
	
	;printregps	xmm1
	mov			r11,r8			;i*4
	mov 		r10,r13			;k*4
	add 		r11,r10			;i*4+k*4
	;mov 		rax,[rbp+label]		
	;mov			rax,[r10]		;label
	;mov			r10,[rbp+j]			;j
	add			r11,dim*6				;(i+k+2)*4
	mov 		[rdx+r11],rcx		
	extractps	rax,xmm4,2		

	;mov			r10,[rbp+min_distance]
	;add			r10,dim*2				;min_distance + (i+k+2)*4
	mov			[rsi+24+r13],rax

if7:
	extractps	rax,xmm3,3
	cmp			rax,0
	jne			fine
	;printregps	xmm1
	mov			r11,r8			;i*4
	mov 		r10,r13			;k*4
	add 		r11,r10			;i*4+k*4
	;mov 		rax,[rbp+label]		
	;mov			rax,[r10]		;label
	;mov			r10,[rbp+j]			;j
	add			r11,dim*7				;(i+k+2)*4
	mov 		[rdx+r11],rcx		
	extractps	rax,xmm4,3	

	;mov			r10,[rbp+min_distance]
	;add			r10,dim*3				;min_distance + (i+k+2)*4
	mov			[rsi+28+r13],rax
	
fine:
	add			r13, 32		; k+=p
	cmp			r13, 128		; (k < dimension) ?
	jb			fork2
	
	popaq						; ripristina i registri generali
	mov		rsp, rbp			; ripristina lo Stack Pointer
	pop		rbp				; ripristina il Base Pointer
	ret						; torna alla funzione C chiamante
