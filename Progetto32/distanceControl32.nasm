%include "sseutils.nasm"

section .data			; Sezione contenente dati inizializzati

align 16
distance		equ		8

align 16
min_distance		equ		12

align 16
label		equ		16

j			equ		20
i			equ		24

dim		equ		4
p		equ		4
UNROLL		equ		4
BLOCKSIZE	equ		16

align 16
inizio:		dd		1.0, 1.0, 1.0, 1.0

section .bss			; Sezione contenente dati non inizializzati

section .text			; Sezione contenente il codice macchina

global	distanceControl32

distanceControl32:

		push		ebp							; salva il Base Pointer
		mov			ebp, esp					; il Base Pointer punta al Record di Attivazione corrente
		push		ebx							; salva i registri da preservare
		push		esi
		push		edi

		mov 	edx,[ebp+i]		;i
		imul	edx,dim			;i*4



		mov		ebx, 0			;z = 0
		;imul	ebx,dim			;z*4
forz:
		mov		ecx, 0			; k = 0
fork:		
		;printregps xmm7
		mov 		esi,[ebp+distance]		;distance
		mov			edi,ebx					;z
		imul		edi,BLOCKSIZE*dim			;z*blocksize
		add			esi,edi					;distance + z*blocksize
		

		movaps		xmm0, [ecx+esi]	; distance[z + k..k+3]
		;printregps  xmm0


		mov 		esi,[ebp+min_distance]		
		;mov			edi,ebx					;z
		;imul		edi,BLOCKSIZE*dim			;z*blocksize
		;add			esi,edi					;distance + z*blocksize

		movaps		xmm1, [ecx+esi]			; min_distance[k..k+3]
		;printregps	xmm1
		movaps		xmm2,xmm0		;sarà la maschera
		cmpnltps	xmm2,xmm1		;>= in modo da farmi dare 0 dove mi serve

		extractps	eax,xmm2,0
		cmp			eax,0
		jne			if1

		mov			esi,edx			;i*4
		mov 		edi,ecx			;k*4
		add 		esi,edi			;i*4+k*4
		mov 		eax,[ebp+label]		
		;mov			eax,[edi]		;label
		mov 		edi,ebx				;z*4
		shr			edi,1
		shr			edi,1
		add			edi,[ebp+j]			;j+z
		
		mov 		[eax+esi],edi		;label[i+k]=j
		extractps	eax,xmm0,0			;

		mov			edi,[ebp+min_distance]		;min_distance
		;mov			esi,ebx					;z
		;imul		esi,BLOCKSIZE*dim			;z*blocksize
		;add			edi,esi					;distance + z*blocksize
		mov			[edi+ecx],eax				;min_distance+4+k*4 

if1:

		extractps	eax,xmm2,1
		cmp			eax,0
		jne			if2
		;printreg	eax	

		;printregps	xmm2
		mov			esi,edx			;i*4
		mov 		edi,ecx			;k*4
		add 		esi,edi			;i*4+k*4
		mov 		eax,[ebp+label]		
		;mov			eax,[edi]		;label
		mov 		edi,ebx				;z*4
		shr			edi,1
		shr			edi,1
		add			edi,[ebp+j]			;j+z


		;printregps	xmm2
		add			esi,dim			;(i+k+1)*4
		mov 		[eax+esi],edi		
		extractps	eax,xmm0,1		
		;printregps	xmm2
			
		mov			edi,[ebp+min_distance]		;min_distance
		add			edi,dim						;min_distance+4
		;mov			esi,ebx					;z
		;imul		esi,BLOCKSIZE*dim			;z*blocksize
		;add			edi,esi					;distance + z*blocksize
		mov			[edi+ecx],eax				;min_distance+4+k*4
		;printregps	xmm2


if2:

		extractps	eax,xmm2,2
		cmp			eax,0
		jne			if3
		
		;printregps	xmm1
		mov			esi,edx			;i*4
		mov 		edi,ecx			;k*4
		add 		esi,edi			;i*4+k*4
		mov 		eax,[ebp+label]		
		;mov			eax,[edi]		;label
		mov 		edi,ebx				;z*4
		shr			edi,1
		shr			edi,1
		add			edi,[ebp+j]			;j+z
		add			esi,dim*2				;(i+k+2)*4
		mov 		[eax+esi],edi		
		extractps	eax,xmm0,2		

		mov			edi,[ebp+min_distance]
		add			edi,dim*2				;min_distance + (i+k+2)*4
		;mov			esi,ebx					;z
		;imul		esi,BLOCKSIZE*dim			;z*blocksize
		;add			edi,esi					;distance + z*blocksize
		mov			[edi+ecx],eax

if3:
		extractps	eax,xmm2,3
		cmp			eax,0
		jne			fine
		;printregps	xmm1
		mov			esi,edx			;i*4
		mov 		edi,ecx			;k*4
		add 		esi,edi			;i*4+k*4
		mov 		eax,[ebp+label]		
		;mov			eax,[edi]		;label
		mov 		edi,ebx				;z*4
		shr			edi,1
		shr			edi,1
		add			edi,[ebp+j]			;j+z
		add			esi,dim*3				;(i+k+2)*4
		mov 		[eax+esi],edi		
		extractps	eax,xmm0,3	

		mov			edi,[ebp+min_distance]
		add			edi,dim*3				;min_distance + (i+k+2)*4
		;mov			esi,ebx					;z
		;imul		esi,BLOCKSIZE*dim			;z*blocksize
		;add			edi,esi					;distance + z*blocksize
		mov			[edi+ecx],eax
		
fine:
		add			ecx, dim*p		; k+=p
		cmp			ecx, p*UNROLL*dim		; (k < dimension) ?
		jb			fork

		add			ebx,dim
		cmp			ebx,BLOCKSIZE*dim
		jb 			forz
		
		pop	edi									; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp							; ripristina lo Stack Pointer
		pop	ebp									; ripristina il Base Pointer
		ret										; torna alla funzione C chiamante