; --------------------------------------------------
; Esempi di utilizzo delle istruzioni SSE
; --------------------------------------------------
; F. Angiulli
;

%include "sseutils.nasm"

section .data			; Sezione contenente dati inizializzati


dataset		equ		8


c		equ		12


distance		equ		16

starti		equ		20
startj		equ		24
dimension		equ		28
n			equ		32
b 			equ		36

dim		equ		4
p		equ		4
UNROLL		equ		4
BLOCKSIZE	equ		32


section .bss			; Sezione contenente dati non inizializzati

section .text			; Sezione contenente il codice macchina

global	colDistance32OptimizedBlock

colDistance32OptimizedBlock:

		push		ebp							; salva il Base Pointer
		mov			ebp, esp					; il Base Pointer punta al Record di Attivazione corrente
		push		ebx							; salva i registri da preservare
		push		esi
		push		edi

		mov		eax, [ebp+starti]	; i = starti
		imul	eax, dim		; i= i*4

		mov		ebx, [ebp+startj]		; j = startj
		imul	ebx,dim		;j=j*4

		mov 	edx,[ebp+dimension]		;d
		imul	edx,dim		;d*4

		mov 	edx,[ebp+b]		;b
		add		edx,BLOCKSIZE		;b+blocksize
		imul	edx,dim		;(b+blocksize)*4

		;movaps		xmm0, [C+eax+edi]	; c0 = C[i..i+p-1][j]  = C[dim*(i+j*n)..dim*(i+j*n+p-1)]		
		mov		ecx, [ebp+b]			; k = b
		imul	ecx,dim					;b*4

		xorps xmm4,xmm4		; xmm2=0;
		xorps xmm5,xmm5		; xmm2=0;
		xorps xmm6,xmm6		; xmm2=0;
		xorps xmm7,xmm7		; xmm2=0;
		;xorps xmm7,xmm7

fork:		
		;printregps xmm7
		mov 		esi,[ebp+c]		;centroids
		mov 		edi,[ebp+dimension]		;d
		imul		edi,ebx		; 4*j*d
		add			esi, edi		; centroids + 4*j*d
		movss		xmm3, [ecx+esi]	; C[j][k] = C[4*k+4*j*d]

		mov 		esi,[ebp+dataset]		;dataset
		mov			edi,[ebp+n]		;n
		imul		edi, ecx		; 4*k*n
		add			esi, edi		;dataset + 4*k*n

		movaps		xmm0, [eax+esi]	; DS[i..i+p-1][k] = DS[4*i+4*k*n..4*(i+*k*n+p-1)]
		movaps		xmm1, [eax+esi+16]	;
		movaps		xmm2, [eax+esi+32]	;
		;movaps		xmm3, [eax+esi+48]	;
		;printregps  xmm0


	
		


		shufps		xmm3, xmm3, 0
		;printregps  xmm1
		subps		xmm0, xmm3
		subps		xmm1, xmm3
		subps		xmm2, xmm3
		;subps		xmm3, xmm3

		mulps		xmm0, xmm0		;tmp[i..i+p-1] rispetto al j-esimo centroide
		mulps		xmm1, xmm1
		mulps		xmm2, xmm2
		;mulps		xmm3, xmm3
		;sqrtps		xmm0, xmm0

		addps		xmm4, xmm0		; distance[i..i+p-1] += tmp
		addps		xmm5, xmm1
		addps		xmm6, xmm2
		;addps		xmm7, xmm3

		movaps		xmm0, [eax+esi+48]	;
		subps		xmm0, xmm3
		mulps		xmm0, xmm0
		addps		xmm7, xmm0

		
		;printregps xmm2

		add			ecx, dim		; k++
		;addps   	xmm7,[inizio]
		cmp			ecx, edx		; (k < dimension) ?
		jb			fork
		
		mov			esi,[ebp+distance]
		movaps		[esi], xmm4	; dist[i..i+p-1]  = distance[i..i+p-1] 
		movaps		[esi+16],xmm5
		movaps		[esi+32],xmm6
		movaps		[esi+48],xmm7


		;printregps		xmm2
fine:
		pop	edi									; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp							; ripristina lo Stack Pointer
		pop	ebp									; ripristina il Base Pointer
		ret										; torna alla funzione C chiamante