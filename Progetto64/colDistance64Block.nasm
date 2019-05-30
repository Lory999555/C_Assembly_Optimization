; --------------------------------------------------
; Esempi di utilizzo delle istruzioni SSE
; --------------------------------------------------
; F. Angiulli
;

%include "sseutils64.nasm"

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
p		equ		8
UNROLL		equ		4
BLOCKSIZE	equ		16


section .bss			; Sezione contenente dati non inizializzati

section .text			; Sezione contenente il codice macchina

global	colDistance32Block


colDistance32Block:

	push		ebp							; salva il Base Pointer
	mov			ebp, esp					; il Base Pointer punta al Record di Attivazione corrente
	push		ebx							; salva i registri da preservare
	push		esi
	push		edi

	mov		eax, [ebp+starti]	; i = starti
	imul	eax, dim		; i= i*4

	mov		ebx, [ebp+startj]		; j = startj
	imul	ebx,dim		;j=j*4

	mov 	edx,[ebp+b]		;b
	add		edx,BLOCKSIZE		;b+blocksize
	imul	edx,dim		;(b+blocksize)*4

	;movaps		xmm0, [C+eax+edi]	; c0 = C[i..i+p-1][j]  = C[dim*(i+j*n)..dim*(i+j*n+p-1)]		
	mov		ecx, [ebp+b]			; k = b
	imul	ecx,dim					;b*4

	xorps xmm2,xmm2		; xmm2=0;
	;xorps xmm7,xmm7

fork:		
	;printregps xmm7
	mov 		esi,[ebp+dataset]		;dataset
	mov			edi,[ebp+n]		;n
	imul		edi, ecx		; 4*k*n
	add			esi, edi		;dataset + 4*k*n

	movaps		xmm0, [eax+esi]	; DS[i..i+p-1][k] = DS[4*i+4*k*n..4*(i+*k*n+p-1)]
	;printregps  xmm0


	mov 		esi,[ebp+c]		;centroids
	mov 		edi,[ebp+dimension]		;d
	imul		edi,ebx		; 4*j*d
	add			esi, edi		; centroids + 4*j*d

	
	movss		xmm1, [ecx+esi]	; C[j][k] = C[4*k+4*j*d]
	


	shufps		xmm1, xmm1, 0
	;printregps  xmm1
	subps		xmm0, xmm1
	mulps		xmm0, xmm0		;tmp[i..i+p-1] rispetto al j-esimo centroide
	;sqrtps		xmm0, xmm0
	addps		xmm2, xmm0		; distance[i..i+p-1] += tmp
	;printregps xmm2

	add			ecx, dim		; k++
	;addps   	xmm7,[inizio]
	cmp			ecx, edx		; (k < (b+blocksize)*4) ?
	jb			fork
	
	mov			esi,[ebp+distance]
	addps		xmm2,[esi]			;sommo quello che c'era con i miei risultati e li scrivo
	movaps		[esi], xmm2	; dist[i..i+p-1]  = distance[i..i+p-1] 
	;printregps		xmm2

	pop	edi									; ripristina i registri da preservare
	pop	esi
	pop	ebx
	mov	esp, ebp							; ripristina lo Stack Pointer
	pop	ebp									; ripristina il Base Pointer
	ret										; torna alla funzione C chiamante
